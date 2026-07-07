# Sketch: one coupled variable, end-to-end, GPU-ready

> Status: **design sketch / proposal** · Companion to
> [`plan-cpp-interface.md`](plan-cpp-interface.md),
> [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md), and
> [`spec-fc-api.md`](spec-fc-api.md) §4a (the `@couple`/`@internal` tier tags and
> the `CppCoupleAccessors` / `FortranDevPtrExport` regions this sketch realizes).
> Shows how a single Tier-A coupling variable would look across *every* layer if
> we (a) keep one source of truth and (b) share device memory with ERF so the
> per-step host bounce + `streamSynchronize` disappear.
>
> Nothing here is implemented yet. This is the target shape to react to.

## 0. What "Tier-A" means here

Recall the three tiers (see `plan-cpp-interface.md`):

| Tier | Example | Needs C++ view? | Device-resident? | ERF-shared (Array4)? |
|------|---------|-----------------|------------------|----------------------|
| **A** ERF-facing coupling (~18) | `SWDOWN`, `U_PHY`, `HFX`, `TSK` | yes | yes | **yes** |
| **B** internal state | soil moisture, snow layers | only if I/O'd | yes | no |
| **C** host-only | tables, namelist | no | no (read-only push) | no |

This sketch is about **Tier A** — the only set that gets the extra ERF
device-sharing glue. Tier B just gets `!$acc enter data`; Tier C is untouched.

The worked example is `SWDOWN` (2-D, direction **in**: ERF → Noah-MP). The
output direction (`HFX`, Noah-MP → ERF) is symmetric — ERF *reads* instead of
*writes*. The 3-D `(i,k,j)` case (`U_PHY`) is in §6.

Ownership model used below: **Option B — Fortran owns the allocation**, makes it
device-resident, and *exports its device pointer* to ERF, which wraps it as an
`amrex::Array4`. This keeps the Noah-MP physics accessing the variable *exactly
as today* (a normal Fortran array); only ERF learns something new.

---

## 1. The single source of truth  ·  `NoahmpIO.H-mc`  (HAND-WRITTEN)

One line in the `@NoahmpMacro:Source m_noahmpio { … }` block. A new
`@couple(dir=…)` annotation marks it Tier-A and states the direction; everything
else is generated from it.

```cpp
// in @NoahmpMacro:Source m_noahmpio { ... }
NoahmpArray2D<noahmp_real> SWDOWN[xstart:xend, ystart:yend]  @couple(dir=in);  // solar down [W/m2]
```

- `@couple(dir=in|out|inout)` is the *only* new syntax vs. today's contract line.
- Absence of `@couple` ⇒ Tier-B: the generator emits the storage + view + the
  `enter data` residency, but **not** the ERF device-sharing accessors of §2–§3.

---

## 2. Generated C++  ·  `NoahmpIO.H`  (GENERATED — do not edit)

```cpp
class NoahmpIO_type {
    // ... @NoahmpMacro:CppStorageFields(m_noahmpio) expands to, among others:

    NoahmpArray2D<noahmp_real> SWDOWN;             // host-side view (init, I/O, restart)

    // ---- @couple(dir=in) extras, generated only for Tier-A -------------------
    noahmp_real*               SWDOWN_devptr();     // GPU address of Noah-MP's storage
    amrex::Array4<noahmp_real> SWDOWN_a4();         // that address, wrapped for ParallelFor
    // -------------------------------------------------------------------------
};
```

## 3. Generated C++ bodies  ·  `NoahmpIO.cpp`  (GENERATED)

```cpp
// Ask the Fortran side for the device address of this block's SWDOWN.
noahmp_real* NoahmpIO_type::SWDOWN_devptr() {
    return NoahmpIO_SWDOWN_devptr_fi(&fptr.level, &fptr.blkid);
}

// Wrap it as an Array4 with the block's global bounds. The slab lives at k=klo;
// z-extent is 1 so (i,j,klo) maps to Fortran (i,j) byte-for-byte.
// NOTE: Array4's `end` is EXCLUSIVE (hi+1).
amrex::Array4<noahmp_real> NoahmpIO_type::SWDOWN_a4() {
    const int klo = this->kds;               // domain z small-end (surface slab)
    return amrex::Array4<noahmp_real>(
        SWDOWN_devptr(),
        amrex::Dim3{xstart,   ystart,   klo  },   // begin (inclusive)
        amrex::Dim3{xend+1,   yend+1,   klo+1},   // end   (exclusive)
        /*ncomp=*/1);
}
```

## 4. Generated Fortran  ·  `NoahmpIO_fi.F90`, `NoahmpIOVarType.F90`, `NoahmpIOVarInitMod.F90`  (GENERATED)

```fortran
! --- NoahmpIOVarType.F90 : storage member (FortranStorageFields region) --------
real(c_kind_noahmp), allocatable :: SWDOWN(:,:)          ! GENERATED

! --- NoahmpIOVarInitMod.F90 : allocate + make device-resident ------------------
if (.not. allocated(NoahmpIO%SWDOWN)) &                  ! GENERATED
    allocate(NoahmpIO%SWDOWN(XSTART:XEND, YSTART:YEND))
!$acc enter data create(NoahmpIO%SWDOWN)                 ! GENERATED (Tier-A/B)

! --- NoahmpIO_fi.F90 : export the device address (bind C) ----------------------
function NoahmpIO_SWDOWN_devptr_fi(level, blkid) bind(C) result(p)   ! GENERATED
    use iso_c_binding, only: c_ptr, c_associated
    use openacc,       only: acc_deviceptr
    type(c_ptr)        :: p
    type(c_ptr), value :: level, blkid
    type(NoahmpIO_type), pointer :: nm
    call resolve_block(level, blkid, nm, require_init=.true.)   ! same guard as every shim
    p = acc_deviceptr(nm%SWDOWN)                                ! GPU address of resident copy
end function
```

That is the *entire* generated surface for one Tier-A variable. Everything below
is hand-written and stays small.

---

## 5. Hand-written usage

### 5a. Once, at init  ·  `ERF_NOAHMP.cpp::Init`  (HAND-WRITTEN, one-time)

Bind Noah-MP's OpenACC queue to AMReX's stream so kernels on both sides are
ordered by the hardware queue — this is what removes the per-step host sync.

```cpp
// after noahmpio_vect.resize(...), once:
acc_set_cuda_stream(NOAHMP_ACC_QUEUE, amrex::Gpu::gpuStream());   // one shared stream
```

### 5b. Every step  ·  `ERF_NOAHMP.cpp::Advance_With_State`  (HAND-WRITTEN)

The pinned buffer + `streamSynchronize` + `LoopOnCpu` transpose all go away. ERF
writes the forcing straight into Noah-MP's device array:

```cpp
NoahmpIO_type* noahmpio = &noahmpio_vect[idb];

Array4<Real>       swdown = noahmpio->SWDOWN_a4();   // device view of Noah-MP memory
const Array4<const Real>& CONS = cons_in.const_array(mfi);

// same stream as DriverMain (5c); no sync between them needed
ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
{
    swdown(i,j,klo) = /* SWDOWN forcing from ERF state */ ;   // writes Noah-MP GPU array
});

noahmpio->itimestep = m_itimestep;
noahmpio->DriverMain();          // offloaded; reads swdown on device, writes fluxes

// output-direction vars (HFX_a4(), TSK_a4(), ...) are READ here, symmetric.
```

Compare to today (`ERF_NOAHMP.cpp:339-366`): the `Gpu::streamSynchronize()`, the
`noahmp_input_tmp` pinned FArrayBox, and the `LoopOnCpu` copy into
`noahmpio->SWDOWN(i,j)` are **all deleted**.

### 5c. In the physics  ·  `NoahmpDriverMainMod.F90`  (HAND-WRITTEN — UNCHANGED indexing)

The driver loop becomes an offload region; the *access to `SWDOWN` is identical to
today* — it is a plain Fortran array that happens to be device-resident.

```fortran
!$acc parallel loop collapse(2) present(NoahmpIO) async(NOAHMP_ACC_QUEUE)
do j = jts, jte
   do i = its, ite
      ! ... pack the column ...
      noahmp%forcing%SWDOWN = NoahmpIO%SWDOWN(i,j)   ! same bytes ERF wrote in 5b
      call NoahmpMain(noahmp)                        ! !$acc routine seq
      ! ... unpack fluxes into NoahmpIO%HFX(i,j) etc. (output-direction Tier-A) ...
   end do
end do
```

---

## 6. The 3-D gotcha  ·  `U_PHY(i, k, j)`  — layer index in the MIDDLE

`U_PHY` is declared `[xstart:xend, kms:kme, ystart:yend]` and used as
`U_PHY(i,1,j)`. A vanilla `Array4(i,j,k)` puts `k` last (stride `nx*ny`); Fortran
puts it in the middle (stride `nx`). Writing `a(i,j,k)` and reading `U_PHY(i,k,j)`
would **silently corrupt**. The contract line already carries the Fortran order,
so the generator emits a **stride-matched accessor** instead of a spatial Array4:

```cpp
// contract line (HAND-WRITTEN):
NoahmpArray3D<noahmp_real> U_PHY[xstart:xend, kms:kme, ystart:yend] @couple(dir=in);

// GENERATED C++ accessor — uses the (i,k,j) offset, NOT Array4's (i,j,k):
struct U_PHY_view {
    noahmp_real* p; int xs, ys, ks, nx, nk;   // nx = xend-xstart+1, nk = kme-kms+1
    AMREX_GPU_HOST_DEVICE noahmp_real& operator()(int i, int k, int j) const {
        return p[(i-xs) + (k-ks)*nx + (j-ys)*nx*nk];   // matches Fortran column-major (i,k,j)
    }
};
U_PHY_view NoahmpIO_type::U_PHY_v();   // p = U_PHY_devptr()
```

ERF then writes in Fortran order and the two agree by construction:

```cpp
auto uphy = noahmpio->U_PHY_v();
ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
    uphy(i, /*k=*/1, j) = /* u forcing */ ;     // same element as Fortran U_PHY(i,1,j)
});
```

> Alternative for §6: since only `k=1` crosses the boundary, declare a 2-D
> coupling alias `U_PHY_cpl[xstart:xend, ystart:yend]` and drop the singleton
> layer — then it reduces to the trivial 2-D case of §2–§5. Pick per-variable.

---

## 7. Option A vs. B (who allocates) — what changes on the Noah-MP side

| | **Option B (this sketch)** | **Option A (ERF/AMReX allocates)** |
|---|---|---|
| Fortran member | `allocatable` (unchanged) | `pointer, contiguous` |
| Noah-MP access | `NoahmpIO%SWDOWN(i,j)` unchanged | same, but bound via `c_f_pointer` first |
| Residency | Fortran `!$acc enter data` | ERF arena; Fortran `deviceptr(...)` clause |
| Extra bind(C) | `*_devptr_fi` (export address) | `*_bind_fi` (receive address, remap bounds) |
| Best for | forcings + most Tier-A | outputs that reuse an existing ERF fab |

Option A binder, for reference:

```fortran
subroutine NoahmpIO_SWDOWN_bind_fi(level, blkid, cptr, ilo, ihi, jlo, jhi) bind(C)
    type(c_ptr), value    :: level, blkid, cptr
    integer(c_int), value :: ilo, ihi, jlo, jhi
    type(NoahmpIO_type), pointer :: nm
    real(c_kind_noahmp), pointer :: tmp(:,:)
    call resolve_block(level, blkid, nm, require_init=.true.)
    call c_f_pointer(cptr, tmp, [ihi-ilo+1, jhi-jlo+1])   ! 1-based...
    nm%SWDOWN(ilo:ihi, jlo:jhi) => tmp                    ! ...remap to global bounds
end subroutine
```

---

## 8. Checklist for adding a Tier-A variable under this scheme

1. Add the contract line with `@couple(dir=…)` in `NoahmpIO.H-mc` (§1).
2. `make codegen` — emits the view, `*_devptr()`, `*_a4()`/`*_v()`, storage,
   allocate, `enter data`, and `*_devptr_fi` (§2–§4).
3. In `ERF_NOAHMP.cpp`: replace the pinned-buffer copy for that variable with a
   direct `*_a4()` read/write inside the shared-stream `ParallelFor` (§5b).
4. In the physics: name it in the offload region's `present(...)` clause; the
   indexing is unchanged (§5c).
5. For I/O/restart of a device-resident var: `!$acc update host(...)` before the
   host NetCDF write (rare, off hot path).

Net vs. today: **one contract line + one `@couple` tag** replaces a pinned
buffer, a `streamSynchronize`, and two `LoopOnCpu` transposes per variable.
