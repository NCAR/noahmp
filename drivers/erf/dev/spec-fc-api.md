# Spec: Fortran ↔ C API & the `@NoahmpMacro` code generator

> Status: living document · Owns: the C++ ↔ Fortran boundary and the generator
> that builds it. · Companion: [`spec-memory-safety.md`](spec-memory-safety.md),
> [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md).

## 1. The contract

Every member of the boundary struct `NoahmpIO_type_fi` is an **opaque pointer** (a
bare memory address with no type attached — see the glossary in
[`spec-overview.md`](spec-overview.md) §7). A correct C ↔ Fortran mapping
therefore requires the *same number* of members, in the *same order*, with the
*same element precision*, and **no padding** (no hidden bytes inserted by the
compiler — the struct must stay a plain flat array of pointers). The struct
exists twice:

- **C++** (`extern "C" struct NoahmpIO_type_fi`, in `NoahmpIO.H-mc`).
- **Fortran** (`type, bind(c) :: NoahmpIO_type_fi`, in `NoahmpIO_fi.F90-mc`).

The two are emitted from **one ordered list** (the `@NoahmpMacro:Source` block),
so count and order match *by construction* — they cannot drift unless the
generator itself is wrong, which `make codegen-check` catches.

### What the boundary carries

| Kind | C++ type | Fortran owner | Who owns the storage |
|------|----------|---------------|----------------------|
| dimension/handle | `int` | `integer(C_INT)` target | C++ scalar, Fortran points back |
| real scalar | `noahmp_real` | `real(c_kind_noahmp)` target | C++ scalar, Fortran points back |
| 2-D / 3-D array | `NoahmpArray{2,3}D<noahmp_real>` | `real(...), allocatable` | **Fortran**, C++ view points in via `C_LOC` |

- **Scalars are C++-owned.** `NoahmpIOScalarInitDefault_fi` calls
  `C_F_POINTER(cptr%X, blk%X)` so the Fortran `pointer` component aliases the
  C++ scalar. The struct is thus *self-referential*: `fptr` stores the addresses
  of the C++ object's own members.
- **Arrays are Fortran-owned.** `NoahmpIOVarInitDefault` allocates them;
  `NoahmpIOVarInitDefault_fi` then `C_LOC`s each allocatable into the matching
  `fptr` handle, and `NoahmpIO_type::VarInitDefault` wraps each handle in a
  `NoahmpArray{2,3}D` view with the right bounds.

## 2. Precision: `noahmp_real` / `c_kind_noahmp`

The build defines `DOUBLE_PREC` iff ERF/AMReX is double precision, so:

```cpp
#ifdef DOUBLE_PREC
using noahmp_real = double;   // == amrex::Real
#else
using noahmp_real = float;
#endif
```

mirrors the Fortran `c_kind_noahmp` parameter (`utility/Machine.F90`,
`C_DOUBLE`/`C_FLOAT`). **Never** write a literal `double` / `C_DOUBLE` on the
boundary. Because the C++ and Fortran sides are *separately-invoked compilers*,
a `DOUBLE_PREC` mismatch is a real hazard that shared source cannot prevent — it
is the one thing the run-time guard checks (§5).

## 3. The C++ API surface (`NoahmpIO_type`)

Public methods (each maps 1:1 to a `bind(C)` `*_fi` shim in `NoahmpIO.cpp-mc`):

```cpp
void ScalarInitDefault();              // wire scalar pointers (must run first)
void VarInitDefault();                 // allocate arrays + build NoahmpArray views
void InitMain();                       // initial physics state
void ReadNamelist();  void ReadTable();
void ReadLandHeader(); void ReadLandMain();
void DriverMain();                     // advance one step
void WriteLand(int filenum);
void WriteRestart(const std::string& dir);
void ReadRestart (const std::string& dir);
```

Plus `NoahmpIO_vector` for the per-level block array (size-once `resize`, then a
read-only/element-access subset of `std::vector`). See
[`spec-memory-safety.md`](spec-memory-safety.md) for the deleted copy/move ops.

### The `bind(C)` entry points (Fortran side)

Each shim (`NoahmpIO_fi.F90-mc`) resolves `(level, blkid)` from the incoming
struct via `resolve_block`, which validates the C pointers and both indices and
aborts with a clear message instead of indexing module state out of bounds. By
default it asserts the block's scalars were already wired
(`ScalarInitDefault` passes `require_init=.false.` since it *is* the wiring).

String arguments (restart dir) cross as `(const char*, int len)` →
`character(kind=C_CHAR)(*)` copied into a `character(len=dir_len)` local; the
C++ side guards the `size_t → int` length cast against overflow
(`checked_dir_len`).

## 4. The code generator (`tools/NoahmpMacro.py`)

### Single source of truth

The member declarations of `class NoahmpIO_type`, wrapped in the
`@NoahmpMacro:Source { ... }` block of `NoahmpIO.H-mc`. **Order in the block =
ABI order.** From this one list the generator regenerates the bodies of every
`@NoahmpMacro:<Region>()` marker across five files:

| Template (`-mc`, hand-edited) | Generated target (gitignored) |
|-------------------------------|-------------------------------|
| `NoahmpIO.H-mc` | `NoahmpIO.H` |
| `NoahmpIO.cpp-mc` | `NoahmpIO.cpp` |
| `NoahmpIO_fi.F90-mc` | `NoahmpIO_fi.F90` |
| `NoahmpIOVarType.F90-mc` | `NoahmpIOVarType.F90` |
| `NoahmpIOVarInitMod.F90-mc` | `NoahmpIOVarInitMod.F90` |

From the Source list it re-derives: the mirrored `fi` struct (both sides), the
constructors and move constructor, `NOAHMP_IO_FI_NUM_MEMBERS`, the Fortran
`bind(C)` type, the coupled members of the storage `NoahmpIO_type`, the
`C_LOC`/`C_F_POINTER` wiring, the `NoahmpArray` extents, **and** the Fortran
`allocate()` of each coupled array. The allocate and the C++ view extent come
from the *same* `bounds_Nd` annotation, so they cannot disagree.

### Annotation grammar

```cpp
@NoahmpMacro:Source {
  int ids, ide, jds, jde;                 // ints: no annotation
  int numrad = 2;                         // int with C++ default (shared with Fortran)

  noahmp_real DTBL;                       // @NoahmpMacro:scalar
  noahmp_real ZLVL = -9999.0;             // @NoahmpMacro:scalar doc="..."

  NoahmpArray2D<noahmp_real> XLAT;        // @NoahmpMacro:bounds_2d(xstart:xend, ystart:yend) doc="latitude [rad]"
  NoahmpArray3D<noahmp_real> U_PHY;       // @NoahmpMacro:bounds_3d(xstart:xend, kms:kme, ystart:yend) doc="..."
}
```

Rules the generator enforces (it raises `SystemExit("NoahmpMacro: …")` otherwise):

- **Kind is inferred from the C++ type** (`int` / `noahmp_real` /
  `NoahmpArray{2,3}D`).
- **Arrays** need a `bounds_{2,3}d(lo:hi, …)` list in **Fortran order**. The rank
  must match the `NoahmpArrayND` type and the number of bound pairs. Every bound
  token must be a **literal** or a **coupled member name** (`xstart`, `xend`,
  `kms`, `kme`, `nsoil`, `numrad`, …) so C++ and Fortran resolve it to the *same*
  value — never a Fortran-only parameter. An optional `doc="…"` is kept as a
  plain comment.
- **Scalars** carry `@NoahmpMacro:scalar` and an optional `doc`.
- **Ints** need no annotation; an optional C++ default (`int numrad = 2;`) is
  shared with Fortran.

### Internal architecture (read `tools/NoahmpMacro.py` top-to-bottom)

> Skip this subsection unless you're modifying the generator itself. *Adding* a
> coupled variable needs only §1 here and
> [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md).

```
parse_source(h_text)   -> ordered [Member(name, kind, rank, begin, end, doc)]
r_*(members)           -> body LINES for one region
H_REGIONS / CPP_REGIONS / F90_REGIONS / VARTYPE_REGIONS / VARINIT_REGIONS  -> region -> renderer
_slice_markers / apply_regions  -> find a region's marker, replace with banner+body+close
render_source_block    -> strip the Source macro down to bare C++ declarations
process(members)       -> {target_path: new_text} for all five targets
main(argv)             -> write targets, or --check (exit 1 + unified diff on drift)
```

Two marker forms in templates:

- **Call form** `@NoahmpMacro:Region();` on its own line — expanded in place to a
  banner-commented body block.
- **Block form** `@NoahmpMacro:Source { … }` — the only one; its wrapper and
  per-line annotations are stripped for the compiled header (`doc` kept as a
  comment).

Region names are PascalCase (matching method naming). Rendering helpers
(`_wrap`, `_decl_lines`, `_runs`) pack several members per line; since targets
regenerate every build, density costs no diff noise.

### Build integration

- **Make**: runs the generator at parse time (`$(info NoahmpMacro codegen: …)`).
- **CMake**: runs it at configure time (`FATAL_ERROR` on failure).
- Convenience targets: `make codegen` / `cmake --build … --target noahmp_codegen`;
  `make codegen-check` / `…noahmp_codegen_check` (exit 1 on drift, for CI).

A fresh checkout needs no manual step — the targets are regenerated before
compiling. The generated targets are gitignored; only the `-mc` templates are
committed.

## 5. ABI guards (what is checked, and why each remains)

Because count/order/types match by construction, only the hazards codegen
**cannot** prevent keep a guard:

1. **Compile-time `static_assert`** (free) in `NoahmpIO.H`:
   `is_standard_layout<NoahmpIO_type_fi>` **and**
   `sizeof == NOAHMP_IO_FI_NUM_MEMBERS * sizeof(void*)`. Pins the C++ side to a
   flat pointer array; the Fortran `bind(C)` side is guaranteed the same layout by
   the Fortran↔C interoperability standard. Catches an incompatible compiler
   pairing / unexpected padding.
2. **Run-time precision check** `NoahmpIO_AssertAbi()`: compares
   `NoahmpRealSize_fi()` (Fortran's `sizeof(real(c_kind_noahmp))`) against
   `sizeof(noahmp_real)`. The struct is all pointers, so its byte size is
   identical for float vs double — *only* a value check can detect a
   `DOUBLE_PREC` mismatch between the two compilers. Runs only once (a
   function-local `static` guards it), from both `NoahmpIO_type`'s constructor and
   `NoahmpIO_vector::resize`, so neither path can skip it.

> **Do not reintroduce** a run-time size or member-order handshake — those only
> re-verified facts codegen already guarantees.

## 6. Acceptance checklist (any boundary or generator change)

- [ ] `make codegen` (or `cmake --build … --target noahmp_codegen`) regenerates
      all five targets without error.
- [ ] `make codegen-check` exits 0; run twice to confirm regeneration is
      idempotent on a clean tree.
- [ ] Project builds → the `static_assert` and `NoahmpIO_AssertAbi()` (the real
      acceptance tests) pass.
- [ ] Precision rule intact (`noahmp_real` / `c_kind_noahmp`, never literal
      `double`/`C_DOUBLE`).
- [ ] No `@NoahmpMacro:` marker survives into any compiled target.
- [ ] [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md) still works
      end to end.
