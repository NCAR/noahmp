# Spec: memory safety & runtime robustness

> Status: living document · Owns: the invariants that keep the self-referential
> coupling struct from corrupting memory, plus the parallel-safe fatal path, plus
> (§7) the device-memory rules for the GPU offload.
> Companion: [`spec-fc-api.md`](spec-fc-api.md),
> [`plan-cpp-interface.md`](plan-cpp-interface.md).

The coupling struct is **self-referential**: `NoahmpIO_type::fptr` stores the
addresses of *this object's own* scalar members (so the Fortran side can alias
them). That single fact drives almost every rule here — any operation that
relocates a `NoahmpIO_type` leaves the Fortran pointers (and the copied handles)
aimed at freed or foreign memory, which is silent corruption in a release build.

## 1. `NoahmpIO_type`: copy banned, move re-points

```cpp
NoahmpIO_type(const NoahmpIO_type&)            = delete;  // copy would alias source's members
NoahmpIO_type& operator=(const NoahmpIO_type&) = delete;
NoahmpIO_type& operator=(NoahmpIO_type&&)      = delete;
NoahmpIO_type(NoahmpIO_type&& o) noexcept;               // allowed, but re-points fptr
```

- **Copy is impossible.** A memberwise copy would leave the new object's `fptr`
  pointing at the *source's* members — dangling once the source dies.
- **Move is required but special.** `std::vector<NoahmpIO_type>::resize` needs the
  type move-constructible, so move cannot simply be deleted. The move ctor copies
  all handles (`fptr(o.fptr)`) — array handles point into Fortran-owned memory and
  are correct as-is — then **re-points the self-referential scalar handles** from
  the source's members to the destination's (`@NoahmpMacro:CppMoveRepoint`).

## 2. `NoahmpIO_vector`: size once, never relocate

`NoahmpIO_vector : private std::vector<NoahmpIO_type>` — it inherits **privately**
from `std::vector`, which hides every vector operation that could move elements in
memory; the compiler simply won't let a caller reach them. Two invariants:

1. The C++ vector size and the Fortran `NoahmpIO_vect(level)%NoahmpIO` extent must
   stay equal — only `resize(size, level)` keeps them in sync.
2. Any grow/reserve/insert/copy relocates elements → dangling `fptr`s.

Re-exported (safe once sized): `operator[]`, `at`, `front`, `back`, `data`,
`begin/end`, `cbegin/cend`, `size`, `empty`. **Everything else** (`push_back`,
`insert`, `erase`, `swap`, `clear`, `reserve`, `assign`, raw `resize`) is a
**compile error** — strictly stronger than a run-time guard.

- `resize(size, level)` is the one sizing entry point. It runs the ABI precision
  guard, refuses a second call (`initialized_`), refuses `size == 0`, sizes the
  C++ vector, then calls `NoahmpIOTypeVectInit_fi`.
- The Fortran `NoahmpIOTypeVectInit_fi` is a **second line of defense**: a repeated
  call with the same size does nothing, a conflicting re-init aborts, it validates
  `level ∈ [NLEVEL_MIN, NLEVEL_MAX]` and `NBlocks >= 1`, and it checks the
  `allocate` status.

## 3. Block resolution: validate before indexing

Every `bind(C)` entry point routes through `resolve_block` (`NoahmpIO_fi.F90-mc`)
which, before dereferencing module-global state, checks:

- `LEVEL` and `BLKID` C pointers are associated (`c_associated`);
- `level ∈ [NLEVEL_MIN, NLEVEL_MAX]`;
- `NoahmpIO_vect(level)%NoahmpIO` is allocated;
- `blkid ∈ [lbound, ubound]`;
- (unless `require_init=.false.`) the block's scalars were wired — `DTBL`
  associated is the proxy, reliable because pointer components default to
  `=> null()`.

Each failure prints a specific message and calls `NoahmpIO_abort()` rather than
indexing out of bounds (silent corruption in release).

## 4. Array bounds checking (`NoahmpArray.H`)

`NoahmpArray{1,2,3}D` are Fortran-layout (column-major, arbitrary lower bound)
views over the Fortran-owned data. Their `operator()` calls
`noahmparray_check(idx, lo, hi, dim)`:

- **Debug** (`NDEBUG` not defined): an out-of-range index prints the offending
  index/range/dimension and calls `NoahmpIO_fatal()`.
- **Release** (`-DNDEBUG`): compiled out entirely — **zero cost**.

`idx` is a `long` so a signed index below a (possibly negative) lower bound — e.g.
the negative-indexed snow-layer arrays — is still caught.

## 5. Fatal-error abstraction (`NoahmpFatal.H` / `NoahmpIO.cpp`)

Noah-MP must not depend on MPI/AMReX, so it **never calls `MPI_Abort`**. Instead:

```cpp
NoahmpIO_set_fatal_handler([](const char* msg){ amrex::Abort(msg); });  // host, once at init
NoahmpIO_fatal(msg);                                                    // anywhere in Noah-MP
```

- `NoahmpIO_fatal()` prints `msg`, calls the host handler (cross-rank propagation), then
  falls back to `std::abort()` if no handler is installed or it returns.
- The Fortran side calls `NoahmpIO_abort()` → the C++ shim `NoahmpIO_fatal_c()` →
  the same handler. A bare Fortran `error stop` would terminate only the calling
  image and deadlock peers in the next collective; the `error stop 1` after the
  shim is the serial fallback.
- **The config/init readers route through this path, not bare `stop`.**
  `NoahmpReadNamelistMod` (missing/malformed `namelist.erf`, and every option
  validation — unset `start_year/month/day`, `NSOIL < 1`, bad timesteps, urban
  levels, …), `NoahmpReadTableMod` (missing `NoahmpTable.TBL`, out-of-range
  `FILOSS`), and `NoahmpReadLandMod` (unsupported level; `ktop`/`kbottom` soil
  interpolation failures) all `call NoahmpIO_abort()` — a serial `stop` in cold
  init would strand peers at the first collective exactly as in the I/O path. The
  diagnostic `write` is rank-0-guarded so only one line prints. (The remaining
  bare `stop` is `CAL_MON_DAY`'s internal day-of-year sanity guard in
  `NoahmpDriverMainMod`, which fires only on a logic error, not on user input.)

This keeps the *entire* parallel-runtime dependency in the host, not in Noah-MP —
a precondition for the GPU port (see [`plan-cpp-interface.md`](plan-cpp-interface.md)).

## 7. Device memory & GPU coupling (forward-looking)

> Applies once the GPU offload ([`plan-cpp-interface.md`](plan-cpp-interface.md))
> lands. The host build ignores all of this — the directives are no-ops without
> offload flags. Recorded here so the safety story is designed, not retrofitted.

The offload keeps the ownership split of §1–§4 and adds a *second* place each
array lives (GPU memory) plus a *second* runtime (the Fortran offload runtime)
touching it alongside AMReX. Three rules keep that safe:

1. **Residency is generator-owned and paired with the allocate.** Each Tier-A/B
   array gets a generated `!$acc enter data create(...)` right after its
   `allocate()` (see [`spec-fc-api.md`](spec-fc-api.md) §4a). It must be matched by
   an `exit data` on teardown, and — crucially — **re-issued if the array is ever
   reallocated**: a reallocate changes the host address, so a stale device mapping
   (and any exported device pointer) dangles. The driver sizes blocks once
   (§2 here), which is what makes this tractable.

2. **Coupling is zero-copy via a shared device pointer, and correctness rests on
   one shared stream — not a host barrier.** For a Tier-A array, Fortran owns the
   allocation and exports its device address (`acc_deviceptr`, Option B); ERF
   wraps that address as an `amrex::Array4` and reads/writes it in a `ParallelFor`
   (see the worked example, [`sketch-couple-variable-gpu.md`](sketch-couple-variable-gpu.md)).
   Two runtimes writing the same GPU bytes **is a data race unless ordered**. The
   ordering is provided by binding Noah-MP's OpenACC queue to AMReX's stream
   (`acc_set_cuda_stream(queue, amrex::Gpu::gpuStream())`) so all coupling kernels
   run in enqueue order on **one** stream. This *replaces* today's per-step
   `Gpu::streamSynchronize()` in `ERF_NOAHMP_Advance.cpp`; do **not** reintroduce a host
   barrier, and do **not** put coupling kernels on a second stream without an
   explicit cross-stream dependency. (Managed/unified memory fixes *coherence* but
   **not** *ordering* — it does not remove this requirement.)

3. **Host and device copies must be reconciled around host-only operations.** I/O,
   restart, cold init, and any host-side inspection see the *host* copy. Before a
   collective NetCDF write, refresh host from device (`!$acc update host`); after a
   restart read, push device from host (`!$acc update device`). See
   [`spec-io-parallel.md`](spec-io-parallel.md) §3.1 and
   [`spec-io-restart.md`](spec-io-restart.md) §4. Skipping this doesn't crash — it
   silently writes or reads stale data, which is worse.

Bounds checking (`NoahmpArray.H`, §4) does not run inside device kernels — but by
design `NoahmpArray` is not on the device hot path (the physics reads the Fortran
allocatables directly; ERF uses the generated `Array4` accessor), so this is
consistent, not a gap. The ABI is unaffected: Tier-B device-resident arrays carry
**no** `fi` slot (see [`spec-fc-api.md`](spec-fc-api.md) §4a), so making state
device-resident never enlarges the self-referential struct or its dangling
surface.

## 8. Invariants to preserve

1. Never make `NoahmpIO_type` copyable, and keep the move ctor's `fptr` re-point.
2. Never widen `NoahmpIO_vector`'s re-exported API to anything that can relocate.
3. Keep every `bind(C)` shim going through `resolve_block`.
4. Keep `noahmparray_check` zero-cost under `NDEBUG`.
5. Never reach for MPI inside Noah-MP — always route fatals through
   `NoahmpIO_fatal` / `NoahmpIO_abort`.
6. (GPU) Keep coupling on the **one shared stream**; never reintroduce a per-step
   host sync or a second stream without an explicit dependency (§7.2).
7. (GPU) Re-issue `enter data` on any reallocation; refresh host↔device around all
   I/O (§7.1, §7.3).
