# Spec: memory safety & runtime robustness

> Status: living document · Owns: the invariants that keep the self-referential
> coupling struct from corrupting memory, plus the parallel-safe fatal path.
> Companion: [`spec-fc-api.md`](spec-fc-api.md).

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
  the source's members to the destination's (`@NoahmpMacro:MoveFptrRepoint`).

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

This keeps the *entire* parallel-runtime dependency in the host, not in Noah-MP —
a precondition for the GPU port (see [`plan-cpp-interface.md`](plan-cpp-interface.md)).

## 6. Invariants to preserve

1. Never make `NoahmpIO_type` copyable, and keep the move ctor's `fptr` re-point.
2. Never widen `NoahmpIO_vector`'s re-exported API to anything that can relocate.
3. Keep every `bind(C)` shim going through `resolve_block`.
4. Keep `noahmparray_check` zero-cost under `NDEBUG`.
5. Never reach for MPI inside Noah-MP — always route fatals through
   `NoahmpIO_fatal` / `NoahmpIO_abort`.
