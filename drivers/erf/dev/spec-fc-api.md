# Spec: Fortran ↔ C API & the `@NoahmpMacro` code generator

> Status: living document · Owns: the C++ ↔ Fortran boundary and the generator
> that builds it. · Companion: [`spec-memory-safety.md`](spec-memory-safety.md),
> [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md),
> [`plan-cpp-interface.md`](plan-cpp-interface.md) (the GPU offload that §4's tier
> tags exist for).

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

### ABI membership is a switch, not a given (tiers) — *planned*

> **Status: planned, not yet implemented.** The tier tags (`@couple` /
> `@internal`) and the per-member ABI switch described here are part of the GPU
> offload ([`plan-cpp-interface.md`](plan-cpp-interface.md)); `NoahmpMacro.py` does
> **not** parse them today. In the current code **every** member of a contract
> block is fully projected — Fortran storage, `allocate()`, `fi` slot, and C++
> view. This design is recorded now so the ABI is built to accept the switch, not
> retrofitted; §4a (also marked planned) gives the region-level detail.

Under the planned design a generated variable is **not** automatically on the ABI:
"generate the Fortran storage" and "put a pointer in `NoahmpIO_type_fi`" become
separate switches, controlled per member by a **tier tag**. This matters because
the GPU port makes hundreds of internal arrays device-resident, and putting each
in `fi` would balloon the flat pointer struct for no reason — ERF never touches
them.

| Tier | Tag (planned) | In `fi` / C++ view? | Generated |
|------|-----|---------------------|-----------|
| **A** ERF-facing coupling | `@couple(dir=in\|out\|inout)` | **yes** | view + `fi` slot + device accessor + `enter data` |
| **B** internal, device-resident | `@internal` | **no** | Fortran storage + `allocate()` + `enter data` only |
| **C** host-only | (not in a contract block) | no | nothing (hand-written) |

Today every contract-block array behaves as Tier A (full projection); the tier
tags are the mechanism that will let new internal state opt out of `fi`. See §4
for how each planned tag maps to region membership, and
[`spec-add-coupled-variable.md`](spec-add-coupled-variable.md) for the author's-eye
workflow.

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

### Names

The public C++ client API is exposed at global scope: `NoahmpIO_type`,
`NoahmpIO_vector`, `NoahmpArray{2,3}D<noahmp_real>`, `noahmp_real`, and
`NoahmpIO_fatal` / `NoahmpIO_set_fatal_handler`. The `extern "C"` Fortran-bridge symbols (the
`*_fi` functions and the `NoahmpIO_type_fi` mirror struct) keep **C language
linkage** — their symbol names stay unmangled and match the Fortran `bind(C)`
side, and the struct is matched by layout, not by name.

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

The field list, written as one **contract block** `@NoahmpMacro:Source m_noahmpio
{ ... }` near the top of `NoahmpIO.H-mc`. The **handle** `m_noahmpio` only names the
binding for reference — the tool derives **no** symbols from it; concrete local
names are passed explicitly at each projection (see the binding grammar below).
**Order in the block = ABI order.** The contract block is *consumed* (replaced by a
breadcrumb comment); the owner class materializes its own fields through an ordinary
`@NoahmpMacro:CppStorageFields(m_noahmpio);` projection, exactly like every other
side of the boundary. From this one list the generator regenerates the bodies of
every `@NoahmpMacro:<Region>(m_noahmpio, …);` projection marker across five files:

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
from the *same* `[lo:hi, …]` bounds clause, so they cannot disagree.

### Annotation grammar

```cpp
@NoahmpMacro:Source m_noahmpio {           // the contract block (near top of header)
  int ids, ide, jds, jde;                  // ints: kind inferred, no annotation
  int numrad = 2;                          // int with C++ default (shared with Fortran)

  noahmp_real DTBL;                        // scalar: kind inferred from the type
  noahmp_real ZLVL = noahmp_zlvl_unset;    // C++ default may be a constexpr, not just a literal

  NoahmpArray2D<noahmp_real> XLAT[xstart:xend, ystart:yend];            // latitude [rad]
  NoahmpArray2D<noahmp_real> SWDOWN[xstart:xend, ystart:yend];          // ERF forcing in
  NoahmpArray2D<noahmp_real> HFX[xstart:xend, ystart:yend];             // flux out to ERF
  NoahmpArray3D<noahmp_real> SMOIS[xstart:xend, nsoil:nsoil, ystart:yend]; // soil moisture
  NoahmpArray3D<noahmp_real> U_PHY[xstart:xend, kms:kme, ystart:yend];  // U wind (i,k,j)
}
```

> **Planned:** the tier tags (`@couple(dir=…)` / `@internal`, §4a) will be appended
> to array lines to classify how far each variable travels. `NoahmpMacro.py` does
> not parse them yet; the grammar below describes what the generator accepts today.

**The annotation rule is uniform: kind is always inferred from the declared type;
you annotate only what the type cannot express.** Concretely, the generator
enforces (raising `SystemExit("NoahmpMacro: …")` otherwise):

- **Kind is inferred from the C++ type** (`int` / `noahmp_real` /
  `NoahmpArray{2,3}D`) — never spelled out. (There is no `@NoahmpMacro:scalar`;
  the `noahmp_real` type already says "scalar", exactly as `int` says "int".) A type
  the vocabulary does not know (`double`, `int64_t`, …) is a hard error.
- **Arrays** carry a trailing `[lo:hi, …]` bounds clause on the name, in **Fortran
  order** — the one fact the type cannot hold. It has no `@NoahmpMacro:` marker; it
  is pure data. The rank is stated **only** by the type (`NoahmpArray{2,3}D`) and the
  number of bound pairs must equal it (the tool's one strictness check). Every bound
  token must be a **literal** or a **coupled member name** (`xstart`, `xend`, `kms`,
  `kme`, `nsoil`, `numrad`, …) so C++ and Fortran resolve it to the *same* value —
  never a Fortran-only parameter.
- **Trailing comments** are free-form doc (no `doc=` form): kept verbatim in the
  owner header, dropped from the Fortran `allocate`, and ignored by the tool.
- **Ints / scalars** need no annotation; an optional C++ default
  (`int numrad = 2;`) is shared with Fortran.
- **Tier tags** (`@couple(dir=…)` / `@internal`) are *planned* (§4a): the one
  annotation the *type* cannot express — how far the variable travels, not what it
  is. `NoahmpMacro.py` does not accept them yet; when it does they will be parsed
  onto the `Member` and consumed by the region filters in §4a.

### 4a. Tiers & the GPU coupling projections — *planned*

> **Status: planned, not yet implemented.** Nothing in this subsection is in
> `NoahmpMacro.py` today — not the `@couple`/`@internal` parsing, nor the
> `CppCoupleAccessors` / `FortranDevPtrExport` regions, nor the `!$acc enter data`
> emitted by `FortranArrayAllocate`. It is the design for the GPU offload
> ([`plan-cpp-interface.md`](plan-cpp-interface.md)); the "Generator
> implementation" note at the end of the subsection is the to-do that would make it
> real. Today `FortranArrayAllocate` emits a plain `allocate()` and every member is
> fully projected.

The tier tags will exist so the **same contract line** that generates the ABI can
also drive GPU device residency and the ERF-facing device accessors — the
single-source-of-truth principle extended to the offload port
([`plan-cpp-interface.md`](plan-cpp-interface.md) §5, principle 3). A tag would
change only **which regions a member flows into**, never how a region renders:

| Region | Tier A (`@couple`) | Tier B (`@internal`) |
|--------|:---:|:---:|
| `CppMirrorFields` / `FortranMirrorFields` (`fi` struct) | ✅ | ❌ |
| `MemberCount` (ABI count / `AssertAbi`) | ✅ | ❌ |
| `CppStorageFields` (C++ owner view) | ✅ | ❌ |
| `FortranStorageFields` (Fortran allocatable) | ✅ | ✅ |
| `FortranArrayAllocate` (`allocate()` + `enter data`) | ✅ | ✅ |
| `CppCoupleAccessors` (device `*_a4()` / `*_v()`) | ✅ (by `dir`) | ❌ |
| `FortranDevPtrExport` (`*_devptr_fi` → `acc_deviceptr`) | ✅ | ❌ |

So a `@internal` member costs **zero** ABI (no `fi` slot, no `MemberCount`
increment, no C++ view) yet still gets generated Fortran storage, `allocate()`,
and its `!$acc enter data` residency directive. A `@couple` member additionally
gets the GPU coupling glue described in
[`sketch-couple-variable-gpu.md`](sketch-couple-variable-gpu.md):

- **`CppCoupleAccessors`** emits, per Tier-A variable, an `amrex::Array4` alias
  (`*_a4()`) over the shared device memory — or, when the array's Fortran axis
  order is not `(i,j,k)` (e.g. `U_PHY(i,k,j)`), a stride-matched device accessor
  (`*_v()`) computed from the `[lo:hi, …]` bounds clause, so C++ and Fortran agree
  by construction. `dir=in`/`out`/`inout` documents intent and can gate a
  `const` view.
- **`FortranDevPtrExport`** emits the `*_devptr_fi` `bind(C)` shim that returns
  `acc_deviceptr(nm%NAME)` (Option B: Fortran owns, exports the device address for
  ERF to wrap). See [`spec-memory-safety.md`](spec-memory-safety.md) §7.
- **`FortranArrayAllocate`** appends `!$acc enter data create(nm%NAME)` after the
  `allocate` for any tiered array — a harmless no-op in a host (offload-off) build.

Generator implementation: add an `internal` flag and a `couple_dir` field to
`Member` (parsed in `parse_members`). The ABI-bearing regions — the Mirror pair,
`MemberCount`, `CppStorageFields`, and the Construct/Wire regions (ctor + move +
`C_LOC`/`C_F_POINTER` plumbing) — filter on `not m.internal`, so a Tier-B member
touches none of them. `CppCoupleAccessors`/`FortranDevPtrExport` filter on
`m.couple_dir is not None`. `FortranStorageFields`/`FortranArrayAllocate` keep
emitting **all** members. `KIND_TRAITS` is untouched — tiers gate *membership*,
not per-kind *rendering*.

### The binding grammar (what you write in the templates)

A **binding** is identified by a **handle** (`m_noahmpio`). The handle is only a
reference — the tool derives **no** identifiers from it. There is ONE marker prefix,
`@NoahmpMacro:`; the forms are distinguished by position, not by a second punctuator:

- **Contract** (one per handle) — a single block, the source of truth, near the top
  of the owning template: `@NoahmpMacro:Source <handle> { …field declarations… }`.
  It is *consumed* (replaced by a breadcrumb); the owner class emits its own fields
  through the `CppStorageFields` projection below.
- **Projection** — names a region, the handle, and any **site-local names** the
  region needs as explicit `key=value` args. Two shapes:
  - **block**: `@NoahmpMacro:<Region>(<handle>[, key=value]…);` on its own line,
    expanded to a banner-commented body block;
  - **value**: `… = @NoahmpMacro:MemberCount(<handle>);` substituted in place as a
    scalar (shows *where* the constant is assigned; the declaration is hand-written).
- **Array bounds clause** — a trailing `[lo:hi, …]` on an array name in the contract
  block (the one thing the C++ type cannot carry). It carries no marker — pure data.

Rigidly enforced: a malformed `Source` line, an unknown region, a block region used
inline (or a value region used as a statement / given params), a marker arg that is
not a bare identifier, and a missing/unknown projection param are each a named error.

Site-local names are **never derived from the handle** — they are passed at the
call so the generated code references the right local. The param names follow a
`<lang>type[_fi]` convention (`f*` = Fortran, `c*` = C++, `*_fi` = the flat interop
struct on that side):

| Param | Meaning | Example |
|-------|---------|---------|
| `ftype` | owning Fortran derived-type instance (`<ftype>%X`) | `blk`, `NoahmpIO` |
| `ftype_fi` | Fortran bind(C) interop struct instance | `NoahmpIO_cptr` |
| `ctype` | C++ owner-object instance | `o` (move source) |
| `ctype_fi` | C++ mirror member holding the interop struct | `fptr` |

`REGION_PARAMS` declares which keys each region requires; the generator errors on a
missing **or** unknown key, so a typo'd local is a named failure, not a miscompile.

Region names are **purpose-based and file-agnostic**, grouped by the four facets a
struct binding must model (type-correspondence / construction / ownership-wiring /
layout). Each is a named registry entry, so `grep <RegionName>` lands on its
definition:

| Facet | Region (required params) | Emits |
|-------|--------------------------|-------|
| **Mirror** | `CppMirrorFields` / `FortranMirrorFields` / `MemberCount` (value) | the flat pointer struct on each side + its element count |
| **Construct** | `CppCtorParams` / `CppCtorInit` / `CppBindAddrs` / `CppMoveInit`(`ctype`) / `CppMoveRepoint`(`ctype_fi`) | ctor + move plumbing for the C++ owner and its mirror |
| **Wire** | `CppArrayViews`(`ctype_fi`) / `FortranScalarWire`(`ftype_fi`,`ftype`) (`C_F_POINTER`) / `FortranArrayWire`(`ftype_fi`,`ftype`) (`C_LOC`) | physically connect the two representations |
| **Storage** | `CppStorageFields` / `FortranStorageFields` / `FortranArrayAllocate`(`ftype`) | each owning side's real fields + the Fortran allocation |
| **Couple/GPU** *(planned, see §4a)* | `CppCoupleAccessors`(`ctype_fi`) / `FortranDevPtrExport`(`ftype`) | Tier-A device accessors (`*_a4()`/`*_v()`) + `acc_deviceptr` export — *not emitted today* |

`CppStorageFields` is the C++ owner class's own fields, emitted **verbatim** from the
contract body (so initializers like `int numrad = 2;`, trailing doc comments, and the
blank-line grouping survive exactly); every other region is reconstructed from the
parsed members.

### Internal architecture (read `tools/NoahmpMacro.py` top-to-bottom)

> Skip this subsection unless you're modifying the generator itself. *Adding* a
> coupled variable needs only §1 here and
> [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md).

The generator is a few declarative layers over a **binding-agnostic** engine.

**Layer 1 — `KIND_TRAITS`: the matching rules, in one place.** A kind→projection
table that *is* the code form of the §1 "What the boundary carries" table. For
each kind it states every projection a region needs — the C++ pointer/ctor type,
the Fortran `bind(C)` field type, the `NoahmpIO_type` storage declaration, and the
wiring (`cfptr` = `C_F_POINTER` for scalar-likes, `cloc` = `C_LOC` for arrays):

| kind | C++ ctor type | Fortran `bind(C)` | `NoahmpIO_type` storage | wiring |
|------|---------------|-------------------|-------------------------|--------|
| `int`    | `int*`         | `type(C_PTR)` | `integer(C_INT), pointer => null()`        | `cfptr` |
| `scalar` | `noahmp_real*` | `type(C_PTR)` | `real(c_kind_noahmp), pointer => null()`   | `cfptr` |
| `array`  | (n/a — `C_LOC`)| `type(C_PTR)` | `real(c_kind_noahmp), allocatable, dim(:…)`| `cloc`  |

Every per-kind type string lives here only; no renderer re-encodes one.

**Layer 2 — `Region(filter, item, layout)`: regions as specs.** Most regions are
"select some members, render one token each from the trait table, lay the tokens
out", so they are declarative `Region` entries. Every region — spec or kept
function — is uniformly `f(binding, params) -> [lines]`, where `params` are the
call-site `key=value` site-locals:

```
Region(SCALAR_LIKE, cpp_ctor_param, L_wrap(...))   # filter, per-member item, layout
```

*Stop-collapsing rule* — a region stays a named `r_*` function only when it (a)
needs run-grouped per-run prefixes (`CppMirrorFields`, `FortranStorageFields` —
keep `_runs`), (b) carries array bounds the trait table doesn't model
(`CppArrayViews`, `FortranArrayAllocate`), or (c) emits the contract body verbatim
(`CppStorageFields`, which reads `Binding.body` rather than the parsed members). The
element count is none of these — it is an inline scalar, so it lives in
`VALUE_REGIONS` (substituted mid-line, not a block). Kept functions still read
`KIND_TRAITS` for type *facts*.

**Layer 3 — `REGION_PARAMS`: explicit call-site site-locals.** There is **no**
name-derivation seam: the handle yields no identifiers. Any region that emits code
referencing a hand-written local (`ftype`/`ftype_fi`/`ctype`/`ctype_fi`, per the
table above) declares those keys in `REGION_PARAMS`, and the projection marker
supplies them as `key=value`. `expand` validates exactly the required keys are
present (missing or unknown → named `SystemExit`); `_validate_registry()` self-checks
the tables at startup (a region cannot be both block and value; `REGION_PARAMS` may
only name block regions). The C++ array-view class is the fixed `ARRAY_TMPL` constant
(`NoahmpArray`), not a site-local.

**Layer 4 — `Binding`: one resolved contract.** the `tag` (handle) + parsed
`members` + the raw `body` (for `CppStorageFields`' verbatim emission) + the
do-not-edit `banner` (built from the file the contract lives in, so it lands only in
generated output).

**Layer 5 — the engine (two passes, binding-agnostic).**

```
parse_members(body, cc, tag) -> ordered [Member(name, kind, rank, begin, end)]
Region / r_*                 -> body LINES for one region: f(binding, params)
REGIONS / VALUE_REGIONS      -> GLOBAL {region -> block-renderer} / {region -> scalar}
collect_bindings(templates)  -> PASS 1: {handle -> Binding} from each Source block
expand(text, cc, bindings)   -> PASS 2: contract block -> breadcrumb; block markers
                                expanded; then ONE strict classifier pass (PROJ_RE)
                                over every remaining marker (value, or precise error)
main(argv)                   -> validate registry, glob *-mc, run both passes, write/--check
```

Templates are discovered by **glob** (`ROOT/*-mc`); each target is the template
with `-mc` stripped, and the comment char is inferred from the extension. There is
no per-file region dict and no `Boundary`/`Target` object — the engine scans each
template for whatever markers it contains.

**Adding a new binding:** add a `@NoahmpMacro:Source <newhandle> { … }` contract
block (in any `*-mc` template) plus projection markers (each supplying its
`REGION_PARAMS` site-locals). No Python edits — the target follows the `*-mc`
convention and `collect_bindings` picks the handle up automatically. (Only a
genuinely new *kind* of field needs a `KIND_TRAITS` entry.)

Rendering helpers pack several members per line; since targets regenerate every
build, density costs no diff noise.

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
- [ ] **Tiers / GPU glue (planned — only once §4a is implemented):** `@internal`
      members produce **no** `fi` slot and leave `NOAHMP_IO_FI_NUM_MEMBERS`
      unchanged; every `@couple` member has a generated accessor and `*_devptr_fi`;
      every tiered array has its `enter data`. *(Not applicable to the current
      generator, which fully projects every member.)*
- [ ] [`spec-add-coupled-variable.md`](spec-add-coupled-variable.md) still works
      end to end.
