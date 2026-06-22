#!/usr/bin/env python3
# =============================================================================
# NoahmpMacro.py -- a tiny "binding contract" projector for C++ <-> Fortran.
# =============================================================================
# THE PROBLEM. Two separately-compiled languages must agree on a shared struct at
# a boundary: same members, same order, same precision, no padding. Keeping that
# agreement by hand rots the instant a member is added. This tool removes the
# possibility of drift: you declare the binding ONCE as a neutral contract, and it
# PROJECTS that contract onto both sides (a C++ class + a Fortran bind(C) type) and
# everywhere the wiring needs it. The two sides cannot disagree unless the
# projector is wrong -- which `make codegen-check` catches in CI.
#
# This is a domain-specific micro-IDL: one contract, many language projections. It
# does NOT parse or translate existing source -- it owns both sides and emits them.
#
# SCOPE: VARIABLES ONLY. A binding describes the DATA members crossing the
# boundary. Methods / API entry points are a different beast (lifetime, dispatch,
# error propagation) and are written by hand, not generated here.
#
# -----------------------------------------------------------------------------
# THE GRAMMAR (what you write in the tracked `*-mc` templates)
# -----------------------------------------------------------------------------
# A binding is identified by a HANDLE (e.g. `m_NoahmpIO`). The handle is ONLY a
# reference -- the tool derives NO identifiers from it; every concrete site-local
# name a projection emits is passed explicitly at the call (see PARAMS below).
# Markers carry the handle, so one tool hosts many bindings. Punctuation splits the
# two surfaces: `:` for STANDALONE projections, `.` for an IN-DECLARATION field
# annotation. Every form ends its statement with `;`.
#
#   * the CONTRACT (exactly one per handle): a single block, the source of truth,
#     near the top of the owning template -- the ordered field list in the owner
#     language's own syntax. It is consumed (the class materializes its own fields
#     through CppStorageFields, like every other projection):
#
#             @NoahmpMacro:Source m_noahmpio {
#                 int ids, ide, jds, jde;             // ints: kind inferred, no annotation
#                 int numrad = 2;                     // int with a shared C++/Fortran default
#                 noahmp_real DTBL;                   // scalar: kind inferred from the type
#                 noahmp_real ZLVL = -9999.0;         // free-form doc comment, optional
#                 NoahmpArray2D<noahmp_real> XLAT(@NoahmpMacro:bounds_2d(xstart:xend, ystart:yend)); // doc
#                 NoahmpArray3D<noahmp_real> U_PHY(@NoahmpMacro:bounds_3d(xstart:xend, kms:kme, ystart:yend)); // doc
#             }
#
#     RULE: kind is ALWAYS inferred from the declared type (int / noahmp_real /
#     NoahmpArrayND). The one fact the type cannot carry -- array `bounds_Nd(...)` --
#     rides in an in-declaration `(@NoahmpMacro:bounds_Nd(...))` annotation (rank
#     stated in both the type and the annotation; they must agree). Anything after
#     the line comment is free-form doc, kept verbatim in the owner header and
#     ignored by the tool (no `doc=` variable).
#
#   * a PROJECTION marker names a region, the handle, and any site-local names the
#     region needs as explicit `key=value` args (never derived from the handle).
#     Two forms:
#
#       - BLOCK: on its own line, expands to a banner-wrapped body block:
#
#             @NoahmpMacro:CppStorageFields(m_noahmpio);
#             @NoahmpMacro:FortranScalarWire(m_noahmpio, ftype_fi=NoahmpIO_cptr, ftype=blk);
#
#       - VALUE: substituted in place as a scalar, showing WHERE it is assigned:
#
#             inline constexpr std::size_t NOAHMP_IO_FI_NUM_MEMBERS = @NoahmpMacro:MemberCount(m_noahmpio);
#
# Region names are PURPOSE-based and file-agnostic (CppMirrorFields, not
# "FiMembers"), grouped into the four facets a struct binding must model:
#   Mirror   -- the flat interop struct on each side + its element count
#   Construct-- ctor/move plumbing that builds the C++ owner and its mirror
#   Wire     -- physically connect the two representations (C_F_POINTER / C_LOC)
#   Storage  -- the owning side's real fields + their allocation
# These four facets are: type-correspondence / construction / ownership-wiring /
# layout -- the things a clean interop contract must capture, not just types.
#
# -----------------------------------------------------------------------------
# HOW A FIELD PROJECTS (the contract, as a table)
# -----------------------------------------------------------------------------
#   kind   | C++ ctor type | Fortran bind(C) | NoahmpIO_type storage          | wiring
#   -------|---------------|-----------------|--------------------------------|-------
#   int    | int*          | type(C_PTR)     | integer(C_INT), pointer=>null()| C_F_POINTER (C++-owned)
#   scalar | noahmp_real*  | type(C_PTR)     | real(c_kind_noahmp),ptr=>null() | C_F_POINTER (C++-owned)
#   array  | (via C_LOC)   | type(C_PTR)     | real(...), allocatable, dim(:..)| C_LOC      (Fortran-owned)
# Every per-kind type fact lives in KIND_TRAITS below and ONLY there.
#
# -----------------------------------------------------------------------------
# ADDING THINGS (the one-line test)
# -----------------------------------------------------------------------------
#   * a FIELD     -> one declaration line in the `@NoahmpMacro:Source <handle> {...}`
#                    contract block.
#   * a KIND      -> one KIND_TRAITS entry (+ a VARTYPE_PREFIX entry if scalar-like).
#   * a BINDING   -> a new handle: a `@NoahmpMacro:Source <handle> { ... }` contract
#                    block plus projection markers (each supplying its REGION_PARAMS
#                    site-locals). No Python edits; targets follow the `*-mc` convention.
#
# INTERNAL ARCHITECTURE (top-to-bottom below):
#   1. KIND_TRAITS   -- kind -> projection traits (the table above), in ONE place.
#   2. Region / r_*  -- a region is (filter, item, layout); regular ones are specs,
#                       the few irregular ones (run-grouped prefix, array bounds)
#                       stay named r_* functions; the element count is an inline
#                       VALUE_REGIONS scalar. All read KIND_TRAITS.
#   3. PARAMS        -- REGION_PARAMS: the site-local names each region needs, passed
#                       explicitly at the call (never derived from the handle) and
#                       validated by expand().
#   4. Binding       -- handle + parsed members + banner (one per contract).
#   5. engine        -- collect_bindings (pass 1: {handle -> Binding}) and expand
#                       (pass 2: replace every marker). Both are binding-agnostic.
#
# Stdlib only (regex). No third-party libraries, any Python 3.
#
# Templates are the ONLY hand-edited files (`*-mc`, discovered by glob in ROOT);
# each generated TARGET is the template with `-mc` stripped, gitignored, rebuilt
# every build. List-regions are packed several members per line; targets
# regenerate every build, so density costs no diff noise.
#
# Usage:
#   python3 tools/NoahmpMacro.py            # (re)generate the targets
#   python3 tools/NoahmpMacro.py --check    # fail (exit 1) if regen would change
# =============================================================================

import os
import re
import sys
import glob
import difflib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)            # drivers/erf

# The macro prefix. ONE punctuator (`@NoahmpMacro:`) for every form; one constant
# hosts all bindings (multiplicity comes from the HANDLE, not the prefix). The forms
# are distinguished by POSITION, not by a second punctuator:
#   * contract block   -- @NoahmpMacro:Source <handle> { ...fields... }
#   * block projection -- @NoahmpMacro:<Region>(<handle>[, key=value]...);   (own line)
#   * value projection -- ... = @NoahmpMacro:MemberCount(<handle>);          (inline)
#   * field annotation -- NoahmpArray2D<noahmp_real> XLAT(@NoahmpMacro:bounds_2d(lo:hi, lo:hi));
MARKER = "@NoahmpMacro:"                       # the one and only marker prefix
SOURCE = "Source"                             # the contract block keyword

# The contract is one block at the top of the owning template -- the single source
# of truth -- and the owner class materializes its own fields through an ordinary
# region projection, exactly like every other side of the boundary:
#     @NoahmpMacro:Source m_noahmpio {     <- the contract block (field list)
#         int ids, ide;  noahmp_real DTBL;  NoahmpArray2D<noahmp_real> XLAT(...);
#     }
#     ...
#     @NoahmpMacro:CppStorageFields(m_noahmpio);   <- emits those fields into the class
# The HANDLE is only a reference: projections supply their own concrete names as
# explicit params (see REGION_PARAMS), so the tool never derives a symbol from it.
# Block opener `@NoahmpMacro:Source <handle> {`; body until the matching bare `}`.
SOURCE_OPEN_RE = re.compile(r"^\s*%s%s\s+(\w+)\s*\{\s*$" % (re.escape(MARKER), SOURCE))


# ---------------------------------------------------------------------------
# A parsed contract field.
# ---------------------------------------------------------------------------
class Member:
    __slots__ = ("name", "kind", "rank", "begin", "end")

    def __init__(self, name, kind, rank=0, begin=None, end=None):
        self.name = name
        self.kind = kind            # 'int' | 'scalar' | 'array'
        self.rank = rank
        self.begin = begin or []
        self.end = end or []


# ---------------------------------------------------------------------------
# 1. Kind -> projection trait table: the RULES for matching a Fortran type to a C
# struct, in ONE place. Every per-kind type fact the renderers need lives here
# (and ONLY here), so the C++ struct and the Fortran bind(C) type are described as
# two projections of the same member rather than re-encoded across renderers.
# ---------------------------------------------------------------------------
class KindTrait:
    __slots__ = ("scalar_like", "cpp_ptr_base", "cpp_ctor_type",
                 "fort_bindc_type", "fort_storage", "wiring")

    def __init__(self, scalar_like, cpp_ptr_base, cpp_ctor_type,
                 fort_bindc_type, fort_storage, wiring):
        self.scalar_like = scalar_like         # passed by pointer & appears in the ctor
        self.cpp_ptr_base = cpp_ptr_base       # C++ fi-struct pointer base type (`*name`)
        self.cpp_ctor_type = cpp_ctor_type     # C++ ctor parameter type, or None (arrays)
        self.fort_bindc_type = fort_bindc_type # Fortran bind(C) field type
        self.fort_storage = fort_storage       # NoahmpIO_type storage token: fn(member) -> str
        self.wiring = wiring                   # 'cfptr' (C_F_POINTER) | 'cloc' (C_LOC)


KIND_TRAITS = {
    "int":    KindTrait(True,  "int",         "int*",
                        "type(C_PTR)", lambda m: "%s=>null()" % m.name, "cfptr"),
    "scalar": KindTrait(True,  "noahmp_real", "noahmp_real*",
                        "type(C_PTR)", lambda m: "%s=>null()" % m.name, "cfptr"),
    "array":  KindTrait(False, "noahmp_real", None,
                        "type(C_PTR)",
                        lambda m: "real(kind=c_kind_noahmp), allocatable, dimension(%s) :: %s"
                                  % (",".join([":"] * m.rank), m.name),
                        "cloc"),
}

# Line-prefix for the grouped NoahmpIO_type storage declaration of scalar-like
# kinds (arrays emit a self-contained line from KIND_TRAITS[...].fort_storage, so
# they have no shared prefix -- this asymmetry is why FortranStorageFields stays a
# named renderer rather than a flat Region spec).
VARTYPE_PREFIX = {
    "int":    "integer(C_INT), pointer :: ",
    "scalar": "real(kind=c_kind_noahmp), pointer :: ",
}


# ---------------------------------------------------------------------------
# Contract parser: read a `<handle> { ... }` definition body into ordered Members.
# ---------------------------------------------------------------------------
def parse_members(body, cc, tag):
    """Parse the definition block BODY lines (everything between the braces) into an
    ordered Member list. `cc` is the owner language's line comment. The type-token
    -> kind recognition (int / noahmp_real / NoahmpArrayND) defines this binding's
    field vocabulary; kind is inferred from the type, never annotated. Array bounds
    are the one fact the type cannot carry, so they ride in an in-declaration
    `(@NoahmpMacro:bounds_Nd(...))` annotation; the trailing comment is free-form
    doc the parser ignores entirely (it survives verbatim into the owner header)."""
    members = []
    seen = set()
    for raw in body:
        if raw.strip().startswith(cc):
            continue                                   # guidance/usage comment line
        # Keep only the declaration; the trailing comment is free-form doc (the
        # bounds annotation now lives in the code, not the comment).
        code = (raw.split(cc, 1)[0] if cc in raw else raw).strip()
        if not code or not code.endswith(";"):
            continue                                   # blank / comment-only line
        code = code[:-1].strip()                        # drop trailing ';'

        # Dispatch on the leading TYPE token, matched at a WORD BOUNDARY so a type
        # the vocabulary does not know (`int64_t`, `double`, `bool`, ...) is a hard
        # error here, never silently mis-parsed as a known kind.
        if re.match(r"NoahmpArray\w*<", code):
            # Rank is stated twice -- in the type and in the bounds annotation -- and
            # the two MUST agree; the annotation also carries the bound list.
            m = re.match(r"NoahmpArray([23])D<noahmp_real>\s+(\w+)\s*\(\s*"
                         r"%sbounds_([23])d\s*\(([^)]*)\)\s*\)$"
                         % re.escape(MARKER), code)
            if not m:
                raise SystemExit("NoahmpMacro: cannot parse array decl (want "
                                 "`NoahmpArray{2,3}D<noahmp_real> NAME(%sbounds_{2,3}d(lo:hi, ...));`): %r"
                                 % (MARKER, raw))
            rank, name, brank = int(m.group(1)), m.group(2), int(m.group(3))
            if brank != rank:
                raise SystemExit("NoahmpMacro: array %s is %d-D but annotated bounds_%dd"
                                 % (name, rank, brank))
            pairs = [p.strip() for p in m.group(4).split(",") if p.strip()]
            if len(pairs) != rank:
                raise SystemExit("NoahmpMacro: array %s is %d-D but has %d bound pairs"
                                 % (name, rank, len(pairs)))
            begin, end = [], []
            for p in pairs:
                if p.count(":") != 1:
                    raise SystemExit("NoahmpMacro: bad bound %r for %s (want lo:hi)"
                                     % (p, name))
                lo, hi = (s.strip() for s in p.split(":"))
                begin.append(lo)
                end.append(hi)
            members.append(Member(name, "array", rank, begin, end))

        elif re.match(r"noahmp_real\b", code):
            m = re.match(r"noahmp_real\s+(\w+)\s*(=\s*\S.*)?$", code)
            if not m:
                raise SystemExit("NoahmpMacro: cannot parse scalar decl (want "
                                 "`noahmp_real NAME[ = default];`): %r" % raw)
            # Kind is inferred from `noahmp_real`; no annotation is required.
            members.append(Member(m.group(1), "scalar"))

        elif re.match(r"int\b", code):
            # Comma-separated identifiers, each optionally `= default`. Every name
            # must be a bare identifier; anything else (a stray annotation, a `(`,
            # an unsupported declarator) is rejected rather than stored verbatim.
            for nm in code[3:].split(","):
                nm = nm.split("=")[0].strip().lstrip("*").strip()   # drop any `= default`
                if not nm:
                    continue
                if not re.match(r"^[A-Za-z_]\w*$", nm):
                    raise SystemExit("NoahmpMacro: bad int member name %r in decl: %r"
                                     % (nm, raw))
                members.append(Member(nm, "int"))
        else:
            raise SystemExit("NoahmpMacro: unrecognised member decl in %s (kind is "
                             "inferred from the type: int / noahmp_real / "
                             "NoahmpArray{2,3}D<noahmp_real>): %r" % (tag, raw))

    for m in members:
        if m.name in seen:
            raise SystemExit("NoahmpMacro: duplicate member name %r in %s definition"
                             % (m.name, tag))
        seen.add(m.name)
    if not members:
        raise SystemExit("NoahmpMacro: %s definition block is empty" % tag)

    # Validate every array bound token now that all members are known. Each token is
    # emitted verbatim into BOTH the C++ NoahmpArray extent and the Fortran
    # allocate(), so it must resolve to the SAME integer on both sides: an integer
    # literal (negative lower bounds allowed) or a declared int member. This turns a
    # typo'd/renamed dimension from an opaque downstream compile error into a named
    # SystemExit here. (It checks token IDENTITY, not dimension semantics or order.)
    int_names = {m.name for m in members if m.kind == "int"}
    for m in members:
        if m.kind != "array":
            continue
        for tok in m.begin + m.end:
            if not (re.match(r"^-?\d+$", tok) or tok in int_names):
                raise SystemExit("NoahmpMacro: array %s bound %r is neither an "
                                 "integer literal nor a declared int member (typo?)"
                                 % (m.name, tok))
    return members


# ---------------------------------------------------------------------------
# Rendering helpers
#
# The generated targets are gitignored and rewritten every build, so packing
# several members per line costs no diff noise -- readability is the only goal.
# Every list-style region is therefore wrapped to a fixed right margin instead of
# exploded one-member-per-line.
# ---------------------------------------------------------------------------
WRAP_COL = 100          # soft right margin for packed lines


def _is_scalar_like(m):      # participates in the constructor / scalar wiring
    return KIND_TRAITS[m.kind].scalar_like


def _commad(items):
    """Append a comma to every entry except the last (for arg/init lists)."""
    return [s + ("," if i < len(items) - 1 else "") for i, s in enumerate(items)]


def _wrap(items, indent):
    """Greedily pack pre-punctuated `items` (each already carrying its own trailing
    ',' / ';' where one is needed) onto `indent`-prefixed lines no wider than
    WRAP_COL. Items sharing a line are space-joined."""
    lines, cur = [], ""
    for it in items:
        cand = it if not cur else cur + " " + it
        if cur and len(indent) + len(cand) > WRAP_COL:
            lines.append(indent + cur)
            cur = it
        else:
            cur = cand
    if cur:
        lines.append(indent + cur)
    return lines


def _wrap_joined(items, indent, sep):
    """Greedily pack raw `items` onto `indent`-prefixed lines no wider than WRAP_COL,
    joining items on a line with `sep` and leaving NO trailing separator. Used for
    Fortran statement lists, where ';' separates statements within a line (it is not
    a continuation), so a trailing ';' would be a redundant empty statement."""
    lines, cur = [], []
    for it in items:
        if cur and len(indent) + len(sep.join(cur + [it])) > WRAP_COL:
            lines.append(indent + sep.join(cur))
            cur = [it]
        else:
            cur.append(it)
    if cur:
        lines.append(indent + sep.join(cur))
    return lines


def _decl_lines(indent, prefix, tokens, end):
    """Render a grouped declaration -- `prefix` once per line, `tokens` joined by
    ', ', closed by `end` -- wrapping so no line exceeds WRAP_COL."""
    lines, cur = [], []
    for t in tokens:
        if cur and len(indent) + len(prefix) + len(", ".join(cur + [t])) + len(end) > WRAP_COL:
            lines.append(indent + prefix + ", ".join(cur) + end)
            cur = [t]
        else:
            cur.append(t)
    if cur:
        lines.append(indent + prefix + ", ".join(cur) + end)
    return lines


def _runs(ms, key):
    """Yield (key_value, [members]) for each maximal run of consecutive members
    sharing key(m), preserving order."""
    i = 0
    while i < len(ms):
        k = key(ms[i])
        j = i + 1
        while j < len(ms) and key(ms[j]) == k:
            j += 1
        yield k, ms[i:j]
        i = j


# ---------------------------------------------------------------------------
# 3. Call-site PARAMS: the concrete local names a projection operates on.
#
# A region that emits code referencing a site-local (the owner-object instance, the
# interop/mirror struct instance, the move source) does NOT guess that name -- and
# crucially never DERIVES it from the binding handle. Each such name is passed
# explicitly at the projection site as a `key=value` argument:
#
#     @NoahmpMacro:FortranScalarWire(m_NoahmpIO, ftype_fi=NoahmpIO_cptr, ftype=blk);
#     @NoahmpMacro:FortranArrayAllocate(m_NoahmpIO, ftype=NoahmpIO);
#     @NoahmpMacro:CppArrayViews(m_NoahmpIO, ctype_fi=fptr);
#     @NoahmpMacro:CppMoveInit(m_NoahmpIO, ctype=o);
#
# Param names follow a `<lang>type[_fi]` convention -- `f*` = a Fortran-side name,
# `c*` = a C++-side name, `*_fi` = the flat interop struct on that side:
#   ftype    -- the owning Fortran derived-type instance (`<ftype>%X`): the local
#               `blk` in the wiring shims, the dummy `NoahmpIO` in the allocate.
#   ftype_fi -- the Fortran bind(C) interop struct instance (the `NoahmpIO_cptr` arg).
#   ctype    -- a C++ owner-object instance (the move source `o` in CppMoveInit).
#   ctype_fi -- the C++ mirror member holding the interop struct (the `fptr` member).
# REGION_PARAMS declares which keys each region requires; expand() validates the
# call site supplies exactly those (missing OR unknown keys are an error), so a
# typo'd local surfaces as a named SystemExit rather than a downstream miscompile.
# ---------------------------------------------------------------------------
ARRAY_TMPL = "NoahmpArray"          # C++ array-view class template (NoahmpArray{2,3}D)

REGION_PARAMS = {
    "CppMoveInit":          ["ctype"],
    "CppMoveRepoint":       ["ctype_fi"],
    "CppArrayViews":        ["ctype_fi"],
    "FortranScalarWire":    ["ftype_fi", "ftype"],
    "FortranArrayWire":     ["ftype_fi", "ftype"],
    "FortranArrayAllocate": ["ftype"],
}


# ---------------------------------------------------------------------------
# 2. Region = (filter, item, layout).
#
# Most regions are "select some members, render one token each, lay the tokens
# out", so they collapse to a declarative Region(filter, item, layout) spec instead
# of a bespoke function. `filter` picks the members, `item(member, params)` projects
# one member to its per-region token (reading KIND_TRAITS for type facts and the
# call-site `params` for site-local names), and `layout` is one of the wrap/decl
# helpers. A Region is callable, so expand() treats specs and the few kept r_*
# functions identically: every region is `f(binding, params) -> [lines]` -- the
# resolved Binding and the call-site params flow through each region.
# ---------------------------------------------------------------------------
class Region:
    __slots__ = ("filter", "item", "layout")

    def __init__(self, filt, item, layout):
        self.filter = filt
        self.item = item
        self.layout = layout

    def __call__(self, b, params):
        return self.layout([self.item(m, params)
                            for m in b.members if self.filter(m)])


# Member filters (named so each region's selection stays greppable).
ALL         = lambda m: True
SCALAR_LIKE = _is_scalar_like
WIRED       = lambda which: (lambda m: KIND_TRAITS[m.kind].wiring == which)


# Layout adapters: close over indent/separator and hand back the EXISTING packing
# helpers unchanged (no new layout behavior is introduced here).
def L_wrap(indent, comma=False):
    return lambda items: _wrap(_commad(items) if comma else items, indent)


def L_wrap_joined(indent, sep):
    return lambda items: _wrap_joined(items, indent, sep)


# Per-member projections (the `item` of a Region). Each takes (member, params);
# type facts come from KIND_TRAITS, site-local names from the call-site `params`.
def cpp_ctor_param(m, p):   return "%s %s" % (KIND_TRAITS[m.kind].cpp_ctor_type, m.name)
def cpp_ctor_init(m, p):    return "%s(%s)" % (m.name, m.name)
def cpp_bind_addr(m, p):    return "&%s" % m.name
def cpp_move_init(m, p):    return "%s(%s.%s)," % (m.name, p["ctype"], m.name)  # all trail a comma
def cpp_move_repoint(m, p): return "%s.%s=&%s;" % (p["ctype_fi"], m.name, m.name)

# Wiring statement, selected by the member's `wiring` trait: scalar-likes are
# aliased into C++-owned storage (C_F_POINTER); arrays are pointed at the
# Fortran-owned storage (C_LOC). Both the interop struct (`ftype_fi`) and the
# derived-type instance (`ftype`) come from the call site, never from the handle.
WIRING_STMT = {
    "cfptr": lambda m, p: "call C_F_POINTER(%s%%%s, %s%%%s)" % (p["ftype_fi"], m.name, p["ftype"], m.name),
    "cloc":  lambda m, p: "%s%%%s = C_LOC(%s%%%s)" % (p["ftype_fi"], m.name, p["ftype"], m.name),
}


# ---------------------------------------------------------------------------
# Kept renderers -- the few regions that do NOT fit the (filter, item, layout)
# shape: run-grouped prefixes or array bounds the trait table does not model. Each
# is `f(binding, params) -> body lines`, reading `.members` / the call-site params
# like every region. (The element count is not a per-member map either, but it is
# an inline VALUE -- a scalar substituted mid-line -- so it lives in VALUE_REGIONS.)
# ---------------------------------------------------------------------------
def r_cpp_mirror_fields(b, p):
    # Grouped declaration per consecutive same-pointer-type run (kept because the
    # line PREFIX varies per run). Base type from KIND_TRAITS; runs are by
    # CONSECUTIVE members, so emitted order follows the contract order.
    out = []
    for ctype, run in _runs(b.members, lambda m: KIND_TRAITS[m.kind].cpp_ptr_base):
        out += _decl_lines("        ", ctype + " ",
                           ["*%s=nullptr" % m.name for m in run], ";")
    return out


def r_cpp_array_views(b, p):
    # Bounds-bearing, one statement per line stays clearest. Element type is the
    # array kind's C++ base type from KIND_TRAITS; view class is the ARRAY_TMPL
    # constant; the C++ mirror member (`ctype_fi`) comes from the call site.
    out = []
    for m in b.members:
        if m.kind != "array":
            continue
        out.append("      %s = %s%dD<%s>(%s.%s, {%s}, {%s});"
                   % (m.name, ARRAY_TMPL, m.rank, KIND_TRAITS[m.kind].cpp_ptr_base,
                      p["ctype_fi"], m.name, ",".join(m.begin), ",".join(m.end)))
    return out


def r_cpp_storage_fields(b, p):
    # The C++ owner class's REAL fields, emitted VERBATIM from the contract body
    # (the single source of truth) so initializers (`int numrad = 2;`), trailing doc
    # comments, and the blank-line grouping survive exactly. Only the in-declaration
    # `(@NoahmpMacro:bounds_Nd(...))` annotation is stripped and the line re-indented
    # to the class body; guidance comment lines are dropped. The contract is written
    # in the owner language (C++), so this region is inherently C++.
    indent = "        "
    out = []
    for raw in b.body:
        core = INLINE_ANNOT_RE.sub("", raw).strip()
        if not core:
            out.append("")                               # blank separator preserved
        elif core.startswith("//"):
            continue                                     # guidance comment -> drop
        else:
            out.append(indent + core)
    return out


def r_fortran_array_allocate(b, p):
    # The Fortran-owned storage the C++ view points into. Generated from the SAME
    # bounds_Nd annotation that drives the C++ extent (r_cpp_array_views), so the
    # two cannot drift. The derived-type instance (`ftype`) comes from the call
    # site; the guard keeps it idempotent.
    out = []
    for m in b.members:
        if m.kind != "array":
            continue
        bounds = ", ".join("%s:%s" % (lo, hi) for lo, hi in zip(m.begin, m.end))
        out.append("    if ( .not. allocated (%s%%%s) ) allocate ( %s%%%s (%s) )"
                   % (p["ftype"], m.name, p["ftype"], m.name, bounds))
    return out


def r_fortran_mirror_fields(b, p):
    # Every member is an opaque pointer on the bind(C) side, so all kinds share one
    # field type -- one grouped declaration. Type from KIND_TRAITS (assert the
    # shared-type invariant rather than hardcode it).
    btypes = {KIND_TRAITS[m.kind].fort_bindc_type for m in b.members}
    assert len(btypes) == 1, "bind(C) field type is not uniform across kinds"
    return _decl_lines("    ", btypes.pop() + " :: ", [m.name for m in b.members], "")


def r_fortran_storage_fields(b, p):
    # The coupled storage in the owning derived type: int dims/handles and real
    # scalars are C-owned pointers (=> null() so association status is well-defined
    # before the C wiring runs); arrays are allocatables handed to C via C_LOC.
    # Scalars/ints pack into grouped declarations; arrays differ in rank, so each
    # gets its own line.
    out = []
    for kind, run in _runs(b.members, lambda m: m.kind):
        if kind in VARTYPE_PREFIX:                      # int / scalar: grouped, shared prefix
            out += _decl_lines("    ", VARTYPE_PREFIX[kind],
                               [KIND_TRAITS[kind].fort_storage(m) for m in run], "")
        else:                                           # array: self-contained line per member
            for m in run:
                out.append("    " + KIND_TRAITS[kind].fort_storage(m))
    return out


# Inline VALUE regions: substituted as a scalar mid-line (e.g.
# `... = @NoahmpMacro:MemberCount(<handle>);`), not expanded to a body block.
# `f(binding) -> str`. The convention shows WHERE the value is assigned -- the
# constant's declaration is hand-written; the macro supplies only the number.
VALUE_REGIONS = {
    "MemberCount": lambda b: str(len(b.members)),
}


# region name -> renderer (a Region spec or a kept function). PURPOSE-based,
# file-agnostic names, grouped by interop facet. Every key stays a named entry so
# `grep <RegionName>` lands on its definition.
REGIONS = {
    # -- Mirror: the flat interop struct on each side (element count -> VALUE_REGIONS).
    "CppMirrorFields":      r_cpp_mirror_fields,                                  # kept: run-grouped prefix
    "FortranMirrorFields":  r_fortran_mirror_fields,                             # kept: uniform-type assert
    # -- Construct: build / move the C++ owner and its mirror.
    "CppCtorParams":  Region(SCALAR_LIKE, cpp_ctor_param,   L_wrap("          ", comma=True)),
    "CppCtorInit":    Region(SCALAR_LIKE, cpp_ctor_init,    L_wrap("          ", comma=True)),
    "CppBindAddrs":   Region(SCALAR_LIKE, cpp_bind_addr,    L_wrap("            ", comma=True)),
    "CppMoveInit":    Region(ALL,         cpp_move_init,    L_wrap("          ")), # item bakes its comma
    "CppMoveRepoint": Region(SCALAR_LIKE, cpp_move_repoint, L_wrap("          ")),
    # -- Wire: physically connect the two representations (ownership/lifetime).
    "CppArrayViews":     r_cpp_array_views,                                       # kept: bounds-bearing
    "FortranScalarWire": Region(WIRED("cfptr"), WIRING_STMT["cfptr"], L_wrap_joined("    ", "; ")),
    "FortranArrayWire":  Region(WIRED("cloc"),  WIRING_STMT["cloc"],  L_wrap_joined("    ", "; ")),
    # -- Storage: each owning side's real fields + the Fortran allocation.
    "CppStorageFields":     r_cpp_storage_fields,                                  # kept: verbatim contract body
    "FortranStorageFields": r_fortran_storage_fields,
    "FortranArrayAllocate": r_fortran_array_allocate,
}


# ---------------------------------------------------------------------------
# 4. Binding: one resolved contract -- the HANDLE, parsed members, and the
# do-not-edit banner (built from the file the contract lives in, so it lands only
# in generated output, never in the tracked template).
#
# The HANDLE (`m_noahmpio`) is what markers reference. It carries NO derived
# identifiers: every site-local name a region emits is passed explicitly at the
# projection (see REGION_PARAMS), so the tool never assumes a variable name from
# the handle. `_display_name` (handle minus a leading `m_`) is used only for the
# human-readable breadcrumb comment, never to generate code.
# ---------------------------------------------------------------------------
def _display_name(handle):
    return handle[2:] if handle.startswith("m_") else handle


class Binding:
    __slots__ = ("tag", "members", "body", "source_file", "banner")

    def __init__(self, tag, members, body, source_file):
        self.tag = tag                        # the marker HANDLE, e.g. m_noahmpio
        self.members = members
        self.body = body                      # raw contract body lines (verbatim owner fields)
        self.source_file = source_file        # basename of the template holding the contract
        self.banner = ("generated by tools/NoahmpMacro.py from %s, do not edit"
                       % source_file)


# ---------------------------------------------------------------------------
# 5. Engine -- binding-agnostic. Two passes: collect_bindings finds every contract
# block (`@NoahmpMacro:Source <handle> { ... }`) as {handle -> Binding}; expand
# replaces every marker in a template.
# ---------------------------------------------------------------------------
def comment_for(path):
    """The target language's line comment, inferred from the template extension."""
    return "!" if path.endswith(".F90-mc") else "//"


def target_of(path):
    """The generated target for a `*-mc` template (strip the `-mc` suffix)."""
    assert path.endswith("-mc"), path
    return path[:-3]


def _close_re(cc):
    return re.compile(r"^\s*(?:%s\s*)?\}\s*$" % re.escape(cc))


def _find_source_blocks(lines, cc):
    """Return [(handle, begin_index, end_index)] for every contract block
    `@NoahmpMacro:Source <handle> { ... }` in `lines` (begin = opener, end = the
    matching bare `}`; field declarations carry no nested braces, so the first bare
    `}` after the opener is the closer). RIGID: a line whose first token is
    `@NoahmpMacro:Source` but which is not a well-formed opener is a hard error
    (catches a missing `{`, an old-style `Source(handle)`, or a stray `;`)."""
    close_re = _close_re(cc)
    blocks = []
    for i, l in enumerate(lines):
        if not SOURCE_LEAD_RE.match(l):
            continue
        m = SOURCE_OPEN_RE.match(l)
        if not m:
            raise SystemExit("NoahmpMacro: malformed %s%s line (want "
                             "`%s%s <handle> {` opening a brace block):\n  %s"
                             % (MARKER, SOURCE, MARKER, SOURCE, l.strip()))
        for j in range(i + 1, len(lines)):
            if close_re.match(lines[j]):
                blocks.append((m.group(1), i, j))
                break
        else:
            raise SystemExit("NoahmpMacro: no closing '}' for %s%s %s {"
                             % (MARKER, SOURCE, m.group(1)))
    return blocks


# A block projection marker on its own line:
#   <indent>@NoahmpMacro:<Region>(<handle>[, key=value]...);
# The arg list is captured whole and split by _parse_marker_args.
CALL_RE = re.compile(r"^(\s*)%s(\w+)\(([^)]*)\)\s*;\s*$" % re.escape(MARKER))

# Any colon-form projection marker `@NoahmpMacro:<Region>(<args>)` anywhere in a
# line. After own-line BLOCK markers are expanded, this is the single authoritative
# classifier for every remaining marker (value substitution, or a precise error for
# an inline block region / unknown region / Source-as-call): see expand().
PROJ_RE = re.compile(r"%s(\w+)\(([^)]*)\)" % re.escape(MARKER))

# A line whose FIRST token is `@NoahmpMacro:Source` MUST be a well-formed block
# opener (SOURCE_OPEN_RE). Prose mentions never lead with it (they sit after a
# comment char), so this cleanly flags a malformed/old-style Source line.
SOURCE_LEAD_RE = re.compile(r"^\s*%s%s\b" % (re.escape(MARKER), SOURCE))

_IDENT_RE = re.compile(r"^[A-Za-z_]\w*$")            # a bare identifier (handle / local name)


def _parse_marker_args(argstr, where):
    """Split a projection marker's argument list into (handle, {key: value}) with
    rigid validation: the first arg is the positional handle (an identifier); each
    remaining arg is `key=value` with identifier key and value; duplicate keys are
    rejected. `where` names the marker for error messages."""
    parts = [a.strip() for a in argstr.split(",")]
    handle = parts[0]
    if not handle:
        raise SystemExit("NoahmpMacro: %s is missing its binding handle" % where)
    if not _IDENT_RE.match(handle):
        raise SystemExit("NoahmpMacro: %s: handle %r is not a bare identifier"
                         % (where, handle))
    params = {}
    for a in parts[1:]:
        if not a:
            raise SystemExit("NoahmpMacro: %s: empty argument (stray comma?)" % where)
        if "=" not in a:
            raise SystemExit("NoahmpMacro: %s: argument %r must be key=value" % (where, a))
        k, v = (s.strip() for s in a.split("=", 1))
        if not _IDENT_RE.match(k):
            raise SystemExit("NoahmpMacro: %s: argument key %r is not a bare identifier"
                             % (where, k))
        if not _IDENT_RE.match(v):
            raise SystemExit("NoahmpMacro: %s: value of %s=%r is not a bare identifier"
                             % (where, k, v))
        if k in params:
            raise SystemExit("NoahmpMacro: %s: duplicate argument %r" % (where, k))
        params[k] = v
    return handle, params


def _check_params(region, params, where):
    """Validate the call site supplies EXACTLY the params `region` requires -- a
    missing key would miscompile, an unknown key is a typo. Both are named errors."""
    required = set(REGION_PARAMS.get(region, []))
    got = set(params)
    if got - required:
        raise SystemExit("NoahmpMacro: %s: unexpected argument(s) %s (expected %s)"
                         % (where, sorted(got - required), sorted(required) or "none"))
    if required - got:
        raise SystemExit("NoahmpMacro: %s: missing argument(s) %s"
                         % (where, sorted(required - got)))


def collect_bindings(templates):
    """Pass 1: scan every template for its contract block
    `@NoahmpMacro:Source <handle> { ... }` and build {handle -> Binding}. A handle
    is declared once across all templates; its body is kept for verbatim emission."""
    bindings = {}
    for path in templates:
        cc = comment_for(path)
        with open(path) as f:
            lines = f.read().splitlines()
        for handle, bi, ei in _find_source_blocks(lines, cc):
            if handle in bindings:
                raise SystemExit("NoahmpMacro: duplicate %s%s %s { ... }"
                                 % (MARKER, SOURCE, handle))
            body = lines[bi + 1:ei]
            members = parse_members(body, cc, handle)
            bindings[handle] = Binding(handle, members, body, os.path.basename(path))
    if not bindings:
        raise SystemExit("NoahmpMacro: no %s%s <handle> { ... } contract block found"
                         % (MARKER, SOURCE))
    return bindings


# An in-declaration annotation group, parens and all: `(@NoahmpMacro:bounds_2d(...))`.
# Stripped from a field declaration (in CppStorageFields) so the owner class sees a
# plain member decl. Single `@NoahmpMacro:` surface -- no separate punctuator.
INLINE_ANNOT_RE = re.compile(r"\(\s*%sbounds_\w+\([^)]*\)\s*\)" % re.escape(MARKER))


def expand(text, cc, bindings):
    """Pass 2 for one template. Line-based edits first: each `Source <handle> {...}`
    contract block (-> a one-line breadcrumb), and every own-line BLOCK projection
    marker (-> banner + rendered body + close). Then ONE authoritative text-wide
    pass over every remaining `@NoahmpMacro:` marker (PROJ_RE): a VALUE region is
    substituted in place; anything else (a block region used inline, a `Source`
    call, an unknown region) is a precise error."""
    lines = text.splitlines()
    edits = []                                           # (begin, end, new_lines)

    # The contract block is the single source of truth, consumed here: the class
    # materializes its fields via @NoahmpMacro:CppStorageFields. Leave a breadcrumb.
    for handle, bi, ei in _find_source_blocks(lines, cc):
        b = bindings[handle]
        indent = re.match(r"\s*", lines[bi]).group(0)
        edits.append((bi, ei, ["%s%s %s binding (%s) -- %s"
                               % (indent, cc, _display_name(handle), handle, b.banner)]))

    for i, line in enumerate(lines):
        m = CALL_RE.match(line)
        if not m:
            continue
        indent, region = m.group(1), m.group(2)
        if region not in REGIONS:
            continue                                     # value/unknown -> classified below
        where = "%s%s(...)" % (MARKER, region)
        handle, params = _parse_marker_args(m.group(3), where)
        if handle not in bindings:
            raise SystemExit("NoahmpMacro: marker %s%s(%s) references an unknown "
                             "binding handle" % (MARKER, region, handle))
        _check_params(region, params, where)
        b = bindings[handle]
        edits.append((i, i,
                      ["%s%s %s -- %s" % (indent, cc, region, b.banner)]
                      + REGIONS[region](b, params)
                      + ["%s%s end %s" % (indent, cc, region)]))

    # Apply bottom-up so earlier indices stay valid.
    for bi, ei, new in sorted(edits, key=lambda e: e[0], reverse=True):
        lines = lines[:bi] + new + lines[ei + 1:]
    out = "\n".join(lines)

    # Authoritative classifier for every marker not consumed as an own-line block.
    def _classify(m):
        region, argstr = m.group(1), m.group(2)
        where = "%s%s(...)" % (MARKER, region)
        handle, params = _parse_marker_args(argstr, where)
        if region in VALUE_REGIONS:
            if params:
                raise SystemExit("NoahmpMacro: %s: value region takes only a handle, "
                                 "got %s" % (where, sorted(params)))
            if handle not in bindings:
                raise SystemExit("NoahmpMacro: %s references an unknown binding "
                                 "handle %r" % (where, handle))
            return VALUE_REGIONS[region](bindings[handle])
        if region in REGIONS:
            raise SystemExit("NoahmpMacro: block region %s%s used inline; it must be "
                             "on its own line as `%s%s(<handle>, ...);`"
                             % (MARKER, region, MARKER, region))
        if region == SOURCE:
            raise SystemExit("NoahmpMacro: %s%s is a contract block, not a call; "
                             "write `%s%s <handle> { ... }`" % (MARKER, SOURCE, MARKER, SOURCE))
        raise SystemExit("NoahmpMacro: unknown region %s%s (known block: %s; value: %s)"
                         % (MARKER, region, sorted(REGIONS), sorted(VALUE_REGIONS)))
    out = PROJ_RE.sub(_classify, out)

    if text.endswith("\n"):
        out += "\n"
    return out


# A live marker is a region/annotation name immediately followed by `(` (call) or
# `{` (block). Matching the GRAMMAR -- not the bare `@NoahmpMacro:` string -- lets
# author prose that survives into generated output (e.g. "generated from
# @NoahmpMacro:Source and...", the static_assert "@NoahmpMacro:Source; run make
# codegen") pass untouched: those are followed by `)`/`,`/`;`/space, never `(`/`{`.
MARKER_RE = re.compile(r"%s\w+\s*[({]" % re.escape(MARKER))


def _assert_no_markers(texts):
    """Fail if any marker survived expansion -- the signature of an unknown/
    misspelled region name, which expand() silently leaves in place. Reports every
    hit as file:line."""
    bad = []
    for path, text in texts.items():
        for n, line in enumerate(text.splitlines(), 1):
            if MARKER_RE.search(line):
                bad.append("  %s:%d: %s" % (os.path.basename(path), n, line.strip()))
    if bad:
        raise SystemExit("NoahmpMacro: unexpanded %s marker(s) left in generated "
                         "output (unknown/misspelled region name?):\n" % MARKER
                         + "\n".join(bad))


def _validate_registry():
    """Self-check the region tables before running -- guards future extension
    mistakes: a region cannot be both block and value, and REGION_PARAMS may only
    name block regions (value regions take only the handle)."""
    both = set(REGIONS) & set(VALUE_REGIONS)
    if both:
        raise SystemExit("NoahmpMacro: region(s) %s are in both REGIONS and "
                         "VALUE_REGIONS (a region is block XOR value)" % sorted(both))
    orphan = set(REGION_PARAMS) - set(REGIONS)
    if orphan:
        raise SystemExit("NoahmpMacro: REGION_PARAMS names non-block region(s) %s "
                         "(only block regions take params)" % sorted(orphan))


def main(argv):
    check = "--check" in argv[1:]
    _validate_registry()
    templates = sorted(glob.glob(os.path.join(ROOT, "*-mc")))
    bindings = collect_bindings(templates)

    new_texts = {}
    for path in templates:
        with open(path) as f:
            text = f.read()
        new_texts[target_of(path)] = expand(text, comment_for(path), bindings)
    _assert_no_markers(new_texts)
    nmembers = sum(len(b.members) for b in bindings.values())

    changed = []
    for path, new in new_texts.items():
        try:
            with open(path) as f:
                old = f.read()
        except FileNotFoundError:
            old = None                                   # target not generated yet
        if old != new:
            changed.append(path)
            if check:
                sys.stdout.writelines(difflib.unified_diff(
                    (old or "").splitlines(keepends=True), new.splitlines(keepends=True),
                    fromfile=path + (" (committed)" if old is not None else " (missing)"),
                    tofile=path + " (regenerated)"))

    if check:
        if changed:
            print("\nNoahmpMacro: %d generated file(s) are out of sync with their "
                  "%sSource contract: %s\nRun `make codegen` and commit the result."
                  % (len(changed), MARKER, ", ".join(os.path.basename(p) for p in changed)),
                  file=sys.stderr)
            return 1
        print("NoahmpMacro: %d coupling members; all generated regions are in sync."
              % nmembers)
        return 0

    for path in changed:
        with open(path, "w") as f:
            f.write(new_texts[path])
    print("NoahmpMacro: %d coupling members; updated %d file(s)%s."
          % (nmembers, len(changed),
             ": " + ", ".join(os.path.basename(p) for p in changed) if changed else " (already in sync)"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
