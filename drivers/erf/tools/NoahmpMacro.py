#!/usr/bin/env python3
# =============================================================================
# NoahmpMacro.py -- expand the Noah-MP <-> ERF coupling member list into every
# place the C<->Fortran ABI needs it.
# =============================================================================
# SINGLE SOURCE OF TRUTH: the member declarations of `class NoahmpIO_type` in
# NoahmpIO.H-mc, wrapped in the `@NoahmpMacro:Source { ... }` marker block.
# Their ORDER is the canonical ABI order. Member kind is inferred from the C++
# type and annotated with a `@NoahmpMacro:` macro in a trailing comment:
#
#     int <names>;                                                       -> int handle(s) (no annotation)
#     noahmp_real NAME [= init];  // @NoahmpMacro:scalar [doc="..."]          -> real scalar
#     NoahmpArray2D<noahmp_real> N; // @NoahmpMacro:bounds_2d(lo:hi, lo:hi) [doc="..."]        -> 2-D array
#     NoahmpArray3D<noahmp_real> N; // @NoahmpMacro:bounds_3d(lo:hi, lo:hi, lo:hi) [doc="..."] -> 3-D array
#
# From that one list this script regenerates the bodies of all `NoahmpMacro:<name>`
# regions in the five generated files -- the mirrored fi struct, the constructors,
# the move constructor, NOAHMP_IO_FI_NUM_MEMBERS, the Fortran bind(C) type, the
# coupled members of the NoahmpIO_type storage type, the C_LOC / C_F_POINTER
# wiring, and the Fortran allocate() of the coupled arrays. Because the C++ struct
# and the Fortran bind(C) type are emitted from the SAME ordered list, their member
# count and order are identical by construction -- so the old runtime size and
# member-order ABI handshakes were redundant and have been removed; only the (free)
# compile-time layout static_assert and the runtime precision check (a build-flag
# hazard codegen cannot prevent) remain.
#
# The Fortran-owned allocate() of each coupled array (in NoahmpIOVarInitMod.F90) is
# emitted from the SAME bounds_Nd annotation that drives the C++ NoahmpArray extent,
# so the two can no longer disagree by construction. (This replaces an older
# read-only audit that compared a hand-written allocate against the annotation.)
#
# Stdlib only (regex). No TOML, no third-party libraries, any Python 3.
#
# Only the tracked `*-mc` TEMPLATES are edited by hand:
#     NoahmpIO.H-mc, NoahmpIO.cpp-mc, NoahmpIO_fi.F90-mc, NoahmpIOVarType.F90-mc,
#     NoahmpIOVarInitMod.F90-mc
# Each carries bare, function-call `@NoahmpMacro:Region();` markers (PascalCase,
# like the codebase's DriverMain / ReadLandMain methods). The lone block marker is
# `@NoahmpMacro:Source { ... }`: an un-commented macro whose body is the real
# member declarations. The generator strips its wrapper and the per-line
# `@NoahmpMacro:` annotations (keeping any "doc" as a plain comment) so no
# marker reaches the compiled header. This script expands every marker into its
# do-not-edit body and writes the five GENERATED, gitignored targets (NoahmpIO.H,
# NoahmpIO.cpp, NoahmpIO_fi.F90, NoahmpIOVarType.F90, NoahmpIOVarInitMod.F90). The
# build runs this before compiling, so the targets are never committed -- only the
# templates are. The
# generated list-regions are packed several members per line (fixed-width wrap);
# the targets are regenerated every build, so density costs no diff noise.
#
# Usage:
#   python3 tools/NoahmpMacro.py            # (re)generate the targets
#   python3 tools/NoahmpMacro.py --check    # fail (exit 1) if regen would change
#
# Adding a coupled variable: add ONE annotated line to the @NoahmpMacro:Source
# block in NoahmpIO.H-mc and rebuild -- for an array the bounds_Nd annotation also
# generates its Fortran allocate(). (Any namelist guard remains manual; see
# dev/spec-add-coupled-variable.md)
# =============================================================================

import os
import re
import sys
import difflib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)            # drivers/erf

# Each generated file is produced from a tracked `*-mc` TEMPLATE into a gitignored
# TARGET that the build (re)creates every time. The templates carry only bare
# `@NoahmpMacro:Region();` markers; the generator expands each into the
# do-not-edit body below. The single source of truth (the @NoahmpMacro:Source
# member list) lives in NoahmpIO.H-mc.
H_MC_FILE      = os.path.join(ROOT, "NoahmpIO.H-mc")
H_FILE         = os.path.join(ROOT, "NoahmpIO.H")
CPP_MC_FILE    = os.path.join(ROOT, "NoahmpIO.cpp-mc")
CPP_FILE       = os.path.join(ROOT, "NoahmpIO.cpp")
F90_MC_FILE    = os.path.join(ROOT, "NoahmpIO_fi.F90-mc")
F90_FILE       = os.path.join(ROOT, "NoahmpIO_fi.F90")
VARTYPE_MC_FILE = os.path.join(ROOT, "NoahmpIOVarType.F90-mc")
VARTYPE_FILE    = os.path.join(ROOT, "NoahmpIOVarType.F90")
VARINIT_MC_FILE = os.path.join(ROOT, "NoahmpIOVarInitMod.F90-mc")
VARINIT_FILE    = os.path.join(ROOT, "NoahmpIOVarInitMod.F90")

# Banner stamped onto every generated region's opening marker in the TARGET. The
# templates carry only the bare function-call marker; this is appended by the
# generator so the warning lands only in the generated output, never in the template.
BANNER = "generated by tools/NoahmpMacro.py from NoahmpIO.H-mc, do not edit"


# ---------------------------------------------------------------------------
# Parse the @NoahmpMacro:Source block of NoahmpIO.H into an ordered member list.
# ---------------------------------------------------------------------------
class Member:
    __slots__ = ("name", "kind", "rank", "begin", "end", "doc")

    def __init__(self, name, kind, rank=0, begin=None, end=None, doc=""):
        self.name = name
        self.kind = kind            # 'int' | 'scalar' | 'array'
        self.rank = rank
        self.begin = begin or []
        self.end = end or []
        self.doc = doc


def _slice_markers(lines, region, cc):
    """Return (begin_index, end_index) for a region's marker in a template. Two
    forms are accepted:

      * function-call (the convention for every expandable region). The template
        is never compiled directly, so the marker carries no comment prefix and
        reads like a call:
            <indent>@NoahmpMacro:Region();
        here begin_index == end_index (the single line is both open and close);
      * block (only @NoahmpMacro:Source, whose body is the real member
        declarations). It is an un-commented macro too -- the generator strips the
        wrapper (and the per-line annotations) so nothing reaches the compiled
        output:
            <indent>@NoahmpMacro:Region { ...
            ...
            <indent>}
        terminated by the next bare `}` line. The optional `cc` prefix is still
        accepted so an older commented block (`cc @...Region { ... cc }`) keeps
        parsing.

    Errors if the marker is missing/duplicated or a block opener has no matching
    closer."""
    call_re  = re.compile(r"^\s*@NoahmpMacro:%s\s*\(\s*\)\s*;\s*$" % re.escape(region))
    block_re = re.compile(r"^\s*(?:%s\s*)?@NoahmpMacro:%s\b.*\{" % (re.escape(cc), re.escape(region)))
    close_re = re.compile(r"^\s*(?:%s\s*)?\}\s*$" % re.escape(cc))
    hits = [(i, "call") for i, l in enumerate(lines) if call_re.match(l)]
    hits += [(i, "block") for i, l in enumerate(lines) if block_re.search(l)]
    if len(hits) != 1:
        raise SystemExit("NoahmpMacro: expected exactly one "
                         "'@NoahmpMacro:%s();' marker (found %d)"
                         % (region, len(hits)))
    b, form = hits[0]
    if form == "call":
        return b, b                                    # function-call -- empty region
    for j in range(b + 1, len(lines)):
        if close_re.match(lines[j]):
            return b, j
    raise SystemExit("NoahmpMacro: no closing '%s }' marker for "
                     "@NoahmpMacro:%s" % (cc, region))


def parse_source(h_text):
    lines = h_text.splitlines()
    bi, ei = _slice_markers(lines, "Source", "//")
    members = []
    seen = set()
    for raw in lines[bi + 1:ei]:
        if raw.strip().startswith("//"):
            continue                                   # comment-only line (usage note, closing brace)
        # Split the declaration from its trailing // comment; the annotation, if
        # any, is a @NoahmpMacro: macro inside that comment.
        if "//" in raw:
            code, annot = raw.split("//", 1)
        else:
            code, annot = raw, ""
        code = code.strip()
        if not code or not code.endswith(";"):
            continue                                   # blank / comment-only line
        code = code[:-1].strip()                        # drop trailing ';'

        if code.startswith("NoahmpArray"):
            m = re.match(r"NoahmpArray([23])D<noahmp_real>\s+(\w+)$", code)
            if not m:
                raise SystemExit("NoahmpMacro: cannot parse array decl: %r" % raw)
            rank, name = int(m.group(1)), m.group(2)
            mb = re.search(r"@NoahmpMacro:bounds_([23])d\s*\(([^)]*)\)", annot)
            if not mb:
                raise SystemExit("NoahmpMacro: array %s needs a "
                                 "@NoahmpMacro:bounds_%dd(lo:hi, ...) annotation" % (name, rank))
            if int(mb.group(1)) != rank:
                raise SystemExit("NoahmpMacro: array %s is %d-D but annotated bounds_%sd"
                                 % (name, rank, mb.group(1)))
            pairs = [p.strip() for p in mb.group(2).split(",") if p.strip()]
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
            members.append(Member(name, "array", rank, begin, end, _doc(annot, name)))

        elif code.startswith("noahmp_real"):
            m = re.match(r"noahmp_real\s+(\w+)\s*(=.*)?$", code)
            if not m:
                raise SystemExit("NoahmpMacro: cannot parse scalar decl: %r" % raw)
            sm = re.search(r"@NoahmpMacro:scalar\b(.*)", annot)
            if not sm:
                raise SystemExit("NoahmpMacro: scalar %s needs a "
                                 "@NoahmpMacro:scalar [doc=\"...\"] annotation"
                                 % m.group(1))
            members.append(Member(m.group(1), "scalar", doc=_doc(annot, m.group(1))))

        elif code.startswith("int"):
            for nm in code[3:].split(","):
                nm = nm.split("=")[0].strip().lstrip("*").strip()   # drop any `= default`
                if nm:
                    members.append(Member(nm, "int"))
        else:
            raise SystemExit("NoahmpMacro: unrecognised member decl in @NoahmpMacro:Source: %r" % raw)

    for m in members:
        if m.name in seen:
            raise SystemExit("NoahmpMacro: duplicate member name %r" % m.name)
        seen.add(m.name)
    if not members:
        raise SystemExit("NoahmpMacro: @NoahmpMacro:Source block is empty")

    # Validate every array bound token now that all members are known (so it is
    # independent of where the int dims sit relative to the arrays). Each token is
    # emitted verbatim into BOTH the C++ NoahmpArray extent and the Fortran
    # allocate(), so it must resolve to the SAME integer on both sides: an integer
    # literal (negative lower bounds allowed) or a declared int member. This turns
    # a typo'd/renamed dimension from an opaque downstream compile error into a
    # named SystemExit here. (It checks token IDENTITY, not whether the int is
    # semantically a dimension, nor the dimension order.)
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


def _doc(annot, where):
    """Extract the optional doc string. STRICT: a doc MUST be written as
    doc="..."; a bare quoted string is rejected so the annotation grammar stays
    unambiguous. `where` names the member for the error message."""
    m = re.search(r'\bdoc\s*=\s*"([^"]*)"', annot)
    doc = m.group(1) if m else ""
    leftover = re.sub(r'\bdoc\s*=\s*"[^"]*"', "", annot)   # strip the valid doc=, keep the rest
    if '"' in leftover:
        raise SystemExit('NoahmpMacro: %s: a doc string must be written as '
                         'doc="..." (a bare "..." is not allowed)' % where)
    return doc


# ---------------------------------------------------------------------------
# Rendering helpers
#
# The generated targets are gitignored and rewritten every build, so packing
# several members per line costs no diff noise -- readability is the only goal.
# Every list-style region is therefore wrapped to a fixed right margin instead
# of exploded one-member-per-line.
# ---------------------------------------------------------------------------
WRAP_COL = 100          # soft right margin for packed lines


def _is_scalar_like(m):      # participates in the constructor / scalar wiring
    return m.kind in ("int", "scalar")


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
# Renderers -- each returns the list of body lines for one region.
# ---------------------------------------------------------------------------
def r_fi_members(ms):
    # Grouped declaration per consecutive same-pointer-type run (every member is a
    # pointer, so int vs noahmp_real is the only split). Docs live in the Source
    # block; the packed mirror drops them.
    out = []
    for ctype, run in _runs(ms, lambda m: "int" if m.kind == "int" else "noahmp_real"):
        out += _decl_lines("        ", ctype + " ",
                           ["*%s=nullptr" % m.name for m in run], ";")
    return out


def r_num_members(ms):
    return ["inline constexpr std::size_t NOAHMP_IO_FI_NUM_MEMBERS = %d;" % len(ms)]


def r_fi_ctor_params(ms):
    items = _commad(["%s %s" % ("int*" if m.kind == "int" else "noahmp_real*", m.name)
                     for m in ms if _is_scalar_like(m)])
    return _wrap(items, "          ")


def r_fi_ctor_init(ms):
    items = _commad(["%s(%s)" % (m.name, m.name) for m in ms if _is_scalar_like(m)])
    return _wrap(items, "          ")


def r_type_ctor_fptr(ms):
    items = _commad(["&%s" % m.name for m in ms if _is_scalar_like(m)])
    return _wrap(items, "            ")


def r_move_init_list(ms):
    # Every member is copied; fptr(o.fptr) follows the region, so all trail a comma.
    items = ["%s(o.%s)," % (m.name, m.name) for m in ms]
    return _wrap(items, "          ")


def r_move_fptr_repoint(ms):
    items = ["fptr.%s=&%s;" % (m.name, m.name) for m in ms if _is_scalar_like(m)]
    return _wrap(items, "          ")


def r_varinit_arrays(ms):
    # Long, individually-bounded statements: one per line stays clearest.
    out = []
    for m in ms:
        if m.kind != "array":
            continue
        out.append("      %s = NoahmpArray%dD<noahmp_real>(fptr.%s, {%s}, {%s});"
                   % (m.name, m.rank, m.name, ",".join(m.begin), ",".join(m.end)))
    return out


def r_varinit_allocate(ms):
    # The Fortran-owned storage the C++ NoahmpArray points into. Generated from the
    # SAME bounds_Nd annotation that drives the C++ extent (r_varinit_arrays), so
    # the two can no longer drift -- this region replaces the old read-only
    # allocate audit. Emitted inside NoahmpIOVarInitDefault's `associate`, so the
    # bound names (xstart, kms, numrad, ...) resolve via its aliases; the guard
    # keeps it idempotent. One statement per line stays clearest.
    out = []
    for m in ms:
        if m.kind != "array":
            continue
        bounds = ", ".join("%s:%s" % (lo, hi) for lo, hi in zip(m.begin, m.end))
        doc = "  ! " + m.doc if m.doc else ""
        out.append("    if ( .not. allocated (NoahmpIO%%%s) ) allocate ( NoahmpIO%%%s (%s) )%s"
                   % (m.name, m.name, bounds, doc))
    return out


def r_fi_type_members(ms):
    return _decl_lines("    ", "type(C_PTR) :: ", [m.name for m in ms], "")


def r_scalarinit_cfptr(ms):
    items = ["call C_F_POINTER(NoahmpIO_cptr%%%s, blk%%%s)" % (m.name, m.name)
             for m in ms if _is_scalar_like(m)]
    return _wrap_joined(items, "    ", "; ")


def r_varinit_cloc(ms):
    items = ["NoahmpIO_cptr%%%s = C_LOC(blk%%%s)" % (m.name, m.name)
             for m in ms if m.kind == "array"]
    return _wrap_joined(items, "    ", "; ")


def r_vartype_members(ms):
    # The coupled storage in the NoahmpIO_type derived type: int dims/handles and
    # real scalars are C-owned pointers (=> null() so their association status is
    # well-defined before the C wiring runs); arrays are allocatables handed to C
    # via C_LOC. Order is irrelevant -- NoahmpIO_type is not bind(C). Scalars and
    # int handles pack into grouped declarations; arrays differ in rank, so each
    # gets its own line.
    out = []
    for kind, run in _runs(ms, lambda m: m.kind):
        if kind == "int":
            out += _decl_lines("    ", "integer(C_INT), pointer :: ",
                               ["%s=>null()" % m.name for m in run], "")
        elif kind == "scalar":
            out += _decl_lines("    ", "real(kind=c_kind_noahmp), pointer :: ",
                               ["%s=>null()" % m.name for m in run], "")
        else:
            for m in run:
                dims = ",".join([":"] * m.rank)
                out.append("    real(kind=c_kind_noahmp), allocatable, dimension(%s) :: %s"
                           % (dims, m.name))
    return out


# region name -> renderer. PascalCase to match the codebase's method naming
# (DriverMain, ReadLandMain, ...). Order is cosmetic.
H_REGIONS = {
    "FiMembers":       r_fi_members,
    "NumMembers":      r_num_members,
    "FiCtorParams":    r_fi_ctor_params,
    "FiCtorInit":      r_fi_ctor_init,
    "TypeCtorFptr":    r_type_ctor_fptr,
    "MoveInitList":    r_move_init_list,
    "MoveFptrRepoint": r_move_fptr_repoint,
}
CPP_REGIONS = {
    "VarInitArrays":   r_varinit_arrays,
}
F90_REGIONS = {
    "FiTypeMembers":   r_fi_type_members,
    "ScalarInitCfptr": r_scalarinit_cfptr,
    "VarInitCloc":     r_varinit_cloc,
}
VARTYPE_REGIONS = {
    "VarTypeMembers":  r_vartype_members,
}
VARINIT_REGIONS = {
    "VarInitAllocate": r_varinit_allocate,
}


# ---------------------------------------------------------------------------
# Region replacement.
# ---------------------------------------------------------------------------
def apply_regions(text, regions, members, cc):
    lines = text.splitlines()
    # Replace from the bottom up so earlier indices stay valid.
    for region in sorted(regions, key=lambda r: _slice_markers(lines, r, cc)[0], reverse=True):
        bi, ei = _slice_markers(lines, region, cc)
        # Expand the `@NoahmpMacro:Region();` template marker into the generated
        # block: an opening comment marker carrying the GENERATED/do-not-edit banner,
        # the rendered body, and a matching closing comment marker -- all at the
        # template marker's own indentation. The banner therefore lands only on the
        # TARGET, never on the tracked template.
        indent = re.match(r"\s*", lines[bi]).group(0)
        open_marker  = "%s%s %s -- %s" % (indent, cc, region, BANNER)
        close_marker = "%s%s end %s" % (indent, cc, region)
        body = regions[region](members)
        lines = lines[:bi] + [open_marker] + body + [close_marker] + lines[ei + 1:]
    out = "\n".join(lines)
    if text.endswith("\n"):
        out += "\n"
    return out


def render_source_block(text):
    """Expand the @NoahmpMacro:Source { ... } block of the H template into the
    bare member declarations for the compiled header. Unlike the call-form regions,
    the body here is the single source of truth, so it is emitted (near-)verbatim --
    only the macro wrapper and the per-line `@NoahmpMacro:` annotations are
    stripped (any "doc" string is kept as a plain comment), leaving valid C++ with
    initializers (e.g. `= -9999.0`) and the blank-line grouping intact. No marker
    survives into NoahmpIO.H."""
    lines = text.splitlines()
    bi, ei = _slice_markers(lines, "Source", "//")
    indent = re.match(r"\s*", lines[bi]).group(0)
    body = []
    for raw in lines[bi + 1:ei]:
        if raw.lstrip().startswith("//"):
            continue                                   # guidance/usage comment -> drop
        if "//" in raw:
            code, annot = raw.split("//", 1)
            doc = _doc(annot, code.strip())
            body.append(code.rstrip() + ("  // " + doc if doc else ""))
        else:
            body.append(raw.rstrip())                   # declaration or blank separator
    lines = lines[:bi] + [indent + "// " + BANNER] + body + lines[ei + 1:]
    out = "\n".join(lines)
    if text.endswith("\n"):
        out += "\n"
    return out


# A live marker is a region name immediately followed by `(` (call form) or `{`
# (block form). Matching the GRAMMAR -- not the bare `@NoahmpMacro:` string -- lets
# the legit author prose that survives into generated output (e.g. "generated from
# @NoahmpMacro:Source) and...", the static_assert "@NoahmpMacro:Source; run make
# codegen") pass untouched: those are followed by `)`/`,`/`;`/space, never `(`/`{`.
MARKER_RE = re.compile(r"@NoahmpMacro:\w+\s*[({]")


def _assert_no_markers(texts):
    """Fail if any @NoahmpMacro marker survived expansion -- the signature of an
    unknown/misspelled region name, which apply_regions silently leaves in place
    (it only iterates the KNOWN region dicts). Reports every hit as file:line."""
    bad = []
    for path, text in texts.items():
        for n, line in enumerate(text.splitlines(), 1):
            if MARKER_RE.search(line):
                bad.append("  %s:%d: %s" % (os.path.basename(path), n, line.strip()))
    if bad:
        raise SystemExit("NoahmpMacro: unexpanded @NoahmpMacro marker(s) left in "
                         "generated output (unknown/misspelled region name?):\n"
                         + "\n".join(bad))


def process(members):
    """Return {target_path: new_text} for every generated file, each expanded from
    its tracked `*-mc` template."""
    with open(H_MC_FILE) as f:
        h = f.read()
    with open(CPP_MC_FILE) as f:
        cpp = f.read()
    with open(F90_MC_FILE) as f:
        f90 = f.read()
    with open(VARTYPE_MC_FILE) as f:
        vartype = f.read()
    with open(VARINIT_MC_FILE) as f:
        varinit = f.read()
    h = render_source_block(h)          # strip the Source macro -> bare declarations
    out = {
        H_FILE:       apply_regions(h,       H_REGIONS,       members, "//"),
        CPP_FILE:     apply_regions(cpp,     CPP_REGIONS,     members, "//"),
        F90_FILE:     apply_regions(f90,     F90_REGIONS,     members, "!"),
        VARTYPE_FILE: apply_regions(vartype, VARTYPE_REGIONS, members, "!"),
        VARINIT_FILE: apply_regions(varinit, VARINIT_REGIONS, members, "!"),
    }
    _assert_no_markers(out)             # catch unknown/misspelled region markers
    return out


def main(argv):
    check = "--check" in argv[1:]
    with open(H_MC_FILE) as f:
        members = parse_source(f.read())

    new_texts = process(members)
    changed = []
    for path, new in new_texts.items():
        try:
            with open(path) as f:
                old = f.read()
        except FileNotFoundError:
            old = None                                  # target not generated yet
        if old != new:
            changed.append(path)
            if check:
                sys.stdout.writelines(difflib.unified_diff(
                    (old or "").splitlines(keepends=True), new.splitlines(keepends=True),
                    fromfile=path + (" (committed)" if old is not None else " (missing)"),
                    tofile=path + " (regenerated)"))
        elif not check:
            pass

    if check:
        if changed:
            print("\nNoahmpMacro: %d file(s) are out of sync with the @NoahmpMacro:Source "
                  "block in NoahmpIO.H: %s\nRun `make codegen` and commit the result."
                  % (len(changed), ", ".join(os.path.basename(p) for p in changed)),
                  file=sys.stderr)
            return 1
        print("NoahmpMacro: %d coupling members; all generated regions are in sync."
              % len(members))
        return 0

    for path in changed:
        with open(path, "w") as f:
            f.write(new_texts[path])
    print("NoahmpMacro: %d coupling members; updated %d file(s)%s."
          % (len(members), len(changed),
             ": " + ", ".join(os.path.basename(p) for p in changed) if changed else " (already in sync)"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
