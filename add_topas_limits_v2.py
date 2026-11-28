"""
add_topas_limits.py

Inserts missing min/max limits into a TOPAS .inp file for:

  * CS / CS_L / CS_G        -> min 32, max 1000   (CS is shorthand for CS_L)
  * LVol_FWHM_CS_G_L(...)   -> its CS_G and CS_L arguments, as above
  * e0_from_Strain(...)     -> its Strain_G and Strain_L arguments, as below
  * PO(...) / Preferred_Orientation(...)   -> min 0.5, max 1.5 (March-Dollase)
  * PO_Two_Directions(...)  -> both March-Dollase values min 0.5, max 1.5;
                               fraction in the first direction min 0, max 1
  * Smart PO (optional, --smart-po): str phases with no preferred-orientation model get a fixed
    (not refined) March-Dollase line in the likely direction, chosen from
    the phase name (po_directions.csv next to the script, plus a built-in
    clay list) or, for unnamed layer silicates, from the cell:
        PO(!po_<phase>, 1 min 0.5 max 1.5, , h k l)
    Value 1 has no effect; change ! to a parameter name or @ to refine.
    The directions are suggestions: the user must confirm each one is correct.
  * Strain_L / Strain_G     -> min 0,  max 0.5
  * Unit cell (a b c al be ga, and the Cubic/Tetragonal/Hexagonal/Trigonal/
    Rhombohedral macros) in every str / hkl_Is phase
                            -> start value -1 % / +1 %

Rules
  * An existing min or max is never changed; only a missing one is inserted.
  * Unit cell: a, b, c are always eligible. Angles are eligible only if the
    phase's space group lets them refine:
        triclinic      al be ga
        monoclinic     the unique-axis angle (the non-90 angle in the file;
                       otherwise from the symbol, default be)
        rhombohedral   al be ga (rhombohedral-axes setting only)
        all others     none (fixed at 90 / 120)
    If a phase has no space_group (or it isn't recognised) an angle is
    treated as eligible only if its start value is not 90 or 120; this is
    reported.
  * Parameters sharing a name (constrained) always get the same limits:
    an existing min (or max) on any instance is copied to instances missing
    it; otherwise the default is used (for the cell, +/-1 % of the first
    instance's start value). Differing existing limits or start values are
    reported and not changed.
  * #ifdef / #ifndef / #if / #elseif / #else / #endif aware: instances of a
    name that sit in mutually exclusive branches (never active together) are
    treated as separate parameters, each with its own limits and checks.
  * Equation values (= ...;) and non-numeric values are skipped.
  * Symmetry check (phases with a recognised space group, not using the
    Cubic/Tetragonal/... macros): reports cell parameters that symmetry
    lets refine but are fixed (no name or @; not reported for hkl_Is phases,
    which are treated as PONKCS with fixed cells), angles that are refined but
    fixed by symmetry, and parameters that must be equal (e.g. a = b) but
    are not constrained by a shared name or an equation.
  * Comments (' ... and /* ... */) and "strings" are ignored.
  * Identifiers created with #define or macro are treated as flags, not
    parameter names.
  * The original file is not modified; output goes to <name>_limits.inp
    (or the path given with -o).

https://github.com/HenryXRD/TOPAS-inp-Limits  (MIT licence)

Usage
  py add_topas_limits_v2.py my_file.inp
  py add_topas_limits_v2.py my_file.inp --smart-po
  py add_topas_limits_v2.py my_file.inp -o new_file.inp
"""

__version__ = "2"

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

# ---- Settings (edit here if needed) -----------------------------------------
LIMITS = {
    "CS": ("32", "1000"),          # TOPAS shorthand for CS_L
    "CS_L": ("32", "1000"),
    "CS_G": ("32", "1000"),
    "Strain_L": ("0", "0.5"),
    "Strain_G": ("0", "0.5"),
    "PO": ("0.5", "1.5"),          # March-Dollase preferred orientation
    "PO_fraction": ("0", "1"),     # fraction of crystals in the first PO_Two_Directions direction
}
CELL_FRACTION = 0.01        # +/- 1 % of start value
CELL_MIN_DECIMALS = 4       # decimals used when writing cell limits
# -----------------------------------------------------------------------------

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
NUM_RE = re.compile(NUM)
# number, optionally followed by an error / limit suffix: `_0.01, `, _0.01, `_LIMIT_MIN_32
VALUE_RE = re.compile(NUM + r"(?:`?_[^\s,();]*|`[^\s,();]*)?")
IDENT_RE = re.compile(r"[!@]?[A-Za-z_]\w*|@")
# macros with size/strain parameters at fixed argument positions:
#   (name arg index, value arg index, limit type)
MULTI_MACROS = {
    # LVol_FWHM_CS_G_L(k, lvol, kf, lvolf, csgc, csgv, cslc, cslv)
    "LVol_FWHM_CS_G_L": [(4, 5, "CS_G"), (6, 7, "CS_L")],
    # e0_from_Strain(e0, sg, sgv, sl, slv)
    "e0_from_Strain": [(1, 2, "Strain_G"), (3, 4, "Strain_L")],
    # March-Dollase: PO(c, v, s, hkl) and Preferred_Orientation(c, v, hkl)
    "PO": [(0, 1, "PO")],
    "Preferred_Orientation": [(0, 1, "PO")],
    # PO_Two_Directions(c1, v1, ang1, hkl1, c2, v2, ang2, hkl2, w1c, w1v)
    "PO_Two_Directions": [(0, 1, "PO", "PO1"), (4, 5, "PO", "PO2"), (8, 9, "PO_fraction", "fraction")],
}
CS_RE = re.compile(r"(?<![\w])(LVol_FWHM_CS_G_L|e0_from_Strain|Preferred_Orientation|PO_Two_Directions|PO|"
                   r"CS_L|CS_G|CS|Strain_L|Strain_G)\(")
# "(" must follow the macro name directly
CELL_KW_RE = re.compile(r"(?<![\w.`])(a|b|c|al|be|ga)(?=\s)")
CELL_MACROS = {
    "Cubic": ("a",),
    "Tetragonal": ("a", "c"),
    "Hexagonal": ("a", "c"),
    "Trigonal": ("a", "c"),
    "Rhombohedral": ("a", "al"),
}
CELL_MACRO_RE = re.compile(r"(?<![\w])(" + "|".join(CELL_MACROS) + r")\(")
PHASE_RE = re.compile(r"(?<![\w])(str|hkl_Is|xo_Is|d_Is|xdd)(?![\w])|(?<![\w])(STR)\(")
ATTR_KW = {"min", "max", "del", "val_on_continue", "update", "stop_when"}
NOT_A_NAME = ATTR_KW | {"a", "b", "c", "al", "be", "ga"}

CUBIC = {s.replace("-", "") for s in """
P23 F23 I23 P213 I213 Pm-3 Pn-3 Fm-3 Fd-3 Im-3 Pa-3 Ia-3 P432 P4232 F432 F4132
I432 P4332 P4132 I4132 P-43m F-43m I-43m P-43n F-43c I-43d Pm-3m Pn-3n Pm-3n
Pn-3m Fm-3m Fm-3c Fd-3m Fd-3c Im-3m Ia-3d""".split()}
R_GROUP_NUMBERS = {146, 148, 155, 160, 161, 166, 167}


# ----------------------------------------------------------------- utilities
def code_mask(text):
    """True for characters that are real code (not comment / string)."""
    mask = [True] * len(text)
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "'":
            j = text.find("\n", i)
            j = n if j == -1 else j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j == -1 else j + 2
        elif c == '"':
            j = text.find('"', i + 1)
            j = n if j == -1 else j + 1
        else:
            i += 1
            continue
        for k in range(i, j):
            mask[k] = False
        i = j
    return mask


def line_start_pos(text, ln):
    pos = 0
    for _ in range(ln - 1):
        pos = text.index("\n", pos) + 1
    return pos


def line_no(text, pos):
    return text.count("\n", 0, pos) + 1


def num(v):
    m = NUM_RE.match(v.strip())
    return float(m.group(0)) if m else None


def plain(v):
    """'=1000;' -> '1000'."""
    v = v.strip()
    if v.startswith("="):
        v = v[1:]
    return v.rstrip(";").strip()


EQ_TOKEN = re.compile(r"\s*(?:((?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)|([-+*/^()]))")


def eval_simple(expr):
    """Evaluate a purely numeric TOPAS expression, e.g. '7.0567 - 7.0567 0.01'
    (a space between operands is multiplication, as in TOPAS). None if it
    contains anything else (parameter names, functions)."""
    toks, i = [], 0
    expr = expr.strip()
    while i < len(expr):
        m = EQ_TOKEN.match(expr, i)
        if not m or m.end() == i:
            if expr[i:].strip() == "":
                break
            return None
        toks.append(("n", float(m.group(1))) if m.group(1) else ("o", m.group(2)))
        i = m.end()
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else (None, None)

    def take():
        pos[0] += 1
        return toks[pos[0] - 1]

    def expr_():
        v = term()
        while peek() in (("o", "+"), ("o", "-")):
            op = take()[1]
            r = term()
            v = v + r if op == "+" else v - r
        return v

    def term():
        v = power()
        while True:
            t = peek()
            if t in (("o", "*"), ("o", "/")):
                op = take()[1]
                r = power()
                v = v * r if op == "*" else v / r
            elif t[0] == "n" or t == ("o", "("):
                v = v * power()                      # implicit multiplication
            else:
                return v

    def power():
        v = unary()
        if peek() == ("o", "^"):
            take()
            v = v ** power()
        return v

    def unary():
        t = peek()
        if t in (("o", "-"), ("o", "+")):
            take()
            v = unary()
            return -v if t[1] == "-" else v
        if t[0] == "n":
            return take()[1]
        if t == ("o", "("):
            take()
            v = expr_()
            if take() != ("o", ")"):
                raise ValueError
            return v
        raise ValueError

    try:
        v = expr_()
        return v if pos[0] == len(toks) else None
    except (ValueError, ZeroDivisionError, IndexError, OverflowError):
        return None


def plain_num(v):
    v = plain(v)
    return float(v) if NUM_RE.fullmatch(v) else eval_simple(v)


def norm(v):
    v = plain(v)
    m = NUM_RE.match(v)
    if m and (m.end() == len(v) or v[m.end()] in "`_"):
        return float(m.group(0))          # ignore error suffix
    e = eval_simple(v)
    if e is not None:
        return round(e, 9)
    return re.sub(r"\s+", "", v)


def skip_ws(code, i, end):
    while i < end and code[i].isspace():
        i += 1
    return i


def read_attr_value(code, i, end):
    """Value after min/max/...: '= ... ;' or a single token. Returns (str, end)."""
    i = skip_ws(code, i, end)
    if i < end and code[i] == "=":
        j = code.find(";", i, end)
        j = end if j == -1 else j + 1
        return code[i:j], j
    m = re.compile(r"[^\s()]+").match(code, i, end)
    return (m.group(0), m.end()) if m else ("", i)


def skip_parens(code, i, end):
    depth = 0
    while i < end:
        if code[i] == "(":
            depth += 1
        elif code[i] == ")":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    return end


def parse_value(code, i, end, flags, max_names=2):
    """Parse '[name tokens] number[`...] [attributes]' starting at i.

    Returns dict(name, value, insert_pos, min, max) or a string reason if
    not editable, or None if nothing parameter-like is there.
    """
    i = skip_ws(code, i, end)
    name = None
    status = "fixed"            # fixed / refined / unknown (flag macro)
    for _ in range(max_names + 1):
        if i < end and code[i] == "=":
            return "equation"
        vm = VALUE_RE.match(code, i, end)
        if vm:
            break
        im = IDENT_RE.match(code, i, end)
        if not im or im.group(0).lstrip("!@") in NOT_A_NAME:
            return None
        tok = im.group(0)
        bare = tok.lstrip("!@")
        if bare and bare not in flags:
            name = bare
        if bare in flags:
            status = "unknown"
        elif status != "unknown":
            status = "fixed" if tok.startswith("!") else "refined"
        i = skip_ws(code, im.end(), end)
    else:
        return None
    if not vm:
        return None

    res = dict(name=name, value=vm.group(0), insert_pos=vm.end(), min=None, max=None,
               status=status)
    # scan attributes following the value
    j = vm.end()
    while True:
        j = skip_ws(code, j, end)
        im = re.compile(r"[A-Za-z_]\w*").match(code, j, end)
        if not im:
            break
        tok = im.group(0)
        if tok in ATTR_KW:
            val, j = read_attr_value(code, im.end(), end)
            if tok in ("min", "max") and res[tok] is None:
                res[tok] = val
                res[tok + "_end"] = j
            continue
        k = skip_ws(code, im.end(), end)
        if k < end and code[k] == "(" and tok not in CELL_MACROS:   # e.g. a user macro on the value
            j = skip_parens(code, k, end)
            continue
        break
    return res


# ----------------------------------------------------------- space groups
def sg_elements(s):
    out, i = [], 0
    while i < len(s):
        el = ""
        if s[i] == "-":
            el, i = "-", i + 1
        if i < len(s) and s[i].isdigit():
            el += s[i]
            i += 1
            if (el[-1] in "2346" and i < len(s) and s[i].isdigit()
                    and "0" < s[i] < el[-1]):
                el += s[i]
                i += 1
            if i < len(s) and s[i] == "/":
                el += "/"
                i += 1
                if i < len(s) and (s[i].isalpha() or s[i].isdigit()):
                    el += s[i]
                    i += 1
        elif i < len(s) and s[i].isalpha() and not el:
            el, i = s[i], i + 1
        else:
            return None
        out.append(el)
    return out


def crystal_system(sg):
    """Returns (system, detail). detail = unique axis for monoclinic,
    'R'/'H'/None for trigonal."""
    raw = sg.strip().strip('"')
    setting = None
    if ":" in raw:
        raw, setting = raw.split(":", 1)
        setting = setting.strip().upper()
    s = raw.replace(" ", "").replace("_", "")
    if not s:
        return None, None
    if s.isdigit():
        n = int(s)
        if 1 <= n <= 2:
            return "triclinic", None
        if n <= 15:
            return "monoclinic", None
        if n <= 74:
            return "orthorhombic", None
        if n <= 142:
            return "tetragonal", None
        if n <= 167:
            return "trigonal", (setting[:1] if setting else None) if n in R_GROUP_NUMBERS else "H"
        if n <= 194:
            return "hexagonal", None
        if n <= 230:
            return "cubic", None
        return None, None
    lattice = s[0].upper()
    rest = s[1:]
    while rest and rest[-1].isupper():          # trailing setting letters, e.g. R-3mH
        setting = setting or rest[-1]
        rest = rest[:-1]
    if lattice not in "PABCIFR" or not rest:
        return None, None
    if (lattice + rest).replace("-", "") in CUBIC:
        return "cubic", None
    first = rest.lstrip("-")[:1]
    if first == "6":
        return "hexagonal", None
    if first == "3":
        if lattice == "R":
            return "trigonal", (setting[:1] if setting in ("R", "H") else None)
        return "trigonal", "H"
    if first == "4":
        return "tetragonal", None
    els = sg_elements(rest)
    if els is None:
        return None, None
    real = [k for k, e in enumerate(els) if e not in ("1", "-1")]
    if not real:
        return "triclinic", None
    if len(real) == 1:
        axis = {0: "al", 1: "be", 2: "ga"}[real[0]] if len(els) == 3 else None
        return "monoclinic", axis
    if len(els) == 3:
        return "orthorhombic", None
    return None, None


def is_special(v):
    return v is not None and (abs(v - 90) < 1e-6 or abs(v - 120) < 1e-6)


def eligible_angles(system, detail, angle_vals):
    """Which of al/be/ga may be refined. angle_vals: dict angle -> start value."""
    if system == "triclinic":
        return {"al", "be", "ga"}, None
    if system == "monoclinic":
        non90 = [k for k in ("al", "be", "ga")
                 if angle_vals.get(k) is not None and not is_special(angle_vals[k])]
        if len(non90) == 1:
            return set(non90), None
        return {detail or "be"}, None
    if system == "trigonal":
        if detail == "R":
            return {"al", "be", "ga"}, None
        if detail == "H":
            return set(), None
        ga, al = angle_vals.get("ga"), angle_vals.get("al")
        if ga is not None and abs(ga - 120) < 1e-6:
            return set(), None
        if al is not None and not is_special(al):
            return {"al", "be", "ga"}, None
        return set(), "rhombohedral setting not stated; assumed hexagonal axes"
    if system in ("orthorhombic", "tetragonal", "hexagonal", "cubic"):
        return set(), None
    # unknown: decide from values
    return {k for k, v in angle_vals.items() if v is not None and not is_special(v)}, \
        "space group not found or not recognised; angles chosen from start values"


# ----------------------------------------------------------- collection
def collect_flags(code):
    flags = set(re.findall(r"#define\s+([A-Za-z_]\w*)", code))
    flags |= set(re.findall(r"(?<![\w])macro\s+[&]?\s*([A-Za-z_]\w*)", code))
    return flags


def split_args(code, open_pos, close):
    """(start, end) of each top-level argument between the brackets."""
    bounds, depth, s0 = [], 0, open_pos + 1
    for k in range(open_pos, close + 1):
        ch = code[k]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                bounds.append((s0, k))
        elif ch == "," and depth == 1:
            bounds.append((s0, k))
            s0 = k + 1
    return bounds


def find_cs(text, code, flags, skipped):
    out = []
    for m in CS_RE.finditer(code):
        macro, open_pos = m.group(1), m.end() - 1
        close = skip_parens(code, open_pos, len(code)) - 1
        ln = line_no(text, m.start())
        if code[close] != ")":
            skipped.append((ln, f"{macro}: no closing bracket"))
            continue
        args = split_args(code, open_pos, close)
        call = text[m.start():close + 1].strip()

        if macro in MULTI_MACROS:
            pairs = [(e[0], e[1], e[2], e[3] if len(e) > 3 else e[2]) for e in MULTI_MACROS[macro]
                     if e[1] < len(args) and code[args[e[1]][0]:args[e[1]][1]].strip()]
            label_fmt = macro + ":{}" if len(MULTI_MACROS[macro]) > 1 else macro
        else:
            pairs = [(0 if len(args) > 1 else None, len(args) - 1, macro, macro)]
            label_fmt = "{}"

        for ni, vi, kind, lbl in pairs:
            val_start, val_end = args[vi]
            name = None
            if ni is not None:
                for tok in code[args[ni][0]:args[ni][1]].split():
                    bare = tok.lstrip("!@")
                    if bare and bare not in flags and IDENT_RE.fullmatch(tok):
                        name = bare
            label = label_fmt.format(lbl)
            p = parse_value(code, val_start, val_end, flags, max_names=0)
            if p is None:
                skipped.append((ln, f"{label}: non-numeric value: {call}"))
                continue
            if p == "equation":
                skipped.append((ln, f"{label}: equation value: {call}"))
                continue
            # min/max anywhere inside the value argument count as present
            arg = code[val_start:val_end]
            for key in ("min", "max"):
                if p[key] is None:
                    mm = re.search(r"(?<![\w])" + key + r"\s+(=[^;]*;|[^\s()]+)", arg)
                    if mm:
                        p[key] = mm.group(1)
                        p[key + "_end"] = val_start + mm.end()
            p.update(kind="cs", label=label, ln=ln, name=name, start=m.start(),
                     end=close + 1, text=call, defaults=LIMITS[kind])
            out.append(p)
    return out


def find_phases(text, code):
    marks = []
    for m in PHASE_RE.finditer(code):
        if m.group(2):                                   # STR(sg ...) macro
            if code[max(0, m.start() - 30):m.start()].split()[-1:] == ["macro"]:
                continue                                 # its definition, not a phase
            marks.append((m.start(), "STR", m.end() - 1))
        else:
            marks.append((m.start(), m.group(1), None))
    phases = []
    for k, (pos, kw, open_pos) in enumerate(marks):
        if kw not in ("str", "hkl_Is", "STR"):
            continue
        end = marks[k + 1][0] if k + 1 < len(marks) else len(code)
        sg = None
        if kw == "STR":
            close = skip_parens(code, open_pos, end) - 1
            args = split_args(code, open_pos, close) if code[close] == ")" else []
            if args:
                sg = text[args[0][0]:args[0][1]].strip().strip('"').strip() or None
            kw = "str"
        sm = None if sg else re.compile(r"(?<![\w])space_group(?![\w])").search(code, pos, end)
        if sm:
            i = skip_ws(text, sm.end(), end)
            if i < end and text[i] == '"':
                j = text.find('"', i + 1, end)
                sg = text[i + 1:j if j != -1 else end]
            else:
                tm = re.compile(r"\S+").match(code, i, end)   # comments blanked
                sg = tm.group(0) if tm else None
        pn_matches = [m for m in
                      re.compile(r"(?<![\w])phase_name(?![\w])").finditer(code, pos, end)
                      if not code[pos:m.start()].rstrip().endswith("(")]   # skip Get(phase_name)
        names = [line_no(text, m.start()) for m in pn_matches]
        pname = None
        if pn_matches:
            i = skip_ws(text, pn_matches[0].end(), end)
            if i < end and text[i] == '"':
                j = text.find('"', i + 1, end)
                pname = text[i + 1:j if j != -1 else end]
            else:
                tm = re.compile(r"\S+").match(code, i, end)
                pname = tm.group(0) if tm else None
        phases.append(dict(kw=kw, start=pos, end=end, ln=line_no(text, pos), sg=sg,
                           names=names, pname=pname,
                           pname_pos=pn_matches[0].start() if pn_matches else None))
    return phases


def find_cell(text, code, flags, phases, skipped, notes, first=None):
    out = []
    if first is None:
        first = phases[0]["start"] if phases else len(code)
    stray = [line_no(text, m.start()) for m in CELL_KW_RE.finditer(code, 0, first)
             if isinstance(parse_value(code, m.end(), first, flags), dict)
             and code[max(0, m.start() - 40):m.start()].split()[-1:] != ["occ"]]
    if stray:
        skipped.append((stray[0], "CHECK: cell parameters before the first recognised phase "
                        f"(lines {', '.join(map(str, stray[:10]))}{' ...' if len(stray) > 10 else ''})"
                        " - a phase may start with a macro the script does not recognise"))
    for ph in phases:
        if len(ph["names"]) > 1:
            skipped.append((ph["ln"], f"CHECK: phase at line {ph['ln']} has {len(ph['names'])} "
                            f"phase_name lines ({', '.join(map(str, ph['names']))}) - a phase "
                            "may start with a macro the script does not recognise (or these are "
                            "#ifdef alternatives); its space group was applied to all of them"))
        found, equations = [], set()
        # keyword form: a @ 5.43 ...
        for m in CELL_KW_RE.finditer(code, ph["start"], ph["end"]):
            before = code[max(ph["start"], m.start() - 40):m.start()].split()
            if before and before[-1] in ("prm", "local", "existing_prm", "occ"):
                continue        # e.g. "occ al 0.5": element symbol, not the cell angle
            p = parse_value(code, m.end(), ph["end"], flags)
            if p is None:
                continue
            ln = line_no(text, m.start())
            if p == "equation":
                equations.add(m.group(1))   # constrained by equation, no limits needed
                continue
            p.update(param=m.group(1), label=m.group(1), ln=ln)
            found.append(p)
        # macro form: Cubic(@ 5.43), Hexagonal(@ 4.9, @ 5.4) ...
        for m in CELL_MACRO_RE.finditer(code, ph["start"], ph["end"]):
            open_pos = m.end() - 1
            close = skip_parens(code, open_pos, ph["end"]) - 1
            bounds, depth, s0 = [], 0, open_pos + 1
            for k in range(open_pos, close + 1):
                ch = code[k]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        bounds.append((s0, k))
                elif ch == "," and depth == 1:
                    bounds.append((s0, k))
                    s0 = k + 1
            for param, (s, e) in zip(CELL_MACROS[m.group(1)], bounds):
                p = parse_value(code, s, e, flags)
                if p is None or p == "equation":
                    continue
                p.update(param=param, label=f"{m.group(1)}:{param}", ln=line_no(text, s),
                         from_macro=True)
                found.append(p)

        angle_vals = {}
        for p in found:
            if p["param"] in ("al", "be", "ga"):
                angle_vals.setdefault(p["param"], num(p["value"].split("`")[0]))
        system, detail = crystal_system(ph["sg"]) if ph["sg"] else (None, None)
        macros_used = {p["label"].split(":")[0] for p in found if p.get("from_macro")}
        if system == "trigonal" and detail is None:
            if macros_used & {"Trigonal", "Hexagonal"}:
                detail = "H"
            elif "Rhombohedral" in macros_used:
                detail = "R"
        angles, warn = eligible_angles(system, detail, angle_vals)
        ph.update(system=system, detail=detail, found=found)
        sg_txt = ph["sg"] if ph["sg"] else "none"
        notes.append(f"{ph['kw']} at line {ph['ln']}: space group {sg_txt} -> "
                     f"{system or 'unknown'}; angles refined: "
                     f"{', '.join(sorted(angles)) if angles else 'none'}"
                     + (f"  ({warn})" if warn else ""))
        if system and not any(p.get("from_macro") for p in found):
            check_symmetry(ph, system, detail, angles, found, equations, skipped)
        for p in found:
            if p["param"] in ("al", "be", "ga") and p["param"] not in angles \
                    and not p.get("from_macro"):
                v = num(p["value"].split("`")[0])
                if not is_special(v):
                    skipped.append((p["ln"], f"{p['param']} = {p['value']} but fixed by "
                                    f"{system} symmetry - check space group; no limits added"))
                continue
            p.update(kind="cell", start=p["insert_pos"], end=p["insert_pos"],
                     text=f"{p['label']} {p['value']}")
            out.append(p)
    return out


def check_symmetry(ph, system, detail, angles, found, equations, out):
    """Check cell parameters are refined/constrained as the space group requires."""
    rhombo = system == "trigonal" and "al" in angles
    if system == "cubic" or rhombo:
        equal = [("a", "b", "c")] + ([("al", "be", "ga")] if rhombo else [])
    elif system in ("tetragonal", "hexagonal", "trigonal"):
        equal = [("a", "b")]
    else:
        equal = []
    dependent = {q for grp in equal for q in grp[1:]}
    independent = [q for q in ("a", "b", "c") if q not in dependent] + \
                  [q for q in ("al", "be", "ga") if q in angles and q not in dependent]
    by_param = defaultdict(list)
    for p in found:
        by_param[p["param"]].append(p)
    tag = f"[{system} {ph['sg']}]"

    ponkcs = ph["kw"] == "hkl_Is"     # hkl_Is phases are PONKCS: cells fixed by design
    for q in ([] if ponkcs else independent):
        for p in by_param.get(q, []):
            if p["status"] == "fixed":
                out.append((p["ln"], f"CHECK {tag}: {q} {p['value']} is fixed (no name or @) "
                            f"but symmetry allows it to refine"))
        if not by_param.get(q) and q not in equations:
            out.append((ph["ln"], f"CHECK {tag}: {q} not found in phase"))

    for q in ("al", "be", "ga"):
        if q in angles:
            continue
        for p in by_param.get(q, []):
            if p["status"] == "refined":
                out.append((p["ln"], f"CHECK {tag}: {q} is refined but fixed by symmetry"))

    for grp in equal:
        lead = by_param.get(grp[0], [])
        if not lead:
            continue
        for q in grp[1:]:
            if q in equations and not by_param.get(q):
                continue                     # e.g. b = Get(a);
            for p in by_param.get(q, []):
                for l in lead:
                    if p["status"] == "unknown" or l["status"] == "unknown":
                        continue
                    if l["status"] == "refined" and (p["name"] is None or p["name"] != l["name"]):
                        out.append((p["ln"], f"CHECK {tag}: {q} must equal {grp[0]} but is "
                                    f"not constrained to it (use the same name as {grp[0]}, "
                                    f"or {q} = Get({grp[0]});)"))
                    elif norm(p["value"]) != norm(l["value"]):
                        out.append((p["ln"], f"CHECK {tag}: {q} {p['value']} differs from "
                                    f"{grp[0]} {l['value']}"))


def cell_defaults(value):
    v = num(value.split("`")[0])
    mant = NUM_RE.match(value).group(0)
    dec = len(mant.split(".")[1].split("e")[0].split("E")[0]) if "." in mant else 0
    d = max(dec, CELL_MIN_DECIMALS)
    lo, hi = sorted((v * (1 - CELL_FRACTION), v * (1 + CELL_FRACTION)))
    return f"{lo:.{d}f}", f"{hi:.{d}f}"


# ----------------------------------------------------------- #ifdef awareness
PP_RE = re.compile(r"#\s*(ifdef|ifndef|if|elseif|elif|else|endif)(?![\w])[ \t]*([A-Za-z_]\w*)?")


def preprocessor_contexts(code):
    """List of (pos, context) where context is a tuple of
    (block_id, branch_index, flag, positive) for every enclosing #if block."""
    states, stack, next_id = [(0, ())], [], 0
    for m in PP_RE.finditer(code):
        kw, flag = m.group(1), m.group(2)
        if kw in ("ifdef", "ifndef", "if"):
            f = flag if kw != "if" else None
            stack.append([next_id, 0, f, kw != "ifndef"])
            next_id += 1
        elif kw in ("elseif", "elif") and stack:
            stack[-1][1] += 1
            stack[-1][2] = None
        elif kw == "else" and stack:
            top = stack[-1]
            top[1] += 1
            # #else of a plain #ifdef/#ifndef is the opposite of that flag
            top[3] = (not top[3]) if (top[2] is not None and top[1] == 1) else top[3]
            if top[1] > 1:
                top[2] = None
        elif kw == "endif" and stack:
            stack.pop()
        states.append((m.end(), tuple(tuple(e) for e in stack)))
    return states


def context_at(states, pos):
    lo, hi = 0, len(states) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if states[mid][0] <= pos:
            lo = mid
        else:
            hi = mid - 1
    return states[lo][1]


def exclusive(c1, c2):
    """True if two #ifdef contexts can never be active together."""
    for b1, i1, f1, p1 in c1:
        for b2, i2, f2, p2 in c2:
            if b1 == b2 and i1 != i2:
                return True
            if f1 is not None and f1 == f2 and p1 != p2:
                return True
    return False


def split_coactive(group):
    """Split instances of one name into sets that can be active together."""
    parent = list(range(len(group)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(group)):
        for j in range(i + 1, len(group)):
            if not exclusive(group[i]["ctx"], group[j]["ctx"]):
                parent[find(i)] = find(j)
    sets = defaultdict(list)
    for i, p in enumerate(group):
        sets[find(i)].append(p)
    return sorted(sets.values(), key=lambda g: g[0]["ln"])



# ----------------------------------------------------------- smart PO
PO_DB_NAME = "po_directions.csv"
PO_EXISTING_RE = re.compile(r"(?<![\w])(?:PO\w*|Preferred_Orientation\w*)\s*\(|"
                            r"(?<![\w])spherical_harmonics_hkl(?![\w])")
# clays: 0 0 1 with the layer-silicate sheet check (a ~5.0-5.5, b ~8.7-9.4 A)
CLAY_NAMES = """kaolinite dickite nacrite halloysite illite smectite montmorillonite
beidellite nontronite saponite hectorite stevensite sauconite vermiculite chlorite
clinochlore chamosite corrensite rectorite talc pyrophyllite clay""".split()
CHECKS = ("", "hex", "a_longest", "b_longest", "c_longest",
          "a_shortest", "b_shortest", "c_shortest", "sheet_ab")


def load_po_db(path, warnings):
    """[(names, (h, k, l), check, source)] from the CSV, then the clay list."""
    import csv
    db = []
    if path and Path(path).is_file():
        with open(path, newline="", encoding="utf-8-sig") as f:
            for n, row in enumerate(csv.DictReader(f), start=2):
                row = {(k or "").strip().lower(): (v or "").strip() for k, v in row.items()}
                if row.get("phase", "").startswith("#"):
                    continue                            # section heading / comment row
                names = [row.get("phase", "")] + [a for a in row.get("aliases", "").split(";")]
                names = [x.strip() for x in names if x.strip()]
                raw = [row.get(k, "").strip().lower() for k in ("h", "k", "l")]
                try:
                    if all(r in ("", "none", "-") for r in raw):
                        hkl = None                      # recognised, but no PO direction
                    else:
                        hkl = tuple(int(r) for r in raw)
                except ValueError:
                    warnings.append(f"{Path(path).name} line {n}: h k l not integers - row ignored")
                    continue
                check = row.get("check", "").lower()
                if check not in CHECKS:
                    warnings.append(f"{Path(path).name} line {n}: unknown check '{check}' - "
                                    f"row ignored (use one of: {', '.join(c for c in CHECKS if c)})")
                    continue
                unsure = row.get("confidence", "").lower() in ("check", "uncertain", "verify")
                if names:
                    db.append((names, hkl, check, Path(path).name, unsure))
    db.append((CLAY_NAMES, (0, 0, 1), "sheet_ab", "built-in clay list", False))
    return db


# chemical prefixes that may be joined to a mineral name, e.g. magnesiohornblende
NAME_PREFIX = (r"(?:(?:magnesio|ferro|ferri|fluoro|chloro|oxo|alumino|mangano|mangani|"
               r"titano|chromio|potassic|sodic|parvo|clino|ortho)-?)*")


def match_po(pname, db, heading=False, strict=False):
    """Longest whole-word name match (chemical prefixes and plurals allowed) -> entry or None.
    heading=True: only group-level matches count (a plural, 'X group', or a generic group
    name such as Amphibole or Feldspar), so a single mineral named in a comment is not
    taken as the heading for the phases after it."""
    low = pname.lower()
    best, hkls = None, set()
    for names, hkl, check, src, unsure in db:
        for nm in names:
            pat = (r"(?<![a-z])" + NAME_PREFIX + "(" + re.escape(nm.lower()) + r")(e?s)?(?![a-z])"
                   r"(\s+(?:super|sub)?(?:group|series|family))?")
            for m in re.finditer(pat, low):
                hkls.add(hkl)
                if heading and not (m.group(2) or m.group(3) or nm.lower() in GENERIC_GROUPS):
                    continue
                if best is None or len(nm) > len(best[0]):
                    best = (nm, hkl, check, src, unsure)
    if (strict or heading) and best and len(hkls) > 1:
        return None          # a comment naming minerals with different directions: ambiguous
    return best


# group names that may stand alone as a section heading
GENERIC_GROUPS = {"amphibole", "clinoamphibole", "orthoamphibole", "pyroxene", "clinopyroxene",
                  "orthopyroxene", "feldspar", "k-feldspar", "plagioclase", "mica", "clay",
                  "olivine", "garnet", "spinel", "chlorite", "smectite"}


def by_system(nm, hkl, ph):
    """Generic 'amphibole' / 'pyroxene': 2 1 0 for an orthorhombic cell."""
    if nm.lower() in ("amphibole", "pyroxene") and ph.get("system") == "orthorhombic":
        return (2, 1, 0)
    return hkl


# ----------------------------------------------------------- comments as hints
HEADING_FAMILIES = [
    # (family, words) - singular or plural, matched as whole words
    ("layer", r"phyllosilicate|sheet silicate|layer silicate|clay mineral|kaolin"),
    ("chain", r"inosilicate|chain silicate|pyroxenoid"),
    ("silicate", r"silicate|nesosilicate|orthosilicate|sorosilicate|cyclosilicate|"
                 r"ring silicate|tectosilicate|framework silicate|zeolite"),
    ("nonsilicate", r"native element|element|sulfide|sulphide|sulfosalt|sulphosalt|halide|"
                    r"oxide|hydroxide|oxyhydroxide|carbonate|nitrate|borate|sulfate|sulphate|"
                    r"chromate|molybdate|tungstate|phosphate|arsenate|vanadate|organic"),
    ("other", r"clinker|cement|hydrate|hydration product|standard|internal standard|spike|"
              r"amorphous|glass|slag|fly ash|ponkcs|unknown"),
]
HEADING_RANK = {"entry": 5, "layer": 4, "chain": 3, "silicate": 2, "nonsilicate": 1, "other": 1}
# Nickel-Strunz codes, e.g. 9.DE.10 or 09.DA
STRUNZ_RE = re.compile(r"(?<![\w.])0?(10|[1-9])\.([A-J])([A-Z])?(?:\.\d{2,3}[a-z]?)?(?![\w])")
STRUNZ_FAMILIES = {"9.DA": "pyroxene", "9.DE": "amphibole", "9.D": "chain", "9.E": "layer"}
CODE_WORDS = set("""str STR hkl_Is xo_Is d_Is phase_name space_group a b c al be ga prm local
scale site occ beq x y z MVW Cubic Tetragonal Hexagonal Trigonal Rhombohedral CS CS_L CS_G
Strain_L Strain_G PO PO_Two_Directions Preferred_Orientation LVol_FWHM_CS_G_L e0_from_Strain
macro Out r_bragg num_posns rigid load lor_fwhm gauss_fwhm spherical_harmonics_hkl
TCHZ_Peak_Type Zero_Error cell_mass cell_volume weight_percent""".split())


def strunz_family(text):
    best = None
    for m in STRUNZ_RE.finditer(text):
        cls, sub, fam = m.group(1), m.group(2), m.group(3) or ""
        if cls != "9":
            f = "nonsilicate"
        else:
            f = STRUNZ_FAMILIES.get(f"9.{sub}{fam}") or STRUNZ_FAMILIES.get(f"9.{sub}") or "silicate"
        if best is None or HEADING_RANK.get(f, 3) > HEADING_RANK.get(best, 3):
            best = f
    return best


def classify_heading(text, db):
    """Section heading -> dict(kind, family, entry) or None."""
    hit = match_po(text, db, heading=True)
    if hit:
        nm = hit[0].lower()
        fam = ("nopo" if hit[1] is None else "layer" if hit[2] == "sheet_ab" else
               "chain" if re.search(r"amphibole|pyroxene", nm) or
               hit[1] in ((1, 1, 0), (2, 1, 0)) and hit[2] == "c_shortest" else "mineral")
        return dict(kind="entry", family=fam, entry=hit)
    sf = strunz_family(text)
    if sf in ("pyroxene", "amphibole"):
        e = match_po(sf, db)
        if e:
            return dict(kind="entry", family="chain", entry=e)
        sf = "chain"
    low = text.lower()
    best = None
    for fam, words in HEADING_FAMILIES:
        if re.search(r"(?<![a-z])(?:" + words + r")s?(?![a-z])", low):
            if best is None or HEADING_RANK[fam] > HEADING_RANK[best]:
                best = fam
    if sf and (best is None or HEADING_RANK.get(sf, 0) > HEADING_RANK[best]):
        best = sf
    return dict(kind="family", family=best, entry=None) if best else None


def comment_text(line):
    t = line.strip()
    t = re.sub(r"^(?:/\*+|\*+/|'+|\*+)", "", t)
    t = re.sub(r"\*+/\s*$", "", t)
    return t.strip(" \t'*=-_#~:.|")


def code_like(t):
    if not t or "smart PO" in t or "@" in t or re.search(r"(?<![\w])(?:phase_name|space_group)(?![\w])", t):
        return True
    first = re.match(r"[A-Za-z_]\w*", t)
    return bool(first and first.group(0) in CODE_WORDS and
                re.match(r"[A-Za-z_]\w*\s*[(@!\d\-\"]|[A-Za-z_]\w*\s+[A-Za-z_!@]\w*\s*[@\d(]", t))


def line_comment(text, mask, a, b):
    """Trailing comment on text[a:b] (outside strings), or ''."""
    in_s = False
    for i in range(a, b):
        ch = text[i]
        if ch == '"':
            in_s = not in_s
        elif not in_s and (ch == "'" or text.startswith("/*", i)):
            return comment_text(text[i:b])
    return ""


def comment_hints(text, code, mask, str_phases, db, first_phase):
    """For each str phase: attached comments (directly above it, or trailing on its str /
    phase_name line) and the nearest section heading above it."""
    import bisect
    starts = [0] + [m.end() for m in re.finditer("\n", text)]
    ends = [s - 1 for s in starts[1:]] + [len(text)]
    is_com = []
    for a, b in zip(starts, ends):
        is_com.append(bool(text[a:b].strip()) and not code[a:b].strip())
    first_ln = text.count("\n", 0, first_phase)            # 0-based line index
    # comments before the first phase only count if they sit directly above it
    lo = first_ln
    while lo > 0 and (is_com[lo - 1] or not text[starts[lo - 1]:ends[lo - 1]].strip()):
        lo -= 1
    headings = []                                          # (line index, text, result)
    for i, (a, b) in enumerate(zip(starts, ends)):
        if is_com[i] and (i >= lo):
            t = comment_text(text[a:b])
            if t and not code_like(t):
                r = classify_heading(t, db)
                if r:
                    headings.append((i, t, r))
    h_lines = [h[0] for h in headings]
    for ph in str_phases:
        li = text.count("\n", 0, ph["start"])
        attached = []
        j = li - 1
        while j >= 0 and (is_com[j] or not text[starts[j]:ends[j]].strip()):
            if is_com[j]:
                t = comment_text(text[starts[j]:ends[j]])
                if t and not code_like(t):
                    attached.append((j + 1, t))
            j -= 1
        for pos in (ph["start"], ph.get("pname_pos")):
            if pos is not None:
                k = text.count("\n", 0, pos)
                t = line_comment(text, mask, starts[k], ends[k])
                if t and not code_like(t):
                    attached.append((k + 1, t))
        ph["attached"] = attached
        k = bisect.bisect_left(h_lines, li) - 1
        ph["heading"] = None if k < 0 else dict(headings[k][2], ln=headings[k][0] + 1,
                                                 text=headings[k][1])


def phase_cell(ph):
    cell = {}
    for p in ph.get("found", []):
        v = num(p["value"].split("`")[0].split("_")[0])
        if v is not None:
            cell.setdefault(p["param"], v)
    a = cell.get("a")
    sysm = ph.get("system")
    if a is not None:
        if sysm in ("cubic", "tetragonal", "hexagonal", "trigonal"):
            cell.setdefault("b", a)
        if sysm == "cubic" or (sysm == "trigonal" and ph.get("detail") == "R"):
            cell.setdefault("c", a)
    if sysm in ("hexagonal",) or (sysm == "trigonal" and ph.get("detail") != "R"):
        cell.setdefault("ga", 120.0)
    return cell


def sheet_signature(cell):
    a, b, c = cell.get("a"), cell.get("b"), cell.get("c")
    return (a is not None and b is not None and c is not None
            and 5.0 <= a <= 5.5 and 8.7 <= b <= 9.4 and abs(b / a - 3 ** 0.5) <= 0.05
            and c >= 7.0)


CHAIN_CELLS = [
    # (label, hkl, a range, b range, c range, beta range or None)
    ("clinoamphibole", (1, 1, 0), (9.5, 10.1), (17.5, 18.6), (5.2, 5.45), (101, 108)),
    ("clinopyroxene", (1, 1, 0), (9.4, 10.0), (8.5, 9.3), (5.1, 5.5), (105, 110)),
    ("orthopyroxene", (2, 1, 0), (18.1, 18.6), (8.6, 9.2), (5.1, 5.35), None),
    ("orthoamphibole", (2, 1, 0), (18.4, 18.8), (17.7, 18.2), (5.2, 5.35), None),
]


def chain_signature(cell):
    a, b, c, be = cell.get("a"), cell.get("b"), cell.get("c"), cell.get("be", 90.0)
    if None in (a, b, c):
        return None
    for label, hkl, ra, rb, rc, rbe in CHAIN_CELLS:
        if ra[0] <= a <= ra[1] and rb[0] <= b <= rb[1] and rc[0] <= c <= rc[1] and \
                (rbe is None and abs(be - 90) < 0.5 or rbe and rbe[0] <= be <= rbe[1]):
            return label, hkl
    return None


def check_cell(check, cell, ph):
    """None if OK, else the reason it failed."""
    a, b, c = cell.get("a"), cell.get("b"), cell.get("c")
    if check == "":
        return None
    if check == "hex":
        if ph.get("system") in ("hexagonal", "trigonal") and ph.get("detail") != "R" \
                and abs(cell.get("ga", 0) - 120) < 0.5:
            return None
        return "not in a hexagonal-axes setting"
    if check == "sheet_ab":
        return None if sheet_signature(cell) else \
            "cell does not have the layer-silicate sheet in a, b (a ~5.0-5.5, b ~8.7-9.4 A)"
    axis = check[0]
    vals = {"a": a, "b": b, "c": c}
    if any(v is None for v in vals.values()):
        return "cell lengths not all found"
    others = [v for k, v in vals.items() if k != axis]
    if check.endswith("_shortest"):
        return None if vals[axis] < min(others) else f"{axis} is not the shortest axis"
    return None if vals[axis] > max(others) else f"{axis} is not the longest axis"


def smart_po(text, code, str_phases, db, edits, skipped, smart_report, no_po=None,
             mask=None, use_comments=True):
    if use_comments and str_phases:
        comment_hints(text, code, mask, str_phases, db, min(p["start"] for p in str_phases))
    idents = set(re.findall(r"[A-Za-z_]\w*", code))
    states = preprocessor_contexts(code)
    nl = "\r\n" if "\r\n" in text else "\n"
    for ph in str_phases:
        if PO_EXISTING_RE.search(code, ph["start"], ph["end"]):
            continue                                   # already has a PO model
        pname = ph.get("pname") or ""
        label_p = pname or "(no phase_name)"
        cell = phase_cell(ph)
        head = ph.get("heading")
        fam = head["family"] if head else None
        h_txt = f"section comment '{head['text'][:40]}' (line {head['ln']})" if head else ""
        hit = match_po(pname, db) if pname else None
        via = ""
        if not hit:                                    # comments directly above / beside it
            for cln, t in ph.get("attached", []):
                hit = match_po(t, db, strict=True)
                if hit:
                    via = f" from comment on line {cln}"
                    break
        if hit and hit[1] is None:
            if no_po is not None:
                no_po.append((ph["ln"], label_p, f"'{hit[0]}'{via}"))
            continue
        unsure = bool(via)
        if hit:
            nm, hkl, check, src, u = hit
            hkl, unsure = by_system(nm, hkl, ph), unsure or u
            why = check_cell(check, cell, ph)
            if why:
                skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}' (matched "
                                f"'{nm}'{via}, {' '.join(map(str, hkl))}): {why} - check the setting"))
                continue
            reason = f"'{nm}' ({src}){via}"
            if via and check == "sheet_ab":
                unsure = u                             # the cell confirms the layer silicate
            if head and head["kind"] == "entry" and head["entry"][1] is not None and \
                    by_system(head["entry"][0], head["entry"][1], ph) != hkl:
                reason += (f"; note: {h_txt} suggests "
                           f"{' '.join(map(str, by_system(head['entry'][0], head['entry'][1], ph)))}")
        elif fam == "nopo":
            if no_po is not None:
                no_po.append((ph["ln"], label_p, f"'{head['entry'][0]}' group, {h_txt}"))
            continue
        elif fam in (None, "silicate", "layer") and sheet_signature(cell):
            hkl, reason = (0, 0, 1), "clay inferred from cell (a, b layer-silicate sheet)"
            if fam == "layer":
                reason += f" and {h_txt}"
        elif head and head["kind"] == "entry" and (
                fam != "chain" and chain_signature(cell) or fam != "layer" and sheet_signature(cell)):
            sig = "a chain silicate" if chain_signature(cell) else "a layer silicate"
            skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}': {h_txt} "
                            f"suggests '{head['entry'][0]}' but the cell looks like {sig} - "
                            "check the phase and the heading"))
            continue
        elif head and head["kind"] == "entry":
            nm, hkl, check, src, u = head["entry"]
            hkl = by_system(nm, hkl, ph)
            why = check_cell(check, cell, ph)
            if why == "cell lengths not all found" or (why and fam != "chain"):
                skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}' ({h_txt} "
                                f"suggests {' '.join(map(str, hkl))}): {why}"))
                continue
            if fam == "chain":
                sig = chain_signature(cell)
                if sig:
                    hkl = sig[1]
                elif ph.get("system") == "monoclinic":
                    hkl = (1, 1, 0)
                elif ph.get("system") == "orthorhombic":
                    hkl = (2, 1, 0)
                else:
                    skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}' "
                                    f"({h_txt}): cell is not monoclinic or orthorhombic"))
                    continue
            reason, unsure = f"'{nm}' group from {h_txt}", True
        elif fam == "layer":
            skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}' ({h_txt} "
                            "suggests 0 0 1): cell does not have the layer-silicate sheet in a, b"))
            continue
        elif fam in (None, "silicate", "chain") and chain_signature(cell):
            label, hkl = chain_signature(cell)
            reason, unsure = f"{label} inferred from cell", True
            if fam == "chain":
                reason += f" and {h_txt}"
        else:
            sig = chain_signature(cell)
            sig = sig[0] if sig else "layer silicate" if sheet_signature(cell) else None
            if sig:
                skipped.append((ph["ln"], f"CHECK: smart PO not added to '{label_p}': the cell "
                                f"looks like a {sig} but it is under {h_txt} - check the phase "
                                "and the heading"))
            continue
        base = "po_" + (re.sub(r"[^A-Za-z0-9]+", "_", pname).strip("_") or f"line{ph['ln']}")
        if not base[3:4].isalpha():
            base = "po_p" + base[3:]
        name, k = base, 2
        while name in idents:
            name, k = f"{base}_{k}", k + 1
        idents.add(name)
        # put it after phase_name if that is in the same #ifdef branch as the phase start,
        # otherwise straight after the str / STR(...) line
        anchor = ph["start"]
        if ph.get("pname_pos") is not None and \
                context_at(states, ph["pname_pos"]) == context_at(states, ph["start"]):
            anchor = ph["pname_pos"]
        eol = text.find("\n", anchor)
        eol = len(text) if eol == -1 else eol + 1
        line_start = text.rfind("\n", 0, anchor) + 1
        indent = re.match(r"[ \t]*", text[line_start:]).group(0)
        if anchor == ph["start"]:
            indent += "   "
        hkl_s = " ".join(map(str, hkl))
        if unsure:
            line = (f"{indent}PO(!{name}, 1 min 0.5 max 1.5, , {hkl_s}) ' smart PO: UNCERTAIN "
                    f"direction (sources disagree or it depends on the setting) from {reason} - "
                    f"CHECK and change {hkl_s} if needed before refining; not refined")
            reason += " - UNCERTAIN, check direction"
        else:
            line = (f"{indent}PO(!{name}, 1 min 0.5 max 1.5, , {hkl_s}) ' smart PO: likely "
                    f"direction from {reason}; not refined - change !{name} to {name} to refine")
        ins = line + nl
        if eol == len(text) and not text.endswith("\n"):
            ins = nl + line
        edits.append((eol, ins))
        smart_report.append((line_no(text, anchor) + 1, pname or "(no phase_name)", hkl_s, reason))


# ----------------------------------------------------------- beq
ELEMENTS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu
Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm
Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np
Pu Am Cm Bk Cf Es Fm Md No Lr D""".split())
BEQ_GUI = "1"                  # value TOPAS GUI assumes when a str / CIF has no beq
# smart beq: H 2; O and similar anion-forming non-metals 1; metals, metalloids and
# positively charged non-metals (S+6 in sulfate, C+4, P+5 ...) 0.5
BEQ_H = {"H", "D"}
BEQ_ANION = {"C", "N", "O", "F", "S", "Cl", "Se", "Br", "I", "He", "Ne", "Ar", "Kr", "Xe", "Rn"}
BEQ_VALUES = {"H": "2", "anion": "1", "other": "0.5"}
SITE_RE = re.compile(r"(?<![\w])site(?![\w])")
SITE_KW = {"x", "y", "z", "occ", "beq", "num_posns", "rand_xyz", "inter", "adps", "ADPs",
           "u11", "u22", "u33", "u12", "u13", "u23", "scale_occ", "in_cartesian"}
BEQ_PRESENT_RE = re.compile(r"(?<![\w])(?:beq|Beq\w*|BEQ\w*|u11|adps|ADPs)(?![\w])")


def element_of(tok):
    """'Si+4' -> ('Si', '+'), 'O-2' -> ('O', '-'), 'al' -> ('Al', ''); None if unknown."""
    m = re.match(r"([A-Za-z])([A-Za-z]?)", tok)
    if not m:
        return None
    two = (m.group(1) + m.group(2)).capitalize()
    el = two if m.group(2) and two in ELEMENTS else m.group(1).upper()
    if el not in ELEMENTS:
        return None
    sign = re.search(r"[+-]", tok)
    return el, sign.group(0) if sign else ""


def smart_beq_value(el, sign):
    if el in BEQ_H:
        return BEQ_VALUES["H"]
    if el in BEQ_ANION and sign != "+":
        return BEQ_VALUES["anion"]
    return BEQ_VALUES["other"]


def consume_value(code, i, end):
    """Skip '[name tokens] value [min .. max .. ...]' or '= ... ;'. Returns (end, number or None)."""
    i = skip_ws(code, i, end)
    if i < end and code[i] == "=":
        j = code.find(";", i, end)
        return (end if j == -1 else j + 1), None
    number = None
    for _ in range(4):
        vm = VALUE_RE.match(code, i, end)
        if vm:
            number = num(vm.group(0).split("`")[0])
            i = vm.end()
            break
        im = IDENT_RE.match(code, i, end)
        if not im or im.group(0).lstrip("!@") in SITE_KW:
            return i, None
        i = skip_ws(code, im.end(), end)
        if i < end and code[i] == "=":                 # name = equation;
            j = code.find(";", i, end)
            return (end if j == -1 else j + 1), None
    while True:
        j = skip_ws(code, i, end)
        im = re.compile(r"[A-Za-z_]\w*").match(code, j, end)
        if not im or im.group(0) not in ATTR_KW:
            return i, number
        _, i = read_attr_value(code, im.end(), end)


def insert_beq(text, code, str_phases, mode, edits, skipped, beq_report):
    for ph in str_phases:
        sites = [m for m in SITE_RE.finditer(code, ph["start"], ph["end"])]
        counts, n_sites = {}, 0
        for k, m in enumerate(sites):
            s_end = sites[k + 1].start() if k + 1 < len(sites) else ph["end"]
            if BEQ_PRESENT_RE.search(code, m.end(), s_end):
                continue                                # beq (or ADPs) already there
            i = skip_ws(code, m.end(), s_end)
            nm = re.compile(r"\S+").match(code, i, s_end)
            if not nm:
                continue
            site_name, i, last = nm.group(0), nm.end(), nm.end()
            species = []
            while True:
                j = skip_ws(code, i, s_end)
                im = re.compile(r"[A-Za-z_]\w*").match(code, j, s_end)
                if not im or im.group(0) not in SITE_KW:
                    break
                i = im.end()
                if im.group(0) == "occ":
                    j = skip_ws(code, i, s_end)
                    em = re.compile(r"[^\s()]+").match(code, j, s_end)
                    if not em:
                        break
                    i, occ = consume_value(code, em.end(), s_end)
                    species.append((em.group(0), occ))
                elif im.group(0) not in ("rand_xyz", "in_cartesian"):
                    i, _ = consume_value(code, i, s_end)
                last = i
            if not species:
                skipped.append((line_no(text, m.start()), f"CHECK: site {site_name} in "
                                f"'{ph.get('pname') or 'str'}': no occ found - beq not inserted"))
                continue
            if mode == "smart":
                parsed = [(element_of(t), o) for t, o in species]
                if any(e is None for e, _ in parsed):
                    bad = [t for (t, _), (e, _) in zip(species, parsed) if e is None]
                    skipped.append((line_no(text, m.start()), f"CHECK: site {site_name} in "
                                    f"'{ph.get('pname') or 'str'}': element '{bad[0]}' not "
                                    "recognised - beq not inserted"))
                    continue
                # the species with the largest numeric occupancy decides (first if unknown)
                best = max(range(len(parsed)),
                           key=lambda q: (parsed[q][1] if parsed[q][1] is not None else -1, -q))
                val = smart_beq_value(*parsed[best][0])
            else:
                val = BEQ_GUI
            edits.append((last, f" beq {val}"))
            counts[val] = counts.get(val, 0) + 1
            n_sites += 1
        if n_sites:
            detail = ", ".join(f"beq {v} x{c}" for v, c in
                               sorted(counts.items(), key=lambda x: -float(x[0])))
            beq_report.append((ph["ln"], ph.get("pname") or "(no phase_name)", n_sites, detail))


# ----------------------------------------------------------- main logic
def pick_limit(group, key, default, warnings, label):
    existing = {}
    for c in group:
        if c[key] is not None:
            existing.setdefault(norm(c[key]), c[key])
    if not existing:
        return default
    if len(existing) == 1:
        return next(iter(existing.values()))
    lines = ", ".join(f"{c['ln']} ({key} {c[key]})" for c in group if c[key] is not None)
    warnings.append(f"{label}: different existing {key} values on lines {lines} "
                    f"- not inserting {key}; please fix by hand")
    return None


def process(text, po_db=None, beq=None):
    if not hasattr(process, "use_comments"):
        process.use_comments = True
    mask = code_mask(text)
    code = "".join(ch if m else " " for ch, m in zip(text, mask))
    flags = collect_flags(code)
    skipped, warnings, notes, changed, edits = [], [], [], [], []

    phases = find_phases(text, code)
    str_phases = [ph for ph in phases if ph["kw"] == "str"]
    n_hkl = len(phases) - len(str_phases)
    if n_hkl:
        notes.append(f"{n_hkl} hkl_Is phase(s) left unchanged (limits are only added in str phases)")
    spans = [(ph["start"], ph["end"]) for ph in str_phases]

    def in_str(pos):
        return any(a <= pos < b for a, b in spans)

    params = [p for p in find_cs(text, code, flags, skipped) if in_str(p["start"])]
    skipped[:] = [x for x in skipped
                  if not x[1].startswith(("CS", "Strain", "LVol", "e0_", "PO", "Preferred"))
                  or in_str(line_start_pos(text, x[0]))]
    for ph in phases:
        if ph["kw"] != "str" and len(ph["names"]) > 1:
            skipped.append((ph["ln"], f"CHECK: {ph['kw']} phase at line {ph['ln']} has "
                            f"{len(ph['names'])} phase_name lines ({', '.join(map(str, ph['names']))})"
                            " - a str phase may start inside it with a macro the script does not "
                            "recognise; those phases were left unchanged"))
    params += find_cell(text, code, flags, str_phases, skipped, notes,
                        first=phases[0]["start"] if phases else len(code))

    states = preprocessor_contexts(code)
    for p in params:
        p["ctx"] = context_at(states, p["insert_pos"])

    def phase_of(pos):
        for k, (a, b) in enumerate(spans):
            if a <= pos < b:
                return k
        return -1

    groups = defaultdict(list)
    for idx, p in enumerate(params):
        if not p["name"]:
            key = f"__unnamed_{idx}"
        else:
            key = (p["name"], phase_of(p["insert_pos"]))     # one parameter per str phase
        p["phase"] = phase_of(p["insert_pos"])
        groups[key].append(p)

    subgroups = []
    for key, group in groups.items():
        parts = split_coactive(group) if group[0]["name"] else [group]
        if len(parts) > 1:
            desc = "; ".join("lines " + ", ".join(str(p["ln"]) for p in part) for part in parts)
            notes.append(f"{group[0]['name']}: in {len(parts)} mutually exclusive #ifdef "
                         f"branches ({desc}) - limits set separately for each")
        for k, part in enumerate(parts):
            subgroups.append((part, k + 1 if len(parts) > 1 else None))

    for group, branch in subgroups:
        first = group[0]
        label = first["name"] or f"line {first['ln']} {first['label']}"
        if branch:
            label += f" [branch from line {first['ln']}]"
        if len({p["kind"] for p in group}) > 1:
            warnings.append(f"{label}: used for both cell and CS/strain - not changed")
            continue
        if first["kind"] == "cs":
            lo_def, hi_def = first["defaults"]
        else:
            lo_def, hi_def = cell_defaults(first["value"])
        if len(group) > 1 and len({norm(p["value"]) for p in group}) > 1:
            lines = ", ".join(f"{p['ln']} ({p['value']})" for p in group)
            warnings.append(f"{label}: start values differ on lines {lines} "
                            f"- not changed; please make them the same")
        lo_c = next((q["min"] for q in group if q["min"] is not None), lo_def)
        hi_c = next((q["max"] for q in group if q["max"] is not None), hi_def)
        lo_v, hi_v = plain_num(lo_c), plain_num(hi_c)
        for p in group:
            v = num(p["value"].split("`")[0])
            if v is not None and ((lo_v is not None and v < lo_v) or
                                  (hi_v is not None and v > hi_v)):
                warnings.append(f"{label}: start value {p['value']} on line {p['ln']} "
                                f"is outside min {plain(lo_c)} / max {plain(hi_c)}")
        lo = pick_limit(group, "min", lo_def, warnings, label)
        hi = pick_limit(group, "max", hi_def, warnings, label)
        for p in group:
            p["final"] = (plain(p["min"]) if p["min"] is not None else lo,
                          plain(p["max"]) if p["max"] is not None else hi)
            add = []
            if p["min"] is None and lo is not None:
                add.append(f"min {lo}")
            if p["max"] is None and hi is not None:
                add.append(f"max {hi}")
            if add:
                ins = " " + " ".join(add)
                pos = p["insert_pos"]
                # only max missing: put it after the existing min, keeping min before max
                if p["min"] is not None and p.get("min_end"):
                    pos = p["min_end"]
                edits.append((pos, ins))
                changed.append((p["ln"], p["label"], p["value"], ins.strip()))

    # names shared between different str phases: limits set per phase, flag for checking
    by_name = defaultdict(list)
    for p in params:
        if p["name"] and "final" in p:
            by_name[p["name"]].append(p)
    for name, inst in by_name.items():
        ph_ids = sorted({p["phase"] for p in inst})
        if len(ph_ids) < 2:
            continue
        lines = ", ".join(str(str_phases[k]["ln"]) for k in ph_ids if k >= 0)
        starts = {norm(p["value"]) for p in inst}
        finals = {tuple(None if x is None else norm(x) for x in p["final"]) for p in inst}
        diff = []
        if len(starts) > 1:
            diff.append("start values differ")
        if len(finals) > 1:
            diff.append("limits differ")
        what = (" - " + " and ".join(diff) + " between them (fine if they are never refined "
                "together; if they are, TOPAS treats them as one parameter)") if diff \
            else " - start values and limits are identical"
        skipped.append((inst[0]["ln"], f"CHECK: {name} is shared by {len(ph_ids)} str phases "
                        f"(phases at lines {lines}){what}; check this is intentional"))

    smart_report, no_po = [], []
    if po_db is not None:
        smart_po(text, code, str_phases, po_db, edits, skipped, smart_report, no_po,
                 mask=mask, use_comments=process.use_comments)
    process.no_po = no_po
    beq_report = []
    if beq:
        insert_beq(text, code, str_phases, beq, edits, skipped, beq_report)
    process.beq_report = beq_report

    seen, uniq = set(), []
    for item in skipped:
        if item not in seen:
            seen.add(item)
            uniq.append(item)
    skipped[:] = uniq

    out = text
    for pos, ins in sorted(edits, reverse=True):
        out = out[:pos] + ins + out[pos:]
    changed.sort()
    return out, changed, skipped, warnings, notes, smart_report


def read_text(path):
    raw = path.read_bytes()
    for enc in ("utf-8", "cp1252"):
        try:
            return raw.decode(enc), enc
        except UnicodeDecodeError:
            pass
    return raw.decode("latin-1"), "latin-1"


def main():
    ap = argparse.ArgumentParser(description="Insert missing min/max limits in a TOPAS .inp")
    ap.add_argument("inp", help="input .inp file")
    ap.add_argument("-o", "--output", help="output file (default: <name>_limits.inp)")
    ap.add_argument("--smart-po", action="store_true",
                    help="add fixed (not refined) March-Dollase PO lines in the likely direction "
                         "to str phases that have no PO model; check every direction yourself")
    ap.add_argument("--po-db", help=f"smart PO database (default: {PO_DB_NAME} next to the script)")
    bq = ap.add_mutually_exclusive_group()
    bq.add_argument("--insert-beq", action="store_true",
                    help=f"give every site with no beq the TOPAS GUI default, beq {BEQ_GUI}")
    bq.add_argument("--smart-beq", action="store_true",
                    help="give every site with no beq a value by element: H 2, O and similar "
                         "non-metals 1, metals 0.5")
    ap.add_argument("--no-comments", action="store_true",
                    help="smart PO: ignore comments (phase names and cells only)")

    args = ap.parse_args()

    src = Path(args.inp)
    if not src.is_file():
        sys.exit(f"File not found: {src}")
    dst = Path(args.output) if args.output else src.with_name(src.stem + "_limits" + src.suffix)
    if dst.resolve() == src.resolve():
        sys.exit("Output would overwrite the input; choose a different -o path.")

    text, enc = read_text(src)
    pre_warn = []
    db = None
    if args.po_db and not args.smart_po:
        sys.exit("--po-db only applies with --smart-po")
    if args.no_comments and not args.smart_po:
        sys.exit("--no-comments only applies with --smart-po")
    process.use_comments = not args.no_comments
    if args.smart_po:
        db_path = Path(args.po_db) if args.po_db else Path(__file__).with_name(PO_DB_NAME)
        if args.po_db and not db_path.is_file():
            sys.exit(f"PO database not found: {db_path}")
        db = load_po_db(db_path, pre_warn)
        db_note = db_path.name if db_path.is_file() else "not found - built-in clay list only"
    beq_mode = "smart" if args.smart_beq else "gui" if args.insert_beq else None
    out, changed, skipped, warnings, notes, smart_report = process(text, db, beq_mode)
    warnings[:0] = pre_warn
    with open(dst, "w", encoding=enc, newline="") as f:
        f.write(out)

    print(f"{Path(__file__).name} v{__version__}")
    print(f"Input : {src}\nOutput: {dst}")
    print("Smart PO: off" if db is None else f"Smart PO database: {db_note}")
    print({None: "Beq: off", "gui": f"Beq: missing beq set to {BEQ_GUI} (TOPAS GUI default)",
           "smart": "Beq: missing beq set by element (H 2, O and similar 1, metals 0.5)"}[beq_mode])

    print(f"\nPhases and notes: {len(notes)}")
    for n in notes:
        print(f"  {n}")
    print(f"\nChanged: {len(changed)}")
    for ln, label, value, ins in changed:
        print(f"  line {ln:>5}  {label} {value}  + {ins}")
    print(f"\nSmart PO added (fixed, not refined): {len(smart_report)}")
    for ln, pname, hkl, reason in smart_report:
        print(f"  line {ln:>5}  {pname}: PO {hkl}  ({reason})")
    beq_report = getattr(process, "beq_report", [])
    if beq_mode:
        print(f"\nBeq inserted: {sum(r[2] for r in beq_report)} sites in {len(beq_report)} phases")
        for ln, pname, n, detail in beq_report:
            print(f"  line {ln:>5}  {pname}: {n} site(s) ({detail})")
    no_po = getattr(process, "no_po", [])
    if db is not None and no_po:
        print(f"\nSmart PO not added - no common PO direction in the database: {len(no_po)}")
        for ln, pname, nm in no_po:
            print(f"  line {ln:>5}  {pname} ({nm if nm.startswith(chr(39)) else repr(nm)})")
    print(f"\nSkipped and symmetry checks: {len(skipped)}")
    for ln, why in skipped:
        print(f"  line {ln:>5}  {why}")
    print(f"\nWarnings: {len(warnings)}")
    for w in warnings:
        print(f"  {w}")


if __name__ == "__main__":
    main()
