"""
add_topas_limits.py

Inserts missing min/max limits into a TOPAS .inp file for:

  * CS / CS_L / CS_G        -> min 32, max 1000   (CS is shorthand for CS_L)
  * LVol_FWHM_CS_G_L(...)   -> its CS_G and CS_L arguments, as above
  * e0_from_Strain(...)     -> its Strain_G and Strain_L arguments, as below
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
  py add_topas_limits.py my_file.inp
  py add_topas_limits.py my_file.inp -o new_file.inp
"""

__version__ = "1.0.5"

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
}
CS_RE = re.compile(r"(?<![\w])(LVol_FWHM_CS_G_L|e0_from_Strain|CS_L|CS_G|CS|Strain_L|Strain_G)\(")
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
            pairs = [(ni, vi, kind) for ni, vi, kind in MULTI_MACROS[macro]
                     if vi < len(args) and code[args[vi][0]:args[vi][1]].strip()]
            label_fmt = macro + ":{}"
        else:
            pairs = [(0 if len(args) > 1 else None, len(args) - 1, macro)]
            label_fmt = "{}"

        for ni, vi, kind in pairs:
            val_start, val_end = args[vi]
            name = None
            if ni is not None:
                for tok in code[args[ni][0]:args[ni][1]].split():
                    bare = tok.lstrip("!@")
                    if bare and bare not in flags and IDENT_RE.fullmatch(tok):
                        name = bare
            label = label_fmt.format(kind)
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
        names = [line_no(text, m.start()) for m in
                 re.compile(r"(?<![\w])phase_name(?![\w])").finditer(code, pos, end)]
        phases.append(dict(kw=kw, start=pos, end=end, ln=line_no(text, pos), sg=sg,
                           names=names))
    return phases


def find_cell(text, code, flags, phases, skipped, notes):
    out = []
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


def process(text):
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
                  if not x[1].startswith(("CS", "Strain", "LVol", "e0_"))
                  or in_str(line_start_pos(text, x[0]))]
    params += find_cell(text, code, flags, str_phases, skipped, notes)

    states = preprocessor_contexts(code)
    for p in params:
        p["ctx"] = context_at(states, p["insert_pos"])

    groups = defaultdict(list)
    for idx, p in enumerate(params):
        groups[p["name"] if p["name"] else f"__unnamed_{idx}"].append(p)

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
    return out, changed, skipped, warnings, notes


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
    args = ap.parse_args()

    src = Path(args.inp)
    if not src.is_file():
        sys.exit(f"File not found: {src}")
    dst = Path(args.output) if args.output else src.with_name(src.stem + "_limits" + src.suffix)
    if dst.resolve() == src.resolve():
        sys.exit("Output would overwrite the input; choose a different -o path.")

    text, enc = read_text(src)
    out, changed, skipped, warnings, notes = process(text)
    with open(dst, "w", encoding=enc, newline="") as f:
        f.write(out)

    print(f"add_topas_limits.py v{__version__}")
    print(f"Input : {src}\nOutput: {dst}")
    print(f"\nPhases and notes: {len(notes)}")
    for n in notes:
        print(f"  {n}")
    print(f"\nChanged: {len(changed)}")
    for ln, label, value, ins in changed:
        print(f"  line {ln:>5}  {label} {value}  + {ins}")
    print(f"\nSkipped and symmetry checks: {len(skipped)}")
    for ln, why in skipped:
        print(f"  line {ln:>5}  {why}")
    print(f"\nWarnings: {len(warnings)}")
    for w in warnings:
        print(f"  {w}")


if __name__ == "__main__":
    main()
