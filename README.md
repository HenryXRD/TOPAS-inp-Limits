# TOPAS-inp-Limits

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24585892.svg)](https://doi.org/10.5281/zenodo.24585892)

Batch insertion of `min`/`max` limits for unit cell, crystallite size, strain and preferred-orientation parameters in TOPAS `.inp` files, plus a check that each phase's unit cell is refined and constrained as its space group requires. Optional switches add fixed preferred-orientation lines in likely directions (smart PO) and fill in missing `beq` values.

Large quantitative phase analysis (QPA) input files can contain hundreds of `str` and `hkl_Is` phases. Missing limits let a refinement wander off into unphysical values (e.g. negative strain, crystallite sizes of millions of nm, cells drifting onto a different phase). This script adds the missing limits without touching anything else in the file.

## Requirements

- Python 3.8 or later (standard library only, nothing to install)
- Windows, macOS or Linux

## Usage

Put `add_topas_limits_v2.py` in the folder with your `.inp` file (and `po_directions.csv` next to the script if you use `--smart-po`), open a command prompt there and run:

```
py add_topas_limits_v2.py my_file.inp
```

(On macOS/Linux use `python3` instead of `py`.)

- If the file name or path contains spaces, put it in double quotes:
  `py add_topas_limits_v2.py "clinker QPA.inp"`
- The original file is **never modified**. The output is written to `my_file_limits.inp` in the same folder.
- Choose a different output name with `-o`:
  `py add_topas_limits_v2.py my_file.inp -o my_file_new.inp`
- Save the report: `py add_topas_limits_v2.py my_file.inp > report.txt`
- Add suggested preferred-orientation lines (see [Smart PO](#smart-po-optional)): `py add_topas_limits_v2.py my_file.inp --smart-po`
- Smart PO from names and cells only, ignoring comments: `py add_topas_limits_v2.py my_file.inp --smart-po --no-comments`
- Give sites with no `beq` the TOPAS GUI default (see [Beq](#beq-optional)): `py add_topas_limits_v2.py my_file.inp --insert-beq`
- Give sites with no `beq` a value by element: `py add_topas_limits_v2.py my_file.inp --smart-beq`

Always open the output in TOPAS and do a test refinement before replacing your working file.

## What it does

### Limits inserted

| Parameter | Where | min | max |
|---|---|---|---|
| Crystallite size | `CS(`, `CS_L(`, `CS_G(`, and the CS arguments of `LVol_FWHM_CS_G_L(` | 32 | 1000 |
| Strain | `Strain_L(`, `Strain_G(`, and the strain arguments of `e0_from_Strain(` | 0 | 0.5 |
| Preferred orientation (March-Dollase) | `PO(`, `Preferred_Orientation(` | 0.5 | 1.5 |
| Two-direction PO: both March-Dollase values | `PO_Two_Directions(` (v1, v2) | 0.5 | 1.5 |
| Two-direction PO: fraction in the first direction | `PO_Two_Directions(` (w1v) | 0 | 1 |
| Unit cell | `a b c al be ga` in every `str` / `STR(...)` phase, and the `Cubic(`, `Tetragonal(`, `Hexagonal(`, `Trigonal(`, `Rhombohedral(` macros | start − 1 % | start + 1 % |

Limits are only added inside `str` / `STR(...)` phases. `hkl_Is` phases (e.g. PONKCS phases, whose small crystallite sizes are deliberate) and anything outside a phase (e.g. an `xdd` block) are left unchanged.

The defaults can be changed in the *Settings* block at the top of the script.

### Rules

- **Existing limits are never changed.** If only one of `min`/`max` is present, only the missing one is inserted, keeping `min` before `max`.
- **Constrained parameters get identical limits.** Within a phase, all instances of the same parameter name that can be active together get the same `min`/`max`: an existing limit on any instance is copied to the others, otherwise the default is used (for the cell, ±1 % of the first instance's start value).
- **Names shared between phases are flagged.** Each `str` phase gets its own limits, and the report lists every parameter name used in more than one `str` phase, saying whether start values or limits differ between them. Check these are intentional: fine for alternative phases that are never refined together; if they are refined together, TOPAS treats them as one parameter, so make the start values and limits identical.
- **`#ifdef` aware.** Instances of a name that sit in mutually exclusive `#ifdef`/`#ifndef`/`#if`/`#elseif`/`#else` branches (never active together) are treated as separate parameters, each with its own limits.
- **Only refinable cell angles get limits.** `a`, `b`, `c` are always eligible. Angles are eligible only if the phase's `space_group` allows them to refine:

  | Crystal system | Angles given limits |
  |---|---|
  | Triclinic | `al be ga` |
  | Monoclinic | unique-axis angle only (the non-90° angle in the file, otherwise from the symbol, default `be`) |
  | Rhombohedral axes (e.g. `R-3c:R`) | `al be ga` |
  | Orthorhombic, tetragonal, hexagonal, trigonal (hexagonal axes), cubic | none |

  Space groups are read in compact (`P21/c`), TOPAS (`P2_1/c`), spaced full (`"P 1 21/c 1"`) or number form, including `:R`/`:H` settings. Phases without a `space_group` fall back to the angle values (an angle that is not 90° or 120° is treated as refinable) and this is reported.
- Equation values (`= ...;`) and non-numeric values (macro arguments, `#define` names) are skipped.
- Comments (`'` and `/* */`) and `"strings"` are ignored, so text such as a phase named "CS" is never touched. Macros are only recognised when `(` follows the name directly, e.g. `CS(`.
- Phases started with the `STR(sg)` macro are recognised, with the space group taken from its first argument.
- Values with an error suffix (`` 5.4310`_0.0002 ``, `5.4310_0.0002`, `` 5.4310` ``) are read as one value; limits go after the suffix.
- Element symbols after `occ` (e.g. `occ al 0.5`, `occ c 1`) are never mistaken for the cell parameters `al`, `c`, `b`, `be`, `ga`.
- Identifiers created with `#define` or `macro` are treated as flags, not parameter names.
- Line endings and text encoding are preserved.

### Example

Before:

```
str
   phase_name "Pigeonite"
   space_group "C2/c"
   a @ 9.7 b @ 8.9 c @ 5.25 al 90 be @ 108.5 ga 90
   LVol_FWHM_CS_G_L( 1, 0, 0.89, 0,,,cs_pigeonite, 1000 min =45;)
   e0_from_Strain( 0,!st_pigeonite, 0.1,,)
```

After:

```
str
   phase_name "Pigeonite"
   space_group "C2/c"
   a @ 9.7 min 9.6030 max 9.7970 b @ 8.9 min 8.8110 max 8.9890 c @ 5.25 min 5.1975 max 5.3025 al 90 be @ 108.5 min 107.4150 max 109.5850 ga 90
   LVol_FWHM_CS_G_L( 1, 0, 0.89, 0,,,cs_pigeonite, 1000 min =45; max 1000)
   e0_from_Strain( 0,!st_pigeonite, 0.1 min 0 max 0.5,,)
```

`al` and `ga` are fixed at 90° by monoclinic symmetry, so they get no limits.

## Smart PO (optional)

With `--smart-po`, each `str` phase that has no preferred-orientation model (`PO`, `Preferred_Orientation`, `PO_Spherical_Harmonics` or any other `PO…` macro) gets a March-Dollase line in a common direction for that mineral. The line is **fixed, not refined**:

```
PO(!po_Calcite, 1 min 0.5 max 1.5, , 1 0 4) ' smart PO: likely direction ... not refined
```

A value of 1 means no preferred orientation, so the refinement is unchanged until you decide to refine it (change `!po_Calcite` to `po_Calcite`, or to `@`).

**The directions are suggestions only. The onus is on the user to confirm that the direction is correct for their sample and for the cell setting in their structure before refining it.**

How the direction is chosen:

1. **`po_directions.csv`**: phase names (and aliases) matched as whole words against `phase_name`, e.g. Calcite → `1 0 4`, Gypsum → `0 1 0`, Muscovite → `0 0 1`. Edit the file to add, change or remove entries; use `--po-db <file>` to use a different one.
2. **Built-in clay list** (kaolinite, illite, smectite, montmorillonite, chlorite, vermiculite, talc, pyrophyllite and others) → `0 0 1`.
3. **Unnamed layer silicates**: phases whose cell has the layer-silicate sheet (a ≈ 5.0–5.5 Å, b ≈ 8.7–9.4 Å, b ≈ √3·a) → `0 0 1`, reported as "clay inferred from cell".
4. **Unnamed chain silicates**: phases whose cell matches a clinoamphibole (a ≈ 9.5–10.1, b ≈ 17.5–18.6, c ≈ 5.2–5.45 Å, β ≈ 101–108°) or clinopyroxene (a ≈ 9.4–10.0, b ≈ 8.5–9.3, c ≈ 5.1–5.5 Å, β ≈ 105–110°) → `1 1 0`, or an orthopyroxene / orthoamphibole → `2 1 0`; always marked UNCERTAIN.

Names are matched as whole words, and chemical prefixes may be joined on, so "Magnesiohornblende", "Ferro-actinolite" and "Fluoro-edenite" all match their root name.

#### Comments in the .inp

If the phase name gives no match, comments are used as hints (switch this off with `--no-comments`):

- **Comments on the phase**: full-line comments directly above the `str` line, and a comment at the end of the `str` or `phase_name` line, are matched like a phase name, e.g. `str ' hedenbergite from AMCSD`.
- **Section headings**: the nearest heading comment above the phase, such as `' 9.DE Amphiboles`, `' Phyllosilicates`, `' Feldspars`, `' Spinel group` or `' 4 OXIDES`. Group names, plurals ("Micas"), "... group" and Nickel-Strunz codes (9.DA pyroxenes, 9.DE amphiboles, 9.D inosilicates, 9.E phyllosilicates; other classes by number) are recognised. A single mineral name is not taken as a heading for the phases below it.

Anything taken from a comment is marked UNCERTAIN, except a layer silicate whose cell also has the sheet signature. A heading never overrides the cell checks:
- a heading that names minerals with different directions (e.g. "Corundum-hematite group") is treated as ambiguous
- under a non-silicate heading (oxides, carbonates, clinker, standards ...) no direction is inferred from the cell
- a cell that contradicts its heading (an amphibole-like cell under "Feldspars") is reported as a CHECK and nothing is added
- a phase name that disagrees with its heading gets a note in the report, which helps catch phases filed in the wrong section

Commented-out code (e.g. `'str phase_name "Calcite"`) and comments before the first phase (unless directly above it) are ignored.

Before inserting, the cell is checked against the expected setting (e.g. hexagonal axes for the carbonates, `c` the stacking axis for platy phases, `b` the longest axis for gypsum). If it does not match, nothing is inserted and a `CHECK` line is reported instead. `hkl_Is` phases are never changed.

`po_directions.csv` columns:

| Column | Meaning |
|---|---|
| `phase` | name matched against `phase_name` |
| `aliases` | other names, separated by `;` |
| `h`, `k`, `l` | PO direction |
| `check` | `hex`, `a_longest`, `b_longest`, `c_longest`, `a_shortest`, `b_shortest`, `c_shortest`, `sheet_ab`, or blank for no check |
| `confidence` | `check` where sources disagree or the direction depends on the cell setting; blank otherwise |
| `notes` | ignored by the script |

Entries marked `check` are still inserted, but their TOPAS comment says the direction is **UNCERTAIN** and must be checked and changed if needed before refining; the report marks them too, e.g.

```
PO(!po_Augite, 1 min 0.5 max 1.5, , 1 1 0) ' smart PO: UNCERTAIN direction (sources disagree or it depends on the setting) ... CHECK and change 1 1 0 if needed before refining; not refined
```

Leave `h`, `k`, `l` blank for minerals with no useful PO direction (e.g. olivine, garnet, pyrite, magnetite, clinker phases): they are recognised but no line is added, and the report lists them. Rows whose `phase` starts with `#` are section headings and are ignored.

The supplied `po_directions.csv` covers about 100 common cement, rock-forming and ore minerals (carbonates, sulfates, hydroxides, oxides, sulfides, halides, phosphates and silicates), with each direction taken from the mineral's cleavage or habit (Mindat, Webmineral, Handbook of Mineralogy) and noted in the `notes` column.

## Beq (optional)

Thermal parameters (`beq`) affect the refined scale factors, which directly influence QPA results. A site with no `beq` is not treated the same way everywhere:

- in TOPAS **launch mode**, a `str` with no `beq` uses **0**
- when the same `str` or CIF is imported in TOPAS **GUI mode**, it uses **1**

This means even the same operator on the same data, using the exact same `str` or CIF files, can get a different QPA result, depending only on how the file was run.

Two optional switches fill in the missing values so the result no longer depends on the mode. Only sites with no `beq` are changed; existing `beq` values, `beq` macros and anisotropic ADPs (`u11`… / `adps`) are left alone. The inserted `beq` is fixed, not refined.

| Switch | Value given to a site with no `beq` |
|---|---|
| `--insert-beq` | `1`, the TOPAS GUI default |
| `--smart-beq` | by element: H (and D) `2`; O and similar anion-forming non-metals (C, N, F, S, Cl, Se, Br, I) `1`; metals and metalloids (incl. Si, B, P, As) `0.5` |

With `--smart-beq`:
- a non-metal given with a positive charge (e.g. `S+6` in sulfate, `C+4`) is treated as a cation and gets `0.5`; `S`, `S-2` or `O-2` get `1`
- element symbols are read in any case (`occ fe`, `occ Fe+3`, `occ FE`)
- on a mixed site the species with the largest occupancy decides (e.g. `occ Fe 0.9 occ Ni 0.1` → `0.5`)
- a site whose element is not recognised is left unchanged and reported as a CHECK

The two switches cannot be used together. The values are starting points, not refined results; check they suit your phases.

## Symmetry check

For `str` phases with a recognised space group (not set up with the `Cubic(`/`Tetragonal(`/... macros, which enforce symmetry themselves), the report lists `CHECK` lines for:

- cell parameters that symmetry allows to refine but are fixed (no parameter name or `@`);
- angles that are refined but fixed by symmetry (e.g. `be @ 90` in a cubic phase);
- parameters that must be equal (cubic a = b = c, tetragonal/hexagonal a = b, rhombohedral a = b = c and al = be = ga) but are not tied by a shared name or an equation such as `b = Get(a);`;
- angles that are not 90°/120° although symmetry fixes them (suggests a wrong space group);
- phases the script may have missed: a block with more than one `phase_name`, or cell parameters before the first recognised phase. This usually means a phase starts with a macro the script does not recognise (only `str`, `STR(` and `hkl_Is` are recognised as phase starts), so the previous phase's space group was applied to it.

Parameters switched by a `#define` flag instead of `@`/`!` are not checked for fixed/refined status, because that depends on the launch mode.

## Report

The script prints:

- **Phases and notes**: each phase's space group, crystal system and refinable angles, plus names found in mutually exclusive `#ifdef` branches;
- **Changed**: every inserted limit with its line number;
- **Beq inserted** (with `--insert-beq` or `--smart-beq`): sites given a `beq`, counted per phase by value;
- **Smart PO added** (with `--smart-po`): every PO line inserted, its direction and why it was chosen;
- **Skipped and symmetry checks**: skipped values and `CHECK` lines;
- **Warnings**: instances of a shared name with different existing limits (nothing inserted for that limit; fix by hand), different start values, or a start value outside its limits. Limits written as purely numeric equations (e.g. `min = 5.086 - 5.086 0.01;`, where TOPAS reads the space as multiplication) are evaluated for this check and for comparing shared limits; equations that use parameter names are not.

## Limitations

- Within a phase, instances that are mutually exclusive by some means other than `#ifdef` (e.g. switched by a macro) cannot be detected, so they are treated as active together.
- Only the macros listed above are recognised. Size or strain parameters inside your own macros are only handled where the standard macro calls appear in the `.inp` itself.

## How to cite

Spratt, H. *TOPAS-inp-Limits* (software). Zenodo. https://doi.org/10.5281/zenodo.24585892

GitHub's **Cite this repository** button (from `CITATION.cff`) gives the same citation in APA and BibTeX formats.

## AI-use disclosure

The script was written by Claude (Anthropic) from Henry Spratt's specifications. Henry Spratt defined the requirements, reviewed and error-checked the output, and tested it on his own TOPAS input files.

## Licence

MIT. See [LICENSE](LICENSE).
