# TOPAS-inp-Limits

Batch insertion of `min`/`max` limits for unit cell, crystallite size and strain parameters in TOPAS `.inp` files, plus a check that each phase's unit cell is refined and constrained as its space group requires.

Large quantitative phase analysis (QPA) input files can contain hundreds of `str` and `hkl_Is` phases. Missing limits let a refinement wander off into unphysical values (e.g. negative strain, crystallite sizes of millions of nm, cells drifting onto a different phase). This script adds the missing limits without touching anything else in the file.

## Requirements

- Python 3.8 or later (standard library only, nothing to install)
- Windows, macOS or Linux

## Usage

Put `add_topas_limits.py` in the folder with your `.inp` file, open a command prompt there and run:

```
py add_topas_limits.py my_file.inp
```

(On macOS/Linux use `python3` instead of `py`.)

- If the file name or path contains spaces, put it in double quotes:
  `py add_topas_limits.py "clinker QPA.inp"`
- The original file is **never modified**. The output is written to `my_file_limits.inp` in the same folder.
- Choose a different output name with `-o`:
  `py add_topas_limits.py my_file.inp -o my_file_new.inp`
- Save the report: `py add_topas_limits.py my_file.inp > report.txt`

Always open the output in TOPAS and do a test refinement before replacing your working file.

## What it does

### Limits inserted

| Parameter | Where | min | max |
|---|---|---|---|
| Crystallite size | `CS(`, `CS_L(`, `CS_G(`, and the CS arguments of `LVol_FWHM_CS_G_L(` | 32 | 1000 |
| Strain | `Strain_L(`, `Strain_G(`, and the strain arguments of `e0_from_Strain(` | 0 | 0.5 |
| Unit cell | `a b c al be ga` in every `str` / `hkl_Is`, and the `Cubic(`, `Tetragonal(`, `Hexagonal(`, `Trigonal(`, `Rhombohedral(` macros | start − 1 % | start + 1 % |

The defaults can be changed in the *Settings* block at the top of the script.

### Rules

- **Existing limits are never changed.** If only one of `min`/`max` is present, only the missing one is inserted, keeping `min` before `max`.
- **Constrained parameters get identical limits.** All instances of the same parameter name that can be active together get the same `min`/`max`: an existing limit on any instance is copied to the others, otherwise the default is used (for the cell, ±1 % of the first instance's start value).
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

## Symmetry check

For phases with a recognised space group (not set up with the `Cubic(`/`Tetragonal(`/... macros, which enforce symmetry themselves), the report lists `CHECK` lines for:

- cell parameters that symmetry allows to refine but are fixed (no parameter name or `@`). Not reported for `hkl_Is` phases, which are assumed to be PONKCS phases with fixed cells;
- angles that are refined but fixed by symmetry (e.g. `be @ 90` in a cubic phase);
- parameters that must be equal (cubic a = b = c, tetragonal/hexagonal a = b, rhombohedral a = b = c and al = be = ga) but are not tied by a shared name or an equation such as `b = Get(a);`;
- angles that are not 90°/120° although symmetry fixes them (suggests a wrong space group).

Parameters switched by a `#define` flag instead of `@`/`!` are not checked for fixed/refined status, because that depends on the launch mode.

## Report

The script prints:

- **Phases and notes**: each phase's space group, crystal system and refinable angles, plus names found in mutually exclusive `#ifdef` branches;
- **Changed**: every inserted limit with its line number;
- **Skipped and symmetry checks**: skipped values and `CHECK` lines;
- **Warnings**: instances of a shared name with different existing limits (nothing inserted for that limit; fix by hand), different start values, or a cell start value outside its limits.

## Limitations

- Instances that are mutually exclusive by some means other than `#ifdef` (e.g. switched by a macro) cannot be detected, so they are treated as active together.
- Only the macros listed above are recognised. Size or strain parameters inside your own macros are only handled where the standard macro calls appear in the `.inp` itself.

## AI-use disclosure

The script was written by Claude (Anthropic) from Henry Spratt's specifications. Henry Spratt defined the requirements, reviewed and error-checked the output, and tested it on his own TOPAS input files.

## Licence

MIT. See [LICENSE](LICENSE).
