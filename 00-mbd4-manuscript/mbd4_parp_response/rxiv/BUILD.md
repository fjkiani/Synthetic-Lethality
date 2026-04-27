# bioRxiv PDF build (recorded commands)

## Where you are (two layouts)

**A. Monorepo (this repo: `crispr-assistant-main`)**  
Bundle path: `publications/00-mbd4-manuscript/mbd4_parp_response/`. All `cd` examples below use this prefix unless you are in the standalone clone (B).

**B. Standalone canonical repo**  
Clone [MBD4-LOF-Dual-Therapeutic-Vulnerability](https://github.com/crispro-ai/MBD4-LOF-Dual-Therapeutic-Vulnerability): the manuscript bundle lives at repo root with `rxiv/` and `artifacts/` siblings. Replace monorepo paths with `rxiv/` and `rxiv/FIGURES/` from that repo root (no `publications/...` prefix).

## Prerequisites

- `pandoc` (tested with Homebrew install on PATH)
- `pandoc-crossref` (`pandoc-crossref` binary on PATH)
- `tectonic` (`--pdf-engine=tectonic`)
- Python 3 + `matplotlib` + `numpy` (figure scripts)
- Optional: `PyPDF2` for automated text-layer checks (`python3 -m pip install PyPDF2`)

## Figure 3 (MSI purge panel title + assets)

**Monorepo:**

```bash
cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES
python3 fig3_stress_tests.py
```

**Standalone:** from clone root:

```bash
cd rxiv/FIGURES
python3 fig3_stress_tests.py
```

Writes `fig3_stress_tests.png` and `fig3_stress_tests.pdf` in that directory.

Panel A title in the script: `A. MSI Purge (DepMap ModelSubtypeFeatures)` (aligned with manuscript stress-test wording and MSI / ModelSubtypeFeatures filter logic).

**Panel C:** p-values are read from `artifacts/canonical_atr_wee1_rerun_20260405/ceralasertib_ln_ic50_leave_one_out_MANUSCRIPT_RECEIPT.json` (`p_value_one_sided_mwu_lof_less_wt` per iteration). Do not paste rounded literals that drift from the frozen max p (`0.045165724128583974`).

## Figure 4 (PARP1 expression + PARPi Z)

**Monorepo:**

```bash
cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES
python3 fig4_parp1_expression.py
```

**Standalone:**

```bash
cd rxiv/FIGURES
python3 fig4_parp1_expression.py
```

- **Repo root:** script resolves `OmicsExpression.csv` and DepMap parquet relative to the repository root (`mbd4_parp_response` → `publications` → `00-mbd4-manuscript` → repo in monorepo; or clone root in standalone). Override with `CRISPR_ASSISTANT_ROOT=/path/to/crispr-assistant-main` if needed.
- **Bundle paths:** `artifacts/axis_c_preclinical/gdsc2_parpi_per_line_receipts.csv` is read relative to `mbd4_parp_response/` (monorepo) or the standalone repo root containing `artifacts/`.
- **Panel A:** two-sided Mann-Whitney (expression hypothesis exploration); **Panel B:** Spearman on merged lines (numbers should align with `parp1_parpi_spearman_MANUSCRIPT_RECEIPT.json` when inputs match).

Requires: `pandas`, `pyarrow` (for parquet), `scipy`, `matplotlib`, `numpy`.

## Manuscript PDF (bibliography + citeproc + figure crossrefs)

**Monorepo:**

```bash
cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv
pandoc manuscript.md \
  -o mbd4_parp_response_biorxiv_submission.pdf \
  --pdf-engine=tectonic \
  --filter pandoc-crossref \
  --citeproc
```

**Standalone:**

```bash
cd rxiv
pandoc manuscript.md \
  -o mbd4_parp_response_biorxiv_submission.pdf \
  --pdf-engine=tectonic \
  --filter pandoc-crossref \
  --citeproc
```

- Bibliography: `references.bib` is set in YAML (`bibliography: references.bib`); `citeproc` resolves `[@key]` citations.
- Crossrefs: `pandoc-crossref` resolves `@fig:...` to numbered “Figure …” in text.
- **Tectonic + `\xmpquote`:** Pandoc’s LaTeX template can emit `\xmpquote` inside `\hypersetup{pdfkeywords=...}`. Tectonic fails with “Undefined control sequence” unless defined. `manuscript.md` includes:

  ```yaml
  header-includes: |
    \providecommand{\xmpquote}[1]{#1}
  ```

- **Author block:** Structured `authors:` does not fill PDF `\author{...}`. Putting `\\` or `\newline` inside YAML `author:` fails because Pandoc escapes backslashes (`\\` → `\textbackslash{}`), so the title page uses `author: []` plus `\AtBeginDocument{\def\@author{...}}` in `header-includes` (raw LaTeX line breaks and `\and`). **Title:** avoid `\mbox{ATR}~...` in YAML `title:` — Pandoc turns `~` into `\textasciitilde` (visible tilde in PDF); use plain “ATR inhibitor” / “PARP1” text.

## Output

- Final PDF: `rxiv/mbd4_parp_response_biorxiv_submission.pdf` (path relative to bundle / standalone `rxiv/`).

## Verification (text layer; optional)

**Monorepo:**

```bash
cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv
python3 << 'PY'
from PyPDF2 import PdfReader
r = PdfReader("mbd4_parp_response_biorxiv_submission.pdf")
t = "".join((p.extract_text() or "") for p in r.pages)
assert "[@" not in t, "raw citation keys"
assert "@fig:" not in t, "raw figure keys"
assert "References" in t
assert "Kiani" in t
print("basic checks ok")
PY
```

**Standalone:** `cd rxiv` then the same Python block.

**Note:** Tectonic may warn about missing Unicode glyphs (e.g. ρ, ≤, ≥, superscripts) in Latin Modern; the PDF still builds. For perfect glyph coverage, use XeLaTeX/LuaLaTeX with a Unicode math font if available.

**Figures 2–3 vs captions:** Panel titles and plot annotations live inside the rasterized PNGs; confirm visually in the PDF. Manuscript figure captions are under `## Figures` in `manuscript.md`.

## Commands run in this verification session

```bash
cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv/FIGURES && python3 fig3_stress_tests.py

cd publications/00-mbd4-manuscript/mbd4_parp_response/rxiv
pandoc manuscript.md -o mbd4_parp_response_biorxiv_submission.pdf \
  --pdf-engine=tectonic --filter pandoc-crossref --citeproc

python3 -m pip install --quiet PyPDF2   # if missing
```
