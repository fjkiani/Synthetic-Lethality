# bioRxiv PDF build — SL platform manuscript

```bash
cd publications/00-mbd4-manuscript/mbd4-framework/rxiv
pandoc manuscript.md \
  -o sl_platform_manuscript_biorxiv_submission.pdf \
  --pdf-engine=tectonic \
  --filter pandoc-crossref \
  --citeproc
```

- Bibliography: `sl_platform_manuscript_references.bib` (do not invent entries)
- Figures: `FIGURES/` (copied from `../suplimentary/sl_platform_manuscript_FIGURES_*.png`)
- Supplementary CSVs: `../suplimentary/Table_S1` / `Table_S2` (condensed in manuscript)
