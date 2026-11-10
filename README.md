# Synthetic Lethality — CrisPRO companion repo

Public companion for MBD4 / synthetic-lethality manuscripts and related analysis code.

## Canon (start here)

| Path | What |
|------|------|
| [`00-mbd4-manuscript/mbd4_parp_response/rxiv/`](00-mbd4-manuscript/mbd4_parp_response/rxiv/) | **MBD4-LOF dual vulnerability preprint** (manuscript, figures, bib, build, submission PDF) |
| [`00-mbd4-manuscript/mbd4-framework/rxiv/`](00-mbd4-manuscript/mbd4-framework/rxiv/) | SL platform manuscript pack |
| [`sl_agent/`](sl_agent/) | Synthetic-lethality agent source |
| [`sl_discovery_pancancer/`](sl_discovery_pancancer/) | Pan-cancer discovery / validation scripts |
| [`docs/archive/`](docs/archive/) | Historical sprint notes (not part of the preprint) |

## Rebuild the MBD4 preprint PDF

```bash
cd 00-mbd4-manuscript/mbd4_parp_response/rxiv
# follow BUILD.md (pandoc + citeproc + tectonic/crossref)
```

Submission draft (field pack): `00-mbd4-manuscript/mbd4_parp_response/rxiv/BIORXIV_SUBMISSION_DRAFT.md`

## Data

Primary analyses use public DepMap / GDSC2 resources:

- https://depmap.org/portal/
- https://www.cancerrxgene.org/

Receipts and derived tables live under `00-mbd4-manuscript/mbd4_parp_response/` (see `SUPPLEMENTAL_DATA_PROVENANCE.md`).

## License / use

Manuscript content is for scholarly preprint distribution (bioRxiv CC BY when posted). Code under `sl_agent/` follows that package’s license/readme.
