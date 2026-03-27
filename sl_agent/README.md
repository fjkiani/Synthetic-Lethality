# Synthetic Lethality Mapping Agent — Backend v2

Production-grade FastAPI backend that executes a 4-step synthetic lethality (SL)
pipeline against DepMap CRISPR data, PRISM drug screens, ChEMBL, and an
**open knowledge-base stack** (CIViC + CGI + JAX + ClinVar).  
No commercial licenses required.  All data sources are CC0 or open-access.

> Part of the **CrisPRO** precision oncology platform (S / P / E framework).  
> Validated use case: **MBD4 + gemcitabine** — npj Precis Oncol 2022 ([PMID 35428381](https://pubmed.ncbi.nlm.nih.gov/35428381/))

---

## Architecture

```
sl_agent/
├── __init__.py             # version export (2.0.0)
├── core/
│   ├── config.py           # Environment-driven settings (Pydantic)
│   ├── models.py           # Typed domain models (SLQueryInput → SLMapResult)
│   ├── sl_engine.py        # Step 1: SL computation (Wilcoxon, FDR, delta-dep)
│   ├── drug_mapper.py      # Step 2–3: PRISM + ChEMBL + open KB + rank scoring
│   ├── pathway_annotator.py # Pathway labels + SLIdR/SLAYER cross-validation
│   └── orchestrator.py     # Wires all layers; singleton data store
├── data/
│   └── depmap_loader.py    # Download/cache DepMap matrices (Figshare 24Q4)
├── api/
│   ├── app.py              # FastAPI factory + lifespan + middleware
│   ├── routes.py           # REST endpoints (/analyze, /genes, /cancer_types)
│   └── kb_routes.py        # KB endpoints (/kb/query, /kb/explain, /kb/sources)
├── kb/
│   ├── models.py           # KBEvidence, DrugRecommendation, EvidenceTier (A/B/C/D)
│   ├── evidence_fuser.py   # Dedup + consensus tier upgrade + confidence scoring
│   ├── kb_engine.py        # Orchestrates CIViC+CGI+JAX in parallel; oncokb_drop_in()
│   └── clients/
│       ├── civic_client.py    # GraphQL + nightly TSV fallback (CC0)
│       ├── cgi_client.py      # Bulk TSV (CC0) — cgi_biomarkers_latest.tsv
│       ├── jax_client.py      # Curated seed: BRCA1/2, EGFR, BRAF, KRAS, MBD4, TP53 …
│       └── clinvar_client.py  # eutils filter layer (oncogenicity weighting)
├── tests/
│   ├── test_sl_engine.py   # 9 tests
│   ├── test_drug_mapper.py # 7 tests
│   └── test_api.py         # 5 tests (mocked)
├── kb/tests/
│   └── test_kb.py          # 28 tests (all offline)
├── utils/                  # Reserved for future shared helpers
├── main.py
├── requirements.txt
├── pyproject.toml          # Build config + pytest settings
└── .env.example            # Annotated environment template
```

---

## Pipeline Steps

### Step 1 — Candidate SL Partner Identification (`sl_engine.py`)
- Loads DepMap Chronos gene effect matrix (24Q4, Figshare)
- Splits cell lines into **mutant vs WT** groups for the query gene
  - Supports: LOF mutation · GOF · homozygous deletion · amplification · any mutation
- Applies **cancer-type-specific** stratification (≥ 30 lines), falls back to pan-cancer
- For each candidate partner gene:
  - **Wilcoxon rank-sum test** (one-sided: mutant more dependent)
  - **BH FDR** multiple-testing correction
  - **Cohen's d** effect size
  - **Co-dependency** (Pearson r of dependency profiles)
- Filters by `fdr_cutoff` and `delta_dep_cutoff`
- **Pan-essential filter**: removes genes where ≥ 90 % of cell lines have Chronos < −0.5
  *(uses dependency frequency — not percentile of medians — to prevent mutant-group
  signal inflation into the blacklist)*

### Step 2 — Drug Mapping (`drug_mapper.py`)
- **PRISM repurposing screen**: Wilcoxon on viability LFC, mutant vs WT
- **ChEMBL REST API**: inhibitor lookup by gene target, MoA, max clinical phase
- **Open KB stack** (replaces OncoKB 1:1, no token required):
  - `oncokb_gene_drugs()` → alias → `open_kb_drugs()` → `oncokb_drop_in()` in `kb_engine.py`
  - Returns identical `max_phase` field so the rank scorer is unchanged
- **GDSC**: Available but not auto-loaded — see [GDSC Enrichment](#gdsc-enrichment) below

### Step 3 — Ranking
Composite `rank_score` (0–1):
```
rank = 0.45 × sl_signal + 0.35 × drug_response_differential + 0.20 × druggability
```
Where:
- `sl_signal        = |Δdependency| / 2 × (1 - FDR)`
- `drug_response    = |Δviability|  / 3 × (1 - drug_FDR)`
- `druggability     = max(phase_score, kb_level_score)`

### Step 4 — Structured JSON Output
Returns `SLMapResult` with:
- `input_context` — gene, cancer type, scope, n lines, frequencies
- `sl_partners[]` — ranked SL genes with stats, pathways, framework support
- `gene_drug_pairs[]` — ranked drug candidates with PRISM + external evidence
- `cross_validation` — SLIdR / SLAYER / literature confirmed vs speculative pairs

---

## Open KB Stack

The agent replaces OncoKB with a fully open, license-free evidence tier system.

### Sources

| Source | License | Access | Notes |
|--------|---------|--------|-------|
| **CIViC** | CC0 | GraphQL API + nightly TSV fallback | `https://civicdb.org` |
| **CGI** | CC0 | Bulk TSV download | `https://cancergenomeinterpreter.org` |
| **JAX CKB** | Curated seed | Local JSON (hand-curated) | BRCA1/2, EGFR, BRAF, KRAS, MBD4, TP53, CCNE1, CDK12, ARID1A, PALB2, RAD51C/D, ATM |
| **ClinVar** | Public domain | NCBI eutils | Oncogenicity filter layer |

### Evidence Tiering (maps 1:1 to OncoKB levels)

| Open KB tier | Equivalent OncoKB level | Sources → Tier |
|---|---|---|
| **Tier A** | LEVEL_1 / R1 | CIViC Level A · CGI "FDA approved" |
| **Tier B** | LEVEL_2 / R2 | CIViC Level B · CGI "Phase III" · 2× Tier C consensus boost |
| **Tier C** | LEVEL_3A / 3B | CIViC Level C/D · CGI "Phase I/II" · JAX seed entries |
| **Tier D** | LEVEL_4 | ClinVar pathogenic, no drug link yet |

`max_phase` mapping: A → 4 · B → 3 · C → 2 · D → 0

### Endpoints

```
GET  /api/v1/kb/query?gene=MBD4              # query open KB for a gene
GET  /api/v1/kb/explain?gene=MBD4&drug=gemcitabine  # explain evidence chain
GET  /api/v1/kb/sources                      # list active KB sources + health
```

---

## Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Copy env template
cp .env.example .env
# (edit .env as needed — all fields have sensible defaults)

# Run server (DepMap matrices download on first startup, ~2–5 min)
uvicorn sl_agent.api.app:app --reload --port 8000
```

---

## Usage

### POST `/api/v1/analyze`
```json
{
  "gene": "MBD4",
  "mutation_type": "loss_of_function",
  "cancer_type": "Colorectal Cancer",
  "top_n_partners": 20,
  "top_n_drugs": 10,
  "fdr_cutoff": 0.25,
  "delta_dep_cutoff": 0.15,
  "include_pathway_context": true,
  "include_codependency": true
}
```

**mutation_type** options:
| Value | Description |
|-------|-------------|
| `any_mutation` | Any non-silent somatic mutation |
| `loss_of_function` | Nonsense, frameshift, splice, truncating |
| `gain_of_function` | Activating missense |
| `homozygous_deletion` | CNA log2 < −1.0 |
| `amplification` | CNA log2 > 1.0 |

### GET `/api/v1/genes?prefix=BRCA`
List all genes in the current CRISPR matrix.

### GET `/api/v1/cancer_types`
List available cancer lineages and primary diseases.

### POST `/api/v1/analyze/async` + GET `/api/v1/result/{job_id}`
Submit async job, poll for result.

---

## Example Output (truncated)
```json
{
  "status": "success",
  "result": {
    "input_context": {
      "query_gene": "MBD4",
      "mutation_type": "loss_of_function",
      "cancer_type": "Colorectal Cancer",
      "scope": "cancer_specific",
      "n_total_lines": 52,
      "n_mut_lines": 11,
      "n_wt_lines": 41,
      "depmap_release": "24Q4"
    },
    "sl_partners": [
      {
        "gene": "POLB",
        "delta_dependency": -0.61,
        "p_value": 3.4e-04,
        "fdr": 0.0091,
        "effect_size_cohend": -1.44,
        "pathway": ["Base_Excision_Repair"],
        "supporting_frameworks": { "SLIdR": true, "SLAYER": false, "literature": true }
      }
    ],
    "gene_drug_pairs": [
      {
        "partner_gene": "POLB",
        "drug_name": "Gemcitabine",
        "drug_class": "Nucleoside analog",
        "external_evidence": {
          "source": "JAX",
          "kb_tier": "Tier_B",
          "max_phase": 3,
          "pmids": ["35428381"]
        },
        "rank_score": 0.74
      }
    ]
  }
}
```

---

## Tests

```bash
# All 47 tests — offline, no external calls
cd /path/to/project   # parent of sl_agent/
python -m pytest sl_agent/tests/ sl_agent/kb/tests/ -v

# Or via pyproject.toml (from sl_agent/ directory)
cd sl_agent && python -m pytest
```

Expected output: **47 passed** (9 engine + 7 drug_mapper + 5 api + 28 kb)

---

## GDSC Enrichment

GDSC IC50 response data is **available but not auto-loaded** (the Excel file is ~150 MB
and most deployments prefer the lighter PRISM-only path).

To enable GDSC enrichment:
1. Call `load_gdsc_viability()` from `data/depmap_loader.py`
2. Add `_gdsc` / `gdsc()` class methods to `DataStore` in `orchestrator.py`
3. Pass the loaded DataFrame to `map_gene_to_drugs(gdsc_meta=…)` in the drug-mapping loop

GDSC download URLs (Sanger release 8.5):
```
GDSC1: https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx
GDSC2: https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx
```

---

## Configuration (`.env`)

See `.env.example` for all options.  Key settings:

```env
DEPMAP_RELEASE=24Q4
FDR_THRESHOLD=0.25
DELTA_DEP_THRESHOLD=0.15
MIN_CELL_LINES_CANCER_SPECIFIC=30
LOG_LEVEL=INFO
```

---

## Extension Points

| What | Where |
|------|--------|
| Enable GDSC IC50 enrichment | `orchestrator.py` + `drug_mapper.gdsc_drugs_for_gene()` |
| Add MSigDB pathway sets | `pathway_annotator.py` → extend `PATHWAY_GENE_MAP` |
| Add SLIdR/SLAYER API calls | `pathway_annotator.py` → `_SLIDR_KNOWN` / `_SLAYER_KNOWN` |
| Add another open KB source | `kb/clients/` → new client + register in `kb_engine.py` |
| Add Redis job queue | `api/routes.py` → replace `_jobs` dict |
| Add TCGA mutation frequencies | `orchestrator.py` → after stratification step |
| Add cGAS-STING pathway scoring | `pathway_annotator.py` → `Immune_Signaling_cGAS_STING` set |
| Add shared utility functions | `utils/` directory (currently reserved) |

---

## Data Sources

| Dataset | Release | URL |
|---------|---------|-----|
| DepMap CRISPR (Chronos) | 24Q4 | [Figshare 47596327](https://figshare.com/ndownloader/files/47596327) |
| DepMap Somatic Mutations | 24Q4 | [Figshare 47596516](https://figshare.com/ndownloader/files/47596516) |
| DepMap CNA | 24Q4 | [Figshare 47596471](https://figshare.com/ndownloader/files/47596471) |
| DepMap Expression | 24Q4 | [Figshare 47596480](https://figshare.com/ndownloader/files/47596480) |
| DepMap Sample Info | 24Q4 | [Figshare 47596612](https://figshare.com/ndownloader/files/47596612) |
| PRISM Repurposing | 2020 | [Figshare 20237786](https://figshare.com/ndownloader/files/20237786) |
| CIViC Evidence | Nightly | [civicdb.org nightly TSV](https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv) |
| CGI Biomarkers | Latest | [cancergenomeinterpreter.org](https://www.cancergenomeinterpreter.org/data/cgi_biomarkers_latest.tsv) |
| ChEMBL | REST API | [ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl/api/data) |

---

## Multi-Modal Evidence Engine (v4)

The core enhancement in v4: the system refuses to hallucinate SL where receipts
are weak, and refuses to discard promising axes just because one modality (CRISPR)
is negative.

### Doctrine

| Principle | Implementation |
|-----------|---------------|
| CRISPR = one axis, not the arbiter | `modality_fuser.py` weights CRISPR = 1 vs clinical = 4 |
| CRISPR negative + pharma unexamined ≠ "no vulnerability" | Guardrail in `matrix_builder.py` |
| Compound sensitivity without CRISPR dep = "synthetic sickness" | `disagree_crispr_neg_pharma_pos` agreement category |
| Cytidine analogs = gold-standard calibration bar | `literature_receipts.py` + `/analyze/calibration` |

### New API Endpoints

```
POST /api/v1/analyze/multimodal    — full multi-modal evidence matrix + narrative
GET  /api/v1/analyze/calibration   — cytidine-analog gold-standard narrative
GET  /api/v1/analyze/evidence_matrix?gene=MBD4   — matrix from frozen receipts only
```

### POST `/api/v1/analyze/multimodal`
```json
{
  "gene": "MBD4",
  "mutation_type": "loss_of_function",
  "cancer_type": "Colorectal Cancer",
  "axes": ["cytidine_analogs", "parp_inhibitors", "atr_wee1", "wrn", "immunotherapy"],
  "include_pharmacologic_stratification": true,
  "include_literature_receipts": true,
  "include_calibration_narrative": true,
  "stratify_by_msi": true
}
```

Returns `MultiModalResult`:
```json
{
  "evidence_matrix": {
    "rows": [
      {
        "axis": "cytidine_analogs",
        "axis_label": "Cytidine Analogs (gemcitabine, cytarabine)",
        "recommendation_tier": "Validated SL therapeutic lever",
        "overall_evidence_level": "High",
        "positive_modalities": ["clinical", "in_vivo", "in_vitro", "expression"],
        "missing_modalities": ["crispr", "prism"]
      },
      {
        "axis": "parp_inhibitors",
        "recommendation_tier": "Mechanistic candidate only",
        "overall_evidence_level": "Mechanistic-only",
        "crispr_pharmacologic_agreement": "insufficient_data"
      }
    ]
  },
  "agreement_report": [...],
  "updated_narrative": "=== MULTI-MODAL EVIDENCE NARRATIVE: MBD4 ===\n...",
  "recommendation_map": {
    "Cytidine Analogs (gemcitabine, cytarabine)": "Validated SL therapeutic lever",
    "PARP Inhibitors (olaparib, niraparib, ...)": "Mechanistic candidate only",
    "ATR / WEE1 Inhibitors": "Mechanistic candidate only",
    "WRN Helicase Inhibitors": "Mechanistic candidate only",
    "Immunotherapy / Checkpoint Inhibitors": "Mechanistic candidate only"
  },
  "guardrails_applied": [
    "[PARP Inhibitors] CRISPR is negative but pharmacologic screens have NOT been examined. Cannot declare 'no vulnerability'..."
  ]
}
```

### Evidence Matrix Structure

Each candidate axis produces a row with 7 modality columns:

| Modality | What it captures |
|----------|-----------------|
| `crispr` | Δdep, FDR, n_mut, n_wt from DepMap Wilcoxon |
| `expression` | RNA/protein co-variation with MBD4 LOF |
| `prism` | PRISM repurposing screen IC50/AUC stratified by MBD4 status |
| `gdsc` | GDSC IC50 stratified by MBD4 status (optional enrichment) |
| `in_vitro` | Published KO/WT/rescue isogenic functional assays |
| `in_vivo` | PDX / xenograft response data |
| `clinical` | Patient-level response, PFS, OS |

Each cell has a `ModalityStatus`: `positive` / `negative` / `mixed` / `missing` / `confounded`

### Recommendation Tiers

| Tier | Criteria |
|------|----------|
| **Validated SL therapeutic lever** | Multiple modalities positive including clinical/in vivo |
| **Strong candidate dependency axis** | CRISPR and/or compound screens + mechanism |
| **Mechanistic candidate only** | Expression + mechanism; CRISPR/compound weak or absent |
| **Not supported / negative** | Flat or conflicting across all modalities examined |

### Guardrail: CRISPR Negative ≠ Dead

If CRISPR is negative but pharmacologic screens have NOT been examined, the system
automatically applies a guardrail and sets tier = "Mechanistic candidate only"
(not "Not supported / negative"). This is enforced by a dedicated test:
`test_crispr_neg_pharma_missing_is_not_negative`.

### New Module Structure

```
sl_agent/multimodal/
├── __init__.py
├── models.py              — EvidenceMatrix, EvidenceRow, ModalityEvidence, MultiModalResult
├── pharmacologic_analyzer.py — Stratified PRISM/GDSC analysis + MSI confound check
├── literature_receipts.py    — Frozen validated receipts (MBD4 gold standard)
├── modality_fuser.py         — Evidence level scoring + recommendation tier logic
├── matrix_builder.py         — End-to-end pipeline orchestration
└── tests/
    └── test_multimodal.py    — 37 tests (all offline)
```

### Cytidine Analog Calibration (Gold Standard)

Access via `GET /api/v1/analyze/calibration?gene=MBD4`:
- Shows all 4 frozen modality receipts (in vitro, in vivo, clinical, expression)
- Defines the bar: "3 modalities + clinical receipt = Validated"
- Everything else is graded against this bar
- Source: npj Precis Oncol 2022, [PMID 35428381](https://pubmed.ncbi.nlm.nih.gov/35428381/)

### Tests

```bash
# All 84 tests — fully offline
cd /path/to/project
python -m pytest sl_agent/tests/ sl_agent/kb/tests/ sl_agent/multimodal/tests/ -v
# Expected: 84 passed
```
