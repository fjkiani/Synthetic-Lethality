"""
Configuration — environment-driven, no hardcoded secrets.
All values can be overridden via environment variables or a .env file.
"""
from functools import lru_cache
from pathlib import Path
from typing import Optional

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        case_sensitive=False,
        extra="ignore",
    )

    # ── App ──────────────────────────────────────────────────────────────────
    app_name: str = "Synthetic Lethality Mapping Agent"
    app_version: str = "1.0.0"
    debug: bool = False
    log_level: str = "INFO"

    # ── DepMap ───────────────────────────────────────────────────────────────
    depmap_base_url: str = "https://depmap.org/portal/api"
    depmap_download_base: str = "https://depmap.org/portal/download/api"
    depmap_release: str = "latest"          # e.g. "25Q2" or "latest"
    depmap_cache_dir: Path = Path(".cache/depmap")
    audit_db_path: Path = Path(".cache/audit_queue.db")

    # Figshare / S3 direct download URLs (updated each DepMap quarter)
    # These are the stable Figshare URLs from the 24Q4 public release.
    depmap_crispr_url: str = (
        "https://figshare.com/ndownloader/files/47596327"  # CRISPRGeneEffect.csv
    )
    depmap_mutation_url: str = (
        "https://figshare.com/ndownloader/files/47596516"  # OmicsSomaticMutations.csv
    )
    depmap_cna_url: str = (
        "https://figshare.com/ndownloader/files/47596471"  # OmicsCNGene.csv
    )
    depmap_expression_url: str = (
        "https://figshare.com/ndownloader/files/47596480"  # OmicsExpressionProteinCodingGenesTPMLogp1.csv
    )
    depmap_sample_info_url: str = (
        "https://figshare.com/ndownloader/files/47596612"  # Model.csv
    )
    depmap_prism_url: str = (
        "https://figshare.com/ndownloader/files/20237786"  # PRISM Repurposing secondary (2020)
    )
    depmap_prism_meta_url: str = (
        "https://figshare.com/ndownloader/files/20237777"  # PRISM Repurposing compound meta
    )

    # ── Analysis thresholds ───────────────────────────────────────────────────
    min_cell_lines_cancer_specific: int = 30   # minimum lines for cancer-type analysis
    fdr_threshold: float = 0.25                 # BH-corrected FDR cutoff
    delta_dep_threshold: float = 0.15           # |Δdependency| minimum to report
    dependency_threshold: float = -0.5          # Chronos/CERES score below this = dependent

    # ── External drug DBs ─────────────────────────────────────────────────────
    # DEPRECATED — OncoKB requires a commercial license and is no longer used.
    # Drug evidence is now provided by the open KB stack (CIViC + CGI + JAX).
    # These fields are retained for backward-compat env parsing only; they are
    # not read by any active code path.  Remove in a future major version bump.
    oncokb_api_token: Optional[str] = Field(default=None, alias="ONCOKB_API_TOKEN")
    oncokb_base_url: str = "https://www.oncokb.org/api/v1"  # DEPRECATED

    chembl_base_url: str = "https://www.ebi.ac.uk/chembl/api/data"

    # GDSC (Sanger) bulk download
    gdsc1_url: str = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC1_fitted_dose_response_27Oct23.xlsx"
    gdsc2_url: str = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx"
    gdsc_cache_dir: Path = Path(".cache/gdsc")

    # ── Cache ────────────────────────────────────────────────────────────────
    cache_ttl_seconds: int = 86_400   # 24 h
    chembl_cache_dir: Path = Path(".cache/chembl")

    @field_validator("depmap_cache_dir", "gdsc_cache_dir", "chembl_cache_dir", mode="before")
    @classmethod
    def _ensure_path(cls, v):
        p = Path(v)
        p.mkdir(parents=True, exist_ok=True)
        return p


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    return Settings()
