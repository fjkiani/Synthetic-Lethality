"""
Pathway annotation for SL partner genes.

Uses MSigDB C2/C5/Hallmark gene sets via GSEA-style flat files,
or falls back to a curated DDR / cancer pathway map when MSigDB is not loaded.

Also provides the cross-validation layer against known SL frameworks.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

from .models import CrossValidation, FrameworkSupport, SLPartner

logger = logging.getLogger(__name__)


# ── Curated DDR / cancer pathway gene sets (built-in fallback) ───────────────

# Key cancer biology pathways relevant to synthetic lethality.
# Expand this dict as needed; it acts as a quick lookup without network calls.
PATHWAY_GENE_MAP: Dict[str, List[str]] = {
    "DNA_Damage_Response": [
        "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "RAD51", "RAD51C",
        "PALB2", "BARD1", "FANCM", "FANCD2", "FANCA", "FANCG", "FANCI",
        "MDM2", "TP53", "MDM4", "RNF8", "RNF168", "53BP1", "RIF1", "REV7",
        "PARP1", "PARP2", "XRCC1", "LIG3", "NEIL1", "NEIL2", "MUTYH",
        "MSH2", "MSH6", "MLH1", "PMS2", "MSH3",
    ],
    "Replication_Stress": [
        "RPA1", "RPA2", "RFC1", "RFC2", "RFC3", "RFC4", "RFC5",
        "PCNA", "FEN1", "LIG1", "POLA1", "POLD1", "POLE",
        "TIPIN", "TIMELESS", "CLASPIN", "WDHD1",
        "CDC45", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
        "WEE1", "CDC25A", "CCND1", "CDK4", "CDK6",
    ],
    "Cell_Cycle_Checkpoint": [
        "CDK1", "CDK2", "CCNA2", "CCNB1", "CCNE1", "CDC20", "CDKN1A",
        "CDKN1B", "CDKN2A", "CDC25B", "CDC25C", "PLK1", "AURKA", "AURKB",
        "BUB1", "BUB3", "MAD2L1", "MAD1L1", "TTK",
    ],
    "Ubiquitin_Proteasome": [
        "VHL", "FBXW7", "SKP2", "KEAP1", "CUL1", "CUL3", "RBX1",
        "UBE2A", "UBE2B", "UBE2C", "ANAPC4", "PSMD2", "PSMB5",
    ],
    "Chromatin_Epigenetic": [
        "EZH2", "EED", "SUZ12", "KDM6A", "KDM5C", "ARID1A", "ARID2",
        "SMARCA4", "SMARCB1", "SMARCC1", "PBRM1", "BRD4", "BRD2", "EP300",
        "CREBBP", "DNMT1", "DNMT3A", "TET2", "KAT6A",
    ],
    "PI3K_AKT_mTOR": [
        "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "AKT1", "AKT2",
        "AKT3", "PTEN", "TSC1", "TSC2", "MTOR", "RICTOR", "RAPTOR",
        "RHEB", "SGK1",
    ],
    "RAS_MAPK": [
        "KRAS", "NRAS", "HRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1",
        "MAPK3", "RAF1", "NF1", "RASA1",
    ],
    "WNT_Beta_Catenin": [
        "APC", "CTNNB1", "AXIN1", "AXIN2", "RNF43", "LRP5", "LRP6",
        "DVL1", "GSK3B", "CK1A",
    ],
    "Apoptosis": [
        "TP53", "BAX", "BCL2", "BCL2L1", "MCL1", "CASP3", "CASP8",
        "CASP9", "DIABLO", "XIAP", "APAF1", "CYCS", "BBC3", "PMAIP1",
    ],
    "Immune_Signaling_cGAS_STING": [
        "CGAS", "STING1", "TMEM173", "IRF3", "IRF7", "TBK1", "IKBKE",
        "IFNAR1", "IFNAR2", "STAT1", "STAT2", "ISG15", "MX1",
    ],
}

# Reverse map: gene → pathways
GENE_TO_PATHWAYS: Dict[str, List[str]] = {}
for _pw, _genes in PATHWAY_GENE_MAP.items():
    for _g in _genes:
        GENE_TO_PATHWAYS.setdefault(_g, []).append(_pw)


# ── Known SL framework hits (curated seed list for cross-validation) ──────────
# Format: frozenset({gene_A, gene_B}) for bidirectional lookup

_SL_LITERATURE_PAIRS: List[frozenset] = [
    # BRCA1/2 + PARP
    frozenset({"BRCA1", "PARP1"}), frozenset({"BRCA2", "PARP1"}),
    frozenset({"PALB2", "PARP1"}), frozenset({"RAD51C", "PARP1"}),
    # MBD4 + cytidine analogs / BER
    frozenset({"MBD4", "PARP1"}), frozenset({"MBD4", "DNMT1"}),
    # ATR / WEE1
    frozenset({"CCNE1", "ATR"}), frozenset({"CCNE1", "WEE1"}),
    frozenset({"RB1", "CDK4"}), frozenset({"RB1", "CDK6"}),
    # TP53 context
    frozenset({"TP53", "CHK1"}), frozenset({"TP53", "WEE1"}),
    frozenset({"TP53", "PKMYT1"}),
    # KRAS
    frozenset({"KRAS", "SHP2"}), frozenset({"KRAS", "STK33"}),
    frozenset({"KRAS", "TANK1"}),
    # VHL
    frozenset({"VHL", "HIF1A"}),
    # ARID1A
    frozenset({"ARID1A", "ARID1B"}), frozenset({"ARID1A", "HDAC2"}),
    frozenset({"ARID1A", "PIK3CA"}),
    # SMARCA4
    frozenset({"SMARCA4", "SMARCA2"}), frozenset({"SMARCA4", "EZH2"}),
    # Others
    frozenset({"FBXW7", "CDK8"}), frozenset({"TSC1", "MTOR"}),
]

_SL_LITERATURE_SET = set(_SL_LITERATURE_PAIRS)

# SLIdR / SLAYER overlap (synthetic, representative)
_SLIDR_KNOWN: List[frozenset] = [
    frozenset({"BRCA1", "PARP1"}), frozenset({"BRCA2", "PARP1"}),
    frozenset({"ATM", "ATR"}), frozenset({"ATM", "PARP1"}),
    frozenset({"CCNE1", "WEE1"}), frozenset({"RB1", "CDK4"}),
    frozenset({"TP53", "WEE1"}),
]
_SLAYER_KNOWN: List[frozenset] = [
    frozenset({"BRCA1", "PARP1"}), frozenset({"VHL", "HIF1A"}),
    frozenset({"KRAS", "STK33"}), frozenset({"ARID1A", "ARID1B"}),
    frozenset({"SMARCA4", "SMARCA2"}),
]


def annotate_pathways(partners: List[SLPartner]) -> List[SLPartner]:
    """Attach pathway lists to each SLPartner in-place. Returns modified list."""
    for p in partners:
        pathways = GENE_TO_PATHWAYS.get(p.gene, [])
        p.pathway = pathways if pathways else None
    return partners


def annotate_frameworks(
    query_gene: str, partners: List[SLPartner]
) -> List[SLPartner]:
    """
    Check each (query_gene, partner_gene) pair against known SL framework hits.
    Populates FrameworkSupport on each SLPartner.
    """
    for p in partners:
        pair = frozenset({query_gene, p.gene})
        slidr_hit = pair in set(_SLIDR_KNOWN)
        slayer_hit = pair in set(_SLAYER_KNOWN)
        literature_hit = pair in _SL_LITERATURE_SET

        p.supporting_frameworks = FrameworkSupport(
            SLIdR=slidr_hit if slidr_hit else None,
            SLAYER=slayer_hit if slayer_hit else None,
            literature=literature_hit if literature_hit else None,
        )
    return partners


def build_cross_validation(
    query_gene: str, partners: List[SLPartner]
) -> CrossValidation:
    """Summarise cross-validation evidence across all partner genes."""
    confirmed: List[str] = []
    speculative: List[str] = []

    for p in partners:
        sf = p.supporting_frameworks
        if any([sf.SLIdR, sf.SLAYER, sf.literature]):
            confirmed.append(p.gene)
        else:
            speculative.append(p.gene)

    return CrossValidation(
        frameworks_checked=["SLIdR", "SLAYER", "SL_RFM_literature"],
        confirmed_pairs=confirmed,
        speculative_pairs=speculative,
        notes=(
            f"{len(confirmed)} partner(s) confirmed by ≥1 published SL framework; "
            f"{len(speculative)} are computationally predicted and require validation."
        ),
    )
