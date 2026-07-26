"""
Dataset resolver + non-overlap guard for biomarker training/validation.

The single most important guard here: a "validation" cohort must not share
any patient/sample with the training cohort. Same-cohort reprocessing does NOT
count as external replication (this is exactly the error the PI3K/mTOR audit
caught: ov_tcga_pan_can_atlas_2018 vs ov_tcga are the SAME patients).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Sequence, Set


class CohortOverlapError(ValueError):
    """Raised when train and validation cohorts share samples/patients."""


@dataclass
class CohortSplit:
    """A recorded, receipt-able train/validation split."""
    train_cohort_id: str
    val_cohort_id: str
    train_ids: List[str] = field(default_factory=list)
    val_ids: List[str] = field(default_factory=list)

    def assert_non_overlap(self) -> "CohortSplit":
        """
        Fail loudly if any id appears in both train and validation.
        Also fail on the trivially-wrong case of identical cohort ids.
        """
        if self.train_cohort_id == self.val_cohort_id:
            raise CohortOverlapError(
                f"train and validation cohort ids are identical ({self.train_cohort_id}); "
                "this is same-cohort reprocessing, not external replication."
            )
        overlap: Set[str] = set(map(_norm, self.train_ids)) & set(map(_norm, self.val_ids))
        if overlap:
            example = sorted(overlap)[:5]
            raise CohortOverlapError(
                f"{len(overlap)} sample/patient id(s) appear in BOTH train ({self.train_cohort_id}) "
                f"and validation ({self.val_cohort_id}); e.g. {example}. "
                "External replication requires non-overlapping cohorts."
            )
        return self

    @property
    def is_independent(self) -> bool:
        try:
            self.assert_non_overlap()
            return True
        except CohortOverlapError:
            return False


def _norm(sample_id: str) -> str:
    """
    Normalise ids for overlap detection. TCGA sample barcodes collapse to the
    patient barcode (first 3 hyphen-fields, e.g. TCGA-13-0720). Everything else
    is compared uppercased/trimmed.
    """
    s = str(sample_id).strip().upper()
    if s.startswith("TCGA-"):
        parts = s.split("-")
        if len(parts) >= 3:
            return "-".join(parts[:3])
    return s


def patients_from_samples(sample_ids: Iterable[str]) -> Set[str]:
    """Collapse a set of sample ids to unique patient ids (TCGA-aware)."""
    return {_norm(s) for s in sample_ids}
