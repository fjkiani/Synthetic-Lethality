"""
Biomarker registry — the single source of truth for which biomarker
capabilities exist and whether they are validated.

Key guarantee: get() returns a VALIDATED model only if it carries a passing
external-replication receipt on an independent cohort. Discovery-only models
are hidden unless explicitly requested.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

from .models import BiomarkerModel, BiomarkerStatus, BiomarkerType


class BiomarkerRegistry:
    def __init__(self) -> None:
        self._models: Dict[str, BiomarkerModel] = {}

    # ── registration ─────────────────────────────────────────────────────────
    def _key(self, cancer: str, btype: BiomarkerType, method_id: str) -> str:
        return f"{cancer}::{btype.value}::{method_id}"

    def register(self, model: BiomarkerModel) -> None:
        # enforce the honesty gate at registration: a model may only claim
        # VALIDATED if it actually carries a passing external-replication receipt
        if model.status == BiomarkerStatus.VALIDATED and model.external_replication() is None:
            raise ValueError(
                f"Refusing to register '{model.method_id}' as VALIDATED without a passing "
                "external-replication receipt on an independent cohort."
            )
        self._models[self._key(model.cancer, model.biomarker_type, model.method_id)] = model

    # ── query ─────────────────────────────────────────────────────────────────
    def get(
        self,
        cancer: str,
        biomarker_type: BiomarkerType,
        include_discovery: bool = False,
    ) -> List[BiomarkerModel]:
        """
        Return biomarker models for (cancer, type). By default ONLY validated
        models are returned. Set include_discovery=True to also see discovery-only
        models (clearly flagged by their .status).
        """
        out = [
            m for m in self._models.values()
            if m.cancer == cancer and m.biomarker_type == biomarker_type
        ]
        if not include_discovery:
            out = [m for m in out if m.is_validated()]
        return out

    def all_models(self) -> List[BiomarkerModel]:
        return list(self._models.values())

    # ── persistence ────────────────────────────────────────────────────────────
    def to_json(self, path: str | Path) -> None:
        payload = {k: json.loads(m.model_dump_json()) for k, m in self._models.items()}
        Path(path).write_text(json.dumps(payload, indent=2))

    @classmethod
    def from_json(cls, path: str | Path) -> "BiomarkerRegistry":
        reg = cls()
        data = json.loads(Path(path).read_text())
        for _, m in data.items():
            reg.register(BiomarkerModel(**m))
        return reg
