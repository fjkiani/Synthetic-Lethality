"""
Grounding guard — CI-failing check that a biomarker interpretation contains no
ungrounded statistic.

Rule: every number that looks like a reported statistic (HR, C-index, AUROC,
p-value, n, CI bound) in a shipped interpretation string must be reproducible
from the model's validation receipts. This is the automated defence against
the exact failure the PI3K/mTOR audit caught — a fabricated 'HR=1.30 / p=0.023'
hardcoded into an interpretation with no receipt behind it.
"""
from __future__ import annotations

import re
from typing import List, Set

from .models import BiomarkerModel

# numbers that look like reported statistics (floats, or n=INT).
# A number must not be preceded by a word char or '.', and must not be FOLLOWED
# by a word char or by '.<digit>' (i.e. a decimal continuation). A trailing
# sentence period (e.g. "HR was 1.30.") must still be caught, so the lookahead
# only rejects '.' when it precedes another digit.
_NUM_RE = re.compile(r"(?<![\w.])\d+(?:\.\d+)?(?!\w|\.\d)")

# tokens we always allow in prose without a receipt (years, panel sizes, etc.)
_ALLOWED_LITERALS = {"0", "1", "2", "3", "4", "5", "0.5", "100", "95"}


def _receipt_numbers(model: BiomarkerModel) -> Set[str]:
    """Collect every numeric value present in the model's receipts + panel size."""
    nums: Set[str] = set()

    def add(x):
        if x is None:
            return
        try:
            f = float(x)
        except (TypeError, ValueError):
            return
        # store multiple string renderings so '1.0'/'1.00'/'1' all match
        nums.add(f"{f:g}")
        nums.add(f"{f:.2f}")
        nums.add(f"{f:.3f}")
        if float(f).is_integer():
            nums.add(str(int(f)))

    for r in model.receipts:
        for v in (r.n, r.n_events, r.metric_value, r.ci95_low, r.ci95_high,
                  r.p_value, r.ph_test_p, r.epv, r.seed):
            add(v)
        for v in r.extra.values():
            add(v)
    add(len(model.gene_panel))
    return nums


def check_interpretation_grounded(model: BiomarkerModel) -> List[str]:
    """
    Return a list of ungrounded numeric tokens found in the interpretation.
    Empty list == fully grounded.
    """
    receipt_nums = _receipt_numbers(model) | _ALLOWED_LITERALS
    violations: List[str] = []
    for m in _NUM_RE.finditer(model.interpretation or ""):
        tok = m.group(0)
        try:
            f = float(tok)
        except ValueError:
            continue
        candidates = {f"{f:g}", f"{f:.2f}", f"{f:.3f}"}
        if float(f).is_integer():
            candidates.add(str(int(f)))
        if not (candidates & receipt_nums):
            violations.append(tok)
    return violations


def assert_grounded(model: BiomarkerModel) -> None:
    v = check_interpretation_grounded(model)
    if v:
        raise AssertionError(
            f"Ungrounded statistic(s) in interpretation of '{model.method_id}': {sorted(set(v))}. "
            "Every reported number must trace to a validation receipt."
        )
