from pathlib import Path
import hashlib
import json
import re

rxiv = Path(__file__).resolve().parent
bundle = rxiv.parent
md_file = rxiv / "manuscript.md"
fig_file = rxiv / "FIGURES/Figure_4.png"
receipt_file = bundle / "artifacts/figure_4_rerender_receipt.json"

print("Starting CrisPRO Forensic Validation Checks...")
text = md_file.read_text()
frontmatter = text.split("---", 2)[1]
title = next(line for line in frontmatter.splitlines() if line.startswith("title:")).lower()

assert not re.search(r"vs\s+42\s+WT", text, re.I), "legacy bowel denominator"
assert "p=0.114" not in text and "p = 0.114" not in text, "legacy bowel p-value"
assert "n=5 MBD4-LOF versus n=41" in text, "canonical bowel denominator missing"
print("PASS: Bowel denominators and p-value are canonical.")

assert not re.search(r"Research\s+Use\s+Only", text, re.I), "RUO remains in manuscript"
print("PASS: Research Use Only text is absent.")

assert "non-resolving for context-specific synthetic lethalities" not in text.lower()
assert "nonresolving for context-specific synthetic lethalities" not in text.lower()
print("PASS: Table S1 poison-pill sentence is absent.")

assert "pan-cancer" in title and "ovarian cancer" not in title
print("PASS: Title is pan-cancer rather than ovarian-restricted.")

assert "95.2%" in text and "single-allele clinical-trial relevance" in text
assert "future targeted sequencing is required to map the exact zygosity" not in text.lower()
print("PASS: Heterozygous LikelyLoF relevance is explicit.")

assert "### Replication Stress Checkpoint Axis: ATR and WEE1 inhibition" in text
assert "CLASS_CONCORDANT_WEE1I_SECONDARY" in text
print("PASS: ATR and WEE1 are consolidated under one checkpoint axis.")

for code_file in bundle.rglob("*.py"):
    if code_file.resolve() == Path(__file__).resolve():
        continue
    code = code_file.read_text()
    assert not re.search(r"['\"]RUO", code)
    assert not re.search(r"Research\s+Use\s+Only", code, re.I)
print("PASS: Figure-generation sources contain no RUO label.")

assert fig_file.exists() and fig_file.stat().st_size > 0
receipt = json.loads(receipt_file.read_text())
assert receipt["bytes_changed"] is True
assert receipt["before_sha256"] != receipt["after_sha256"]
assert hashlib.sha256(fig_file.read_bytes()).hexdigest() == receipt["after_sha256"]
print("PASS: Figure 4 was physically re-rendered and its bytes changed.")

print("[VERIFICATION STATUS: 100% CANONICAL DEPLOYMENT SUCCESS]")
