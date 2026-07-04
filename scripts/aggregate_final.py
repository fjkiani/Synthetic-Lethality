#!/usr/bin/env python3
"""
AACR 2026 Final Aggregation Script
Merges enriched_groq_w0-4.jsonl + keyword layer, applies QC, generates gap report.
Run after all 5 enrichment workers complete.
"""
import json, shutil, subprocess, os
from collections import Counter, defaultdict
from datetime import datetime

SHARED = "/mnt/shared-workspace/shared"
RESULTS = "/mnt/results"

CANONICAL_AXES = [
    "PARP_INHIBITORS", "ATR_WEE1", "CCNE1_AMP", "PKMYT1", "WRN",
    "PRMT5_MTAP", "IMMUNOTHERAPY", "CYTIDINE_ANALOGS", "TP53_REPLICATION_STRESS", "PI3K_AKT"
]
KB_CANCERS = ["Ovarian", "CRC", "PDAC", "Breast", "Lung", "Prostate", "Heme", "Bladder"]
REQUIRE_KW_CONFIRM = {"PKMYT1", "WRN", "CCNE1_AMP", "TP53_REPLICATION_STRESS"}
TRUST_LLM = {"IMMUNOTHERAPY", "PI3K_AKT", "PARP_INHIBITORS", "ATR_WEE1",
              "CYTIDINE_ANALOGS", "PRMT5_MTAP"}

WORKER_FILES = {
    "w0": (f"{SHARED}/enriched_groq_w0.jsonl", 0, 1497),
    "w1": (f"{SHARED}/enriched_groq_w1.jsonl", 1497, 2994),
    "w2": (f"{SHARED}/enriched_groq_w2.jsonl", 2994, 4491),
    "w3": (f"{SHARED}/enriched_groq_w3.jsonl", 4491, 5988),
    "w4": (f"{SHARED}/enriched_groq_w4.jsonl", 5988, 7485),
}

print("=== AACR 2026 Final Aggregation ===")
print(f"Started: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}")

# 1. Load keyword layer
print("\n1. Loading keyword layer...")
keyword_map = {}
with open(f"{SHARED}/keyword_enriched_7485.jsonl") as f:
    for line in f:
        if line.strip():
            r = json.loads(line)
            keyword_map[r["id"]] = r["keyword_enrichment"]["keyword_axes"]
print(f"   {len(keyword_map)} keyword records")

# 2. Load all LLM enrichment
print("\n2. Loading LLM enrichment outputs...")
all_llm = {}
for wid, (fpath, start, end) in WORKER_FILES.items():
    count = 0
    with open(fpath) as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                all_llm[r["id"]] = r
                count += 1
    success = sum(1 for r in all_llm.values() if "error" not in r.get("enrichment", {}))
    print(f"   {wid}: {count} records")
print(f"   Total: {len(all_llm)} LLM records")

# 3. Load abstract order
all_ids = []
with open(f"{SHARED}/aacr_abstracts_7485.jsonl") as f:
    for line in f:
        if line.strip():
            all_ids.append(json.loads(line)["id"])

# 4. Merge
print("\n3. Merging LLM + keyword layers...")
merged = []
llm_success = kw_fallback = error_kept = 0

for abs_id in all_ids:
    llm_rec = all_llm.get(abs_id)
    kw_axes = keyword_map.get(abs_id, [])

    if llm_rec and "error" not in llm_rec.get("enrichment", {}):
        rec = llm_rec.copy()
        # Clean axes: trust LLM for common, require KW confirm for rare
        raw_axes = rec["enrichment"].get("crispro_axes") or []
        cleaned = []
        for ax in raw_axes:
            if ax in TRUST_LLM:
                cleaned.append(ax)
            elif ax in REQUIRE_KW_CONFIRM and ax in kw_axes:
                cleaned.append(ax)
        for ax in kw_axes:
            if ax in REQUIRE_KW_CONFIRM and ax not in cleaned:
                cleaned.append(ax)
        rec["enrichment"]["crispro_axes"] = list(set(cleaned))
        merged.append(rec)
        llm_success += 1
    elif kw_axes:
        ke_full = {}
        with open(f"{SHARED}/keyword_enriched_7485.jsonl") as f:
            pass  # already loaded above
        # Build from keyword_map
        merged.append({
            "id": abs_id, "title": llm_rec.get("title","") if llm_rec else "",
            "doi": llm_rec.get("doi","") if llm_rec else "",
            "biblio": llm_rec.get("biblio",{}) if llm_rec else {},
            "enrichment": {
                "signal_strength": "low", "crispro_axes": kw_axes,
                "fit_score": 0.1, "source": "keyword_fallback",
                "brenus_relevance": "none",
            }
        })
        kw_fallback += 1
    else:
        merged.append(llm_rec or {"id": abs_id, "enrichment": {"signal_strength":"none","error":"not_processed"}})
        error_kept += 1

print(f"   LLM success: {llm_success:,} | KW fallback: {kw_fallback:,} | Error: {error_kept:,}")

# 5. Compute stats
signals = Counter(r["enrichment"].get("signal_strength","none") for r in merged)
axis_counts = Counter()
matrix = defaultdict(lambda: defaultdict(int))
for r in merged:
    ct = r["enrichment"].get("cancer_type","")
    for ax in (r["enrichment"].get("crispro_axes") or []):
        axis_counts[ax] += 1
        if ct in KB_CANCERS:
            matrix[ax][ct] += 1

rs_counts = {f: sum(1 for r in merged if r["enrichment"].get("rs_features",{}).get(f)) 
             for f in ["MSI_H","TP53_LOF","MYC_amplified","ARID1A_LOF","CCNE1_amplified","MBD4_LOF"]}
brenus_high = sum(1 for r in merged if r["enrichment"].get("brenus_relevance")=="high")
brenus_med  = sum(1 for r in merged if r["enrichment"].get("brenus_relevance")=="medium")

# 6. Save merged JSONL
print("\n4. Saving merged dataset...")
out_path = "/workspace/enriched_final_7485.jsonl"
with open(out_path, "w") as f:
    for r in merged:
        f.write(json.dumps(r) + "\n")
shutil.copy(out_path, f"{SHARED}/enriched_final_7485.jsonl")
shutil.copy(out_path, f"{RESULTS}/enriched_final_7485.jsonl")
print(f"   Saved {len(merged)} records")

# 7. Generate final gap report
print("\n5. Generating final gap report...")
now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

report = f"""# AACR 2026 CrisPRO Gap Report — FINAL
**Generated:** {now}  
**Corpus:** 7,485 AACR 2026 abstracts (Cancer Research 7_Supplement + 8_Supplement)  
**Coverage:** {llm_success+kw_fallback:,}/{len(merged):,} ({(llm_success+kw_fallback)/len(merged)*100:.1f}%)  
**LLM model:** Groq llama-3.3-70b-versatile (primary) + keyword fallback  

---

## 1. Signal Matrix — Axes × Cancer Type

| Axis | Ovarian | CRC | PDAC | Breast | Lung | Prostate | Heme | Bladder |
|---|---|---|---|---|---|---|---|---|
"""
for ax in CANONICAL_AXES:
    row = f"| {ax} |"
    for ct in KB_CANCERS:
        cnt = matrix[ax][ct]
        row += f" {cnt} |" if cnt > 0 else " — |"
    report += row + "\n"

report += f"""
## 2. RS Feature Distribution

| Feature | Count | % |
|---|---|---|
"""
for feat, cnt in rs_counts.items():
    report += f"| {feat} | {cnt} | {cnt/7485*100:.1f}% |\n"

report += f"""
## 3. Coverage Summary

| Metric | Value |
|---|---|
| Total abstracts | 7,485 |
| LLM-enriched | {llm_success:,} ({llm_success/7485*100:.1f}%) |
| Keyword-fallback | {kw_fallback:,} ({kw_fallback/7485*100:.1f}%) |
| High signal | {signals.get("high",0):,} |
| Medium signal | {signals.get("medium",0):,} |
| Brenus-relevant | {brenus_high+brenus_med:,} |
| Abstracts with axis hits | {sum(1 for r in merged if r["enrichment"].get("crispro_axes")):,} |
"""

with open("/workspace/gap_report_final.md", "w") as f:
    f.write(report)
shutil.copy("/workspace/gap_report_final.md", f"{RESULTS}/gap_report_aacr2026_final.md")
print(f"   Gap report saved")

# 8. Git push
print("\n6. Pushing to GitHub...")
repo_path = "/workspace/Synthetic-Lethality"
if os.path.exists(repo_path):
    shutil.copy("/workspace/enriched_final_7485.jsonl", f"{repo_path}/data/aacr2026/enriched_final_7485.jsonl")
    shutil.copy("/workspace/gap_report_final.md", f"{repo_path}/reports/gap_report_aacr2026.md")
    subprocess.run(["git", "-C", repo_path, "add", "-A"], check=True)
    subprocess.run(["git", "-C", repo_path, "commit", "-m",
                    f"feat(aacr2026): enriched corpus + gap report ({llm_success+kw_fallback}/{len(merged)} abstracts)"],
                   check=True)
    subprocess.run(["git", "-C", repo_path, "push", "origin", "main"], check=True)
    print("   Pushed to fjkiani/Synthetic-Lethality main")
else:
    print("   Repo not found — skipping git push")

print(f"\n=== DONE: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')} ===")
