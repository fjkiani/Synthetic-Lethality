"""
Curated literature-anchored KNOWN synthetic-lethality panel (pan-cancer).

No local SL reference DB (SynLethDB/ISLE/SLKB) exists, so this panel is hand-curated
from primary literature. Inclusion criterion (locked in PLAN.md):
  each pair needs (1) published primary evidence AND (2) a targetable node.

Each entry:
  driver         : the altered gene (loss or activation)
  mode           : 'LOF' (loss) or 'ACT' (activating) -> how to gate mutant lines
  partner        : the dependency / SL partner gene (CRISPR target)
  drug_axis      : canonical drug axis for pharmacologic confirmation (or None)
  falsification  : agents that should NOT show the effect (mechanistic control) or None
  cn_event       : True if driver is a copy-number event (detect via CN/expression, not point mut)
  block          : which worker block runs it (A=DDR/checkpoint, B=chromatin/paralog/metabolic, C=oncogene-context)
  citation       : primary reference (short)
"""

KNOWN_SL_PANEL = [
    # ---- BLOCK A: DDR + PARP + replication checkpoint ----
    dict(driver="BRCA1", mode="LOF", partner="PARP1", drug_axis="parp_inhibitors",
         falsification=None, cn_event=False, block="A",
         citation="Farmer 2005 Nature 10.1038/nature03445; Bryant 2005 Nature 10.1038/nature03443"),
    dict(driver="BRCA2", mode="LOF", partner="PARP1", drug_axis="parp_inhibitors",
         falsification=None, cn_event=False, block="A",
         citation="Bryant 2005 Nature 10.1038/nature03443"),
    dict(driver="PALB2", mode="LOF", partner="PARP1", drug_axis="parp_inhibitors",
         falsification=None, cn_event=False, block="A",
         citation="Buisson 2010 Nat Struct Mol Biol 10.1038/nsmb.1899"),
    dict(driver="RAD51C", mode="LOF", partner="PARP1", drug_axis="parp_inhibitors",
         falsification=None, cn_event=False, block="A",
         citation="Min 2013 Mol Cancer Ther 10.1158/1535-7163.MCT-13-0290"),
    dict(driver="ATM", mode="LOF", partner="ATR", drug_axis="atr_wee1",
         falsification=None, cn_event=False, block="A",
         citation="Reaper 2011 Nat Chem Biol 10.1038/nchembio.573"),
    dict(driver="ATM", mode="LOF", partner="PARP1", drug_axis="parp_inhibitors",
         falsification=None, cn_event=False, block="A",
         citation="Weston 2010 Blood 10.1182/blood-2010-05-284984"),
    dict(driver="CCNE1", mode="ACT", partner="WEE1", drug_axis="atr_wee1",
         falsification=None, cn_event=True, block="A",
         citation="Kok 2020 Oncogene; Chen 2018 Cancer Cell (CCNE1-amp WEE1/ATR)"),
    dict(driver="CCNE1", mode="ACT", partner="ATR", drug_axis="atr_wee1",
         falsification=None, cn_event=True, block="A",
         citation="Toledo 2011 Nat Struct Mol Biol 10.1038/nsmb.2189"),
    dict(driver="CCNE1", mode="ACT", partner="PKMYT1", drug_axis="pkmyt1",
         falsification=None, cn_event=True, block="A",
         citation="Gallo 2022 Nature 10.1038/s41586-022-04638-9"),
    dict(driver="TP53", mode="LOF", partner="WEE1", drug_axis="atr_wee1",
         falsification=None, cn_event=False, block="A",
         citation="Hirai 2009 Mol Cancer Ther 10.1158/1535-7163.MCT-09-0463"),
    dict(driver="TP53", mode="LOF", partner="PKMYT1", drug_axis="pkmyt1",
         falsification=None, cn_event=False, block="A",
         citation="Gallo 2022 Nature 10.1038/s41586-022-04638-9"),
    dict(driver="TP53", mode="LOF", partner="CHEK1", drug_axis="atr_wee1",
         falsification=None, cn_event=False, block="A",
         citation="Ma 2012 Clin Cancer Res (CHK1i in p53-deficient)"),

    # ---- BLOCK B: chromatin remodeler paralogs + metabolic ----
    dict(driver="ARID1A", mode="LOF", partner="ARID1B", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Helming 2014 Nat Med 10.1038/nm.3480"),
    dict(driver="ARID1A", mode="LOF", partner="EZH2", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Bitler 2015 Nat Med 10.1038/nm.3799"),
    dict(driver="SMARCA4", mode="LOF", partner="SMARCA2", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Hoffman 2014 PNAS 10.1073/pnas.1316793111; Oike 2013 Cancer Res"),
    dict(driver="SMARCA4", mode="LOF", partner="EZH2", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Kim 2015 Nat Med 10.1038/nm.3968"),
    dict(driver="MTAP", mode="LOF", partner="PRMT5", drug_axis=None,
         falsification=None, cn_event=True, block="B",
         citation="Mavrakis 2016 Science 10.1126/science.aad5944; Kryukov 2016 Science 10.1126/science.aad5214"),
    dict(driver="MTAP", mode="LOF", partner="MAT2A", drug_axis=None,
         falsification=None, cn_event=True, block="B",
         citation="Marjon 2016 Cell Rep 10.1016/j.celrep.2016.09.004"),
    dict(driver="KEAP1", mode="LOF", partner="GLS", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Romero 2017 Nat Med 10.1038/nm.4407 (KEAP1/NFE2L2 -> glutaminolysis)"),
    dict(driver="STK11", mode="LOF", partner="MTOR", drug_axis=None,
         falsification=None, cn_event=False, block="B",
         citation="Shackelford 2013 Cancer Cell 10.1016/j.ccr.2012.12.008"),

    # ---- BLOCK C: oncogene-context + tumor-suppressor targetable node ----
    dict(driver="KRAS", mode="ACT", partner="PTPN11", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Mainardi 2018 Nat Med 10.1038/s41591-018-0023-9 (SHP2/PTPN11)"),
    dict(driver="KRAS", mode="ACT", partner="STK33", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Scholl 2009 Cell 10.1016/j.cell.2009.03.017"),
    dict(driver="RB1", mode="LOF", partner="CDK4", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Fry 2004 Mol Cancer Ther (RB status & CDK4/6i)"),
    dict(driver="RB1", mode="LOF", partner="SKP2", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Zhao 2013 (RB1-loss SKP2 dependency)"),
    dict(driver="PTEN", mode="LOF", partner="PIK3CB", drug_axis="pi3k_beta",
         falsification=["taselisib", "alpelisib"], cn_event=False, block="C",
         citation="Wee 2008 PNAS 10.1073/pnas.0802655105; Jia 2008 Nature 10.1038/nature07091"),
    dict(driver="NF1", mode="LOF", partner="MAP2K1", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Jessen 2013 (NF1-loss MEK dependency)"),
    dict(driver="VHL", mode="LOF", partner="EPAS1", drug_axis=None,
         falsification=None, cn_event=False, block="C",
         citation="Kaelin 2008 Nat Rev Cancer (VHL/HIF2A/EPAS1)"),
]

# Novel-discovery anchors (literature drivers whose NEW partners we will discover)
NOVEL_ANCHORS = [
    dict(driver="PTEN", mode="LOF"), dict(driver="ARID1A", mode="LOF"),
    dict(driver="KEAP1", mode="LOF"), dict(driver="STK11", mode="LOF"),
    dict(driver="SMARCA4", mode="LOF"), dict(driver="VHL", mode="LOF"),
    dict(driver="RB1", mode="LOF"), dict(driver="KRAS", mode="ACT"),
    dict(driver="NF1", mode="LOF"), dict(driver="TP53", mode="LOF"),
]

if __name__ == "__main__":
    from collections import Counter
    print(f"KNOWN_SL_PANEL: {len(KNOWN_SL_PANEL)} pairs")
    print("by block:", dict(Counter(e["block"] for e in KNOWN_SL_PANEL)))
    print("with drug axis:", sum(1 for e in KNOWN_SL_PANEL if e["drug_axis"]))
    print("CN events:", [e["driver"] for e in KNOWN_SL_PANEL if e["cn_event"]])
    print(f"NOVEL_ANCHORS: {len(NOVEL_ANCHORS)}")
