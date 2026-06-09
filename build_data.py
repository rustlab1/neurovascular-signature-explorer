#!/usr/bin/env python3
"""
build_data.py -- bake the cross-disease neurovascular single-cell DEG analysis
into a single static data file (data.js) consumed by index.html.

Mirrors the bbb-proteome-explorer pattern: preprocessing happens here, the web
app loads pre-computed constants and renders client-side (no server, no build).

Source data: scRNAseq_human/Data/Analysis (MAST DEGs, GO enrichment, Venn
overlaps) + the AD/HD/DM (=FTD-GRN) sheets of DGE_allClusters.xlsx.

Disease code mapping: the data uses "DM" for the FTD-GRN cohort; the paper and
this explorer label it FTD-GRN. AD = Alzheimer's, HD = Huntington's.
"""
import os, json, math, csv, glob

# Point NV_DATA_DIR at the source analysis directory (scRNAseq_human/Data).
DATA = os.environ.get("NV_DATA_DIR", "source_data")
OUT  = os.path.dirname(os.path.abspath(__file__))  # write data.js next to this script

DGE   = os.path.join(DATA, "Analysis/2_MAST_DGE_all/DGE_allClusters.xlsx")
ENDO  = os.path.join(DATA, "Analysis/3_Endothelium/3_MAST_DGE/DEG_endothel.xlsx")
GO_DIR   = os.path.join(DATA, "Analysis/2_MAST_DGE_all/GO_outputs")
VENN_DIR = os.path.join(DATA, "Analysis/2_MAST_DGE_all/Venn_outputs_filtered")
PERI_DIR = os.path.join(DATA, "Analysis/4_Pericytes")

# sheet code (in source) -> app disease code
SHEET2APP = {"AD": "AD", "HD": "HD", "DM": "FTD"}
DISEASES = [
    {"code": "AD",  "label": "Alzheimer's (AD)",   "short": "AD",  "region": "Cortex",      "color": "#d97706"},
    {"code": "HD",  "label": "Huntington's (HD)",  "short": "HD",  "region": "Striatum",    "color": "#2563eb"},
    {"code": "FTD", "label": "FTD-GRN",            "short": "FTD", "region": "Frontal/temporal cortex", "color": "#9333ea"},
]
DISEASE_ORDER = ["AD", "HD", "FTD"]

# lineage display order + labels + functional group
LIN = [
    ("Endothelial", "Endothelial",            "Vascular"),
    ("Pericyte",    "Pericyte",               "Vascular"),
    ("vSMC",        "Smooth muscle (vSMC)",   "Vascular"),
    ("M_Fibro",     "Meningeal fibroblast",   "Vascular"),
    ("P_Fibro",     "Perivascular fibroblast","Vascular"),
    ("P_Mac",       "Perivascular macrophage","Immune"),
    ("Microglia",   "Microglia",              "Immune"),
    ("Astro",       "Astrocyte",              "Glia"),
    ("Oligo",       "Oligodendrocyte",        "Glia"),
    ("OPC",         "OPC",                    "Glia"),
    ("Neuron_exc",  "Excitatory neuron",      "Neuron"),
    ("Neuron_inh",  "Inhibitory neuron",      "Neuron"),
]
LIN_CODES  = [x[0] for x in LIN]
LIN_LABEL  = {x[0]: x[1] for x in LIN}
LIN_GROUP  = {x[0]: x[2] for x in LIN}
LIN_IDX    = {x[0]: i for i, x in enumerate(LIN)}

def r(x, n=3):
    """round, return None for NaN/inf"""
    try:
        if x is None: return None
        x = float(x)
        if math.isnan(x) or math.isinf(x): return None
        return round(x, n)
    except Exception:
        return None

def sig3(x):
    """compact significant-figure representation for p-values"""
    try:
        if x is None: return None
        x = float(x)
        if math.isnan(x): return None
        if x <= 0: return 0.0
        return float(f"{x:.3g}")
    except Exception:
        return None

# ---------------------------------------------------------------- core DEGs
import pandas as pd

def build_core():
    rows = []           # [gene, dIdx, lIdx, fc, p1, p2, q]
    genes = set()
    stats = {}
    xl = pd.ExcelFile(DGE)
    for sheet in ["AD", "HD", "DM"]:
        app = SHEET2APP[sheet]
        dIdx = DISEASE_ORDER.index(app)
        df = pd.read_excel(DGE, sheet_name=sheet)
        kept = 0
        for _, row in df.iterrows():
            lin = str(row.get("lineage", "")).strip()
            if lin not in LIN_IDX:
                continue
            gene = str(row.get("gene", "")).strip()
            if not gene or gene == "nan":
                continue
            fc = row.get("avg_log2FC")
            q  = row.get("p_adj_fdr")
            try:
                fcf = float(fc); qf = float(q)
            except Exception:
                continue
            if math.isnan(fcf):
                continue
            # keep significant OR strong-effect rows (trims NS low-effect noise)
            if not ((not math.isnan(qf) and qf <= 0.05) or abs(fcf) >= 0.5):
                continue
            rows.append([
                gene, dIdx, LIN_IDX[lin],
                r(fcf, 3), r(row.get("pct.1"), 3), r(row.get("pct.2"), 3), sig3(qf)
            ])
            genes.add(gene)
            kept += 1
        stats[app] = {"raw": int(df.shape[0]), "kept": kept}
    return rows, sorted(genes), stats

# ---------------------------------------------------------- common (shared) genes
def build_common():
    df = pd.read_excel(DGE, sheet_name="Common Genes")
    out = []
    for _, row in df.iterrows():
        lin = str(row.get("lineage", "")).strip()
        if lin not in LIN_IDX:
            continue
        def split(v):
            if v is None or (isinstance(v, float) and math.isnan(v)): return []
            return [g.strip() for g in str(v).split(",") if g.strip() and g.strip() != "nan"]
        up = split(row.get("common_genes_up"))
        dn = split(row.get("common_genes_down"))
        out.append({
            "lin": LIN_IDX[lin],
            "nUp": int(row.get("n_common_up") or 0),
            "nDown": int(row.get("n_common_down") or 0),
            "up": up[:80],
            "down": dn[:80],
        })
    out.sort(key=lambda d: -(d["nUp"] + d["nDown"]))
    return out

# ----------------------------------------------------------- endothelial zonation
def build_endo():
    rows = []
    segset = []
    try:
        xl = pd.ExcelFile(ENDO)
    except Exception as e:
        print("  endo: skip", e); return [], []
    sheet2app = {"AD_endo": "AD", "HD_endo": "HD", "DM_endo": "FTD"}
    # collect segments first
    for sh in xl.sheet_names:
        df = pd.read_excel(ENDO, sheet_name=sh)
        for s in df.get("endo_segment", pd.Series([])).dropna().unique():
            if str(s) not in segset:
                segset.append(str(s))
    seg_order = [s for s in ["Arterial", "Capillary", "Venous"] if s in segset] + \
                [s for s in segset if s not in ["Arterial", "Capillary", "Venous"]]
    seg_idx = {s: i for i, s in enumerate(seg_order)}
    for sh in xl.sheet_names:
        app = sheet2app.get(sh)
        if not app: continue
        dIdx = DISEASE_ORDER.index(app)
        df = pd.read_excel(ENDO, sheet_name=sh)
        for _, row in df.iterrows():
            seg = str(row.get("endo_segment", "")).strip()
            if seg not in seg_idx: continue
            gene = str(row.get("gene", "")).strip()
            if not gene or gene == "nan": continue
            try:
                fcf = float(row.get("avg_log2FC")); qf = float(row.get("p_adj_fdr"))
            except Exception:
                continue
            if math.isnan(fcf): continue
            if not ((not math.isnan(qf) and qf <= 0.05) or abs(fcf) >= 0.5): continue
            rows.append([gene, dIdx, seg_idx[seg], r(fcf,3), r(row.get("pct.1"),3), r(row.get("pct.2"),3), sig3(qf)])
    return rows, seg_order

# ------------------------------------------------------------------- pathways (GO)
def read_go(path, topn=12):
    if not os.path.exists(path): return []
    try:
        df = pd.read_csv(path)
    except Exception:
        return []
    if df.empty or "p.adjust" not in df.columns: return []
    df = df.sort_values("p.adjust").head(topn)
    out = []
    for _, row in df.iterrows():
        genes = str(row.get("geneID", "")).split("/")
        out.append({
            "id": str(row.get("ID","")),
            "desc": str(row.get("Description","")),
            "fold": r(row.get("FoldEnrichment"), 2),
            "q": sig3(row.get("p.adjust")),
            "count": int(row.get("Count") or 0),
            "genes": [g for g in genes if g][:10],
        })
    return out

def build_pathways():
    # NOTE: in GO_DIR the _All_/_Endothelial_/_Pericyte_ files are byte-identical
    # duplicates (only the disease varies), so we use _All_ as the disease-level
    # enrichment and pull genuine pericyte-specific enrichment from 4_Pericytes/GO.
    out = {}
    onts = ["BP", "CC", "MF"]
    for sheet, app in SHEET2APP.items():
        out[app] = {"All": {}, "Pericyte": {}}
        for ont in onts:
            out[app]["All"][ont] = read_go(os.path.join(GO_DIR, f"{sheet}_All_{ont}_GO.csv"))
            out[app]["Pericyte"][ont] = read_go(os.path.join(
                PERI_DIR, "GO", f"{sheet}_Pericyte_ANY_{ont}_lfc0p5_fdr0p05_pcteither0p1_GO.csv"))
    return out

# ----------------------------------------------------------------- venn overlaps
def count_genes(path):
    if not os.path.exists(path): return 0
    n = 0
    with open(path, newline="") as f:
        rd = csv.reader(f)
        for i, line in enumerate(rd):
            if not line: continue
            v = line[0].strip()
            if i == 0 and v.lower() in ("gene", "x", ""):  # header
                continue
            if v:
                n += 1
    return n

def venn_for(group):
    base = f"Venn_{group}_ANY_lfc1_fdr0p05_pcteither0p1_"
    def c(suffix): return count_genes(os.path.join(VENN_DIR, base + suffix + ".csv"))
    AD, DM, HD = c("AD"), c("DM"), c("HD")
    AD_DM, AD_HD, DM_HD = c("AD_DM_intersection"), c("AD_HD_intersection"), c("DM_HD_intersection")
    trip = c("AD_DM_HD_intersection")
    # disjoint regions (FTD == DM)
    return {
        "AD_only":  max(AD - AD_DM - AD_HD + trip, 0),
        "HD_only":  max(HD - AD_HD - DM_HD + trip, 0),
        "FTD_only": max(DM - AD_DM - DM_HD + trip, 0),
        "AD_HD":    max(AD_HD - trip, 0),
        "AD_FTD":   max(AD_DM - trip, 0),
        "HD_FTD":   max(DM_HD - trip, 0),
        "all3":     trip,
        "totals":   {"AD": AD, "HD": HD, "FTD": DM},
    }

def build_venn():
    return {g: venn_for(g) for g in ["All", "Endothelial", "Pericyte"]}

# -------------------------------------------------------------- pericyte focus
def read_gene_list(path):
    if not os.path.exists(path): return []
    out = []
    with open(path, newline="") as f:
        for i, line in enumerate(csv.reader(f)):
            if not line: continue
            v = line[0].strip()
            if i == 0 and v.lower() in ("gene", "x"): continue
            if v: out.append(v)
    return out

def read_overlap_counts(path):
    d = {}
    if not os.path.exists(path): return d
    with open(path, newline="") as f:
        for i, line in enumerate(csv.reader(f)):
            if i == 0 or len(line) < 2: continue
            key = line[0].replace("DM", "FTD")
            try: d[key] = int(line[1])
            except Exception: pass
    return d

def build_pericyte():
    convergent = read_gene_list(os.path.join(PERI_DIR, "DEGs/Common_Pericyte_DEGs_AD_DM_HD.csv"))
    overlap = read_overlap_counts(os.path.join(PERI_DIR, "DEGs/Pericyte_DEG_overlap_counts.csv"))
    # M-peri / T-peri subtype proportions (matrix vs transport pericytes), from the paper.
    # Selective depletion of matrix-associated M-peri is the convergent finding.
    subtype = {
        "AD":  {"M_ctrl": 27.0, "M_dis": 6.7},
        "FTD": {"M_ctrl": 25.9, "M_dis": 9.8},
        "HD":  {"M_ctrl": 54.9, "M_dis": 33.5},
    }
    markers = {
        "M": ["COL4A1", "COL4A2", "ADAMTS1", "ADAMTS9", "LAMA4"],
        "T": ["SLC6A1", "SLC6A13", "APOD", "SLC20A2"],
    }
    return {
        "convergent": convergent,
        "overlap": overlap,
        "subtype": subtype,
        "markers": markers,
    }

# ----------------------------------------- ligand-receptor (curated from paper)
# No CellChat interaction table is exported in the source data; these pairs are
# curated from the paper's pericyte<->endothelial communication analysis.
# change codes: 1 = enhanced/up in disease, -1 = reduced/down, 0 = not prominent.
def build_lr():
    pairs = [
        # pathway, sender, receiver, {AD,HD,FTD}, shared, note
        ("LAMININ",  "Pericyte",    "Endothelial", {"AD":-1,"HD":-1,"FTD":-1}, True,  "BBB basement membrane; reduced across all three diseases"),
        ("COLLAGEN", "Pericyte",    "Endothelial", {"AD":-1,"HD":-1,"FTD":-1}, True,  "Matrix support from M-pericytes; broadly downregulated"),
        ("FN1",      "Pericyte",    "Endothelial", {"AD":-1,"HD":-1,"FTD":-1}, True,  "Fibronectin matrix signaling; shared reduction"),
        ("NCAM",     "Endothelial", "Pericyte",    {"AD":-1,"HD":-1,"FTD":-1}, True,  "Cell-adhesion signaling; reduced across diseases"),
        ("NOTCH",    "Endothelial", "Pericyte",    {"AD":-1,"HD":0,"FTD":-1},  False, "Mural-cell maintenance; attenuated mainly in AD/FTD"),
        ("VEGF",     "Pericyte",    "Endothelial", {"AD":-1,"HD":-1,"FTD":0},  False, "Angiogenic support; most reduced in HD"),
        ("ANGPTL",   "Pericyte",    "Endothelial", {"AD":-1,"HD":0,"FTD":-1},  False, "Vascular stability; reduced in AD/FTD"),
        ("JAM",      "Endothelial", "Endothelial", {"AD":-1,"HD":-1,"FTD":-1}, True,  "Junctional adhesion; barrier integrity"),
        ("PECAM1",   "Endothelial", "Endothelial", {"AD":-1,"HD":0,"FTD":-1},  False, "Endothelial junctions; reduced in AD/FTD"),
        ("BMP",      "Endothelial", "Pericyte",    {"AD":0,"HD":0,"FTD":-1},   False, "Trophic signaling; reduced in FTD-GRN"),
        ("VISFATIN", "Pericyte",    "Endothelial", {"AD":1,"HD":1,"FTD":1},    True,  "NAMPT metabolic signaling; broadly enhanced"),
        ("EPHB",     "Endothelial", "Pericyte",    {"AD":0,"HD":-1,"FTD":-1},  False, "Contact guidance; reduced in HD/FTD"),
        ("NECTIN",   "Endothelial", "Pericyte",    {"AD":1,"HD":-1,"FTD":0},   False, "Adhesion; heterogeneous, attenuated in HD"),
        ("NRG",      "Endothelial", "Pericyte",    {"AD":1,"HD":0,"FTD":0},    False, "Trophic; enhanced pericyte->endothelial routing in AD"),
    ]
    out = []
    for pw, snd, rcv, chg, shared, note in pairs:
        out.append({"pathway": pw, "sender": snd, "receiver": rcv,
                    "chg": chg, "shared": shared, "note": note})
    return out

# ----------------------------------------------------------------------- assemble
def main():
    print("Reading core DEGs ...")
    deg, genes, stats = build_core()
    print("  kept rows:", len(deg), "unique genes:", len(genes), stats)
    print("Reading common (shared) genes ...")
    common = build_common()
    print("  lineages with shared genes:", sum(1 for c in common if c["nUp"]+c["nDown"]>0))
    print("Reading endothelial zonation ...")
    endo, segs = build_endo()
    print("  endo rows:", len(endo), "segments:", segs)
    print("Reading GO pathways ...")
    pathways = build_pathways()
    print("  pathway terms (AD/All/BP):", len(pathways.get("AD",{}).get("All",{}).get("BP",[])),
          " (AD/Pericyte/BP):", len(pathways.get("AD",{}).get("Pericyte",{}).get("BP",[])))
    print("Reading pericyte focus ...")
    peri = build_pericyte()
    print("  convergent pericyte genes:", len(peri["convergent"]), "overlap:", peri["overlap"])
    lr = build_lr()

    data = {
        "meta": {
            "title": "Neurovascular Molecular Signatures in Neurodegeneration",
            "nDeg": len(deg), "nGenes": len(genes),
            "stats": stats,
            "note": "Cross-disease single-cell atlas of brain vascular & perivascular cells (AD, HD, FTD-GRN).",
        },
        "diseases": DISEASES,
        "diseaseOrder": DISEASE_ORDER,
        "lineages": [{"code": c, "label": LIN_LABEL[c], "group": LIN_GROUP[c]} for c in LIN_CODES],
        "genes": genes,
        "deg": deg,
        "common": common,
        "endo": endo, "segs": segs,
        "pathways": pathways,
        "pericyte": peri,
        "lr": lr,
    }

    os.makedirs(OUT, exist_ok=True)
    out_path = os.path.join(OUT, "data.js")
    payload = json.dumps(data, separators=(",", ":"))
    with open(out_path, "w") as f:
        f.write("window.NV = " + payload + ";\n")
    size = os.path.getsize(out_path)
    print(f"\nWrote {out_path}  ({size/1e6:.2f} MB)")

if __name__ == "__main__":
    main()
