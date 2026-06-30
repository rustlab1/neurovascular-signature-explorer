# Neurovascular Molecular Signatures in Neurodegeneration

Interactive companion to a comparative single-cell / single-nucleus RNA-seq study of the
human brain neurovascular unit across three neurodegenerative diseases: Alzheimer's disease
(AD), Huntington's disease (HD) and GRN-related frontotemporal dementia (FTD-GRN).

The explorer shows disease-versus-control differential expression per cell type, the shared
cross-disease signature, the pericyte and endothelial focus of the study, GO pathway
enrichment, and pericyte-endothelial signalling.

## Sections

- **Cell atlas** - harmonized single-nucleus UMAP, downsampled per disease, coloured by cell type or condition.
- **Gene search** - any detected gene, two views: disease-vs-control log2 fold-change, or absolute mean expression in control and disease per cell type (confirm markers).
- **Cell-type DEGs** - ranked up/down genes for a chosen disease and cell type.
- **Cross-disease** - Venn overlap of significant DEGs plus the concordant genes shared by all three diseases.
- **Pericytes** - M-peri / T-peri proportion shift, the convergent pericyte gene set, and subtype markers.
- **Endothelium** - DEGs resolved by vascular segment (arterial / capillary / venous).
- **Pathways** - GO enrichment (BP / CC / MF) by disease.
- **Signaling** - curated pericyte-endothelial ligand-receptor changes from the study's CellChat analysis.

## Files

| File | Purpose |
|------|---------|
| `index.html` | Page structure and styles |
| `app.js` | Client logic and rendering |
| `data.js` | Pre-computed DEG / pathway / overlap data (`window.NV`) |
| `umap.js` | Downsampled atlas embedding (`window.NV_UMAP`) |
| `expr.js` | Mean expression per gene × cell type × condition (`window.NV_EXPR`) |
| `build_data.py` | Builds `data.js` from the source analysis tables |
| `build_umap.R` | Builds `umap.js` from the Seurat objects |
| `build_expr.R` | Builds `expr.js` from the Seurat objects |

## Rebuilding the data

`build_data.py` reads the MAST DEG tables, GO enrichment, and overlap files and writes `data.js`.
`build_umap.R` reads the per-disease Seurat objects and writes the downsampled atlas (`umap.js`).
Both take the source path from the `NV_DATA_DIR` environment variable (the `scRNAseq_human/Data` directory):

```bash
export NV_DATA_DIR=/path/to/scRNAseq_human/Data
python3 build_data.py     # needs pandas, openpyxl
Rscript build_umap.R      # needs Seurat
Rscript build_expr.R      # needs Seurat
```

| File | Generated from |
|------|----------------|
| `data.js` | `2_MAST_DGE_all/DGE_allClusters.xlsx`, `3_Endothelium`, `4_Pericytes`, GO outputs |
| `umap.js` | `1_Obj_Clusterin_QC/Objects/*_reclustered.rds` (5,000 nuclei/disease) |
| `expr.js` | `1_Obj_Clusterin_QC/Objects/*_reclustered.rds` (mean expression, detected genes) |

## Data notes

- The source data labels the FTD-GRN cohort `DM`; it is shown here as **FTD-GRN**.
- Differential expression is per cell type, disease vs. matched control, computed with MAST.
  A positive log2 fold-change means higher expression in disease.
- The signalling table is curated from the study's CellChat findings; the raw interaction
  tables can be added to make it fully data-driven.

## Deployment

Static site. The included GitHub Actions workflow publishes the repository root to GitHub Pages
on push to `main`.
