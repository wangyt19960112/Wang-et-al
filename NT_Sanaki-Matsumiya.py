"""
RNA velocity analysis of the Bar1 dataset projected onto
the PCW3 reference UMAP space.

Required input files:
    data/
    ├── expr_Bar1_new.csv
    ├── meta_Bar1_new.csv
    ├── NT_5uM.loom
    ├── subset_pcw3_umap.csv
    └── Bar1_umap_new.csv

Output:
    results/rna_velocity/Bar1_RNA_velocity_PCW3.svg
"""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvelo as scv


# ============================================================
# 1. File paths
# ============================================================

DATA_DIR = Path("data")
RESULTS_DIR = Path("results") / "rna_velocity"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

EXPRESSION_FILE = DATA_DIR / "expr_Bar1_new.csv"
METADATA_FILE = DATA_DIR / "meta_Bar1_new.csv"
LOOM_FILE = DATA_DIR / "NT_5uM.loom"
PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
BAR1_UMAP_FILE = DATA_DIR / "Bar1_umap_new.csv"


# ============================================================
# 2. Construct AnnData object
# ============================================================

X = pd.read_csv(
    EXPRESSION_FILE,
    index_col=0
).T

# Standardize barcode format
X.index = X.index.str.replace(
    ".",
    "-",
    regex=False
)

Bar1_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)

Bar1 = ad.AnnData(
    X=X,
    obs=Bar1_obs
)

# Remove "-1" suffix from cell barcodes
Bar1.obs_names = Bar1.obs_names.str.replace(
    r"-1$",
    "",
    regex=True
)

Bar1.obs_names_make_unique()

print(
    f"Bar1 AnnData: "
    f"{Bar1.n_obs} cells × {Bar1.n_vars} genes"
)


# ============================================================
# 3. Load spliced/unspliced counts
# ============================================================

Bar1_loom = sc.read_loom(
    LOOM_FILE
)

Bar1_loom.var_names_make_unique()


# Convert loom cell names to the barcode format used by Bar1
Bar1_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in Bar1_loom.obs_names
]

Bar1_loom.obs_names_make_unique()


# ============================================================
# 4. Match cells between expression data and loom file
# ============================================================

common_barcodes = Bar1.obs_names.intersection(
    Bar1_loom.obs_names
)

print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)

if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Bar1 expression matrix and loom file."
    )


Bar1 = Bar1[
    common_barcodes
].copy()

Bar1_loom = Bar1_loom[
    common_barcodes
].copy()


# ============================================================
# 5. Match genes and add RNA velocity layers
# ============================================================

common_genes = Bar1.var_names.intersection(
    Bar1_loom.var_names
)

print(
    f"Matched genes: {len(common_genes)}"
)

if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Bar1 expression matrix and loom file."
    )


Bar1_sub = Bar1[
    :,
    common_genes
].copy()

Bar1_loom_sub = Bar1_loom[
    :,
    common_genes
].copy()


Bar1_sub.layers["spliced"] = (
    Bar1_loom_sub.layers["spliced"].copy()
)

Bar1_sub.layers["unspliced"] = (
    Bar1_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Bar1 object: "
    f"{Bar1_sub.n_obs} cells × {Bar1_sub.n_vars} genes"
)


# ============================================================
# 6. Load PCW3 reference UMAP
# ============================================================

pcw3 = pd.read_csv(
    PCW3_UMAP_FILE
)

pcw3.index = pcw3["cell"]

pcw3.index = pcw3.index.str.replace(
    r"-1$",
    "",
    regex=True
)

pcw3 = pcw3[
    ["UMAP_1", "UMAP_2"]
]


# ============================================================
# 7. Load Bar1 coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_Bar1 = pd.read_csv(
    BAR1_UMAP_FILE
)

ref_umap_Bar1.index = ref_umap_Bar1["cell"]

ref_umap_Bar1.index = (
    ref_umap_Bar1.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)

ref_umap_Bar1 = ref_umap_Bar1[
    ["UMAP_1", "UMAP_2"]
]


# Remove duplicated cell IDs
ref_umap_Bar1 = ref_umap_Bar1[
    ~ref_umap_Bar1.index.duplicated(
        keep="first"
    )
]


# Match Bar1 cells with reference UMAP coordinates
common_barcodes = (
    Bar1_sub.obs_names.intersection(
        ref_umap_Bar1.index
    )
)

print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)

if len(common_barcodes) == 0:
    raise ValueError(
        "No Bar1 cells could be matched to "
        "the reference UMAP coordinates."
    )


Bar1_sub = Bar1_sub[
    common_barcodes
].copy()

ref_umap_Bar1 = ref_umap_Bar1.loc[
    common_barcodes
].copy()


# Add PCW3 reference UMAP coordinates
Bar1_sub.obsm["X_umap"] = (
    ref_umap_Bar1.to_numpy()
)


# ============================================================
# 8. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    Bar1_sub
)

scv.pp.moments(
    Bar1_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 9. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    Bar1_sub,
    mode="stochastic"
)

scv.tl.velocity_graph(
    Bar1_sub
)


# ============================================================
# 10. Plot RNA velocity in PCW3 reference UMAP space
# ============================================================

# Keep text editable in SVG
plt.rcParams["svg.fonttype"] = "none"

fig, ax = plt.subplots(
    figsize=(8, 6)
)


# PCW3 reference cells shown in the background
ax.scatter(
    pcw3["UMAP_1"],
    pcw3["UMAP_2"],
    color="lightgrey",
    alpha=0.3,
    s=5,
    label="PCW3 background"
)


# Cell-type colors
celltype_colors = {
    "NMP": "#e74c3c",
    "Cervical-Thoracic NT": "#27ae60",
    "Hindbrain": "#3498db",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    Bar1_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="predicted.celltype",
    palette=celltype_colors,
    density=2.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="Bar1 RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 11. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Bar1_RNA_velocity_PCW3.svg"
)

fig.savefig(
    OUTPUT_FILE,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


print(
    f"RNA velocity analysis completed.\n"
    f"Figure saved to: {OUTPUT_FILE}"
)