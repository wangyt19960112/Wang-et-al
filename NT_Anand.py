"""
RNA velocity analysis of the HNT dataset projected onto
the PCW3 reference UMAP space.

Required input files:
    data/
    ├── expr_HNT_new.csv
    ├── meta_HNT_new.csv
    ├── HNT.loom
    ├── subset_pcw3_umap.csv
    └── HNT_umap_new.csv

Output:
    results/rna_velocity/HNT_new.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_HNT_new.csv"
METADATA_FILE = DATA_DIR / "meta_HNT_new.csv"
LOOM_FILE = DATA_DIR / "HNT.loom"
PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
HNT_UMAP_FILE = DATA_DIR / "HNT_umap_new.csv"


# ============================================================
# 2. Construct AnnData object from expression and metadata
# ============================================================

X = pd.read_csv(
    EXPRESSION_FILE,
    index_col=0
).T

# Match barcode format
X.index = X.index.str.replace(
    ".",
    "-",
    regex=False
)

HNT_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)

HNT = ad.AnnData(
    X=X,
    obs=HNT_obs
)

# Remove the "-1" suffix from cell barcodes
HNT.obs_names = HNT.obs_names.str.replace(
    r"-1$",
    "",
    regex=True
)

HNT.obs_names_make_unique()

print(
    f"AnnData object: "
    f"{HNT.n_obs} cells × {HNT.n_vars} genes"
)


# ============================================================
# 3. Load spliced/unspliced counts from loom file
# ============================================================

HNT_loom = sc.read_loom(
    LOOM_FILE
)

HNT_loom.var_names_make_unique()


# Convert loom cell names to barcode format used in HNT
HNT_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in HNT_loom.obs_names
]

HNT_loom.obs_names_make_unique()


# ============================================================
# 4. Match cells between expression data and loom file
# ============================================================

common_barcodes = HNT.obs_names.intersection(
    HNT_loom.obs_names
)

print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)

if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the expression matrix and loom file."
    )


HNT = HNT[
    common_barcodes
].copy()

HNT_loom = HNT_loom[
    common_barcodes
].copy()


# ============================================================
# 5. Match genes and add spliced/unspliced layers
# ============================================================

common_genes = HNT.var_names.intersection(
    HNT_loom.var_names
)

print(
    f"Matched genes: {len(common_genes)}"
)

if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the expression matrix and loom file."
    )


HNT_sub = HNT[
    :,
    common_genes
].copy()

HNT_loom_sub = HNT_loom[
    :,
    common_genes
].copy()


HNT_sub.layers["spliced"] = (
    HNT_loom_sub.layers["spliced"].copy()
)

HNT_sub.layers["unspliced"] = (
    HNT_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final HNT object: "
    f"{HNT_sub.n_obs} cells × {HNT_sub.n_vars} genes"
)


# ============================================================
# 6. Load PCW3 reference UMAP coordinates
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
# 7. Load HNT coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_HNT = pd.read_csv(
    HNT_UMAP_FILE
)

ref_umap_HNT.index = ref_umap_HNT["cell"]

ref_umap_HNT.index = (
    ref_umap_HNT.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)

ref_umap_HNT = ref_umap_HNT[
    ["UMAP_1", "UMAP_2"]
]


# Remove duplicated cell IDs if present
ref_umap_HNT = ref_umap_HNT[
    ~ref_umap_HNT.index.duplicated(
        keep="first"
    )
]


# Match HNT cells with reference UMAP coordinates
common_barcodes = (
    HNT_sub.obs_names.intersection(
        ref_umap_HNT.index
    )
)

print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)

if len(common_barcodes) == 0:
    raise ValueError(
        "No HNT cells could be matched to "
        "the reference UMAP coordinates."
    )


HNT_sub = HNT_sub[
    common_barcodes
].copy()

ref_umap_HNT = ref_umap_HNT.loc[
    common_barcodes
].copy()


# Add PCW3 reference UMAP coordinates
HNT_sub.obsm["X_umap"] = (
    ref_umap_HNT.to_numpy()
)


# ============================================================
# 8. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    HNT_sub
)

scv.pp.moments(
    HNT_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 9. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    HNT_sub,
    mode="stochastic"
)

scv.tl.velocity_graph(
    HNT_sub
)


# ============================================================
# 10. Plot velocity in PCW3 reference UMAP space
# ============================================================

plt.rcParams["svg.fonttype"] = "none"

fig, ax = plt.subplots(
    figsize=(8, 6)
)


# PCW3 reference cells shown as background
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
    "Forebrain-Midbrain": "#F5761A",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    HNT_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="predicted.celltype",
    palette=celltype_colors,
    density=3.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="HNT RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 11. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "HNT_RNA_velocity_PCW3.svg"
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