"""
RNA velocity analysis of the Bar2 (TK_San) axial organoid dataset
projected onto the PCW3 reference UMAP space.

Required input files:
    data/
    ├── expr_Bar2_new.csv
    ├── meta_Bar2_new.csv
    ├── Trunk_7uM.loom
    ├── subset_pcw3_umap.csv
    └── Bar2_umap_new.csv

Output:
    results/rna_velocity/Bar2_RNA_velocity_PCW3.svg
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

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)

EXPRESSION_FILE = DATA_DIR / "expr_Bar2_new.csv"
METADATA_FILE = DATA_DIR / "meta_Bar2_new.csv"

LOOM_FILE = DATA_DIR / "Trunk_7uM.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
BAR2_UMAP_FILE = DATA_DIR / "Bar2_umap_new.csv"


# ============================================================
# 2. Construct Bar2 AnnData object
# ============================================================

X = pd.read_csv(
    EXPRESSION_FILE,
    index_col=0
).T


# Standardize cell barcode format
X.index = X.index.str.replace(
    ".",
    "-",
    regex=False
)


Bar2_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)


Bar2 = ad.AnnData(
    X=X,
    obs=Bar2_obs
)


# Remove "-1" suffix from cell barcodes
Bar2.obs_names = (
    Bar2.obs_names.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


Bar2.obs_names_make_unique()


print(
    f"Bar2 AnnData: "
    f"{Bar2.n_obs} cells × "
    f"{Bar2.n_vars} genes"
)


# ============================================================
# 3. Load spliced/unspliced counts
# ============================================================

Bar2_loom = sc.read_loom(
    LOOM_FILE
)


Bar2_loom.var_names_make_unique()


# Standardize loom cell barcodes
Bar2_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in Bar2_loom.obs_names
]


Bar2_loom.obs_names_make_unique()


print(
    f"Bar2 loom: "
    f"{Bar2_loom.n_obs} cells × "
    f"{Bar2_loom.n_vars} genes"
)


# ============================================================
# 4. Match cells between expression and loom data
# ============================================================

common_barcodes = (
    Bar2.obs_names.intersection(
        Bar2_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Bar2 expression matrix and loom file."
    )


Bar2 = Bar2[
    common_barcodes
].copy()


Bar2_loom = Bar2_loom[
    common_barcodes
].copy()


# ============================================================
# 5. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    Bar2.var_names.intersection(
        Bar2_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Bar2 expression matrix and loom file."
    )


Bar2_sub = Bar2[
    :,
    common_genes
].copy()


Bar2_loom_sub = Bar2_loom[
    :,
    common_genes
].copy()


Bar2_sub.layers["spliced"] = (
    Bar2_loom_sub.layers["spliced"].copy()
)

Bar2_sub.layers["unspliced"] = (
    Bar2_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Bar2 object: "
    f"{Bar2_sub.n_obs} cells × "
    f"{Bar2_sub.n_vars} genes"
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
# 7. Load Bar2 coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_Bar2 = pd.read_csv(
    BAR2_UMAP_FILE
)


ref_umap_Bar2.index = (
    ref_umap_Bar2["cell"]
)


ref_umap_Bar2.index = (
    ref_umap_Bar2.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


ref_umap_Bar2 = (
    ref_umap_Bar2[
        ["UMAP_1", "UMAP_2"]
    ]
)


# Remove duplicated cell IDs
ref_umap_Bar2 = (
    ref_umap_Bar2[
        ~ref_umap_Bar2.index.duplicated(
            keep="first"
        )
    ]
)


# ============================================================
# 8. Match cells with PCW3 reference UMAP coordinates
# ============================================================

common_barcodes = (
    Bar2_sub.obs_names.intersection(
        ref_umap_Bar2.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No Bar2 cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


Bar2_sub = Bar2_sub[
    common_barcodes
].copy()


ref_umap_Bar2 = (
    ref_umap_Bar2.loc[
        common_barcodes
    ].copy()
)


# Add PCW3 reference UMAP coordinates
Bar2_sub.obsm["X_umap"] = (
    ref_umap_Bar2.to_numpy()
)


# ============================================================
# 9. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    Bar2_sub
)


scv.pp.moments(
    Bar2_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 10. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    Bar2_sub,
    mode="stochastic"
)


scv.tl.velocity_graph(
    Bar2_sub
)


# ============================================================
# 11. Plot RNA velocity in PCW3 reference UMAP space
# ============================================================

# Keep text editable in SVG
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
    "Somite": "#9b59b6",
    "Anterior PSM": "#F5C26B",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    Bar2_sub,
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
    title="TK_San RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 12. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Bar2_RNA_velocity_PCW3.svg"
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