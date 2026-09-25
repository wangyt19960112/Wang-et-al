"""
RNA velocity analysis of the trunk organoid dataset projected
onto the PCW3 reference UMAP space.

Required input files:
    data/
    ├── expr_trunk_new.csv
    ├── meta_trunk_new.csv
    ├── trunk.loom
    ├── subset_pcw3_umap.csv
    └── trunk_umap_new.csv

Output:
    results/rna_velocity/Trunk_RNA_velocity_PCW3.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_trunk_new.csv"
METADATA_FILE = DATA_DIR / "meta_trunk_new.csv"
LOOM_FILE = DATA_DIR / "trunk.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
TRUNK_UMAP_FILE = DATA_DIR / "trunk_umap_new.csv"


# ============================================================
# 2. Construct trunk AnnData object
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


trunk_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)


trunk = ad.AnnData(
    X=X,
    obs=trunk_obs
)


# Remove "-1" suffix from cell barcodes
trunk.obs_names = (
    trunk.obs_names.str.replace(
        r"-1$",
        "",
        regex=True
    )
)

trunk.obs_names_make_unique()


print(
    f"Trunk AnnData: "
    f"{trunk.n_obs} cells × "
    f"{trunk.n_vars} genes"
)


# ============================================================
# 3. Load spliced/unspliced counts
# ============================================================

trunk_loom = sc.read_loom(
    LOOM_FILE
)


trunk_loom.var_names_make_unique()


# Standardize loom cell barcodes
trunk_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in trunk_loom.obs_names
]


trunk_loom.obs_names_make_unique()


print(
    f"Trunk loom: "
    f"{trunk_loom.n_obs} cells × "
    f"{trunk_loom.n_vars} genes"
)


# ============================================================
# 4. Match cells between expression and loom data
# ============================================================

common_barcodes = (
    trunk.obs_names.intersection(
        trunk_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the trunk expression matrix and loom file."
    )


trunk = trunk[
    common_barcodes
].copy()


trunk_loom = trunk_loom[
    common_barcodes
].copy()


# ============================================================
# 5. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    trunk.var_names.intersection(
        trunk_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the trunk expression matrix and loom file."
    )


trunk_sub = trunk[
    :,
    common_genes
].copy()


trunk_loom_sub = trunk_loom[
    :,
    common_genes
].copy()


trunk_sub.layers["spliced"] = (
    trunk_loom_sub.layers["spliced"].copy()
)

trunk_sub.layers["unspliced"] = (
    trunk_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final trunk object: "
    f"{trunk_sub.n_obs} cells × "
    f"{trunk_sub.n_vars} genes"
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
# 7. Load trunk coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_trunk = pd.read_csv(
    TRUNK_UMAP_FILE
)


ref_umap_trunk.index = (
    ref_umap_trunk["cell"]
)


ref_umap_trunk.index = (
    ref_umap_trunk.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


ref_umap_trunk = (
    ref_umap_trunk[
        ["UMAP_1", "UMAP_2"]
    ]
)


# Remove duplicated cell IDs
ref_umap_trunk = (
    ref_umap_trunk[
        ~ref_umap_trunk.index.duplicated(
            keep="first"
        )
    ]
)


# ============================================================
# 8. Match cells with PCW3 reference UMAP coordinates
# ============================================================

common_barcodes = (
    trunk_sub.obs_names.intersection(
        ref_umap_trunk.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No trunk cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


trunk_sub = trunk_sub[
    common_barcodes
].copy()


ref_umap_trunk = (
    ref_umap_trunk.loc[
        common_barcodes
    ].copy()
)


# Add PCW3 reference UMAP coordinates
trunk_sub.obsm["X_umap"] = (
    ref_umap_trunk.to_numpy()
)


# ============================================================
# 9. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    trunk_sub
)


scv.pp.moments(
    trunk_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 10. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    trunk_sub,
    mode="stochastic"
)


scv.tl.velocity_graph(
    trunk_sub
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
    "Somite": "#9b59b6",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    trunk_sub,
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
    title="Trunk RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 12. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Trunk_RNA_velocity_PCW3.svg"
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