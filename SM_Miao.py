"""
RNA velocity analysis of the Segmentoid dataset projected onto
the PCW3 reference UMAP space.

The Segmentoid dataset contains four time points:
24 h, 48 h, 72 h, and 98 h.

NMP cells are further separated by time point for visualization.

Required input files:
    data/
    ├── expr_miao_new.csv
    ├── meta_miao_new.csv
    ├── Segmentoid_24h.loom
    ├── Segmentoid_48h.loom
    ├── Segmentoid_72h.loom
    ├── Segmentoid_98h.loom
    ├── subset_pcw3_umap.csv
    └── miao_umap_new.csv

Output:
    results/rna_velocity/Segmentoid_RNA_velocity_PCW3.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_miao_new.csv"
METADATA_FILE = DATA_DIR / "meta_miao_new.csv"

LOOM_24H_FILE = DATA_DIR / "Segmentoid_24h.loom"
LOOM_48H_FILE = DATA_DIR / "Segmentoid_48h.loom"
LOOM_72H_FILE = DATA_DIR / "Segmentoid_72h.loom"
LOOM_98H_FILE = DATA_DIR / "Segmentoid_98h.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
SEGMENTOID_UMAP_FILE = DATA_DIR / "miao_umap_new.csv"


# ============================================================
# 2. Construct Segmentoid AnnData object
# ============================================================

X = pd.read_csv(
    EXPRESSION_FILE,
    index_col=0
).T


# Extract cell barcode after the final "."
X.index = X.index.str.split(".").str[-1]


Segmentoid_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)


# Extract cell barcode after ":"
Segmentoid_obs.index = (
    Segmentoid_obs.index
    .str.split(":")
    .str[-1]
)


Segmentoid = ad.AnnData(
    X=X,
    obs=Segmentoid_obs
)


Segmentoid.obs_names_make_unique()


print(
    f"Segmentoid AnnData: "
    f"{Segmentoid.n_obs} cells × "
    f"{Segmentoid.n_vars} genes"
)


# ============================================================
# 3. Load loom files
# ============================================================

Segmentoid_24h_loom = sc.read_loom(
    LOOM_24H_FILE
)

Segmentoid_48h_loom = sc.read_loom(
    LOOM_48H_FILE
)

Segmentoid_72h_loom = sc.read_loom(
    LOOM_72H_FILE
)

Segmentoid_98h_loom = sc.read_loom(
    LOOM_98H_FILE
)


# Ensure unique gene names
for loom_obj in [
    Segmentoid_24h_loom,
    Segmentoid_48h_loom,
    Segmentoid_72h_loom,
    Segmentoid_98h_loom,
]:
    loom_obj.var_names_make_unique()


# ============================================================
# 4. Combine loom files
# ============================================================

Segmentoid_loom = ad.concat(
    [
        Segmentoid_24h_loom,
        Segmentoid_48h_loom,
        Segmentoid_72h_loom,
        Segmentoid_98h_loom,
    ],
    join="outer",
    label="time",
    keys=[
        "24h",
        "48h",
        "72h",
        "98h",
    ]
)


# Standardize loom cell barcodes
Segmentoid_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in Segmentoid_loom.obs_names
]


Segmentoid_loom.obs_names_make_unique()


print(
    f"Combined loom object: "
    f"{Segmentoid_loom.n_obs} cells × "
    f"{Segmentoid_loom.n_vars} genes"
)


# ============================================================
# 5. Match cells between expression and loom data
# ============================================================

common_barcodes = (
    Segmentoid.obs_names.intersection(
        Segmentoid_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Segmentoid expression matrix and loom files."
    )


Segmentoid = Segmentoid[
    common_barcodes
].copy()

Segmentoid_loom = Segmentoid_loom[
    common_barcodes
].copy()


# ============================================================
# 6. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    Segmentoid.var_names.intersection(
        Segmentoid_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Segmentoid expression matrix and loom files."
    )


Segmentoid_sub = Segmentoid[
    :,
    common_genes
].copy()

Segmentoid_loom_sub = Segmentoid_loom[
    :,
    common_genes
].copy()


Segmentoid_sub.layers["spliced"] = (
    Segmentoid_loom_sub.layers["spliced"].copy()
)

Segmentoid_sub.layers["unspliced"] = (
    Segmentoid_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Segmentoid object: "
    f"{Segmentoid_sub.n_obs} cells × "
    f"{Segmentoid_sub.n_vars} genes"
)


# ============================================================
# 7. Load PCW3 reference UMAP coordinates
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
# 8. Load Segmentoid coordinates in PCW3 UMAP space
# ============================================================

ref_umap_Segmentoid = pd.read_csv(
    SEGMENTOID_UMAP_FILE
)


ref_umap_Segmentoid.index = (
    ref_umap_Segmentoid["cell"]
    .str.split(":")
    .str[-1]
)


ref_umap_Segmentoid = (
    ref_umap_Segmentoid[
        ["UMAP_1", "UMAP_2"]
    ]
)


# Remove duplicated cell IDs
ref_umap_Segmentoid = (
    ref_umap_Segmentoid[
        ~ref_umap_Segmentoid.index.duplicated(
            keep="first"
        )
    ]
)


# ============================================================
# 9. Match cells with PCW3 reference UMAP coordinates
# ============================================================

common_barcodes = (
    Segmentoid_sub.obs_names.intersection(
        ref_umap_Segmentoid.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No Segmentoid cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


Segmentoid_sub = Segmentoid_sub[
    common_barcodes
].copy()

ref_umap_Segmentoid = (
    ref_umap_Segmentoid.loc[
        common_barcodes
    ].copy()
)


# ============================================================
# 10. Prepare neighborhood graph
# ============================================================

sc.pp.highly_variable_genes(
    Segmentoid_sub
)

sc.pp.scale(
    Segmentoid_sub
)

sc.tl.pca(
    Segmentoid_sub,
    n_comps=30
)

sc.pp.neighbors(
    Segmentoid_sub,
    n_neighbors=30,
    n_pcs=30
)


# ============================================================
# 11. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    Segmentoid_sub,
    mode="stochastic"
)

scv.tl.velocity_graph(
    Segmentoid_sub
)


# ============================================================
# 12. Add PCW3 reference UMAP coordinates
# ============================================================

Segmentoid_sub.obsm["X_umap"] = (
    ref_umap_Segmentoid.to_numpy()
)


# ============================================================
# 13. Separate NMP cells by time point
# ============================================================

Segmentoid_sub.obs["celltype_time"] = (
    Segmentoid_sub.obs[
        "predicted.celltype"
    ].astype(str)
)


mask_nmp = (
    Segmentoid_sub.obs[
        "celltype_time"
    ] == "NMP"
)


Segmentoid_sub.obs.loc[
    mask_nmp,
    "celltype_time"
] = (
    Segmentoid_sub.obs.loc[
        mask_nmp,
        "celltype_time"
    ]
    + "_"
    + Segmentoid_sub.obs.loc[
        mask_nmp,
        "timepoint"
    ].astype(str)
)


print(
    "Cell types used for plotting:"
)

print(
    Segmentoid_sub.obs[
        "celltype_time"
    ].value_counts()
)


# ============================================================
# 14. Plot RNA velocity in PCW3 reference UMAP space
# ============================================================

# Keep text editable in SVG files
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

    "NMP_24h": "#F19890",
    "NMP_48h": "#FDEC8A",
    "NMP_72h": "#E79B2C",
    "NMP_98h": "#008EBE",

    "Cervical-Thoracic NT": "#27ae60",
    "Forebrain-Midbrain": "#F5761A",

    "Somite": "#9b59b6",
    "Anterior PSM": "#F5C26B",
    "Posterior PSM2": "#CBD6E2",
    "Posterior PSM1": "#74B72E",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    Segmentoid_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="celltype_time",
    palette=celltype_colors,
    density=2.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="Segmentoid RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 15. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Segmentoid_RNA_velocity_PCW3.svg"
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