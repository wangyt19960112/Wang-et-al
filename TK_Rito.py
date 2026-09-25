"""
RNA velocity analysis of the Tiago trunk organoid dataset
projected onto the PCW3 reference UMAP space.

The Tiago dataset contains three time points:
day 3, day 5, and day 7.

NMP cells are separated by time point for visualization.

Required input files:
    data/
    ├── expr_tiago_new.csv
    ├── meta_tiago_new.csv
    ├── tiago_day3.loom
    ├── tiago_day5.loom
    ├── tiago_day7.loom
    ├── subset_pcw3_umap.csv
    └── tiago_umap_new.csv

Output:
    results/rna_velocity/Tiago_RNA_velocity_PCW3.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_tiago_new.csv"
METADATA_FILE = DATA_DIR / "meta_tiago_new.csv"

LOOM_DAY3_FILE = DATA_DIR / "tiago_day3.loom"
LOOM_DAY5_FILE = DATA_DIR / "tiago_day5.loom"
LOOM_DAY7_FILE = DATA_DIR / "tiago_day7.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
TIAGO_UMAP_FILE = DATA_DIR / "tiago_umap_new.csv"


# ============================================================
# 2. Construct Tiago AnnData object
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

# Remove leading "X" introduced during data export
X.index = X.index.str.lstrip("X")


tiago_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)


tiago = ad.AnnData(
    X=X,
    obs=tiago_obs
)


# Remove "-1" suffix
tiago.obs_names = (
    tiago.obs_names.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove day prefixes such as "day3_", "day5_", and "day7_"
tiago.obs_names = (
    tiago.obs_names.str.replace(
        r"^day\d+_",
        "",
        regex=True
    )
)


tiago.obs_names_make_unique()


print(
    f"Tiago AnnData: "
    f"{tiago.n_obs} cells × "
    f"{tiago.n_vars} genes"
)


# ============================================================
# 3. Load loom files
# ============================================================

tiago_day3_loom = sc.read_loom(
    LOOM_DAY3_FILE
)

tiago_day5_loom = sc.read_loom(
    LOOM_DAY5_FILE
)

tiago_day7_loom = sc.read_loom(
    LOOM_DAY7_FILE
)


# Ensure unique gene names
for loom_obj in [
    tiago_day3_loom,
    tiago_day5_loom,
    tiago_day7_loom,
]:
    loom_obj.var_names_make_unique()


# ============================================================
# 4. Combine loom files
# ============================================================

tiago_loom = ad.concat(
    [
        tiago_day3_loom,
        tiago_day5_loom,
        tiago_day7_loom,
    ],
    join="outer",
    label="day",
    keys=[
        "day3",
        "day5",
        "day7",
    ]
)


# Standardize loom cell barcodes
tiago_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in tiago_loom.obs_names
]


tiago_loom.obs_names_make_unique()


print(
    f"Combined Tiago loom object: "
    f"{tiago_loom.n_obs} cells × "
    f"{tiago_loom.n_vars} genes"
)


# ============================================================
# 5. Match cells between expression and loom data
# ============================================================

common_barcodes = (
    tiago.obs_names.intersection(
        tiago_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Tiago expression matrix and loom files."
    )


tiago = tiago[
    common_barcodes
].copy()


tiago_loom = tiago_loom[
    common_barcodes
].copy()


# ============================================================
# 6. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    tiago.var_names.intersection(
        tiago_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Tiago expression matrix and loom files."
    )


tiago_sub = tiago[
    :,
    common_genes
].copy()


tiago_loom_sub = tiago_loom[
    :,
    common_genes
].copy()


tiago_sub.layers["spliced"] = (
    tiago_loom_sub.layers["spliced"].copy()
)

tiago_sub.layers["unspliced"] = (
    tiago_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Tiago object: "
    f"{tiago_sub.n_obs} cells × "
    f"{tiago_sub.n_vars} genes"
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
# 8. Load Tiago coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_tiago = pd.read_csv(
    TIAGO_UMAP_FILE
)


ref_umap_tiago.index = (
    ref_umap_tiago["cell"]
)


# Remove "-1" suffix
ref_umap_tiago.index = (
    ref_umap_tiago.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove day prefixes
ref_umap_tiago.index = (
    ref_umap_tiago.index.str.replace(
        r"^day\d+_",
        "",
        regex=True
    )
)


ref_umap_tiago = (
    ref_umap_tiago[
        ["UMAP_1", "UMAP_2"]
    ]
)


# Remove duplicated cell IDs
ref_umap_tiago = (
    ref_umap_tiago[
        ~ref_umap_tiago.index.duplicated(
            keep="first"
        )
    ]
)


# ============================================================
# 9. Match cells with PCW3 reference UMAP coordinates
# ============================================================

common_barcodes = (
    tiago_sub.obs_names.intersection(
        ref_umap_tiago.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No Tiago cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


tiago_sub = tiago_sub[
    common_barcodes
].copy()


ref_umap_tiago = (
    ref_umap_tiago.loc[
        common_barcodes
    ].copy()
)


# Add PCW3 reference UMAP coordinates
tiago_sub.obsm["X_umap"] = (
    ref_umap_tiago.to_numpy()
)


# ============================================================
# 10. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    tiago_sub
)


scv.pp.moments(
    tiago_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 11. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    tiago_sub,
    mode="stochastic"
)


scv.tl.velocity_graph(
    tiago_sub
)


# ============================================================
# 12. Separate NMP cells by time point
# ============================================================

if "predicted.celltype" not in tiago_sub.obs.columns:
    raise KeyError(
        "'predicted.celltype' was not found in Tiago metadata."
    )


if "timepoint" not in tiago_sub.obs.columns:
    raise KeyError(
        "'timepoint' was not found in Tiago metadata."
    )


# Start with the original predicted cell type
tiago_sub.obs["celltype_time"] = (
    tiago_sub.obs[
        "predicted.celltype"
    ].astype(str)
)


# Separate NMP cells according to time point
mask_nmp = (
    tiago_sub.obs[
        "celltype_time"
    ] == "NMP"
)


tiago_sub.obs.loc[
    mask_nmp,
    "celltype_time"
] = (
    tiago_sub.obs.loc[
        mask_nmp,
        "celltype_time"
    ]
    + "_"
    + tiago_sub.obs.loc[
        mask_nmp,
        "timepoint"
    ].astype(str)
)


print(
    "Cell types used for plotting:"
)

print(
    tiago_sub.obs[
        "celltype_time"
    ].value_counts()
)


# ============================================================
# 13. Plot RNA velocity in PCW3 reference UMAP space
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

    "NMP_day3_24h": "#e74c3c",
    "NMP_day5_24h": "#eee313",
    "NMP_day7_24h": "#A35308",

    "Cervical-Thoracic NT": "#27ae60",
    "Hindbrain": "#3498db",
    "Forebrain-Midbrain": "#F5761A",

    "Somite": "#9b59b6",
    "Notochord": "#67032F",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    tiago_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="celltype_time",
    palette=celltype_colors,
    density=3.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="Tiago RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 14. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Tiago_RNA_velocity_PCW3.svg"
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