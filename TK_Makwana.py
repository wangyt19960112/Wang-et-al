"""
RNA velocity analysis of the Naomi axial organoid dataset
projected onto the PCW3 reference UMAP space.

The Naomi dataset contains three time points:
72 h, 96 h, and 120 h.

NMP cells are separated by time point for visualization.

Required input files:
    data/
    ├── expr_naomi_new.csv
    ├── meta_naomi_new.csv
    ├── naomi_72h.loom
    ├── naomi_96h.loom
    ├── naomi_120h.loom
    ├── subset_pcw3_umap.csv
    └── naomi_umap_new.csv

Output:
    results/rna_velocity/Naomi_RNA_velocity_PCW3.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_naomi_new.csv"
METADATA_FILE = DATA_DIR / "meta_naomi_new.csv"

LOOM_72H_FILE = DATA_DIR / "naomi_72h.loom"
LOOM_96H_FILE = DATA_DIR / "naomi_96h.loom"
LOOM_120H_FILE = DATA_DIR / "naomi_120h.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
NAOMI_UMAP_FILE = DATA_DIR / "naomi_umap_new.csv"


# ============================================================
# 2. Construct Naomi AnnData object
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


naomi_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)


naomi = ad.AnnData(
    X=X,
    obs=naomi_obs
)


# Remove "-1" suffix
naomi.obs_names = (
    naomi.obs_names.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove time-point prefixes such as "72h_", "96h_", etc.
naomi.obs_names = (
    naomi.obs_names.str.replace(
        r"^\d+h_",
        "",
        regex=True
    )
)


naomi.obs_names_make_unique()


print(
    f"Naomi AnnData: "
    f"{naomi.n_obs} cells × "
    f"{naomi.n_vars} genes"
)


# ============================================================
# 3. Load loom files
# ============================================================

naomi_72h_loom = sc.read_loom(
    LOOM_72H_FILE
)

naomi_96h_loom = sc.read_loom(
    LOOM_96H_FILE
)

naomi_120h_loom = sc.read_loom(
    LOOM_120H_FILE
)


# Ensure unique gene names
for loom_obj in [
    naomi_72h_loom,
    naomi_96h_loom,
    naomi_120h_loom,
]:
    loom_obj.var_names_make_unique()


# ============================================================
# 4. Combine loom files
# ============================================================

naomi_loom = ad.concat(
    [
        naomi_72h_loom,
        naomi_96h_loom,
        naomi_120h_loom,
    ],
    join="outer",
    label="time",
    keys=[
        "72h",
        "96h",
        "120h",
    ]
)


# Standardize loom cell barcodes
naomi_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in naomi_loom.obs_names
]


naomi_loom.obs_names_make_unique()


print(
    f"Combined Naomi loom object: "
    f"{naomi_loom.n_obs} cells × "
    f"{naomi_loom.n_vars} genes"
)


# ============================================================
# 5. Match cells between expression and loom data
# ============================================================

common_barcodes = (
    naomi.obs_names.intersection(
        naomi_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Naomi expression matrix and loom files."
    )


naomi = naomi[
    common_barcodes
].copy()


naomi_loom = naomi_loom[
    common_barcodes
].copy()


# ============================================================
# 6. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    naomi.var_names.intersection(
        naomi_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Naomi expression matrix and loom files."
    )


naomi_sub = naomi[
    :,
    common_genes
].copy()


naomi_loom_sub = naomi_loom[
    :,
    common_genes
].copy()


naomi_sub.layers["spliced"] = (
    naomi_loom_sub.layers["spliced"].copy()
)

naomi_sub.layers["unspliced"] = (
    naomi_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Naomi object: "
    f"{naomi_sub.n_obs} cells × "
    f"{naomi_sub.n_vars} genes"
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
# 8. Load Naomi coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_naomi = pd.read_csv(
    NAOMI_UMAP_FILE
)


ref_umap_naomi.index = (
    ref_umap_naomi["cell"]
)


# Remove "-1" suffix
ref_umap_naomi.index = (
    ref_umap_naomi.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove time-point prefixes
ref_umap_naomi.index = (
    ref_umap_naomi.index.str.replace(
        r"^\d+h_",
        "",
        regex=True
    )
)


ref_umap_naomi = (
    ref_umap_naomi[
        ["UMAP_1", "UMAP_2"]
    ]
)


# Remove duplicated cell IDs
ref_umap_naomi = (
    ref_umap_naomi[
        ~ref_umap_naomi.index.duplicated(
            keep="first"
        )
    ]
)


# ============================================================
# 9. Match cells with PCW3 reference UMAP coordinates
# ============================================================

common_barcodes = (
    naomi_sub.obs_names.intersection(
        ref_umap_naomi.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No Naomi cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


naomi_sub = naomi_sub[
    common_barcodes
].copy()


ref_umap_naomi = (
    ref_umap_naomi.loc[
        common_barcodes
    ].copy()
)


# Add PCW3 reference UMAP coordinates
naomi_sub.obsm["X_umap"] = (
    ref_umap_naomi.to_numpy()
)


# ============================================================
# 10. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    naomi_sub
)


scv.pp.moments(
    naomi_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 11. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    naomi_sub,
    mode="stochastic"
)


scv.tl.velocity_graph(
    naomi_sub
)


# ============================================================
# 12. Separate NMP cells by time point
# ============================================================

if "predicted.celltype" not in naomi_sub.obs.columns:
    raise KeyError(
        "'predicted.celltype' was not found in Naomi metadata."
    )


if "timepoint" not in naomi_sub.obs.columns:
    raise KeyError(
        "'timepoint' was not found in Naomi metadata."
    )


# Start with the original predicted cell type
naomi_sub.obs["celltype_time"] = (
    naomi_sub.obs[
        "predicted.celltype"
    ].astype(str)
)


# Separate NMP cells according to time point
mask_nmp = (
    naomi_sub.obs[
        "celltype_time"
    ] == "NMP"
)


naomi_sub.obs.loc[
    mask_nmp,
    "celltype_time"
] = (
    naomi_sub.obs.loc[
        mask_nmp,
        "celltype_time"
    ]
    + "_"
    + naomi_sub.obs.loc[
        mask_nmp,
        "timepoint"
    ].astype(str)
)


print(
    "Cell types used for plotting:"
)

print(
    naomi_sub.obs[
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

    "NMP_72h": "#F19890",
    "NMP_96h": "#E79B2C",
    "NMP_120h": "#008EBE",

    "Cervical-Thoracic NT": "#27ae60",

    "Somite": "#9b59b6",
    "Anterior PSM": "#F5C26B",
    "Posterior PSM2": "#CBD6E2",
    "Posterior PSM1": "#74B72E",
}


# RNA velocity stream
scv.pl.velocity_embedding_stream(
    naomi_sub,
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
    title="Naomi RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 14. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Naomi_RNA_velocity_PCW3.svg"
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