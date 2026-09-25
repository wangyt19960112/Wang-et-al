"""
RNA velocity analysis of the Jianping dataset projected onto
the PCW3 reference UMAP space.

The Jianping dataset contains Day 4 and Day 9 samples.
The corresponding loom files are combined before RNA velocity analysis.

Required input files:
    data/
    ├── expr_Jianping_new.csv
    ├── meta_Jianping_new.csv
    ├── Jianping.loom
    ├── Jianping_day9.loom
    ├── subset_pcw3_umap.csv
    └── Jianping_umap_new.csv

Output:
    results/rna_velocity/Jianping_RNA_velocity_PCW3.svg
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

EXPRESSION_FILE = DATA_DIR / "expr_Jianping_new.csv"
METADATA_FILE = DATA_DIR / "meta_Jianping_new.csv"

LOOM_DAY4_FILE = DATA_DIR / "Jianping.loom"
LOOM_DAY9_FILE = DATA_DIR / "Jianping_day9.loom"

PCW3_UMAP_FILE = DATA_DIR / "subset_pcw3_umap.csv"
JIANPING_UMAP_FILE = DATA_DIR / "Jianping_umap_new.csv"


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

Jianping_obs = pd.read_csv(
    METADATA_FILE,
    index_col=0
)

Jianping = ad.AnnData(
    X=X,
    obs=Jianping_obs
)


# Remove "-1" suffix
Jianping.obs_names = (
    Jianping.obs_names.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove prefixes such as "day4_" and "day9_"
Jianping.obs_names = (
    Jianping.obs_names.str.replace(
        r"^day\d+_",
        "",
        regex=True
    )
)

Jianping.obs_names_make_unique()


print(
    f"Jianping AnnData: "
    f"{Jianping.n_obs} cells × {Jianping.n_vars} genes"
)


# ============================================================
# 3. Load Day 4 and Day 9 loom files
# ============================================================

Jianping_loom_day4 = sc.read_loom(
    LOOM_DAY4_FILE
)

Jianping_loom_day9 = sc.read_loom(
    LOOM_DAY9_FILE
)


Jianping_loom_day4.var_names_make_unique()
Jianping_loom_day9.var_names_make_unique()


# ============================================================
# 4. Combine loom files
# ============================================================

Jianping_loom = ad.concat(
    [
        Jianping_loom_day4,
        Jianping_loom_day9
    ],
    join="outer",
    label="day",
    keys=[
        "day4",
        "day9"
    ]
)


# Convert loom cell names to barcode format used in Jianping
Jianping_loom.obs_names = [
    name.split(":")[-1].rstrip("x")
    for name in Jianping_loom.obs_names
]


Jianping_loom.obs_names_make_unique()


print(
    f"Combined loom object: "
    f"{Jianping_loom.n_obs} cells × "
    f"{Jianping_loom.n_vars} genes"
)


# ============================================================
# 5. Match cells between expression data and loom data
# ============================================================

common_barcodes = (
    Jianping.obs_names.intersection(
        Jianping_loom.obs_names
    )
)


print(
    f"Matched cells between expression and loom data: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No matching cell barcodes were found between "
        "the Jianping expression matrix and loom files."
    )


Jianping = Jianping[
    common_barcodes
].copy()

Jianping_loom = Jianping_loom[
    common_barcodes
].copy()


# ============================================================
# 6. Match genes and add RNA velocity layers
# ============================================================

common_genes = (
    Jianping.var_names.intersection(
        Jianping_loom.var_names
    )
)


print(
    f"Matched genes: {len(common_genes)}"
)


if len(common_genes) == 0:
    raise ValueError(
        "No matching genes were found between "
        "the Jianping expression matrix and loom files."
    )


Jianping_sub = Jianping[
    :,
    common_genes
].copy()

Jianping_loom_sub = Jianping_loom[
    :,
    common_genes
].copy()


Jianping_sub.layers["spliced"] = (
    Jianping_loom_sub.layers["spliced"].copy()
)

Jianping_sub.layers["unspliced"] = (
    Jianping_loom_sub.layers["unspliced"].copy()
)


print(
    f"Final Jianping object: "
    f"{Jianping_sub.n_obs} cells × "
    f"{Jianping_sub.n_vars} genes"
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
# 8. Load Jianping coordinates in PCW3 reference UMAP space
# ============================================================

ref_umap_Jianping = pd.read_csv(
    JIANPING_UMAP_FILE
)

ref_umap_Jianping.index = (
    ref_umap_Jianping["cell"]
)


# Remove "-1" suffix
ref_umap_Jianping.index = (
    ref_umap_Jianping.index.str.replace(
        r"-1$",
        "",
        regex=True
    )
)


# Remove day4_ / day9_ prefixes
ref_umap_Jianping.index = (
    ref_umap_Jianping.index.str.replace(
        r"^day\d+_",
        "",
        regex=True
    )
)


ref_umap_Jianping = ref_umap_Jianping[
    ["UMAP_1", "UMAP_2"]
]


# Remove duplicated cell IDs
ref_umap_Jianping = ref_umap_Jianping[
    ~ref_umap_Jianping.index.duplicated(
        keep="first"
    )
]


# ============================================================
# 9. Match cells with reference UMAP coordinates
# ============================================================

common_barcodes = (
    Jianping_sub.obs_names.intersection(
        ref_umap_Jianping.index
    )
)


print(
    f"Cells with reference UMAP coordinates: "
    f"{len(common_barcodes)}"
)


if len(common_barcodes) == 0:
    raise ValueError(
        "No Jianping cells could be matched to "
        "the PCW3 reference UMAP coordinates."
    )


Jianping_sub = Jianping_sub[
    common_barcodes
].copy()

ref_umap_Jianping = ref_umap_Jianping.loc[
    common_barcodes
].copy()


# Add PCW3 reference UMAP coordinates
Jianping_sub.obsm["X_umap"] = (
    ref_umap_Jianping.to_numpy()
)


# ============================================================
# 10. RNA velocity preprocessing
# ============================================================

scv.pp.filter_and_normalize(
    Jianping_sub
)

scv.pp.moments(
    Jianping_sub,
    n_pcs=30,
    n_neighbors=30
)


# ============================================================
# 11. Estimate RNA velocity
# ============================================================

# Stochastic RNA velocity model
scv.tl.velocity(
    Jianping_sub,
    mode="stochastic"
)

scv.tl.velocity_graph(
    Jianping_sub
)


# ============================================================
# 12. Plot RNA velocity in PCW3 reference UMAP space
# ============================================================

# Keep SVG text editable
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
    Jianping_sub,
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
    title="Jianping RNA velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)


fig.tight_layout()


# ============================================================
# 13. Save figure
# ============================================================

OUTPUT_FILE = (
    RESULTS_DIR
    / "Jianping_RNA_velocity_PCW3.svg"
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