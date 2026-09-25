import pandas as pd
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad 
import loompy
print(loompy.__version__)
# 读取数据
X = pd.read_csv("expr_Barcelona_new.csv", index_col=0).T
X.index = X.index.str.replace('.', '-', regex=False)
Barcelona_obs = pd.read_csv("meta_Barcelona_new.csv", index_col=0)
Barcelona = ad.AnnData(X=X, obs=Barcelona_obs)
print(X.index)    # 细胞名
print(Barcelona_obs.index)  # 元数据中的细胞名
Barcelona_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/chapter2/Barcelona_somite.loom')
Barcelona_loom.var_names_make_unique()
print(Barcelona_loom.obs_names)
Barcelona_loom_barcodes = [x.split(':')[-1].rstrip('x') for x in Barcelona_loom.obs_names]
Barcelona_loom.obs_names = Barcelona_loom_barcodes
Barcelona.obs_names = Barcelona_obs.index.str.replace(r'-1$', '', regex=True)

print(Barcelona.obs_names)

if Barcelona_loom.obs_names.has_duplicates:
   Barcelona_loom.obs_names_make_unique()
Barcelona.obs_names_make_unique()
common_barcodes_Barcelona = Barcelona.obs_names.intersection(Barcelona_loom.obs_names)
print(Barcelona.obs_names)
print(common_barcodes_Barcelona)
Barcelona = Barcelona[common_barcodes_Barcelona].copy()
Barcelona_loom = Barcelona_loom[common_barcodes_Barcelona].copy()

Barcelona.obs_names_make_unique()
Barcelona_loom.obs_names_make_unique()
common_genes_Barcelona = Barcelona.var_names.intersection(Barcelona_loom.var_names)
Barcelona_sub = Barcelona[:, common_genes_Barcelona].copy()
Barcelona_loom_sub = Barcelona_loom[:, common_genes_Barcelona].copy()
Barcelona_sub.layers['spliced'] = Barcelona_loom_sub.layers['spliced']
Barcelona_sub.layers['unspliced'] = Barcelona_loom_sub.layers['unspliced']
print(Barcelona_sub.shape)
pcw3 = pd.read_csv("subset_pcw3_umap.csv")
pcw3.index = pcw3['cell']
pcw3.index = pcw3.index.str.replace("-1$", "", regex=True)
pcw3 = pcw3[['UMAP_1', 'UMAP_2']]
# 首先对 ref_umap_Jianping 保留所需列
# 重新读取 ref_umap_Jianping 并清洗
ref_umap_Barcelona = pd.read_csv("Barcelona_umap_new.csv")
ref_umap_Barcelona.index = ref_umap_Barcelona['cell']
ref_umap_Barcelona.index = ref_umap_Barcelona.index.str.replace("-1$", "", regex=True)

ref_umap_Barcelona = ref_umap_Barcelona[['UMAP_1', 'UMAP_2']]
print(ref_umap_Barcelona.index)  # 检查索引
ref_umap_Barcelona = ref_umap_Barcelona[~ref_umap_Barcelona.index.duplicated(keep='first')]
common_barcodes = Barcelona_sub.obs_names.intersection(ref_umap_Barcelona.index)

Barcelona_sub = Barcelona_sub[common_barcodes].copy()
ref_umap_Barcelona = ref_umap_Barcelona.loc[common_barcodes].copy()


# 添加 UMAP 坐标
Barcelona_sub.obsm["X_umap"] = ref_umap_Barcelona.to_numpy()

# scvelo preprocessing
scv.pp.filter_and_normalize(Barcelona_sub)

# compute moments using scvelo
scv.pp.moments(Barcelona_sub, n_pcs=30, n_neighbors=30)


# Step 3: 计算 RNA velocity（基于 dynamical model）

scv.tl.velocity(Barcelona_sub, mode="stochastic")
scv.tl.velocity_graph(Barcelona_sub)

fig, ax = plt.subplots(figsize=(8,6))

# 先画PCW3背景细胞
ax.scatter(
    pcw3["UMAP_1"], pcw3["UMAP_2"],
    color='lightgrey', alpha=0.3, s=5,
    label='PCW3 background'
)

my_colors = {
    #'NMP_day3_24h': '#e74c3c',
    #'NMP_day5_24h': "#eee313",
    #'NMP_day7_24h': "#A35308",
    'NMP': '#e74c3c',
    #'Cervical-Thoracic NT': '#27ae60',
    #'Hindbrain': '#3498db',
    #'Forebrain-Midbrain': '#F5761A'
    'Somite': '#9b59b6',
    'Anterior PSM': '#F5C26B'
    #'Posterior PSM2': '#CBD6E2',
    #'Posterior PSM1': '#74B72E',
    #'Notochord': '#67032F'
}

# 再画trunk velocity流
scv.pl.velocity_embedding_stream(
    Barcelona_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="predicted.celltype", 
    palette=my_colors,
    density=2.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="SM_San RNA Velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)

plt.tight_layout()

# 保存图片，格式可以是png、pdf、svg等
plt.savefig("SM_San.svg", dpi=300)