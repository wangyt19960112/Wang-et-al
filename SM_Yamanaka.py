import pandas as pd
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad 
import loompy
X = pd.read_csv("expr_cantas_new.csv", index_col=0).T
X.index = X.index.str.replace('.', '-', regex=False)
X.index = X.index.str.lstrip('X') 
cantas_obs = pd.read_csv("meta_cantas_new.csv", index_col=0)
cantas = ad.AnnData(X=X, obs=cantas_obs)
print(X.index)    # 细胞名
print(cantas_obs.index) 
print(cantas.obs.columns)  # 元数据中的细胞名
cantas_obs.index = cantas_obs.index.str.replace(r'^\d+h_', '', regex=True)
cantas_obs.index = cantas_obs.index.str.replace(r'-1$', '', regex=True)

cantas.obs_names = cantas_obs.index

print(cantas.obs_names)
cantas_48h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/hAxioloid_48h/Cantas_48h_loom/Cantas_48h.loom')
cantas_72h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/hAxioloid_72h/Cantas_72h_loom/Cantas_72h.loom')
cantas_96h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/hAxioloid_96h/Cantas_96h_loom/Cantas_96h.loom')
cantas_120h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/hAxioloid_120h/Cantas_120h_loom/Cantas_120h.loom')

cantas_48h_loom.var_names_make_unique()
cantas_72h_loom.var_names_make_unique()
cantas_96h_loom.var_names_make_unique()
cantas_120h_loom.var_names_make_unique()

cantas_loom = ad.concat(
    [cantas_48h_loom, cantas_72h_loom, cantas_96h_loom, cantas_120h_loom],
    join='outer',
    label='time',
    keys=['48h', '72h', '96h', '120h']
)
print(cantas_loom.obs_names)
cantas_loom_barcodes = [x.split(':')[-1].rstrip('x') for x in cantas_loom.obs_names]
cantas_loom.obs_names = cantas_loom_barcodes

if cantas_loom.obs_names.has_duplicates:
   cantas_loom.obs_names_make_unique()
cantas.obs_names_make_unique()
common_barcodes_cantas = cantas.obs_names.intersection(cantas_loom.obs_names)
print(common_barcodes_cantas)
cantas = cantas[common_barcodes_cantas].copy()
cantas_loom = cantas_loom[common_barcodes_cantas].copy()


cantas_loom.obs_names_make_unique()
common_genes_cantas = cantas.var_names.intersection(cantas_loom.var_names)
cantas_sub = cantas[:, common_genes_cantas].copy()
cantas_loom_sub = cantas_loom[:, common_genes_cantas].copy()
cantas_sub.layers['spliced'] = cantas_loom_sub.layers['spliced']
cantas_sub.layers['unspliced'] = cantas_loom_sub.layers['unspliced']
print(cantas_sub.shape)
pcw3 = pd.read_csv("subset_pcw3_umap.csv")
pcw3.index = pcw3['cell']
pcw3.index = pcw3.index.str.replace("-1$", "", regex=True)
pcw3 = pcw3[['UMAP_1', 'UMAP_2']]
# 首先对 ref_umap_Jianping 保留所需列
# 重新读取 ref_umap_Jianping 并清洗
ref_umap_cantas = pd.read_csv("cantas_umap_new.csv")
ref_umap_cantas.index = ref_umap_cantas['cell']
ref_umap_cantas.index = ref_umap_cantas.index.str.replace("-1$", "", regex=True)
ref_umap_cantas.index = ref_umap_cantas.index.str.replace(r'^\d+h_', '', regex=True)
ref_umap_cantas = ref_umap_cantas[['UMAP_1', 'UMAP_2']]
print(ref_umap_cantas.index)  # 检查索引
ref_umap_cantas = ref_umap_cantas[~ref_umap_cantas.index.duplicated(keep='first')]
common_barcodes = cantas_sub.obs_names.intersection(ref_umap_cantas.index)

cantas_sub = cantas_sub[common_barcodes].copy()
ref_umap_cantas = ref_umap_cantas.loc[common_barcodes].copy()

sc.pp.highly_variable_genes(cantas_sub, n_top_genes=2000, flavor='seurat_v3', inplace=True)
cantas_sub._inplace_subset_var(cantas_sub.var['highly_variable'])
sc.pp.scale(cantas_sub)
sc.tl.pca(cantas_sub, n_comps=30)

# -------------------------------
# 6) Neighbors and RNA velocity
# -------------------------------
# neighbors on PCA
sc.pp.neighbors(cantas_sub, n_neighbors=30, n_pcs=30)
scv.pp.moments(cantas_sub, n_pcs=30, n_neighbors=30)
scv.tl.velocity(cantas_sub, mode="stochastic")
scv.tl.velocity_graph(cantas_sub)

cantas_sub.obsm["X_umap"] = ref_umap_cantas.to_numpy()

cantas_sub.obs["celltype_time"] = cantas_sub.obs["predicted.celltype"].astype(str)
print(cantas_sub.obs.columns)


# 细分 NMP
mask_nmp = cantas_sub.obs["celltype_time"] == "NMP"
cantas_sub.obs.loc[mask_nmp, "celltype_time"] = (
    cantas_sub.obs.loc[mask_nmp, "celltype_time"] + "_" + cantas_sub.obs.loc[mask_nmp, "timepoint"].astype(str)
)

fig, ax = plt.subplots(figsize=(8,6))

# 先画PCW3背景细胞
ax.scatter(
    pcw3["UMAP_1"], pcw3["UMAP_2"],
    color='lightgrey', alpha=0.3, s=5,
    label='PCW3 background'
)
#time ={'48h': "#e41a1c", '72h': "#eee313", '96h': "#A35308", '120h': "#377eb8"}
my_colors = {
    'NMP_48h': '#F19890',
    'NMP_72h': "#FDEC8A",
    'NMP_96h': "#E79B2C",
    'NMP_120h': "#008EBE",
    #'Cervical-Thoracic NT': '#27ae60',
    #'Hindbrain': '#3498db',
    #'Forebrain-Midbrain': '#F5761A',
    'Somite': '#9b59b6',
    'Anterior PSM': '#F5C26B',
    'Posterior PSM2': '#CBD6E2',
    'Posterior PSM1': '#74B72E'
    #'Notochord': '#67032F'
}

# 再画trunk velocity流
scv.pl.velocity_embedding_stream(
    cantas_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="celltype_time", 
    palette=my_colors,
    density=3.0,
    arrow_size=1,
    linewidth=0.3,
    smooth=0.8,
    min_mass=3,
    title="RNA Velocity in PCW3 UMAP space",
    legend_loc='right',  # 有些版本 prefer None
    show=False
)


plt.tight_layout()


plt.savefig("cantas_new2.svg", dpi=300)
