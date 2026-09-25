import pandas as pd
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
import anndata as ad 
import loompy
print(loompy.__version__)

X = pd.read_csv("expr_RA_10000.csv", index_col=0).T
print(X.index[:5])
print(X.columns[:5])
# 去掉开头的时间前缀（例如 96h_、24h_ 等）
X.index = X.index.str.split('_', n=1).str[-1]

# 去掉结尾的 -1 或类似数字后缀
X.index = X.index.str.replace(r'.\d+$', '', regex=True)
print(X.index[:5])

RA_obs = pd.read_csv("meta_RA_10000.csv", index_col=0)
print(RA_obs.index[:5])
RA_obs.index = RA_obs.index.str.split('_', n=1).str[-1]
RA_obs.index = RA_obs.index.str.replace(r'-\d+$', '', regex=True)
print(RA_obs.index[:5])
# 把 X 的分隔符 - 改成 :
RA = ad.AnnData(X=X.astype(np.float32), obs=RA_obs)
print(RA)
print(X.index)    # 细胞名
print(RA_obs.index)   # 元数据中的细胞名
RA_96h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/RA_gastruloid_96h/RA_96h_loom/possorted_genome_bam_BRB6H.loom')
RA_120h_loom = sc.read_loom('/home/wangyt/data_pi-drukkerm/wangyt/ncbi/RA_gastruloid_120h/RA_120h_loom/possorted_genome_bam_TOPF7.loom')
RA_96h_loom.var_names_make_unique()
RA_120h_loom.var_names_make_unique()

RA_loom = ad.concat(
    [RA_96h_loom, RA_120h_loom],
    join='outer',
    label='time',
    keys=['96h', '120h']
)
print(RA_loom.obs_names)


RA_loom.obs_names = [x.split(':')[-1].rstrip('x') for x in RA_loom.obs_names]
RA.obs_names = [x.split(':')[-1].rstrip('x') for x in RA.obs_names]

print(RA.obs_names)
print(RA_loom.obs_names)


if RA_loom.obs_names.has_duplicates:
   RA_loom.obs_names_make_unique()
RA.obs_names_make_unique()
common_barcodes_RA = RA.obs_names.intersection(RA_loom.obs_names)
print(RA.obs_names)
print(common_barcodes_RA)
RA = RA[common_barcodes_RA].copy()
RA_loom = RA_loom[common_barcodes_RA].copy()

RA.obs_names_make_unique()
RA_loom.obs_names_make_unique()
common_genes_RA = RA.var_names.intersection(RA_loom.var_names)
RA_sub = RA[:, common_genes_RA].copy()
RA_loom_sub = RA_loom[:, common_genes_RA].copy()
RA_sub.layers['spliced'] = RA_loom_sub.layers['spliced']
RA_sub.layers['unspliced'] = RA_loom_sub.layers['unspliced']
print(RA_sub.shape)
pcw3 = pd.read_csv("subset_pcw3_umap.csv")
pcw3.index = pcw3['cell']
pcw3.index = pcw3.index.str.replace("-1$", "", regex=True)
pcw3 = pcw3[['UMAP_1', 'UMAP_2']]
# 首先对 ref_umap_Jianping 保留所需列
# 重新读取 ref_umap_Jianping 并清洗
ref_umap_RA = pd.read_csv("RA_umap_10000.csv")
ref_umap_RA.index = ref_umap_RA['cell']
ref_umap_RA.index = [x.split('_')[1].replace("-1", "") for x in ref_umap_RA.index]

ref_umap_RA = ref_umap_RA[['UMAP_1', 'UMAP_2']]
print(ref_umap_RA.index)  # 检查索引
ref_umap_RA = ref_umap_RA[~ref_umap_RA.index.duplicated(keep='first')]
common_barcodes = RA_sub.obs_names.intersection(ref_umap_RA.index)

RA_sub = RA_sub[common_barcodes].copy()
ref_umap_RA = ref_umap_RA.loc[common_barcodes].copy()

print(RA_sub)
# 添加 UMAP 坐标
RA_sub.obsm["X_umap"] = ref_umap_RA.to_numpy()

# scvelo preprocessing
scv.pp.filter_and_normalize(RA_sub)

# compute moments using scvelo
scv.pp.moments(RA_sub, n_pcs=30, n_neighbors=30)


# Step 3: 计算 RNA velocity（基于 dynamical model）

scv.tl.velocity(RA_sub, mode="stochastic")
scv.tl.velocity_graph(RA_sub)
RA_sub.write("RA_sub_velocity.h5ad")
RA_sub = scv.read("RA_sub_velocity.h5ad")
print(RA_sub.obs['predicted.celltype'].unique())
RA_sub.obs['predicted.celltype.score']

# 先把 celltype_time 赋值为原始 predicted.celltype
RA_sub.obs['celltype_time'] = RA_sub.obs['predicted.celltype'].astype(str)

mask_nmp = RA_sub.obs['predicted.celltype'] == "NMP"
RA_sub.obs.loc[mask_nmp, 'celltype_time'] = (
    RA_sub.obs.loc[mask_nmp, 'predicted.celltype'].astype(str) + "_" +
    RA_sub.obs.loc[mask_nmp, 'timepoint'].astype(str)
)
print(RA_sub.obs['celltype_time'].unique())

fig, ax = plt.subplots(figsize=(8,6))

# 先画PCW3背景细胞
ax.scatter(
    pcw3["UMAP_1"], pcw3["UMAP_2"],
    color='lightgrey', alpha=0.3, s=5,
    label='PCW3 background'
)

my_colors = {
    'NMP_96h': '#e74c3c',
    'NMP_120h': "#eee313",
    'Cervical-Thoracic NT': '#27ae60',
    #'Hindbrain': '#3498db',
    #'Forebrain-Midbrain': '#F5761A',
    'Somite': '#9b59b6',
    #'Anterior PSM': '#F5C26B',
    #'Posterior PSM2': '#CBD6E2',
    #'Posterior PSM1': '#74B72E',
    #'Notochord': '#67032F'
}

# 再画trunk velocity流
scv.pl.velocity_embedding_stream(
    RA_sub,
    basis="umap",
    ax=ax,
    size=200,
    color="celltype_time", 
    palette=my_colors,
    density=3.0,
    arrow_size=1,
    linewidth=2.0,
    smooth=0.8,
    min_mass=3,
    title="RNA Velocity in PCW3 UMAP space",
    legend_loc="right margin",
    show=False
)

plt.tight_layout()

# 保存图片，格式可以是png、pdf、svg等
plt.savefig("RA_new.png", dpi=300)

#-----------------------------------------------
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# -----------------------------
# 1️⃣ 读取 AnnData
# -----------------------------
RA_sub = scv.read("RA_sub_velocity.h5ad")

# -----------------------------
# 2️⃣ 定义颜色映射
# -----------------------------
my_colors = {
    'NMP': '#e74c3c',
    'Somite': '#9b59b6',
    'Cervical-Thoracic NT': '#27ae60'
}

# -----------------------------
# 3️⃣ 绘图：基础 UMAP 背景
# -----------------------------
fig, ax = plt.subplots(figsize=(8,6))

# 所有细胞灰色背景
ax.scatter(
    RA_sub.obsm['X_umap'][:,0],
    RA_sub.obsm['X_umap'][:,1],
    color='lightgrey', alpha=0.3, s=5,
    label='All cells'
)

# -----------------------------
# 4️⃣ 根据 predicted.celltype 绘制不同颜色
# -----------------------------
for ct, color in my_colors.items():
    mask = RA_sub.obs['predicted.celltype'] == ct
    ax.scatter(
        RA_sub.obsm['X_umap'][mask,0],
        RA_sub.obsm['X_umap'][mask,1],
        color=color,
        s=20,
        label=ct,
        alpha=0.8
    )

# -----------------------------
# 5️⃣ 用 prediction.score 标记高低可靠性
# -----------------------------
# 这里假设 score > 0.8 是高置信度
high_conf_mask = RA_sub.obs['predicted.celltype.score'] > 0.8
ax.scatter(
    RA_sub.obsm['X_umap'][high_conf_mask,0],
    RA_sub.obsm['X_umap'][high_conf_mask,1],
    facecolors='none', edgecolors='black', s=50, linewidths=0.8,
    label='high confidence'
)

# -----------------------------
# 6️⃣ 绘制 RNA velocity
# -----------------------------
scv.pl.velocity_embedding_stream(
    RA_sub,
    basis='umap',
    ax=ax,
    color=None,            # 已用手动颜色
    density=2.0,
    arrow_size=1,
    linewidth=1.5,
    smooth=0.6,
    min_mass=3,
    show=False
)

# -----------------------------
# 7️⃣ 美化图
# -----------------------------
ax.legend(markerscale=1.5, fontsize=10)
ax.set_xlabel('UMAP1')
ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.savefig("RA_label_transfer_check.png", dpi=300)
plt.show()
