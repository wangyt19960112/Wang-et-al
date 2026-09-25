# Axial organoid scRNA-seq reference mapping
# GitHub-ready version
#
# Before running:
#   1. Put/download the required datasets under DATA_DIR using the folder
#      structure documented in your repository README, OR
#   2. Set environment variable AXIAL_DATA_DIR to the directory containing
#      the data.
#
# Large raw scRNA-seq files and generated results should normally not be
# committed to GitHub.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(svglite)
})

DATA_DIR <- Sys.getenv("AXIAL_DATA_DIR", unset = "data")
RESULTS_DIR <- Sys.getenv("AXIAL_RESULTS_DIR", unset = "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

PCW3_seurat <- readRDS(file.path(DATA_DIR, "chapter2", "PCW3_0802.rds"))
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#3182bd", "#31a354", "#756bb1", "#636363", "#e6550d",
  "#9ecae1", "#74c476", "#fd8d3c", "#969696", "#6baed6",
  "#e7ba52"
)
P <- DimPlot(
  PCW3_seurat,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
) +
  scale_color_manual(
    values = my_colors[1:length(unique(PCW3_seurat$cell_type))]
  ) +

  scale_x_continuous(breaks = NULL, name = NULL) +
  scale_y_continuous(breaks = NULL, name = NULL) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),      
    axis.line = element_blank(),       
    axis.ticks = element_blank(),      
    axis.text = element_blank(),       
    axis.title = element_blank(),      
    legend.position = "none"
  )

ggsave(file.path(RESULTS_DIR, "PCW3_cell_types_0802.svg"), plot = P, width = 15, height = 10)





subset_PCW3 <- readRDS(file.path(DATA_DIR, "chapter2", "subset_PCW3_0802.rds"))
Day4_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping", "filtered", "matrix.mtx.gz"))
Day4_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping", "filtered", "features.tsv.gz"), header = FALSE)
Day4_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping", "filtered", "barcodes.tsv.gz"), header = FALSE)


if (nrow(Day4_genes) != nrow(Day4_mat)) stop("Day4: gene count mismatch")
if (nrow(Day4_barcodes) != ncol(Day4_mat)) stop("Day4: barcode count mismatch")
rownames(Day4_mat) <- make.unique(Day4_genes$V2)
colnames(Day4_mat) <- Day4_barcodes$V1


day4 <- CreateSeuratObject(counts = Day4_mat, project = "day4")
day4[["percent.mt"]] <- PercentageFeatureSet(day4, pattern = "^MT-")
day4 <- subset(day4, subset = nFeature_RNA > 4000 & nFeature_RNA < 12000 & percent.mt < 15)
ribosomal_genes_day4 <- grep("^RPS|^RPL", rownames(day4), value = TRUE)
day4 <- subset(day4, features = setdiff(rownames(day4), ribosomal_genes_day4))
day4$timepoint <- "day4"
saveRDS(day4, file.path(RESULTS_DIR, "day4_Jianping_seurat.rds"))
# -------- Day9 --------
Day9_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping_day9", "filtered", "matrix.mtx.gz"))
Day9_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping_day9", "filtered", "features.tsv.gz"), header = FALSE)
Day9_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_Jianping_day9", "filtered", "barcodes.tsv.gz"), header = FALSE)


if (nrow(Day9_genes) != nrow(Day9_mat)) stop("Day9: gene count mismatch")
if (nrow(Day9_barcodes) != ncol(Day9_mat)) stop("Day9: barcode count mismatch")
rownames(Day9_mat) <- make.unique(Day9_genes$V2)
colnames(Day9_mat) <- Day9_barcodes$V1


day9 <- CreateSeuratObject(counts = Day9_mat, project = "day9")
day9[["percent.mt"]] <- PercentageFeatureSet(day9, pattern = "^MT-")
day9 <- subset(day9, subset = nFeature_RNA > 4000 & nFeature_RNA < 12000 & percent.mt < 15)
ribosomal_genes_day9 <- grep("^RPS|^RPL", rownames(day9), value = TRUE)
day9 <- subset(day9, features = setdiff(rownames(day9), ribosomal_genes_day9))
day9$timepoint <- "day9"
saveRDS(day9, file.path(RESULTS_DIR, "day9_Jianping_seurat.rds"))


combined_Jianping <- merge(day4, y = day9, add.cell.ids = c("day4", "day9"), project = "combined_Jianping")
cell_names_Jianping <- colnames(combined_Jianping)

combined_Jianping <- NormalizeData(combined_Jianping)
combined_Jianping <- FindVariableFeatures(combined_Jianping)
combined_Jianping <- ScaleData(combined_Jianping)
combined_Jianping <- RunPCA(combined_Jianping)

day9 <- NormalizeData(day9)
day9 <- FindVariableFeatures(day9)
day9 <- ScaleData(day9)
day9 <- RunPCA(day9)

anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = day9,
  reference.reduction = "pca",
  dims = 1:30
)


day9 <- MapQuery(
  anchorset = anchors, 
  reference = PCW3_seurat, 
  query = day9, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)


pred_table <- table(combined_Jianping$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
combined_Jianping <- subset(combined_Jianping, subset = predicted.celltype %in% valid_labels)
table(combined_Jianping$predicted.celltype)

PCW3_umap <- Embeddings(PCW3_seurat, "umap")
PCW3_umap <- as.data.frame(PCW3_umap)
colnames(PCW3_umap) <- c("UMAP_1", "UMAP_2") 

combined_umap <- Embeddings(combined_Jianping, "ref.umap") |> as.data.frame()
colnames(combined_umap) <- c("UMAP_1", "UMAP_2")
combined_umap$cell <- rownames(combined_umap)
combined_umap$celltype <- combined_Jianping$predicted.celltype
combined_umap$time <- combined_Jianping$timepoint
write.csv(combined_umap, file.path(RESULTS_DIR, "Jianping_umap.csv"), row.names = FALSE)

combined_Jianping <- JoinLayers(combined_Jianping)
str(combined_Jianping)
write.csv(GetAssayData(combined_Jianping, slot = "data"), "expr_Jianping.csv")


write.csv(combined_Jianping@meta.data, file.path(RESULTS_DIR, "meta_Jianping.csv"))

combined_umap$timepoint <- gsub("_.*", "", combined_umap$cell)
table(combined_umap$timepoint)
time_colors <- setNames(
  c("#ff0019ff", "#87CEEB"),  # 玫红 & 天蓝
  c("day4", "day9")
)

# ---------------------------------------------------------

library(ggplot2)
p <- ggplot() +
  geom_point(data = PCW3_umap, aes(x = UMAP_1, y = UMAP_2), 
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = combined_umap, 
             aes(x = UMAP_1, y = UMAP_2, color = timepoint), 
             alpha = 0.5, size = 0.4) +
  scale_color_manual(values = time_colors, name = "Timepoint") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  ) +
  ggtitle("Jianping projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Jianping_projection_on_PCW3_by_timepoint.svg"), plot = p, width = 8, height = 6)


day9_umap <- subset(combined_umap, timepoint == "day9")

# 画图
p_day9 <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = day9_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   
  ) +
  ggtitle("Day9 projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Day9_projection_on_PCW3.png"), plot = p_day9, width = 8, height = 6)
#------------------------------------------------------------------------------
Barcelona_mat <- readMM(file.path(DATA_DIR, "Barcelona_somite", "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz"))
Barcelona_genes <- read.table(file.path(DATA_DIR, "Barcelona_somite", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
Barcelona_barcodes <- read.table(file.path(DATA_DIR, "Barcelona_somite", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(Barcelona_genes) != nrow(Barcelona_mat)) stop("not matched")
if (nrow(Barcelona_barcodes) != ncol(Barcelona_mat)) stop("not matched")


rownames(Barcelona_mat) <- make.unique(Barcelona_genes$V2)
colnames(Barcelona_mat) <- Barcelona_barcodes$V1
Barcelona_seurat <- CreateSeuratObject(counts = Barcelona_mat, project = "Barcelona")
Barcelona_seurat[["percent.mt"]] <- PercentageFeatureSet(Barcelona_seurat, pattern = "^MT-")
Barcelona_seurat <- subset(Barcelona_seurat, subset = nFeature_RNA > 4000 & percent.mt < 15)

ribosomal_genes_Barcelona <- grep("^RPS|^RPL", rownames(Barcelona_seurat), value = TRUE)
Barcelona_seurat <- subset(Barcelona_seurat, features = setdiff(rownames(Barcelona_seurat), ribosomal_genes_Barcelona))
saveRDS(Barcelona_seurat, file.path(RESULTS_DIR, "Barcelona_seurat_2.rds"))


Barcelona_seurat$sample <- "Barcelona"
subset_PCW3$sample <- "PCW3"
dim(Barcelona_seurat)

Barcelona_seurat <- NormalizeData(Barcelona_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
genes_use_Barcelona <- rownames(Barcelona_seurat)[Matrix::rowSums(Barcelona_seurat@assays$RNA$counts > 0) >= 3]
Barcelona_seurat <- subset(Barcelona_seurat, features = genes_use_Barcelona)
Barcelona_seurat <- FindVariableFeatures(Barcelona_seurat, selection.method = "vst", nfeatures = 2000)
Barcelona_seurat <- ScaleData(Barcelona_seurat)
Barcelona_seurat <- RunPCA(Barcelona_seurat, npcs = 30)
Barcelona_anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = Barcelona_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
Barcelona_seurat <- MapQuery(
  anchorset = Barcelona_anchors, 
  reference = PCW3_seurat, 
  query = Barcelona_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)




Barcelona_umap <- Embeddings(Barcelona_seurat, "ref.umap")
Barcelona_umap <- as.data.frame(Barcelona_umap)
colnames(Barcelona_umap) <- c("UMAP_1", "UMAP_2") 

p_Barcelona <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = Barcelona_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"  
  ) +
  ggtitle("Barcelona projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Barcelona_projection_on_PCW3.png"), plot = p_Barcelona, width = 8, height = 6)
#--------------------------------------------------------------------------------
HNT_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_HNT", "filtered", "matrix.mtx.gz"))
HNT_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_HNT", "filtered", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
HNT_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_HNT", "filtered", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(HNT_genes) != nrow(HNT_mat)) stop("not matched")
if (nrow(HNT_barcodes) != ncol(HNT_mat)) stop("not matched")


rownames(HNT_mat) <- make.unique(HNT_genes$V2)
colnames(HNT_mat) <- HNT_barcodes$V1
HNT_seurat <- CreateSeuratObject(counts = HNT_mat, project = "Barcelona")
HNT_seurat[["percent.mt"]] <- PercentageFeatureSet(HNT_seurat, pattern = "^MT-")
library(ggplot2)
P <- VlnPlot(HNT_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(RESULTS_DIR, "HNT_QC_violin.png"), plot = P, width = 12)
HNT_seurat <- subset(HNT_seurat, subset = nFeature_RNA > 2500 & percent.mt < 20)

ribosomal_genes_HNT <- grep("^RPS|^RPL", rownames(HNT_seurat), value = TRUE)
HNT_seurat <- subset(HNT_seurat, features = setdiff(rownames(HNT_seurat), ribosomal_genes_HNT))
saveRDS(HNT_seurat, file.path(RESULTS_DIR, "HNT_seurat_2.rds"))

HNT_seurat <- NormalizeData(HNT_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
genes_use_HNT <- rownames(HNT_seurat)[Matrix::rowSums(HNT_seurat@assays$RNA$counts > 0) >= 3]
HNT_seurat <- subset(HNT_seurat, features = genes_use_HNT)
HNT_seurat <- FindVariableFeatures(HNT_seurat, selection.method = "vst", nfeatures = 2000)
HNT_seurat <- ScaleData(HNT_seurat)
HNT_seurat <- RunPCA(HNT_seurat, npcs = 30)
HNT_anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = HNT_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
HNT_seurat <- MapQuery(
  anchorset = HNT_anchors, 
  reference = PCW3_seurat, 
  query = HNT_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)


HNT_umap <- Embeddings(HNT_seurat, "ref.umap")
HNT_umap <- as.data.frame(HNT_umap)
colnames(HNT_umap) <- c("UMAP_1", "UMAP_2") 

p_HNT <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = HNT_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  ) +
  ggtitle("HNT projection on PCW3")

ggsave(file.path(RESULTS_DIR, "HNT_projection_on_PCW3.png"), plot = p_HNT, width = 8, height = 6)
#---------------------------------------------------------------------------------
Bar1_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar1", "matrix.mtx.gz"))
Bar1_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar1", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
Bar1_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar1", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(Bar1_genes) != nrow(Bar1_mat)) stop("not matched")
if (nrow(Bar1_barcodes) != ncol(Bar1_mat)) stop("not matched")


rownames(Bar1_mat) <- make.unique(Bar1_genes$V1)
colnames(Bar1_mat) <- Bar1_barcodes$V1
Bar1_seurat <- CreateSeuratObject(counts = Bar1_mat, project = "NT")
Bar1_seurat[["percent.mt"]] <- PercentageFeatureSet(Bar1_seurat, pattern = "^MT-")
Bar1_seurat <- subset(Bar1_seurat, subset = nFeature_RNA > 2500 & percent.mt < 20)


ribosomal_genes_Bar1 <- grep("^RPS|^RPL", rownames(Bar1_seurat), value = TRUE)
Bar1_seurat <- subset(Bar1_seurat, features = setdiff(rownames(Bar1_seurat), ribosomal_genes_Bar1))
saveRDS(Bar1_seurat, file.path(RESULTS_DIR, "Bar1_seurat_2.rds"))

# TODO: In the original analysis script, Bar2_seurat is used here before
# it is created below. Move this mapping block after the Bar2 input/QC block
# once you confirm the intended analysis order.
Bar2_seurat$sample <- "Bar2"
subset_PCW3$sample <- "PCW3"


Bar2_seurat <- NormalizeData(Bar2_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
Bar2_seurat <- FindVariableFeatures(Bar2_seurat, selection.method = "vst", nfeatures = 2000)
Bar2_seurat <- ScaleData(Bar2_seurat)
Bar2_seurat <- RunPCA(Bar2_seurat, npcs = 30)
anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = Bar2_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
Bar2_seurat <- MapQuery(
  anchorset = anchors, 
  reference = PCW3_seurat, 
  query = Bar2_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)

Bar2_umap <- Embeddings(Bar2_seurat, "ref.umap")
Bar2_umap <- as.data.frame(Bar2_umap)
colnames(Bar2_umap) <- c("UMAP_1", "UMAP_2") 

p_Bar2 <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = Bar2_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   
  ) +
  ggtitle("Bar2 projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Bar2_projection_on_PCW3.png"), plot = p_Bar2, width = 8, height = 6)
#--------------------------------------------------------------------------------
Bar2_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar2", "matrix.mtx.gz"))
Bar2_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar2", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
Bar2_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "filtered_feature_bc_matrix_S03", "Bar2", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(Bar2_genes) != nrow(Bar2_mat)) stop("not matched")
if (nrow(Bar2_barcodes) != ncol(Bar2_mat)) stop("not matched")


rownames(Bar2_mat) <- make.unique(Bar2_genes$V1)
colnames(Bar2_mat) <- Bar2_barcodes$V1
Bar2_seurat <- CreateSeuratObject(counts = Bar2_mat, project = "NT")
Bar2_seurat[["percent.mt"]] <- PercentageFeatureSet(Bar2_seurat, pattern = "^MT-")
Bar2_seurat <- subset(Bar2_seurat, subset = nFeature_RNA > 2500 & percent.mt < 20)
saveRDS(Bar2_seurat, file.path(RESULTS_DIR, "Bar2_seurat_2.rds"))


ribosomal_genes_Bar2 <- grep("^RPS|^RPL", rownames(Bar2_seurat), value = TRUE)
Bar2_seurat <- subset(Bar2_seurat, features = setdiff(rownames(Bar2_seurat), ribosomal_genes_Bar2))
saveRDS(Bar2_seurat, file.path(RESULTS_DIR, "Bar2_seurat_2.rds"))


# TODO: Bar3_seurat is not created or loaded anywhere in the supplied script.
# Add the Bar3 input/QC code (or readRDS call) before running this section.
Bar3_seurat <- NormalizeData(Bar3_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
Bar3_seurat <- FindVariableFeatures(Bar3_seurat, selection.method = "vst", nfeatures = 2000)
Bar3_seurat <- ScaleData(Bar3_seurat)
Bar3_seurat <- RunPCA(Bar3_seurat, npcs = 30)
anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = Bar3_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
Bar3_seurat <- MapQuery(
  anchorset = anchors, 
  reference = PCW3_seurat, 
  query = Bar3_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)

Bar3_umap <- Embeddings(Bar3_seurat, "ref.umap")
Bar3_umap <- as.data.frame(Bar3_umap)
colnames(Bar3_umap) <- c("UMAP_1", "UMAP_2") 

p_Bar3 <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = Bar3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   
  ) +
  ggtitle("Bar3 projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Bar3_projection_on_PCW3.png"), plot = p_Bar3, width = 8, height = 6)
#----------------------------------------------------------------------------------
h48_mat <- readMM(file.path(DATA_DIR, "hAxioloid_48h", "outs", "filtered_feature_bc_matrix_48h", "matrix.mtx.gz"))
h48_genes <- read.table(file.path(DATA_DIR, "hAxioloid_48h", "outs", "filtered_feature_bc_matrix_48h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
h48_barcodes <- read.table(file.path(DATA_DIR, "hAxioloid_48h", "outs", "filtered_feature_bc_matrix_48h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(h48_genes) != nrow(h48_mat)) {
  stop("not matched")
}
if (nrow(h48_barcodes) != ncol(h48_mat)) {
  stop("not matched")
}

rownames(h48_mat) <- h48_genes$V2
colnames(h48_mat) <- h48_barcodes$V1
rownames(h48_mat) <- make.unique(rownames(h48_mat))

h48_cantas <- CreateSeuratObject(counts = h48_mat, project = "Cantas")
h48_cantas[["percent.mt"]] <- PercentageFeatureSet(h48_cantas, pattern = "^MT-")
h48_cantas <- subset(h48_cantas, subset = nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt > 2 & percent.mt < 15)
saveRDS(h48_cantas, file.path(RESULTS_DIR, "day2_cantas_seurat.rds"))

h72_mat <- readMM(file.path(DATA_DIR, "hAxioloid_72h", "outs", "filtered_feature_bc_matrix_72h", "matrix.mtx.gz"))
h72_genes <- read.table(file.path(DATA_DIR, "hAxioloid_72h", "outs", "filtered_feature_bc_matrix_72h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
h72_barcodes <- read.table(file.path(DATA_DIR, "hAxioloid_72h", "outs", "filtered_feature_bc_matrix_72h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(h72_genes) != nrow(h72_mat)) {
  stop("not matched")
}
if (nrow(h72_barcodes) != ncol(h72_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(h72_mat) <- h72_genes$V2
colnames(h72_mat) <- h72_barcodes$V1
rownames(h72_mat) <- make.unique(rownames(h72_mat))
# 创建 Seurat 对象
h72_cantas <- CreateSeuratObject(counts = h72_mat, project = "Cantas")
h72_cantas[["percent.mt"]] <- PercentageFeatureSet(h72_cantas, pattern = "^MT-")
h72_cantas <- subset(h72_cantas, subset = nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt > 2 & percent.mt < 15)
saveRDS(h72_cantas, file.path(RESULTS_DIR, "day3_cantas_seurat.rds"))
h96_mat <- readMM(file.path(DATA_DIR, "hAxioloid_96h", "outs", "filtered_feature_bc_matrix_96h", "matrix.mtx.gz"))
h96_genes <- read.table(file.path(DATA_DIR, "hAxioloid_96h", "outs", "filtered_feature_bc_matrix_96h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
h96_barcodes <- read.table(file.path(DATA_DIR, "hAxioloid_96h", "outs", "filtered_feature_bc_matrix_96h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 确保维度匹配
if (nrow(h96_genes) != nrow(h96_mat)) {
  stop("基因数量与矩阵行数不匹配")
}
if (nrow(h96_barcodes) != ncol(h96_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(h96_mat) <- h96_genes$V2
colnames(h96_mat) <- h96_barcodes$V1
rownames(h96_mat) <- make.unique(rownames(h96_mat))
# 创建 Seurat 对象
h96_cantas <- CreateSeuratObject(counts = h96_mat, project = "Cantas")
h96_cantas[["percent.mt"]] <- PercentageFeatureSet(h96_cantas, pattern = "^MT-")
h96_cantas <- subset(h96_cantas, subset = nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt > 2 & percent.mt < 15)
saveRDS(h96_cantas, file.path(RESULTS_DIR, "day4_cantas_seurat.rds"))
h120_mat <- readMM(file.path(DATA_DIR, "hAxioloid_120h", "outs", "filtered_feature_bc_matrix_120h", "matrix.mtx.gz"))
h120_genes <- read.table(file.path(DATA_DIR, "hAxioloid_120h", "outs", "filtered_feature_bc_matrix_120h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
h120_barcodes <- read.table(file.path(DATA_DIR, "hAxioloid_120h", "outs", "filtered_feature_bc_matrix_120h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 确保维度匹配
if (nrow(h120_genes) != nrow(h120_mat)) {
  stop("基因数量与矩阵行数不匹配")
}
if (nrow(h120_barcodes) != ncol(h120_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(h120_mat) <- h120_genes$V2
colnames(h120_mat) <- h120_barcodes$V1
rownames(h120_mat) <- make.unique(rownames(h120_mat))
# 创建 Seurat 对象
h120_cantas <- CreateSeuratObject(counts = h120_mat, project = "Cantas")
h120_cantas[["percent.mt"]] <- PercentageFeatureSet(h120_cantas, pattern = "^MT-")
h120_cantas <- subset(h120_cantas, subset = nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt > 2 & percent.mt < 15)
saveRDS(h120_cantas, file.path(RESULTS_DIR, "day5_cantas_seurat.rds"))
h48_cantas$timepoint <- "48h"
h72_cantas$timepoint <- "72h"
h96_cantas$timepoint <- "96h"
h120_cantas$timepoint <- "120h"

# -------- 合并 + QC + 降维 --------
combined_cantas <- merge(
  h48_cantas,
  y = list(h72_cantas, h96_cantas, h120_cantas),
  add.cell.ids = c("48h", "72h", "96h", "120h"),
  project = "combined_cantas"
)
dim(combined_cantas)
# 可视化用降维（不影响 RNA velocity）
combined_cantas <- NormalizeData(combined_cantas)
combined_cantas <- FindVariableFeatures(combined_cantas)
combined_cantas <- ScaleData(combined_cantas)
combined_cantas <- RunPCA(combined_cantas)
# 4. 计算数据集迁移的锚点
cantas_anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = combined_cantas,
  reference.reduction = "pca",
  dims = 1:30
)
# 5. 在 PCW3 的 UMAP 坐标上投射 CFS 细胞
combined_cantas <- MapQuery(
  anchorset = cantas_anchors, 
  reference = PCW3_seurat, 
  query = combined_cantas, 
  refdata = list(celltype = PCW3_seurat$cell_type), # 迁移细胞类型
  reference.reduction = "pca", 
  reduction.model = "umap"
)
# 过滤预测标签占比 < 5% 的细胞
pred_table <- table(combined_cantas$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
combined_cantas <- subset(combined_cantas, subset = predicted.celltype %in% valid_labels)
# 获取 PCW3 的 UMAP 坐标
# PCW3 灰色背景

# 获取细胞投射后的 UMAP 坐标
combined_cantas_umap <- Embeddings(combined_cantas, "ref.umap")
combined_cantas_umap <- as.data.frame(combined_cantas_umap)
colnames(combined_cantas_umap) <- c("UMAP_1", "UMAP_2") 
combined_cantas_umap$celltype <- combined_cantas$predicted.celltype
table(combined_cantas$predicted.celltype)
combined_cantas_umap$time <- combined_cantas$timepoint
combined_cantas_umap$cell <- rownames(combined_cantas_umap)
combined_cantas_umap$timepoint <- gsub("_.*", "", combined_cantas_umap$cell)
table(combined_cantas_umap$timepoint)

time_colors <- setNames(
  c("#32CD32", "#9370DB", "#ff0019ff", "#87CEEB"),  # 玫红 & 天蓝
  c("48h", "72h", "96h", "120h")
)
# ---------------------------------------------------------
# 绘图
library(ggplot2)
p <- ggplot() +
  geom_point(data = PCW3_umap, aes(x = UMAP_1, y = UMAP_2), 
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = combined_cantas_umap, 
             aes(x = UMAP_1, y = UMAP_2, color = timepoint), 
             alpha = 0.5, size = 0.4) +
  scale_color_manual(values = time_colors, name = "Timepoint") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  ) +
  ggtitle("Cantas projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Cantas_projection_on_PCW3_by_timepoint.png"), plot = p, width = 8, height = 6)
cantas_120h_umap <- subset(combined_cantas_umap, timepoint == "120h")

# 画图
p_120h <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = naomi_120h_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   # 不需要图例
  ) +
  ggtitle("Cantas_120h projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Cantas_120h_projection_on_PCW3.png"), plot = p_120h, width = 8, height = 6)
#-----------------------------------------------------------------------------------
day3_mat <- readMM(file.path(DATA_DIR, "day3_24h_analysis", "outs", "filtered_feature_bc_matrix_day3_24h", "matrix.mtx.gz"))
genes_day3 <- read.table(file.path(DATA_DIR, "day3_24h_analysis", "outs", "filtered_feature_bc_matrix_day3_24h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcodes_day3 <- read.table(file.path(DATA_DIR, "day3_24h_analysis", "outs", "filtered_feature_bc_matrix_day3_24h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(genes_day3) != nrow(day3_mat)) stop("Day3: gene行数不匹配")
if (nrow(barcodes_day3) != ncol(day3_mat)) stop("Day3: barcode列数不匹配")

rownames(day3_mat) <- make.unique(genes_day3$V2)
colnames(day3_mat) <- barcodes_day3$V1

day3_seurat <- CreateSeuratObject(counts = day3_mat, project = "Tiago")
day3_seurat$timepoint <- "day3_24h"
saveRDS(day3_seurat, file = file.path(RESULTS_DIR, "day3_tiago_seurat.rds"))
# --- Day5_24h ---
day5_mat <- readMM(file.path(DATA_DIR, "day5_24h_analysis", "outs", "filtered_feature_bc_matrix_day5_24h", "matrix.mtx.gz"))
genes_day5 <- read.table(file.path(DATA_DIR, "day5_24h_analysis", "outs", "filtered_feature_bc_matrix_day5_24h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcodes_day5 <- read.table(file.path(DATA_DIR, "day5_24h_analysis", "outs", "filtered_feature_bc_matrix_day5_24h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(genes_day5) != nrow(day5_mat)) stop("Day5: gene行数不匹配")
if (nrow(barcodes_day5) != ncol(day5_mat)) stop("Day5: barcode列数不匹配")

rownames(day5_mat) <- make.unique(genes_day5$V2)
colnames(day5_mat) <- barcodes_day5$V1

day5_seurat <- CreateSeuratObject(counts = day5_mat, project = "Tiago")
day5_seurat$timepoint <- "day5_24h"
saveRDS(day5_seurat, file = file.path(RESULTS_DIR, "day5_tiago_seurat.rds"))
# --- Day7_24h ---
day7_mat <- readMM(file.path(DATA_DIR, "day7_24h_analysis", "outs", "filtered_feature_bc_matrix_day7_24h", "matrix.mtx.gz"))
genes_day7 <- read.table(file.path(DATA_DIR, "day7_24h_analysis", "outs", "filtered_feature_bc_matrix_day7_24h", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
barcodes_day7 <- read.table(file.path(DATA_DIR, "day7_24h_analysis", "outs", "filtered_feature_bc_matrix_day7_24h", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)

if (nrow(genes_day7) != nrow(day7_mat)) stop("Day7: gene行数不匹配")
if (nrow(barcodes_day7) != ncol(day7_mat)) stop("Day7: barcode列数不匹配")

rownames(day7_mat) <- make.unique(genes_day7$V2)
colnames(day7_mat) <- barcodes_day7$V1

day7_seurat <- CreateSeuratObject(counts = day7_mat, project = "Tiago")
day7_seurat$timepoint <- "day7_24h"
saveRDS(day7_seurat, file = file.path(RESULTS_DIR, "day7_tiago_seurat.rds"))
combined_tiago <- merge(
  day3_seurat,
  y = list(day5_seurat, day7_seurat),
  add.cell.ids = c("day3", "day5", "day7"),
  project = "combined_tiago"
)
# -------------------- QC + 降维 --------------------
combined_tiago[["percent.mt"]] <- PercentageFeatureSet(combined_tiago, pattern = "^MT-")
VlnPlot(
  combined_tiago,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

combined_tiago <- subset(combined_tiago, subset = nFeature_RNA > 4000 & nFeature_RNA < 8000 & percent.mt < 15)

combined_tiago <- NormalizeData(combined_tiago)
combined_tiago <- FindVariableFeatures(combined_tiago, selection.method = "vst", nfeatures = 2000)


combined_tiago <- ScaleData(combined_tiago)
combined_tiago <- RunPCA(combined_tiago)


dim(combined_tiago)
# -------------------- 投射到 PCW3 UMAP --------------------
anchors_tiago <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = combined_tiago,
  reference.reduction = "pca",
  dims = 1:30
)

combined_tiago <- MapQuery(
  anchorset = anchors_tiago,
  reference = PCW3_seurat,
  query = combined_tiago,
  refdata = list(celltype = PCW3_seurat$cell_type),
  reference.reduction = "pca",
  reduction.model = "umap"
)

# -------------------- 过滤预测标签占比 < 3% 的细胞 --------------------
pred_table <- table(combined_tiago$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
combined_tiago <- subset(combined_tiago, subset = predicted.celltype %in% valid_labels)



# 投射后的 Tiago 细胞
# 获取细胞投射后的 UMAP 坐标
tiago_umap <- Embeddings(combined_tiago, "ref.umap")
tiago_umap <- as.data.frame(tiago_umap)
colnames(tiago_umap) <- c("UMAP_1", "UMAP_2") 
tiago_umap$celltype <- combined_tiago$predicted.celltype
table(combined_tiago$predicted.celltype)
tiago_umap$time <- combined_tiago$timepoint
tiago_umap$cell <- rownames(tiago_umap)
tiago_umap$timepoint <- gsub("_.*", "", tiago_umap$cell)
table(tiago_umap$timepoint)

time_colors <- setNames(
  c("#9370DB", "#ff0019ff", "#87CEEB"),  # 玫红 & 天蓝
  c("day3", "day5", "day7")
)
# ---------------------------------------------------------
# 绘图
library(ggplot2)
p <- ggplot() +
  geom_point(data = PCW3_umap, aes(x = UMAP_1, y = UMAP_2), 
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = tiago_umap, 
             aes(x = UMAP_1, y = UMAP_2, color = timepoint), 
             alpha = 0.5, size = 0.4) +
  scale_color_manual(values = time_colors, name = "Timepoint") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  ) +
  ggtitle("Tiago projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Tiago_projection_on_PCW3_by_timepoint.png"), plot = p, width = 8, height = 6)

tiago_day7_umap <- subset(tiago_umap, timepoint == "day7")
p_day7 <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = tiago_day7_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   # 不需要图例
  ) +
  ggtitle("Tiago_day7_projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Tiago_day7_projection_on_PCW3.png"), plot = p_day7, width = 8, height = 6)
#-----------------------------------------------------------------------------------
trunk_mat <- readMM(file.path(DATA_DIR, "trunk", "outs", "filtered_feature_bc_matrix_trunk", "matrix.mtx.gz"))
trunk_genes <- read.table(file.path(DATA_DIR, "trunk", "outs", "filtered_feature_bc_matrix_trunk", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
trunk_barcodes <- read.table(file.path(DATA_DIR, "trunk", "outs", "filtered_feature_bc_matrix_trunk", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 维度检查
if (nrow(trunk_genes) != nrow(trunk_mat)) stop("基因数量与矩阵行数不匹配")
if (nrow(trunk_barcodes) != ncol(trunk_mat)) stop("细胞数量与矩阵列数不匹配")

# 设置行名列名
rownames(trunk_mat) <- make.unique(trunk_genes$V2)
colnames(trunk_mat) <- trunk_barcodes$V1
trunk_seurat <- CreateSeuratObject(counts = trunk_mat, project = "trunk")
dim(trunk_seurat)
trunk_seurat[["percent.mt"]] <- PercentageFeatureSet(trunk_seurat, pattern = "^MT-")
trunk_seurat <- subset(trunk_seurat, subset = nFeature_RNA > 2500 & nFeature_RNA < 14000 & percent.mt < 15)

# 去除核糖体基因
ribosomal_genes_trunk <- grep("^RPS|^RPL", rownames(trunk_seurat), value = TRUE)
trunk_seurat <- subset(trunk_seurat, features = setdiff(rownames(trunk_seurat), ribosomal_genes_trunk))
saveRDS(trunk_seurat, file.path(RESULTS_DIR, "trunk_seurat_2.rds"))
# 预处理与 PCA
trunk_seurat <- NormalizeData(trunk_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
trunk_seurat <- FindVariableFeatures(trunk_seurat, selection.method = "vst", nfeatures = 2000)
trunk_seurat <- ScaleData(trunk_seurat)
trunk_seurat <- RunPCA(trunk_seurat, npcs = 30)
anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = trunk_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
trunk_seurat <- MapQuery(
  anchorset = anchors, 
  reference = PCW3_seurat, 
  query = trunk_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)
pred_table <- table(trunk_seurat$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
trunk_seurat <- subset(trunk_seurat, subset = predicted.celltype %in% valid_labels)


# 获取细胞投射后的 UMAP 坐标
trunk_umap <- Embeddings(trunk_seurat, "ref.umap")
trunk_umap <- as.data.frame(trunk_umap)
colnames(trunk_umap) <- c("UMAP_1", "UMAP_2") 
trunk_umap$celltype <- trunk_seurat$predicted.celltype
table(trunk_seurat$predicted.celltype)
p_trunk <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = trunk_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   # 不需要图例
  ) +
  ggtitle("trunk projection on PCW3")

ggsave(file.path(RESULTS_DIR, "trunk_projection_on_PCW3.png"), plot = p_trunk, width = 8, height = 6)
#---------------------------------------------------------------------------------
hb72_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "trunk2", "72h_trunk2", "filtered", "matrix.mtx.gz"))
hb72_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "72h_trunk2", "filtered", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
hb72_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "72h_trunk2", "filtered", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 确保维度匹配
if (nrow(hb72_genes) != nrow(hb72_mat)) {
  stop("基因数量与矩阵行数不匹配")
}
if (nrow(hb72_barcodes) != ncol(hb72_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(hb72_mat) <- hb72_genes$V2
colnames(hb72_mat) <- hb72_barcodes$V1
rownames(hb72_mat) <- make.unique(rownames(hb72_mat))
# 创建 Seurat 对象
hb72_naomi <- CreateSeuratObject(counts = hb72_mat, project = "Naomi")
hb72_naomi[["percent.mt"]] <- PercentageFeatureSet(hb72_naomi, pattern = "^MT-")
hb72_naomi <- subset(hb72_naomi, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt > 2 & percent.mt < 15)
saveRDS(hb72_naomi, file.path(RESULTS_DIR, "day3_naomi_seurat.rds"))

hb96_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "trunk2", "96h_trunk2", "filtered", "matrix.mtx.gz"))
hb96_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "96h_trunk2", "filtered", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
hb96_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "96h_trunk2", "filtered", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 确保维度匹配
if (nrow(hb96_genes) != nrow(hb96_mat)) {
  stop("基因数量与矩阵行数不匹配")
}
if (nrow(hb96_barcodes) != ncol(hb96_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(hb96_mat) <- hb96_genes$V2
colnames(hb96_mat) <- hb96_barcodes$V1
rownames(hb96_mat) <- make.unique(rownames(hb96_mat))
# 创建 Seurat 对象
hb96_naomi <- CreateSeuratObject(counts = hb96_mat, project = "Naomi")
hb96_naomi[["percent.mt"]] <- PercentageFeatureSet(hb96_naomi, pattern = "^MT-")
hb96_naomi <- subset(hb96_naomi, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt > 2 & percent.mt < 15)
saveRDS(hb96_naomi, file.path(RESULTS_DIR, "day4_naomi_seurat.rds"))
hb120_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "trunk2", "120h_trunk2", "filtered", "matrix.mtx.gz"))
hb120_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "120h_trunk2", "filtered", "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
hb120_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "trunk2", "120h_trunk2", "filtered", "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 确保维度匹配
if (nrow(hb120_genes) != nrow(hb120_mat)) {
  stop("基因数量与矩阵行数不匹配")
}
if (nrow(hb120_barcodes) != ncol(hb120_mat)) {
  stop("细胞数量与矩阵列数不匹配")
}
# 设置行名和列名
rownames(hb120_mat) <- hb120_genes$V2
colnames(hb120_mat) <- hb120_barcodes$V1
rownames(hb120_mat) <- make.unique(rownames(hb120_mat))
# 创建 Seurat 对象
hb120_naomi <- CreateSeuratObject(counts = hb120_mat, project = "Naomi")
hb120_naomi[["percent.mt"]] <- PercentageFeatureSet(hb120_naomi, pattern = "^MT-")
hb120_naomi <- subset(hb120_naomi, subset = nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt > 2 & percent.mt < 15)
saveRDS(hb120_naomi, file.path(RESULTS_DIR, "day5_naomi_seurat.rds"))

hb72_naomi$timepoint <- "72h"
hb96_naomi$timepoint <- "96h"
hb120_naomi$timepoint <- "120h"

# -------- 合并 + QC + 降维 --------
combined_naomi <- merge(
  hb72_naomi,
  y = list(hb96_naomi, hb120_naomi),
  add.cell.ids = c("72h", "96h", "120h"),
  project = "combined_naomi"
)

# 可视化用降维（不影响 RNA velocity）
combined_naomi <- NormalizeData(combined_naomi)
combined_naomi <- FindVariableFeatures(combined_naomi)
combined_naomi <- ScaleData(combined_naomi)
combined_naomi <- RunPCA(combined_naomi)
# 4. 计算数据集迁移的锚点
naomi_anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = combined_naomi,
  reference.reduction = "pca",
  dims = 1:30
)
# 5. 在 PCW3 的 UMAP 坐标上投射 CFS 细胞
combined_naomi <- MapQuery(
  anchorset = naomi_anchors, 
  reference = PCW3_seurat, 
  query = combined_naomi, 
  refdata = list(celltype = PCW3_seurat$cell_type), # 迁移细胞类型
  reference.reduction = "pca", 
  reduction.model = "umap"
)
# 过滤预测标签占比 < 5% 的细胞
pred_table <- table(combined_naomi$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
combined_naomi <- subset(combined_naomi, subset = predicted.celltype %in% valid_labels)
# 获取 PCW3 的 UMAP 坐标
# PCW3 灰色背景

# 获取细胞投射后的 UMAP 坐标
combined_naomi_umap <- Embeddings(combined_naomi, "ref.umap")
combined_naomi_umap <- as.data.frame(combined_naomi_umap)
colnames(combined_naomi_umap) <- c("UMAP_1", "UMAP_2") 
combined_naomi_umap$celltype <- combined_naomi$predicted.celltype
table(combined_naomi$predicted.celltype)
combined_naomi_umap$time <- combined_naomi$timepoint
combined_naomi_umap$cell <- rownames(combined_naomi_umap)
combined_naomi_umap$timepoint <- gsub("_.*", "", combined_naomi_umap$cell)
table(combined_naomi_umap$timepoint)

time_colors <- setNames(
  c("#9370DB", "#ff0019ff", "#87CEEB"),  # 玫红 & 天蓝
  c("72h", "96h", "120h")
)
# ---------------------------------------------------------
# 绘图
library(ggplot2)
p <- ggplot() +
  geom_point(data = PCW3_umap, aes(x = UMAP_1, y = UMAP_2), 
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = combined_naomi_umap, 
             aes(x = UMAP_1, y = UMAP_2, color = timepoint), 
             alpha = 0.5, size = 0.4) +
  scale_color_manual(values = time_colors, name = "Timepoint") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  ) +
  ggtitle("Naomi projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Naomi_projection_on_PCW3_by_timepoint.png"), plot = p, width = 8, height = 6)
naomi_120h_umap <- subset(combined_naomi_umap, timepoint == "120h")

# 画图
p_120h <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = naomi_120h_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   # 不需要图例
  ) +
  ggtitle("Naomi_120h projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Naomi_120h_projection_on_PCW3.png"), plot = p_120h, width = 8, height = 6)
#-----------------------------------------------------------------------------------
Yam_mat <- readMM(file.path(DATA_DIR, "all_cellranger", "Olivier", "velocyto_output", "GSM6806916_matrix.mtx.gz"))
Yam_genes <- read.table(file.path(DATA_DIR, "all_cellranger", "Olivier", "velocyto_output", "GSM6806916_features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
Yam_barcodes <- read.table(file.path(DATA_DIR, "all_cellranger", "Olivier", "velocyto_output", "GSM6806916_barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
# 维度检查
if (nrow(Yam_genes) != nrow(Yam_mat)) stop("基因数量与矩阵行数不匹配")
if (nrow(Yam_barcodes) != ncol(Yam_mat)) stop("细胞数量与矩阵列数不匹配")

# 设置行名列名
rownames(Yam_mat) <- make.unique(Yam_genes$V2)
colnames(Yam_mat) <- Yam_barcodes$V1
Yam_seurat <- CreateSeuratObject(counts = Yam_mat, project = "Yam")
Yam_seurat[["percent.mt"]] <- PercentageFeatureSet(Yam_seurat, pattern = "^MT-")
Yam_seurat <- subset(Yam_seurat, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt < 25)
# 去除核糖体基因
ribosomal_genes_Yam <- grep("^RPS|^RPL", rownames(Yam_seurat), value = TRUE)
Yam_seurat <- subset(Yam_seurat, features = setdiff(rownames(Yam_seurat), ribosomal_genes_Yam))
saveRDS(Yam_seurat, file.path(RESULTS_DIR, "Yam_seurat.rds"))



# 预处理与 PCA
Yam_seurat <- NormalizeData(Yam_seurat, normalization.method = "LogNormalize", scale.factor = 1e4)
genes_use_Yam <- rownames(Yam_seurat)[Matrix::rowSums(Yam_seurat@assays$RNA$counts > 0) >= 3]
Yam_seurat <- subset(Yam_seurat, features = genes_use_Yam)
Yam_seurat <- FindVariableFeatures(Yam_seurat, selection.method = "vst", nfeatures = 2000)
Yam_seurat <- ScaleData(Yam_seurat)
Yam_seurat <- RunPCA(Yam_seurat, npcs = 30)
Yam_anchors <- FindTransferAnchors(
  reference = PCW3_seurat,
  query = Yam_seurat,
  reference.reduction = "pca",
  dims = 1:30
)
Yam_seurat <- MapQuery(
  anchorset = Yam_anchors, 
  reference = PCW3_seurat, 
  query = Yam_seurat, 
  refdata = list(celltype = PCW3_seurat$cell_type), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)
# 提取投射后的 UMAP
PCW3_umap <- Embeddings(PCW3_seurat, "umap")
PCW3_umap <- as.data.frame(PCW3_umap)
colnames(PCW3_umap) <- c("UMAP_1", "UMAP_2") 


# 获取细胞投射后的 UMAP 坐标
Yam_umap <- Embeddings(Yam_seurat, "ref.umap")
Yam_umap <- as.data.frame(Yam_umap)
colnames(Yam_umap) <- c("UMAP_1", "UMAP_2") 

p_Yam <- ggplot() +
  geom_point(data = PCW3_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey80", alpha = 0.3, size = 0.3) +
  geom_point(data = Yam_umap, 
             aes(x = UMAP_1, y = UMAP_2),
             color = "#87CEEB", alpha = 0.6, size = 0.4) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none"   # 不需要图例
  ) +
  ggtitle("Yam projection on PCW3")

ggsave(file.path(RESULTS_DIR, "Yam_projection_on_PCW3.svg"), plot = p_Yam, width = 8, height = 6)
ggsave(file.path(RESULTS_DIR, "Yam_projection_on_PCW3.png"), plot = p_Yam, width = 8, height = 6)
pred_table <- table(Yam_seurat$predicted.celltype)
pred_freq <- pred_table / sum(pred_table)
valid_labels <- names(pred_freq[pred_freq >= 0.03])
Yam_seurat <- subset(Yam_seurat, subset = predicted.celltype %in% valid_labels)

table(Yam_seurat$predicted.celltype)
# 计算每个细胞类型的数量
celltype_counts <- table(Yam_seurat$predicted.celltype)

# 计算百分比
celltype_percent <- prop.table(celltype_counts) * 100

# 转成数据框方便查看或导出
celltype_df <- data.frame(
  celltype = names(celltype_percent),
  count = as.numeric(celltype_counts),
  percent = round(as.numeric(celltype_percent), 2)
)

# 按占比从高到低排序
celltype_df <- celltype_df[order(-celltype_df$percent), ]

# 查看结果
print(celltype_df)
