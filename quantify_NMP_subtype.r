# ============================================================
# NMP TBXT/SOX2 expression analysis across axial organoid models
#
# This script:
#   1. Transfers PCW3 reference cell-type labels to organoid datasets.
#   2. Selects cells predicted as NMPs.
#   3. Calculates log1p(TBXT) - log1p(SOX2).
#   4. Compares NMP expression states using ridge plots.
#
# Required input:
#   data/reference/subset_PCW3_1030.rds
#   data/seurat_objects/*.rds
#
# Output:
#   results/NMP_TBXT_SOX2_ridgeplots/*.pdf
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(purrr)
library(tibble)
library(future)

# Run sequentially for reproducibility and memory control
plan(sequential)


# ============================================================
# 2. File paths
# ============================================================

DATA_DIR <- "data"
REFERENCE_DIR <- file.path(DATA_DIR, "reference")
SEURAT_DIR <- file.path(DATA_DIR, "seurat_objects")

RESULTS_DIR <- file.path(
  "results",
  "NMP_TBXT_SOX2_ridgeplots"
)

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Load PCW3 reference
# ============================================================

reference_file <- file.path(
  REFERENCE_DIR,
  "subset_PCW3_1030.rds"
)

subset_PCW3 <- readRDS(reference_file)

Idents(subset_PCW3) <- "cell_type"
DefaultAssay(subset_PCW3) <- "RNA"


message(
  "Reference loaded: ",
  ncol(subset_PCW3),
  " cells"
)


# ============================================================
# 4. Load query Seurat objects
# ============================================================

query_list <- list(

  # Neural tube models
  NT_Xue_D4 = readRDS(
    file.path(SEURAT_DIR, "day4_Jianping_seurat.rds")
  ),

  NT_Xue_D9 = readRDS(
    file.path(SEURAT_DIR, "day9_Jianping_seurat.rds")
  ),


  # Somite-containing models
  SM_YaA_D2 = readRDS(
    file.path(SEURAT_DIR, "day2_cantas_seurat.rds")
  ),

  SM_YaA_D3 = readRDS(
    file.path(SEURAT_DIR, "day3_cantas_seurat.rds")
  ),

  SM_YaA_D4 = readRDS(
    file.path(SEURAT_DIR, "day4_cantas_seurat.rds")
  ),

  SM_YaA_D5 = readRDS(
    file.path(SEURAT_DIR, "day5_cantas_seurat.rds")
  ),


  SM_Mia_D1 = readRDS(
    file.path(SEURAT_DIR, "day1_Segmentoid_seurat.rds")
  ),

  SM_Mia_D2 = readRDS(
    file.path(SEURAT_DIR, "day2_Segmentoid_seurat.rds")
  ),

  SM_Mia_D3 = readRDS(
    file.path(SEURAT_DIR, "day3_Segmentoid_seurat.rds")
  ),

  SM_Mia_D4 = readRDS(
    file.path(SEURAT_DIR, "day4_Segmentoid_seurat.rds")
  ),


  # Trunk models
  TK_Rit_D3 = readRDS(
    file.path(SEURAT_DIR, "day3_tiago_seurat.rds")
  ),

  TK_Rit_D5 = readRDS(
    file.path(SEURAT_DIR, "day5_tiago_seurat.rds")
  ),

  TK_Rit_D7 = readRDS(
    file.path(SEURAT_DIR, "day7_tiago_seurat.rds")
  ),


  TK_Mak_D3 = readRDS(
    file.path(SEURAT_DIR, "day3_naomi_seurat.rds")
  ),

  TK_Mak_D4 = readRDS(
    file.path(SEURAT_DIR, "day4_naomi_seurat.rds")
  ),

  TK_Mak_D5 = readRDS(
    file.path(SEURAT_DIR, "day5_naomi_seurat.rds")
  ),


  TK_Ham_D4 = readRDS(
    file.path(SEURAT_DIR, "D4_RA_merged_seurat.rds")
  ),

  TK_Ham_D5 = readRDS(
    file.path(SEURAT_DIR, "D5_RA_merged_seurat.rds")
  )
)


message(
  "Loaded ",
  length(query_list),
  " query datasets."
)


# ============================================================
# 5. Define plotting groups
# ============================================================

plot_groups <- list(

  NT_Xue = c(
    "Reference_PCW3",
    "NT_Xue_D4",
    "NT_Xue_D9"
  ),

  SM_YaA = c(
    "Reference_PCW3",
    "SM_YaA_D5",
    "SM_YaA_D4",
    "SM_YaA_D3",
    "SM_YaA_D2"
  ),

  SM_Mia = c(
    "Reference_PCW3",
    "SM_Mia_D4",
    "SM_Mia_D3",
    "SM_Mia_D2",
    "SM_Mia_D1"
  ),

  TK_Rit = c(
    "Reference_PCW3",
    "TK_Rit_D7",
    "TK_Rit_D5",
    "TK_Rit_D3"
  ),

  TK_Mak = c(
    "Reference_PCW3",
    "TK_Mak_D5",
    "TK_Mak_D4",
    "TK_Mak_D3"
  ),

  TK_Ham = c(
    "Reference_PCW3",
    "TK_Ham_D5",
    "TK_Ham_D4"
  )
)


# ============================================================
# 6. Label transfer function
# ============================================================

transfer_to_PCW3 <- function(
    query_obj,
    query_name,
    reference_obj
) {

  message(
    "Processing query: ",
    query_name
  )

  DefaultAssay(query_obj) <- "RNA"


  # Standard preprocessing
  query_obj <- NormalizeData(
    query_obj,
    verbose = FALSE
  )

  query_obj <- FindVariableFeatures(
    query_obj,
    verbose = FALSE
  )

  query_obj <- ScaleData(
    query_obj,
    verbose = FALSE
  )

  query_obj <- RunPCA(
    query_obj,
    verbose = FALSE
  )


  # Identify transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference_obj,
    query = query_obj,
    reference.reduction = "pca",
    dims = 1:30
  )


  # Transfer PCW3 cell-type labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_obj$cell_type,
    dims = 1:30
  )


  query_obj <- AddMetaData(
    query_obj,
    metadata = predictions
  )


  query_obj$protocol <- query_name

  return(query_obj)
}


# ============================================================
# 7. Perform label transfer for all query datasets
# ============================================================

transferred_queries <- imap(
  query_list,
  ~ transfer_to_PCW3(
    query_obj = .x,
    query_name = .y,
    reference_obj = subset_PCW3
  )
)


# ============================================================
# 8. Function to calculate TBXT/SOX2 expression difference
# ============================================================

get_log_ratio <- function(
    seurat_obj,
    protocol_name,
    cell_subset = NULL
) {

  DefaultAssay(seurat_obj) <- "RNA"


  if (!is.null(cell_subset)) {

    seurat_obj <- subset(
      seurat_obj,
      cells = cell_subset
    )
  }


  genes <- c(
    "TBXT",
    "SOX2"
  )


  # Check whether required genes are present
  missing_genes <- setdiff(
    genes,
    rownames(seurat_obj)
  )

  if (length(missing_genes) > 0) {

    stop(
      "Missing gene(s) in ",
      protocol_name,
      ": ",
      paste(
        missing_genes,
        collapse = ", "
      )
    )
  }


  # Retrieve normalized expression
  if (
    "data" %in%
    SeuratObject::Layers(
      seurat_obj[["RNA"]]
    )
  ) {

    expr <- FetchData(
      seurat_obj,
      vars = genes,
      layer = "data"
    )

  } else {

    expr <- FetchData(
      seurat_obj,
      vars = genes
    )
  }


  # Expression difference
  expr$log_ratio <- (
    log1p(expr$TBXT)
    - log1p(expr$SOX2)
  )


  expr$cell <- rownames(expr)
  expr$protocol <- protocol_name


  return(
    expr[
      ,
      c(
        "cell",
        "protocol",
        "log_ratio"
      )
    ]
  )
}


# ============================================================
# 9. Calculate expression difference in reference NMPs
# ============================================================

ref_nmp_cells <- WhichCells(
  subset_PCW3,
  idents = "NMP"
)


message(
  "Reference NMP cells: ",
  length(ref_nmp_cells)
)


ref_logratio <- get_log_ratio(
  subset_PCW3,
  protocol_name = "Reference_PCW3",
  cell_subset = ref_nmp_cells
)

ref_logratio$cell_type <- "NMP"


# ============================================================
# 10. Calculate expression difference in predicted query NMPs
# ============================================================

query_logratio_list <- map(
  transferred_queries,
  function(q) {

    nmp_cells <- rownames(
      q@meta.data[
        q@meta.data$predicted.id == "NMP",
        ,
        drop = FALSE
      ]
    )


    if (length(nmp_cells) == 0) {

      message(
        "No predicted NMP cells in ",
        q$protocol[1],
        "; skipping."
      )

      return(NULL)
    }


    message(
      q$protocol[1],
      ": ",
      length(nmp_cells),
      " predicted NMP cells"
    )


    get_log_ratio(
      q,
      protocol_name = q$protocol[1],
      cell_subset = nmp_cells
    ) %>%
      mutate(
        predicted.celltype = "NMP"
      )
  }
)


query_logratio <- bind_rows(
  query_logratio_list
)


# ============================================================
# 11. Combine reference and query NMPs
# ============================================================

log_ratio_all <- bind_rows(
  ref_logratio,
  query_logratio
)


message(
  "Total NMP cells used for ridge plots: ",
  nrow(log_ratio_all)
)


# ============================================================
# 12. Generate ridge plots
# ============================================================

for (group_name in names(plot_groups)) {

  group_protocols <- plot_groups[
    [group_name]
  ]


  df_group <- log_ratio_all %>%

    filter(
      protocol %in% group_protocols
    ) %>%

    mutate(
      protocol = factor(
        protocol,
        levels = group_protocols
      )
    )


  if (nrow(df_group) == 0) {

    message(
      "No cells available for group: ",
      group_name
    )

    next
  }


  output_file <- file.path(
    RESULTS_DIR,
    paste0(
      "Ridge_",
      group_name,
      ".pdf"
    )
  )


  pdf(
    file = output_file,
    width = 6.5,
    height = 4 +
      length(group_protocols) * 0.5
  )


  p <- ggplot(
    df_group,
    aes(
      x = log_ratio,
      y = protocol,
      fill = protocol
    )
  ) +

    geom_density_ridges(
      alpha = 0.7,
      scale = 1.1,
      color = "gray30",
      bandwidth = 0.3
    ) +

    labs(
      title = paste0(
        group_name,
        " — NMP TBXT/SOX2 expression"
      ),
      x = "log1p(TBXT) - log1p(SOX2)",
      y = NULL
    ) +

    coord_cartesian(
      xlim = c(-4, 4)
    ) +

    theme_ridges() +

    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_text(
        size = 11
      ),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 14,
        hjust = 0.5
      )
    )


  print(p)

  dev.off()


  message(
    "Ridge plot saved: ",
    output_file
  )
}


message(
  "All grouped ridge plots saved to: ",
  RESULTS_DIR
)