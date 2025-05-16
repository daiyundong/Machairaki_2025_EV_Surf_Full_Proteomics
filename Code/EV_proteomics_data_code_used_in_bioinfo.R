# -------------------------------------------------------------------
# EV Surface vs Full Proteome – Analytical Workflow (R)
# -------------------------------------------------------------------
# This script imports batch‑1 and batch‑2 proteomics exports, attaches
# clinical metadata, creates unified sample labels, merges the two
# batches (with and without batch‑effect correction), filters proteins
# by missing‑value threshold, performs log2 + median normalisation and a
# differential‑expression test contrasting organoid (Org) versus
# extracellular‑vesicle (EV) samples. Downstream visualisations include a
# Volcano plot and annotation‑rich heatmaps.
# -------------------------------------------------------------------

# ===================================================================
# 1  Load libraries & set working directory
# ===================================================================
library(readxl)
library(ggplot2)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(VennDiagram)
library(patchwork)
library(HPAanalyze)
library(randomForest)
library(caret)
library(pROC)
library(ROSE)
library(ggrepel)
setwd("./Data_used_in_Bioinfo")


# ===================================================================
# 2  Import raw datasets and annotation tables
# ===================================================================
surf.data.list <- list()
full.data.list <- list()

surf.data.list[["raw"]] <- read_xlsx("112524-JHMI-Vasso-EV-surface-A-O-quant.xlsx")
full.data.list[["raw"]] <- read_xlsx("112124-JHMI-Vasso-EV-full-A-O-quant.xlsx")

colnames(surf.data.list[["raw"]])[2] <- 'Gene'
colnames(full.data.list[["raw"]])[2] <- 'Gene'

exocarta_anno <- read.csv("exocarta_top.csv")
exocarta_anno$Gene.Symbol <- toupper(exocarta_anno$Gene.Symbol)
multidata_anno <- read.csv("multi_data_anno.csv")
neuroPro_anno <- read.csv("neuroPro.csv")

braincell_ntpm <- read.csv("genes_nTPM_brain_cell.csv")
cluster_to_celltype <- read.csv("rna_single_nuclei_cluster_type_cluster_types.csv")

# ===================================================================
# 3  Check ExoCarta coverage
# ===================================================================

temp_exocarta <- intersect(full.data.list[["raw"]]$Gene, exocarta_anno$Gene.Symbol)
paste0("full-exocarta: ", length(temp_exocarta))

temp_exocarta <- intersect(surf.data.list[["raw"]]$Gene, exocarta_anno$Gene.Symbol)
paste0("surf-exocarta: ", length(temp_exocarta))

# ===================================================================
# 4  Column pruning & renaming
# ===================================================================

surf.data.list[["clean"]] <- surf.data.list[["raw"]][,c(1,2,13:58)]
full.data.list[["clean"]] <- full.data.list[["raw"]][,c(1,2,13:58)]

temp <- rowSums(is.na(surf.data.list[["clean"]][, 4:ncol(surf.data.list[["clean"]])])) == (ncol(surf.data.list[["clean"]]) - 3)
surf.data.list[["clean"]] <- surf.data.list[["clean"]][!temp, ]

colnames(surf.data.list[["clean"]])[2] <- 'Gene'
colnames(surf.data.list[["clean"]])[4:48] <- c(
  paste0("AD_", rep(LETTERS[1:8], each = 3), 1:3),
  paste0("CL_", rep(LETTERS[9:15], each = 3), 1:3))

temp <- rowSums(is.na(full.data.list[["clean"]][, 4:ncol(full.data.list[["clean"]])])) == (ncol(full.data.list[["clean"]]) - 3)
full.data.list[["clean"]] <- full.data.list[["clean"]][!temp, ]

colnames(full.data.list[["clean"]])[2] <- 'Gene'
colnames(full.data.list[["clean"]])[4:48] <- c(
  paste0("AD_", rep(LETTERS[1:8], each = 3), 1:3),
  paste0("CL_", rep(LETTERS[9:15], each = 3), 1:3))

# ===================================================================
# 5  Missing‑value bar chart helper
# ===================================================================

miss_value_bar_plot <- function(my.clean.data, note, mytitle){
  missing_percent <- rowMeans(is.na(my.clean.data[, 4:ncol(my.clean.data)])) * 100
  plot_data <- data.frame(
    missing = cut(missing_percent,
                  breaks = seq(0, 100, by = 2),
                  include.lowest = TRUE,
                  labels = seq(0, 98, by = 2), 
                  right = FALSE)                
  )

  plot_data$missing <- factor(plot_data$missing, 
                              levels = seq(0, 98, by = 2))

  custom_labels <- function(x) {
    ifelse(x == 98, "100", x)}

  ggplot(plot_data, aes(x = missing)) +
    geom_bar(aes(y = after_stat(count)/sum(after_stat(count))*100),
             fill = "#4DBBD5",
             width = 0.8) +
    labs(x = "Percent of Missing Value", 
         y = "Percent of Total Quantified Proteins",
         title = mytitle) +
    scale_x_discrete(
      breaks = seq(0, 98, by = 20),  
      labels = custom_labels         
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      panel.grid.major.y = element_line(color = "grey90")
    )
  ggsave(paste0(note,'.tiff'),width=6, height = 4)
}

miss_value_bar_plot(surf.data.list[["clean"]], note = 'Suppfigure1a_surf_miss_value_bar', mytitle = "EV Surface")
miss_value_bar_plot(full.data.list[["clean"]], note = 'Suppfigure1b_full_miss_value_bar', mytitle = "Full EV")

# ===================================================================
# 6  Outlier masking (|Z| > 3)
# ===================================================================

removing_outlier_3zscore <- function(my.clean.data){
  expr_data <- my.clean.data[, 4:ncol(my.clean.data)]
  zscore_matrix <- t(scale(t(expr_data)))
  outlier_logical <- abs(zscore_matrix) > 3
  expr_data[outlier_logical] <- NA
  nooutlier.data <- cbind(my.clean.data[,1:3],expr_data)
  
  return(nooutlier.data)
}
surf.data.list[["nooutlier"]] <- removing_outlier_3zscore(surf.data.list[["clean"]])
full.data.list[["nooutlier"]] <- removing_outlier_3zscore(full.data.list[["clean"]])

# ===================================================================
# 7  Missing‑value row filter (≥ 80 % in AD or CL)
# ===================================================================

filter_na_rows_improve <- function(my.nooutlier.data, na_threshold = 0.8) {
  ad_cols <- grep("^AD_", colnames(my.nooutlier.data))  
  cl_cols <- grep("^CL_", colnames(my.nooutlier.data)) 
  
  frac_ad_nonzero <- apply(my.nooutlier.data[, ad_cols, drop = FALSE], 1, function(x) {
    mean(x != 0 & !is.na(x))
  })
  
  frac_cl_nonzero <- apply(my.nooutlier.data[, cl_cols, drop = FALSE], 1, function(x) {
    mean(x != 0 & !is.na(x))
  })
  
  keep_rows <- frac_ad_nonzero > na_threshold | frac_cl_nonzero > na_threshold
  filtered_data <- my.nooutlier.data[keep_rows, , drop = FALSE]
  
  return(filtered_data)
}
surf.data.list[["filtermiss"]] <- filter_na_rows_improve(surf.data.list[["nooutlier"]], na_threshold = 0.8)
full.data.list[["filtermiss"]] <- filter_na_rows_improve(full.data.list[["nooutlier"]], na_threshold = 0.8)


# ===================================================================
# 8  Log2 transform + median normalisation
# ===================================================================
norm_transform <- function(my.filtermiss.data){
  metadata_cols <- my.filtermiss.data[, 1:3]
  quant_cols <- my.filtermiss.data[, -c(1:3)]  
  log2_data <- log2(quant_cols + 1)
  
  col_medians <- apply(log2_data, 2, median, na.rm=TRUE)
  norm_data <- sweep(log2_data, 2, col_medians, FUN="-")

  final_data <- cbind(metadata_cols, norm_data)
  
  return(final_data)
}
surf.data.list[["norm"]] <- norm_transform(surf.data.list[["filtermiss"]])
full.data.list[["norm"]] <- norm_transform(full.data.list[["filtermiss"]])

# ===================================================================
# 9  Enrichment analysis (GO/KEGG)
# ===================================================================

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

GO_KEGG_analysis <- function(genelist){
  
  GO_gene <- bitr(genelist,fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
  
  GO <-enrichGO( GO_gene$ENTREZID, 
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05, 
                 readable = T)
  
  KEGG <-enrichKEGG(GO_gene$ENTREZID,
                    organism = KEGG_database,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
  
  
  return(list(gene=GO_gene, GO=GO, KEGG=KEGG))
}

surf.data.list[["ori_GO_KEGG"]] <- GO_KEGG_analysis(surf.data.list[["raw"]]$Accession)
full.data.list[["ori_GO_KEGG"]] <- GO_KEGG_analysis(full.data.list[["raw"]]$Accession)

# ===================================================================
# 10  Venn diagram (full vs surface)
# ===================================================================

dev.off()
raw_data_Venn_plot <- function(){
  full_accession <- full.data.list[["raw"]]$Accession
  surf_accession <- surf.data.list[["raw"]]$Accession
  
  full_set <- unique(full_accession)
  surf_set <- unique(surf_accession)
  
  area_full <- length(full_set)
  area_surf <- length(surf_set)
  cross_area <- length(intersect(full_set, surf_set))
  
  venn.plot <- draw.pairwise.venn(
    area1      = area_full,          
    area2      = area_surf,          
    cross.area = cross_area,         
    category   = c("full EV", "EV surface"),
    fill       = c("#b1d3d2", "#f1e2af"),   
    alpha      = c(0.6, 0.7),              
    cat.col    = c("black", "black"),  
    cat.pos    = c(-20, 20),             
    cat.dist   = c(0.05, 0.05)        
  )
  
  tiff("Figure1c_surf_full_EV_Venn.tiff", width=2, height=2, units="in", res=300)
  grid.draw(venn.plot)
  dev.off()
}
raw_data_Venn_plot()

raw_full_vs_surf_GO_analysis <- function(){
  full_accession <- full.data.list[["raw"]]$Accession
  surf_accession <- surf.data.list[["raw"]]$Accession
  
  full_set <- unique(full_accession)
  surf_set <- unique(surf_accession)
  
  full_uniq <- setdiff(full_set, surf_set)
  surf_uniq <- setdiff(surf_set, full_set)
  full_surf_over <- intersect(full_set, surf_set)
  
  full_uniq_GO <- GO_KEGG_analysis(full_uniq)
  surf_uniq_GO <- GO_KEGG_analysis(surf_uniq)
  full_surf_over_GO <- GO_KEGG_analysis(full_surf_over)
  
  return(list(full_uniq_GO = full_uniq_GO , surf_uniq_GO =surf_uniq_GO , full_surf_over_GO=full_surf_over_GO))
}
full.surf.comp.list <- raw_full_vs_surf_GO_analysis()

raw_data_Venn_GO_plot <- function(){
  temp_surf <- full.surf.comp.list[["surf_uniq_GO"]][["GO"]]@result
  temp_surf <- temp_surf[order(temp_surf$p.adjust), ][c("GO:0098742","GO:0007156","GO:1902495","GO:0031901","GO:0097060"), ]
  temp_surf$negLog10Padj <- -log10(temp_surf$p.adjust)
  
  p_surf <- ggplot(temp_surf, aes(x = reorder(Description, negLog10Padj), y = negLog10Padj)) +
    geom_bar(stat = "identity", fill = "#f6ebc8") +
    scale_x_discrete(drop = FALSE) +
    labs(x = "EV Surface Uniqe", y = "-log10(p.adjust)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = ggplot2::margin(10, 10, 10, 120)
    )
  
  temp_full <- full.surf.comp.list[["full_uniq_GO"]][["GO"]]@result
  temp_full <- temp_full[order(temp_full$p.adjust), ][c("GO:0099111","GO:0048193","GO:0044782","GO:0006858","GO:0004674"), ]
  temp_full$negLog10Padj <- -log10(temp_full$p.adjust)
  
  p_full <- ggplot(temp_full, aes(x = reorder(Description, negLog10Padj), y = negLog10Padj)) +
    geom_bar(stat = "identity", fill = "#d3e5e4") +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Full EV Uniqe", y = "-log10(p.adjust)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = ggplot2::margin(10, 10, 10, 30)
    )
  
  temp_over <- full.surf.comp.list[["full_surf_over_GO"]][["GO"]]@result
  temp_over <- temp_over[order(temp_over$p.adjust), ][c("GO:0061564","GO:0099177","GO:0099003", "GO:0031983","GO:0034774"), ]
  temp_over$negLog10Padj <- -log10(temp_over$p.adjust)
  
  p_over <- ggplot(temp_over, aes(x = reorder(Description, negLog10Padj), y = negLog10Padj)) +
    geom_bar(stat = "identity", fill = "#cedcce") +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Overlap", y = "-log10(p.adjust)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = ggplot2::margin(10, 10, 10, 30)
    )
  
  combined_plot <- p_surf +  p_over + p_full +
    plot_layout(ncol = 3) 
  
  return(combined_plot)
}
raw_data_Venn_GO_plot()
ggsave("Figure1d_raw_data_GO_plot_3panels.tiff", width = 8, height = 6, dpi = 300)

# ===================================================================
# 11  Neuronal Specific Marker
# ===================================================================

EV_neuron_marker <- c("ATP1A3","NEFM","GAP43","NEFL","CNTNAP2","STX1B","SNAP25","NSG2","GABBR2","SYN1")

neuron_expression_boxplot <- function(mydata = surf.data.list[["norm"]], mygene = EV_neuron_marker){
  data <- mydata
  genes_of_interest <- mygene
  subset_data <- data[data$Gene %in% genes_of_interest, ]
  long_data <- melt(subset_data, id.vars = c("Gene", "Accession", "Modifications"), variable.name = "Sample", value.name = "Expression")
  
  gene_order <- long_data %>%
    group_by(Gene) %>%
    summarize(Median = median(Expression, na.rm = TRUE)) %>%
    arrange(desc(Median)) %>%
    pull(Gene)
  long_data$Gene <- factor(long_data$Gene, levels = gene_order)
  
  ggplot(long_data, aes(x = Gene, y = Expression)) +
    geom_boxplot(fill = NA, color = "black") +  # 无填充色
    theme_minimal() +
    labs(x = "Gene", y = "Log2 Expression Level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
neuron_expression_boxplot()
ggsave("Figure2a_neuron_gene_expression.tiff",width=5,height=2.5)

tissue <- c("cerebral cortex","cerebellum","hippocampus","caudate","thyroid gland",
            "parathyroid gland","adrenal gland","nasopharynx","bronchus","lung",
            "oral mucosa","salivary gland","esophagus","stomach","duodenum",
            "small intestine","colon","rectum","liver","gallbladder","pancreas",
            "kidney","urinary bladder","testis","epididymis","seminal vesicle",
            "prostate","vagina","ovary","fallopian tube","endometrium","cervix",
            "placenta","breast","heart muscle","smooth muscle","skeletal muscle",
            "soft tissue","adipose tissue","skin","appendix","spleen","lymph node",
            "tonsil","bone marrow")

make_combination_frame <- function(my.tissue = tissue, my.gene=EV_neuron_marker){
  temp <- hpaVisTissue(targetGene = my.gene)
  temp <- temp[["data"]][,c(2,3,5)]
  temp2 <- temp %>%
    mutate(tissue = gsub(" \\d+$", "", tissue))
  temp2 <- temp2 %>%
    mutate(level = factor(level, levels = c("High", "Medium", "Low", "Not detected"), ordered = TRUE)) %>%
    group_by(gene, tissue) %>%
    slice_min(order_by = level) %>%
    ungroup()
  temp3 <- temp2 %>%
    distinct(gene, tissue, .keep_all = TRUE)
  
  gene_tissue_combinations <- expand.grid(gene = my.gene, tissue = my.tissue, stringsAsFactors = FALSE)
  gene_tissue_combinations <- gene_tissue_combinations %>%
    left_join(temp3, by = c("gene", "tissue")) 
  gene_tissue_combinations <- gene_tissue_combinations %>%
    pivot_wider(names_from = gene, values_from = level)
  
  gene_tissue_combinations <- as.data.frame(gene_tissue_combinations)
  rownames(gene_tissue_combinations) <- gene_tissue_combinations[,1]
  gene_tissue_combinations <- gene_tissue_combinations[,-1]
  
  return(gene_tissue_combinations)
}

neuronal_marker_tissue <- make_combination_frame()
neuronal_marker_tissue[is.na(neuronal_marker_tissue)] <- "Not detected"
neuronal_marker_tissue <- neuronal_marker_tissue[
  rowSums(neuronal_marker_tissue == "Not detected") < ncol(neuronal_marker_tissue),
]
neuronal_marker_tissue["other tissues",] <- "Not detected"

plot_neuronal_marker_tissue_heatmap <- function(df) {
  tissue_order <- rownames(df)
  gene_order <- colnames(df)
  
  df$tissue <- rownames(df)
  df_long <- melt(df, id.vars = "tissue", 
                  variable.name = "gene", 
                  value.name = "expression")
  
  df_long$tissue <- factor(df_long$tissue, levels = tissue_order)
  df_long$gene <- factor(df_long$gene, levels = gene_order)
  
  color_map <- c("High"         = "#990000",
                 "Medium"       = "#FF6666",
                 "Low"          = "#FFCCCC",
                 "Not detected" = "#FFFFFF")
  
  ggplot(df_long, aes(x = gene, y = tissue, fill = expression)) +
    geom_tile(color = NA) +  
    scale_fill_manual(values = color_map) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(x = "Tissue", y = "Gene", fill = "Expression")
}
plot_neuronal_marker_tissue_heatmap(neuronal_marker_tissue)
ggsave("Figure2b_neuronal_marker_tissue_specific_heatmap.tiff",width=8,height=4)

plot_braincell_ntpm_heatmap_by_celltype <- function(braincell_ntpm, 
                                                    cluster_to_celltype, 
                                                    EV_neuron_marker) {

  cell_type_order <- c("Blood & immune cells",
                       "Mesenchymal cells",
                       "Muscle cells",
                       "Endothelial cells",
                       "Glial cells",
                       "Neuronal cells")
  
  gene_order <- EV_neuron_marker
  
  df <- braincell_ntpm %>%
    filter(Gene %in% gene_order) %>%
    left_join(cluster_to_celltype, by = "Cluster_type") %>%
    filter(Cell_type %in% cell_type_order)
  
  df_summarized <- df %>%
    group_by(Cell_type, Gene) %>%
    summarise(mean_nTPM = mean(nTPM, na.rm = TRUE), .groups = "drop")
  
  df_neuronal <- df_summarized %>%
    dplyr::filter(Cell_type == "Neuronal cells") %>%
    dplyr::select(Gene, mean_nTPM) %>%
    dplyr::rename(neuronal_mean_nTPM = mean_nTPM)
  
  df_summarized <- df_summarized %>%
    left_join(df_neuronal, by = "Gene") %>%
    mutate(comparative_expr = mean_nTPM / neuronal_mean_nTPM)
  
  df_summarized$Cell_type <- factor(df_summarized$Cell_type, levels = cell_type_order)
  df_summarized$Gene      <- factor(df_summarized$Gene,      levels = gene_order)
  
  p <- ggplot(df_summarized, aes(x = Gene, y = Cell_type, fill = comparative_expr)) +
    geom_tile(color = NA) +  # 不加边框
    scale_fill_gradient(low = "#FFEBCD", high = "#FF8C00") +  
    scale_y_discrete(limits = rev(cell_type_order)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    ) +
    labs(x   = "Gene",
         y   = "Cell Type",
         fill = "Comparative expression")
  
  return(p)
}

plot_braincell_ntpm_heatmap_by_celltype(braincell_ntpm, 
                                        cluster_to_celltype, 
                                        EV_neuron_marker)
ggsave("Figure2c_neuronal_marker_cell_specific_heatmap.tiff",width=8,height=3)

# ===================================================================
# 12  Differential expression (t‑test AD vs CL)
# ===================================================================

surf.deg.list <- list()
full.deg.list <- list()
perform_t_test <- function(data_matrix) {
  if (!is.data.frame(data_matrix) && !is.matrix(data_matrix)) {
    warning("Input data is not a data frame or matrix. Returning NULL.")
    return(NULL)
  }
  
  data_matrix <- as.data.frame(data_matrix)
  required_cols <- c("Gene", "Accession", "Modifications")
  if (!all(required_cols %in% colnames(data_matrix))) {
    warning("Missing required columns: Gene, Accession, or Modifications. Returning NULL.")
    return(NULL)
  }
  
  sample_names <- colnames(data_matrix)[-(1:3)] 
  group_labels <- ifelse(grepl("^AD", sample_names), "AD", "CL")
  if (sum(group_labels == "AD") == 0 || sum(group_labels == "CL") == 0) {
    warning("Data does not contain both AD and CL groups. Returning NULL.")
    return(NULL)
  }
  
  results <- data.frame(
    Gene = data_matrix$Gene,
    Accession = data_matrix$Accession,
    Modifications = data_matrix$Modifications,
    p_value = NA,
    padj = NA,
    log2FC = NA
  )
  
  group1_idx <- which(group_labels == "AD")
  group2_idx <- which(group_labels == "CL")
  
  for (i in 1:nrow(data_matrix)) {
    # 提取数值并转换为数值型（确保数据可计算）
    group1_values <- suppressWarnings(na.omit(as.numeric(data_matrix[i, group1_idx + 3])))
    group2_values <- suppressWarnings(na.omit(as.numeric(data_matrix[i, group2_idx + 3])))
    
    if (length(group1_values) > 1 && length(group2_values) > 1) {
      t_test_result <- tryCatch(
        {
          t.test(group1_values, group2_values, var.equal = FALSE)
        },
        error = function(e) {
          warning(paste("t.test failed for row", i, "- setting p_value to NA"))
          return(list(p.value = NA))
        }
      )
      
      log2FC <- mean(group1_values) - mean(group2_values)
      
      results$p_value[i] <- t_test_result$p.value
      results$log2FC[i] <- log2FC
    } else {
      results$p_value[i] <- NA
      results$log2FC[i] <- NA
    }
  }
  
  results$padj <- p.adjust(results$p_value, method = "BH")
  final_results <- cbind(results[, c("p_value", "padj", "log2FC")], data_matrix)
  
  return(final_results)
}
surf.deg.list[["result"]] <- perform_t_test(surf.data.list[["norm"]])
full.deg.list[["result"]] <- perform_t_test(full.data.list[["norm"]])

surf.deg.list[["p005"]] <- surf.deg.list[["result"]][which(surf.deg.list[["result"]]$padj<=0.05 
                                                           & abs(surf.deg.list[["result"]]$log2FC)>=1),]
surf.deg.list[["p001"]] <- surf.deg.list[["result"]][which(surf.deg.list[["result"]]$padj<=0.01 
                                                           & abs(surf.deg.list[["result"]]$log2FC)>=1),]
full.deg.list[["p005"]] <- full.deg.list[["result"]][which(full.deg.list[["result"]]$padj<=0.05 
                                                           & abs(full.deg.list[["result"]]$log2FC)>=1),]
full.deg.list[["p001"]] <- full.deg.list[["result"]][which(full.deg.list[["result"]]$padj<=0.01 
                                                           & abs(full.deg.list[["result"]]$log2FC)>=1),]

surf.deg.list[["p005_GO_KEGG"]] <- GO_KEGG_analysis(surf.deg.list[["p005"]]$Accession)
surf.deg.list[["p001_GO_KEGG"]] <- GO_KEGG_analysis(surf.deg.list[["p001"]]$Accession)
full.deg.list[["p005_GO_KEGG"]] <- GO_KEGG_analysis(full.deg.list[["p005"]]$Accession)
full.deg.list[["p001_GO_KEGG"]] <- GO_KEGG_analysis(full.deg.list[["p001"]]$Accession)

surf.deg.list[["p005_multidata_anno"]] <- merge(surf.deg.list[["p005"]][,1:6], multidata_anno, by='Gene',all=F)
surf.deg.list[["p005_multidata_anno"]] <- subset(surf.deg.list[["p005_multidata_anno"]], 
                                                 (log2FC > 0 & Increased.in.AD == TRUE) | 
                                                   (log2FC < 0 & Decreased.in.AD == TRUE))
full.deg.list[["p005_multidata_anno"]] <- merge(full.deg.list[["p005"]][,1:6], multidata_anno, by='Gene',all=F)
full.deg.list[["p005_multidata_anno"]] <- subset(full.deg.list[["p005_multidata_anno"]], 
                        (log2FC > 0 & Increased.in.AD == TRUE) | 
                          (log2FC < 0 & Decreased.in.AD == TRUE))

full.deg.list[["p005_multidata_GO"]] <- GO_KEGG_analysis(full.deg.list[["p005_multidata_anno"]]$Accession)
surf.deg.list[["p005_multidata_GO"]] <- GO_KEGG_analysis(surf.deg.list[["p005_multidata_anno"]]$Accession)

temp_neuroPro <- intersect(full.deg.list[["p005"]]$Gene, neuroPro_anno$Gene)
paste0("full-neuroPro: ", length(temp_neuroPro))

temp_neuroPro <- intersect(surf.deg.list[["p005"]]$Gene, neuroPro_anno$Gene)
paste0("surf-neuroPro: ", length(temp_neuroPro))

# ===================================================================
# 13  Volcano plot 
# ===================================================================

volcano_plot <- function(deg_data = surf.deg.list[["result"]], 
                         select_data = surf.deg.list[["p005_multidata_anno"]]){
  temp <- deg_data[,1:5]
  temp_select <- select_data$Accession
  temp$group <- ifelse(temp$log2FC > 1 & temp$padj < 0.05, "Upregulated", ifelse(temp$log2FC < -1 & temp$padj < 0.05,"Downregulated","No Difference"))
  temp$select <- ifelse(temp$Accession %in% temp_select,temp$Gene,NA)
  
  ggplot(
    temp, aes(x = log2FC, y = -log10(padj), color = factor(group))
  ) + 
    geom_point(size=1.5) +
    geom_label_repel(
      color = "white",
      arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
      point.padding = NA,
      box.padding = 0.1,
      aes(label = select, fill = group), 
      size = 2,
      fontface = "bold",
      max.overlaps = 40
    ) +
    scale_fill_manual(
      values = c("Upregulated" = "tomato", "Downregulated" = "skyblue", "No Difference" = "grey"),labels=NULL
    ) +
    scale_color_manual(
      values = c("Upregulated" = "tomato", "Downregulated" = "skyblue", "No Difference" = "grey")
    ) +
    theme_bw() +
    labs(x="log2(fold change)", y="-log10 (p-value)") +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      legend.title = element_blank()
    ) +
    geom_hline(yintercept = -log10(0.05),linetype=2,cex=0.5,color = "grey")+ 
    geom_vline(xintercept = c(-1,1),linetype=2,cex=0.5,color = "grey")
}
volcano_plot(surf.deg.list[["result"]], surf.deg.list[["p005_multidata_anno"]])
ggsave("Figure4b_surf_deg_volcano.tiff",width=10, height = 9)
volcano_plot(full.deg.list[["result"]], full.deg.list[["p005_multidata_anno"]])
ggsave("Figure3b_full_deg_volcano.tiff",width=7, height = 7)


# ===================================================================
# 14  Random‑Forest biomarker pipeline (all proteins & top‑5)
# ===================================================================

random_forest_model_auc <- function(mydata, cv_n=5, predictor_class = "AD", ntree = 200, mytitle, myfeatures) {
  mydata[["status"]] <- as.factor(mydata[["status"]])
  
  set.seed(123)
  folds <- createFolds(mydata$status, k = cv_n, list = TRUE, returnTrain = TRUE)
  auc_values <- numeric(cv_n)
  rocs <- list()
  
  for(i in seq_along(folds)) {
    train_indices <- folds[[i]]
    test_indices <- setdiff(seq_len(nrow(mydata)), train_indices)
    
    train_data <- mydata[train_indices, ]
    test_data <- mydata[test_indices, ]
    train_control <- trainControl(
      method = "none",  
      classProbs = TRUE,
      summaryFunction = twoClassSummary
    )
    
    set.seed(123)
    rfModel <- train(
      as.formula(paste("status", "~" , paste(myfeatures, collapse = " + "))),
      data = train_data,
      method = "rf",
      ntree = ntree,
      trControl = train_control,
      metric = "ROC",
      importance = TRUE
    )
    
    pred_probs <- predict(rfModel, newdata = test_data, type = "prob")
    roc_obj <- roc(response = test_data$status, predictor = pred_probs[[predictor_class]])
    rocs[[i]] <- roc_obj
    auc_values[i] <- roc_obj$auc
  }
  
  mean_auc <- mean(auc_values)
  se_auc <- sd(auc_values) / sqrt(length(auc_values))
  alpha <- 0.05
  t_value <- qt(1 - alpha/2, df = length(auc_values) - 1)
  confidence_interval <- mean_auc + c(-1, 1) * t_value * se_auc
  
  ggroc_plot <- ggroc(rocs, aes = c("color")) +
    labs(title = mytitle,
         x = "False Positive Rate",
         y = "True Positive Rate") +
    theme_minimal() +
    annotate("text", x = 0.65, y = 0.07, 
             label = sprintf("Average AUC:\n%.2f", 
                             mean_auc), 
             size = 3, color = "black", hjust = 0, vjust = 0) +
    theme(legend.position = "none")
  
  print(ggroc_plot)
  
  return(list(
    folds = folds,
    model = rfModel,
    auc = list(mean_auc = mean_auc, confidence_interval = confidence_interval),
    for_curve = list(train_data = train_data, test_data = test_data)
  ))
}

rf_predict_data <- function(mygenedata = surf.deg.list[["p005_multidata_anno"]],
                            quant_data = surf.data.list[["norm"]],
                            EV_note = "surf",
                            mytitle = "EV Surface",
                            myfillcolor){
  temp <- mygenedata
  temp <- merge(temp[,3:5], quant_data, by='Accession',all=F)
  rownames(temp) <- temp$Gene
  temp <- temp[,6:ncol(temp)]
  temp <- t(temp)
  temp <- t(apply(temp, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    return(row)
  }))
  temp <- as.data.frame(temp)
  temp$sample <- rownames(temp)
  temp$status <- ifelse(grepl("AD",temp$sample),"AD","CL")
  
  png(paste0(EV_note,"_rf_AUC_all.tiff"),width=1500, height=1500, res = 500)
  all_rfresult <- random_forest_model_auc(temp, mytitle = paste0(mytitle,"- All Proteins"), myfeatures = colnames(temp)[1:(ncol(temp)-2)])
  dev.off()
  
  set.seed(123)
  rf.biop <- randomForest(status ~ ., data = all_rfresult[["for_curve"]][["train_data"]],ntree = 500,importance = TRUE )
  png(paste0(EV_note,"_rf_select_n_all.tiff"),width=3000, height=3000, res = 500)
  plot(rf.biop)
  dev.off()

  var_importance_df <- as.data.frame(all_rfresult[["model"]][["finalModel"]][["importance"]])
  var_importance_df$Variable <- rownames(var_importance_df)
  var_importance_df <- var_importance_df[,-c(1,2)]
  colnames(var_importance_df) <- c("MeanDecreaseAccuracy", "MeanDecreaseGini", "Variable")
  
  ggplot(var_importance_df, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity", fill = myfillcolor) +
    coord_flip() +  
    labs(title = paste0(mytitle,"- All Proteins"),
         y = "Mean Decrease in Gini") +
    theme_minimal()
  ggsave(paste0(EV_note,"_rf_var_all_impor.tiff"),width=5, height=6)
  
  png(paste0(EV_note,"_rf_AUC_top5.tiff"),width=1500, height=1500, res = 500)
  top5_rfresult <- random_forest_model_auc(temp, 
                                           mytitle = paste0(mytitle," - Top5 Proteins"),
                                           myfeatures = var_importance_df[order(var_importance_df$MeanDecreaseGini,decreasing = T),3][1:5])
  dev.off()

  set.seed(123)
  rf.biop <- randomForest(status ~ ., data = top5_rfresult[["for_curve"]][["train_data"]],ntree = 500,importance = TRUE )
  png(paste0(EV_note,"_rf_select_n_top5.tiff"),width=3000, height=3000, res = 500)
  plot(rf.biop)
  dev.off()
  
  var_importance_df <- as.data.frame(top5_rfresult[["model"]][["finalModel"]][["importance"]])
  var_importance_df$Variable <- rownames(var_importance_df)
  var_importance_df <- var_importance_df[,-c(1,2)]
  colnames(var_importance_df) <- c("MeanDecreaseAccuracy", "MeanDecreaseGini", "Variable")
  ggplot(var_importance_df, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity", fill = myfillcolor) +
    coord_flip() + 
    labs(title = paste0(mytitle," - Top5 Proteins"),
         y = "Mean Decrease in Gini") +
    theme_minimal()
  ggsave(paste0(EV_note,"_rf_var_top5_impor.tiff"),width=2.5, height=3)
  
  return(list(all_rfresult = all_rfresult,
              top5_rfresult = top5_rfresult))
}

surf_rf <- rf_predict_data(mygenedata = surf.deg.list[["p005_multidata_anno"]],
                           quant_data = surf.data.list[["norm"]],
                           EV_note = "Figure4_surf",
                           mytitle = "EV Surface",
                           myfillcolor = "#ffe599")
full_rf <- rf_predict_data(mygenedata = full.deg.list[["p005_multidata_anno"]],
                           quant_data = full.data.list[["norm"]],
                           EV_note = "Figure3_full",
                           mytitle = "Full EV",
                           myfillcolor = "#9dc8d9")

rf_top5_boxplot <- function(myclean.data, my.gene, mynote){
  temp <- myclean.data[myclean.data$Gene %in% my.gene,]
  
  temp <- temp %>% 
    dplyr::select(Gene, matches("^(AD|CL)_")) %>%  
    pivot_longer(cols = matches("^(AD|CL)_"), 
                 names_to = "Sample", 
                 values_to = "Expression") %>% 
    arrange(Gene)
  temp <- na.omit(temp)
  temp$Expression <- log2(temp$Expression)
  
  temp <- temp %>%
    mutate(Group = if_else(grepl("^AD_", Sample), "AD", "CL"))
  temp <- temp %>%
    mutate(Gene = factor(Gene, levels = my.gene ))
  
  ggplot(temp, aes(x = Gene, y = log2(Expression), fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("AD" = "pink", "CL" = "skyblue")) +
    labs(
      x = "Gene",
      y = "Expression Level (Log2)",
      title = "AD vs CL Boxplot by Gene"
    ) +
    theme_classic()
  
  ggsave(paste0(mynote,"_rf_top5","_boxplot.tiff"),width=6, height=3)
}

rf_top5_boxplot(myclean.data = surf.data.list[["clean"]], 
                my.gene = c("PRDX5","TAGLN3","MACROH2A1","PGLS","HPX"), 
                mynote = "Figure4g_surf")

rf_top5_boxplot(myclean.data = full.data.list[["clean"]], 
                my.gene = c("CTNNA2","FGB","PDCD6IP","STX1B","ESD"), 
                mynote = "Figure3g_full")

# ===================================================================
# 15  EV's role in AD Pathology
# ===================================================================

GO_DEG_boxplot <- function(p005_data, GO_result, term){
  temp <- GO_result[term,"geneID"]
  temp <- unlist(strsplit(temp, "/"))
  temp_df <- p005_data[p005_data$Gene %in% temp, ]
  temp <- temp_df
  
  samples <- grep("^AD_|^CL_", colnames(temp), value = TRUE)
  temp2 <- cbind(temp[,1:6],2^temp[7:ncol(temp)])
  
  for(i in seq_len(nrow(temp2))) {
    cl_cols <- grep("^CL_", samples, value = TRUE)
    median_cl <- median(unlist(temp2[i, cl_cols]), na.rm = TRUE)
    if (!is.na(median_cl) && median_cl != 1) {
      temp2[i, samples] <- temp2[i, samples] / median_cl
    }
  }
  
  temp <- temp2 %>%
    dplyr::select(Gene, log2FC, all_of(samples)) %>%
    pivot_longer(
      cols      = all_of(samples),
      names_to  = "Sample",
      values_to = "Expression"
    )
  
  temp <- temp %>%
    mutate(Group = if_else(grepl("^AD_", Sample), "AD", "CL"))
  
  temp <- temp %>%
    mutate(Gene = reorder(Gene, log2FC))
  
  ggplot(temp, aes(x = log2(Expression), y = Gene, fill = Group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("AD" = "pink", "CL" = "skyblue")) +
    labs(
      x = "Relative Expression (Log2 Fold Change)",
      y = "Gene",
      title = "AD vs CL Boxplot by Gene"
    ) +
    theme_bw()
}

GO_DEG_boxplot(full.deg.list[["p005"]], 
               full.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
               c("GO:0050821","GO:0031647","GO:0097110","GO:0140662"))
ggsave("Figure5c_full_GO_protein.tiff",width=6,height=7)

GO_DEG_boxplot(full.deg.list[["p005"]], 
               full.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
               c("GO:0001540"))
ggsave("Figure5e_full_GO_Abeta.tiff",width=6,height=7)

GO_DEG_boxplot(surf.deg.list[["p005"]], 
               surf.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
               c("GO:0007409","GO:0061564","GO:0010975"))
ggsave("Figure5h_surf_synapse2.tiff",width=6,height=7)

GO_DEG_boxplot(surf.deg.list[["p005"]], 
               surf.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
               c("GO:0045296"))
ggsave("Figure5j_surf_cadhesion.tiff",width=6,height=7)

GO_DEG_boxplot(full.deg.list[["p005"]], 
               full.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
               c("GO:0097242"))
ggsave("SuppFigure4_full_GO_Abeta_clear.tiff",width=6,height=7)

GO_plot_BP_CC_MF <- function(GO_result = full.data.list[["ori_GO_KEGG"]][["GO"]]@result,
                             my_color_plate){
  
  temp_full <- GO_result
  
  temp_BP <- temp_full[temp_full$ONTOLOGY == "BP", ]
  temp_BP <- temp_BP[order(temp_BP$p.adjust), ][1:5, ]
  temp_BP$Group <- "BP"
  temp_BP$negLog10Padj <- -log10(temp_BP$p.adjust)
  temp_BP$PlotTerm <- temp_BP$Description
  
  temp_CC <- temp_full[temp_full$ONTOLOGY == "CC", ]
  temp_CC <- temp_CC[order(temp_CC$p.adjust), ][1:5, ]
  temp_CC$Group <- "CC"
  temp_CC$negLog10Padj <- -log10(temp_CC$p.adjust)
  temp_CC$PlotTerm <- temp_CC$Description
  
  temp_MF <- temp_full[temp_full$ONTOLOGY == "MF", ]
  temp_MF <- temp_MF[order(temp_MF$p.adjust), ][1:5, ]
  temp_MF$Group <- "MF"
  temp_MF$negLog10Padj <- -log10(temp_MF$p.adjust)
  temp_MF$PlotTerm <- temp_MF$Description
  
  BP_levels <- temp_BP$PlotTerm[order(temp_BP$p.adjust)]
  CC_levels <- temp_CC$PlotTerm[order(temp_CC$p.adjust)]
  MF_levels <- temp_MF$PlotTerm[order(temp_MF$p.adjust)]
  
  gap1 <- "__gap1__"
  gap2 <- "__gap2__"
  all_levels <- c(BP_levels, gap1, CC_levels, gap2, MF_levels)
  
  temp_BP$PlotTerm <- factor(temp_BP$PlotTerm, levels = all_levels)
  temp_CC$PlotTerm <- factor(temp_CC$PlotTerm, levels = all_levels)
  temp_MF$PlotTerm <- factor(temp_MF$PlotTerm, levels = all_levels)
  temp_combined <- rbind(temp_BP, temp_CC, temp_MF)
  
  ggplot(temp_combined, aes(x = PlotTerm, y = negLog10Padj, fill = Group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("BP" = my_color_plate[1], "CC" = my_color_plate[2], "MF"=my_color_plate[3])) +
    scale_x_discrete(drop = FALSE, labels = function(x) {
      x[x %in% c("__gap1__", "__gap2__")] <- "" 
      return(x)
    }) +
    labs(x = NULL, y = "-log10(p.adjust)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = ggplot2::margin(10, 10, 10, 90)
    )
}

GO_plot_BP_CC_MF(full.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
                 c("#9dc8d9","#6cadc7","#03729d"))
ggsave("Figure5a_EV_full_GO_plot.tiff",width=13, height = 4)

GO_plot_BP_CC_MF(surf.deg.list[["p005_GO_KEGG"]][["GO"]]@result, 
                 c("#ffe599","#ffda66","#ffc000"))
ggsave("Figure5f_EV_surf_GO_plot.tiff",width=13, height = 5)

# ===================================================================
# 16  Surface DEP vs Full DEP
# ===================================================================

data_df <- data.frame(
  Category = factor(rep(c("Significant\nDEP","NeuroPro\nReported","Consistently\nAltered"), each=4),
                    levels = c("Significant\nDEP","NeuroPro\nReported","Consistently\nAltered")),
  Group     = rep(c("Full EV","Full EV","EV Surface","EV Surface"), 3),
  Regulation= rep(c("Up","Down","Up","Down"), 3),
  Value     = c( 38, -77, 128, -154,   
                 31, -60,  99, -121,  
                 4, -7,  17, -23))

data_up   <- subset(data_df, Regulation == "Up")
data_down <- subset(data_df, Regulation == "Down")
data_down$Value <- abs(data_down$Value)

label_df <- unique(data.frame(Category = data_df$Category))

p_label <- ggplot(label_df, aes(x = Category, y = 0)) +
  geom_text(aes(label = Category), size = 5) +
  theme_void() +
  scale_x_discrete(limits = levels(data_df$Category))

p_up <- ggplot(data_up, aes(x = Category, y = Value, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = Value), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 4) + 
  scale_fill_manual(values = c("Full EV"="#0a76a0","EV Surface"="#ffc000")) +
  labs(y = "Up-regulated DEPs") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank()
  )

p_down <- ggplot(data_down, aes(x = Category, y = Value, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = Value), 
            position = position_dodge(width = 0.8), 
            vjust = 1.2, size = 4) + 
  scale_fill_manual(values = c("Full EV"="#0a76a0","EV Surface"="#ffc000")) +
  scale_y_reverse() +
  labs(y = "Down-regulated DEPs") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank()
  )

p_up / p_label / p_down + plot_layout(heights = c(5, 1, 5))
ggsave("Suppfigure3a_full_surf_DEPs_number.tiff",width=7, height = 9)


full_deg_data <- full.deg.list[["result"]] %>% mutate(Group = "Full EV")
surf_deg_data <- surf.deg.list[["result"]] %>% mutate(Group = "EV Surface")
plot_data <- bind_rows(full_deg_data, surf_deg_data)

plot_data <- plot_data %>%
  group_by(Group) %>%
  mutate(rank_desc = rank(-log2FC, ties.method = "first"),
         rank_asc  = rank(log2FC,  ties.method = "first"),
         Label = if_else(rank_desc <= 3 | rank_asc <= 3, Gene, NA_character_)) %>%
  ungroup()

ggplot(plot_data, aes(x = Group, y = log2FC)) +
  geom_hline(yintercept = c(1, -1), linetype = "dashed", color = "gray50") +
  geom_jitter(aes(color = log2FC > 0), width = 0.2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "skyblue")) +
  geom_text_repel(
    data = subset(plot_data, !is.na(Label)),
    aes(label = Label),
    size = 3,
    color = "black",
    max.overlaps = Inf
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "log2 Fold Change",
    title = "Comparison of log2FC between Full EV and EV Surface"
  )

ggsave("Suppfigure3b_full_surf_DEPs_log2FC.tiff",width=6, height = 6)

