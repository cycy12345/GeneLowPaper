#' Single-Cell Data Visualization
#'
#' Visualizes single-cell RNA-seq data including cell type annotation UMAP/t-SNE,
#' dot plot, and feature plot for a target gene.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param reduction Character. Dimensionality reduction method: \code{"umap"} or \code{"tsne"}.
#'   Default is \code{"umap"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width for UMAP/feature plots. Default is 8.
#' @param plotHeight Numeric. Plot height for UMAP/feature plots. Default is 6.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves plots to disk.
#' @details
#' Output directory: \code{outdir/6_Fig/}.
#'
#' @examples
#' \dontrun{
#' sc_basicPlot(scPath = "NSCLC_GSE148071_Single.Rds", gene = "IL36RN")
#' }
#' @export
sc_basicPlot <- function(scPath = NULL,
                         gene = "TP53",
                         reduction = "umap",
                         outdir = "./",
                         plotWidth = 8,
                         plotHeight = 6,
                         color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  sce <- readRDS(scPath)

  p <- SCNT::scPlot(
    sce,
    group_by = "celltype",
    reduction = reduction,
    point_size = 1,
    label_size = 5,
    colors = NULL,
    legend_cols = 1,
    raster_dpi = 1000
  ) +
    labs(color = "CellType") +
    theme_test(base_size = 1.5) +
    mytheme

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/sc_Anno"),
    style = "g", width = plotWidth, height = plotHeight
  )

  p <- SCNT::scDot(
    sce,
    gene,
    celltype_var = "celltype",
    assay_layer = "data",
    standardize = TRUE
  ) +
    theme_test(base_rect_size = 1.5) +
    mytheme

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/", gene, "_Dotplot"),
    style = "g", width = 6, height = 5
  )

  p <- SeuratExtend::DimPlot2(
    sce, pt.size = 1,
    features = gene,
    theme = theme_test(base_size = 1.5) +
      mytheme +
      theme(plot.title = element_text(hjust = 0.5))
  )

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/", gene, "_Featureplot"),
    style = "g", width = plotWidth, height = plotHeight
  )
}


#' Single-Gene Cell Communication Network Plot
#'
#' Performs CellChat analysis on a Seurat object after splitting a target cell type
#' into gene-positive and gene-negative groups. Draws overall, target-centric,
#' and source-centric communication circle networks.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param Tar_celltype Character. Cell type to split by gene expression.
#'   Default is \code{"Malignant"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 8.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves network plots to disk.
#' @details
#' Output directory: \code{outdir/6_Fig/cellchat/}.
#'
#' @examples
#' \dontrun{
#' gene_Cellchat1(scPath = "NSCLC_GSE148071_Single.Rds", gene = "IL36RN",
#'                Tar_celltype = "Malignant")
#' }
#' @export
gene_Cellchat1 <- function(scPath = NULL,
                           gene = "TP53",
                           Tar_celltype = "Malignant",
                           outdir = "./",
                           plotWidth = 8,
                           plotHeight = 6,
                           color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/cellchat/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  sce <- readRDS(scPath)
  exp <- Seurat::FetchData(object = sce, vars = gene)
  exp$celltype <- sce$celltype
  Cname <- stringr::str_sub(Tar_celltype, 1, 3)
  id1 <- paste0(Cname, "_", gene, "+")
  id2 <- paste0(Cname, "_", gene, "-")
  exp$Group <- ifelse(exp$celltype == Tar_celltype,
    ifelse(exp[, 1] > 0, id1, id2),
    exp$celltype
  )

  sce$Group <- exp$Group
  data.input <- sce@assays$RNA$data
  meta.data <- sce@meta.data
  gc()

  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta.data,
    group.by = "Group"
  )
  cellchat <- CellChat::addMeta(cellchat, meta = meta.data)
  cellchatDB <- CellChat::CellChatDB.human
  cellchat@DB <- cellchatDB

  cellchat <- CellChat::subsetData(cellchat)
  gc()
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  gc()

  cellchat <- CellChat::computeCommunProb(cellchat, population.size = FALSE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)

  df.net <- CellChat::subsetCommunication(cellchat)
  df.pathway <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  groupSize <- table(cellchat@idents) %>% as.numeric()

  # Overall network
  graphics::par(mfrow = c(1, 2), xpd = TRUE)
  CellChat::netVisual_circle(cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE,
    title.name = "number of Interaction"
  )
  CellChat::netVisual_circle(cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE,
    title.name = "Interaction Weight"
  )
  plot <- grDevices::recordPlot()

  grDevices::pdf(paste0(outpath, "/通讯网络.pdf"), onefile = TRUE,
    width = plotWidth, height = plotHeight)
  grDevices::replayPlot(plot)
  grDevices::dev.off()

  grDevices::png(paste0(outpath, "/通讯网络.png"),
    width = plotWidth, height = plotHeight, res = 300, units = "in"
  )
  grDevices::replayPlot(plot)
  grDevices::dev.off()

  # Target network
  graphics::par(mfrow = c(1, 2), xpd = TRUE)
  CellChat::netVisual_circle(cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE, targets.use = c(id1, id2),
    title.name = "number of Interaction"
  )
  CellChat::netVisual_circle(cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE, targets.use = c(id1, id2),
    title.name = "Interaction Weight"
  )
  plot <- grDevices::recordPlot()

  grDevices::pdf(paste0(outpath, "/", Cname, "target_通讯网络.pdf"),
    onefile = TRUE, width = plotWidth, height = plotHeight)
  grDevices::replayPlot(plot)
  grDevices::dev.off()

  grDevices::png(paste0(outpath, "/", Cname, "target_通讯网络.png"),
    width = 8, height = 6, res = 300, units = "in"
  )
  grDevices::replayPlot(plot)
  grDevices::dev.off()

  # Source network
  graphics::par(mfrow = c(1, 2), xpd = TRUE)
  CellChat::netVisual_circle(cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE, sources.use = c(id1, id2),
    title.name = "number of Interaction"
  )
  CellChat::netVisual_circle(cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE, sources.use = c(id1, id2),
    title.name = "Interaction Weight"
  )
  plot <- grDevices::recordPlot()

  grDevices::pdf(paste0(outpath, "/", Cname, "source_通讯网络.pdf"),
    onefile = TRUE, width = plotWidth, height = plotHeight)
  grDevices::replayPlot(plot)
  grDevices::dev.off()

  grDevices::png(paste0(outpath, "/", Cname, "source_通讯网络.png"),
    width = 8, height = 6, res = 300, units = "in"
  )
  grDevices::replayPlot(plot)
  grDevices::dev.off()
}


#' Single-Gene Differential Cell Communication (Ligand-Receptor)
#'
#' Compares cell-cell communication between gene-high and gene-low groups
#' within a target cell type. Outputs comparison bar plots, scatter plots of
#' signaling roles, and ranked signaling pathway differences.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param Tar_celltype Character. Cell type to split. Default is \code{"Malignant"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 8.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves comparison plots to disk.
#' @details
#' Output directory: \code{outdir/6_Fig/cellchat/}.
#' Temporary RDS files \code{tmp_Low.Rds} and \code{tmp_High.Rds} are created
#' in the working directory.
#'
#' @examples
#' \dontrun{
#' gene_Cellchat2(scPath = "NSCLC_GSE148071_Single.Rds", gene = "IL36RN",
#'                Tar_celltype = "Malignant")
#' }
#' @export
gene_Cellchat2 <- function(scPath = NULL,
                           gene = "TP53",
                           Tar_celltype = "Malignant",
                           outdir = "./",
                           plotWidth = 8,
                           plotHeight = 6,
                           color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/cellchat/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  sce <- readRDS(scPath)
  exp <- Seurat::FetchData(object = sce, vars = gene)
  exp$celltype <- sce$celltype
  Cname <- stringr::str_sub(Tar_celltype, 1, 3)
  id1 <- paste0(Cname, "_", gene, "+")
  id2 <- paste0(Cname, "_", gene, "-")
  exp$Group <- ifelse(exp$celltype == Tar_celltype,
    ifelse(exp[, 1] > 0, id1, id2),
    exp$celltype
  )

  sce$Group <- exp$Group

  # Low group
  tmp <- subset(sce, Group != id1)
  data.input <- tmp@assays$RNA$data
  meta.data <- tmp@meta.data

  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta.data,
    group.by = "celltype"
  )
  cellchat <- CellChat::addMeta(cellchat, meta = meta.data)
  gc()

  cellchatDB <- CellChat::CellChatDB.human
  cellchat@DB <- cellchatDB

  cellchat <- CellChat::subsetData(cellchat)
  gc()
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::updateCellChat(cellchat)
  gc()
  cellchat <- CellChat::computeCommunProb(cellchat, population.size = FALSE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)

  df.net <- CellChat::subsetCommunication(cellchat)
  df.pathway <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file = "tmp_Low.Rds")

  # High group
  tmp <- subset(sce, Group != id2)
  data.input <- tmp@assays$RNA$data
  meta.data <- tmp@meta.data

  cellchat <- CellChat::createCellChat(
    object = data.input,
    meta = meta.data,
    group.by = "celltype"
  )
  cellchat <- CellChat::addMeta(cellchat, meta = meta.data)
  gc()

  cellchatDB <- CellChat::CellChatDB.human
  cellchat@DB <- cellchatDB

  cellchat <- CellChat::subsetData(cellchat)
  gc()
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::updateCellChat(cellchat)
  gc()
  cellchat <- CellChat::computeCommunProb(cellchat, population.size = FALSE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)

  df.net <- CellChat::subsetCommunication(cellchat)
  df.pathway <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  saveRDS(cellchat, file = "tmp_High.Rds")

  High <- readRDS("tmp_High.Rds")
  Low <- readRDS("tmp_Low.Rds")

  colist <- list(Low = Low, High = High)

  object.list <- lapply(colist, function(x) {
    x <- CellChat::netAnalysis_computeCentrality(x)
  })

  cellchat <- CellChat::mergeCellChat(colist,
    add.names = c(id2, id1),
    cell.prefix = TRUE
  )

  gg1 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2)) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = color_dis[2:1])
  gg2 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2),
    measure = "weight") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = color_dis[2:1])
  p <- gg1 + gg2
  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/相互作用总数与强度差异"),
    style = "g", width = plotWidth, height = plotHeight
  )

  num.link <- sapply(colist, function(x) {
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
  })
  weight.MinMax <- c(min(num.link), max(num.link))
  gg <- list()
  id <- c(id2, id1)
  for (i in seq_along(colist)) {
    gg[[i]] <- CellChat::netAnalysis_signalingRole_scatter(object.list[[i]],
      title = id[i],
      weight.MinMax = weight.MinMax,
      label.size = 6, font.size.title = 20
    ) +
      mytheme+theme(legend.position = "top")
  }

  p <- patchwork::wrap_plots(plots = gg)
  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/传入与传出细胞差异"),
    style = "g", width = 10, height = 6
  )

  gg1 <- CellChat::rankNet(cellchat,
    mode = "comparison", stacked = TRUE,
    do.stat = TRUE, do.flip = FALSE,
    comparison = c(1, 2), color.use = color_dis[2:1]
  ) +
    mytheme +
    theme(legend.position = "top")

  CyDataPro::save_plot(gg1,
    filename = paste0(outpath, "/信号通路信息流差异"),
    style = "g", height = 4, width = 16
  )
}


#' Single-Cell Gene Expression Group Enrichment Analysis
#'
#' Subsets a target cell type by gene expression (positive vs negative),
#' performs differential expression analysis, and visualizes volcano plots
#' along with GO and KEGG enrichment for Up/Down genes.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param Tar_celltype Character. Cell type to subset. Default is \code{"Malignant"}.
#' @param pvlCut Numeric. P-value threshold for DEGs. Default is 0.05.
#' @param avgFC_cut Numeric. Average log2 fold change threshold. Default is 0.3.
#' @param organism Character. Species for enrichment: \code{"hsa"} (human) or \code{"mmu"} (mouse).
#' @param kegg_data_path Character. Local KEGG GSON data file (.Rds) path.
#' @param showterms Numeric. Number of enrichment terms to show. Default is 6.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 6.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves volcano plot, GO/KEGG enrichment plots, and CSV results.
#' @details
#' Output directory: \code{outdir/6_Fig/Gene_Enrichr/}.
#'
#' @export
sc_GeneEnrichr <- function(scPath = NULL,
                           gene = "TP53",
                           Tar_celltype = "Malignant",
                           pvlCut = 0.05,
                           avgFC_cut = 0.3,
                           organism = "hsa",
                           kegg_data_path = "F:/BioMed/预分析/data/KEGG/KEGG_hsa_data.Rds",
                           showterms = 6,
                           outdir = "./",
                           plotWidth = 6,
                           plotHeight = 6,
                           color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/Gene_Enrichr/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  sce <- readRDS(scPath)
  sub <- subset(sce, celltype == Tar_celltype)
  exp <- Seurat::FetchData(object = sub, vars = gene)
  rm(sce)
  gc()
  id1 <- paste0(gene, "+")
  id2 <- paste0(gene, "-")
  exp$Group <- ifelse(exp[, 1] > 0, id1, id2)
  sub$Group <- exp$Group
  message("Differential analysis...")
  DEG <- Seurat::FindMarkers(sub, group.by = "Group", ident.1 = id1, ident.2 = id2, logfc.threshold = 0)
  DEG$Gene <- rownames(DEG)
  DEG <- DEG %>% dplyr::filter(Gene != gene)
  K1 <- (DEG$p_val < pvlCut) & (DEG$avg_log2FC < -avgFC_cut)
  K2 <- (DEG$p_val < pvlCut) & (DEG$avg_log2FC > avgFC_cut)
  DEG$change <- ifelse(K1, "Down", ifelse(K2, "Up", "not"))
  DEG$Difference <- DEG$pct.1 - DEG$pct.2

  p <- ggplot(DEG, aes(x = Difference, y = avg_log2FC, color = change)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Down" = color_dis[2], "Up" = color_dis[1], "not" = color_dis[3])) +
    theme_test(base_rect_size = 1.5) + mytheme +
    labs(color = "", title = paste0(id1, " vs ", id2)) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
  CyDataPro::save_plot(p, filename = paste0(outpath, "/", gene, "_DEG_Volcanoplot"),
    style = "g", width = plotWidth, height = plotHeight)
  data.table::fwrite(DEG, file = paste0(outpath, "/", gene, "_DEG_result.csv"))
  message("Volcano plot done")

  message("Enrichment analysis...")
  if (organism == "hsa") {
    ORDB <- org.Hs.eg.db::org.Hs.eg.db
  } else if (organism == "mmu") {
    ORDB <- org.Mm.eg.db::org.Mm.eg.db
  }
  gene_df <- clusterProfiler::bitr(DEG$Gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = ORDB
  )

  sig <- DEG %>% dplyr::filter(change == "Up")
  keygene <- gene_df %>% dplyr::filter(SYMBOL %in% sig$Gene)
  gc()
  GO_Up <- clusterProfiler::enrichGO(
    gene = keygene$ENTREZID,
    OrgDb = ORDB,
    universe = gene_df$ENTREZID,
    keyType = "ENTREZID",
    ont = "BP",
    readable = FALSE,
    minGSSize = 10
  )
  gc()
  GO_Up <- clusterProfiler::simplify(GO_Up)
  GO_Up <- clusterProfiler::setReadable(GO_Up, OrgDb = ORDB, keyType = "ENTREZID")
  GO_Up_result <- GO_Up@result
  data.table::fwrite(GO_Up_result, file = paste0(outpath, "/", gene, "_Up_GO.csv"))
  message("KEGG analysis...")
  message(paste0("Local KEGG database: ", kegg_data_path))
  kk_gson <- readRDS(kegg_data_path)
  gc()
  KEGG_Up <- clusterProfiler::enrichKEGG(
    gene = keygene$ENTREZID,
    organism = kk_gson,
    universe = gene_df$ENTREZID,
    use_internal_data = FALSE
  )
  KEGG_Up <- clusterProfiler::setReadable(KEGG_Up, OrgDb = ORDB, keyType = "ENTREZID")
  KEGG_Up_result <- KEGG_Up@result
  KEGG_Up_result$ratio <- DOSE::parse_ratio(KEGG_Up_result$GeneRatio)
  KEGG_Up_result <- KEGG_Up_result %>% dplyr::filter(category != "Human Diseases")
  data.table::fwrite(KEGG_Up_result, file = paste0(outpath, "/", gene, "_Up_KEGG.csv"))

  GO_Up_result$ratio <- DOSE::parse_ratio(GO_Up_result$GeneRatio)
  pd_Up <- GO_Up_result %>% dplyr::top_n(-showterms, pvalue)
  pd_Up$Group <- id1

  sig <- DEG %>% dplyr::filter(change == "Down")
  keygene <- gene_df %>% dplyr::filter(SYMBOL %in% sig$Gene)
  gc()
  GO_Down <- clusterProfiler::enrichGO(
    gene = keygene$ENTREZID,
    OrgDb = ORDB,
    universe = gene_df$ENTREZID,
    keyType = "ENTREZID",
    ont = "BP",
    readable = FALSE,
    minGSSize = 10
  )
  gc()
  GO_Down <- clusterProfiler::simplify(GO_Down)
  GO_Down <- clusterProfiler::setReadable(GO_Down, OrgDb = ORDB, keyType = "ENTREZID")
  GO_Down_result <- GO_Down@result
  data.table::fwrite(GO_Down_result, file = paste0(outpath, "/", gene, "_Down_GO.csv"))
  gc()
  KEGG_Down <- clusterProfiler::enrichKEGG(
    gene = keygene$ENTREZID,
    organism = kk_gson,
    universe = gene_df$ENTREZID,
    use_internal_data = FALSE
  )
  KEGG_Down <- clusterProfiler::setReadable(KEGG_Down, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  KEGG_Down_result <- KEGG_Down@result
  KEGG_Down_result$ratio <- DOSE::parse_ratio(KEGG_Down_result$GeneRatio)
  KEGG_Down_result <- KEGG_Down_result %>% dplyr::filter(category != "Human Diseases")
  data.table::fwrite(KEGG_Down_result, file = paste0(outpath, "/", gene, "_Down_KEGG.csv"))

  GO_Down_result$ratio <- DOSE::parse_ratio(GO_Down_result$GeneRatio)
  pd_Down <- GO_Down_result %>% dplyr::top_n(showterms, -pvalue)
  pd_Down$Group <- id2

  GO_pd <- rbind(pd_Up, pd_Down)
  labt <- stringr::str_wrap(GO_pd$Description, width = 40)
  p <- ggplot(GO_pd, aes(x = -log(pvalue), y = FoldEnrichment, color = Group)) +
    geom_point(aes(size = ratio)) +
    ggrepel::geom_text_repel(aes(label = labt), max.overlaps = Inf, size = 5,
      fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = color_dis[2:1]) +
    theme_test(base_rect_size = 1.5) + mytheme +
    labs(x = "-Log(pvalue)", color = "", title = "GO Enrichment(BP)") +
    theme(plot.title = element_text(hjust = 0.5))
  CyDataPro::save_plot(p, filename = paste0(outpath, "/", gene, "_DEG_GO"),
    style = "g", width = plotWidth + 4, height = plotHeight + 2)

  pk_up <- KEGG_Up_result %>% dplyr::top_n(showterms, -pvalue)
  pk_up$Group <- id1

  pk_down <- KEGG_Down_result %>% dplyr::top_n(showterms, -pvalue)
  pk_down$Group <- id2

  KEGG_pd <- rbind(pk_up, pk_down)
  labt <- stringr::str_wrap(KEGG_pd$Description, width = 40)
  p <- ggplot(KEGG_pd, aes(x = -log(pvalue), y = FoldEnrichment, color = Group)) +
    geom_point(aes(size = ratio)) +
    ggrepel::geom_text_repel(aes(label = labt), hjust = 0.5, vjust = 1, size = 5,
      fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = color_dis[2:1]) +
    theme_test(base_rect_size = 1.5) + mytheme +
    labs(x = "-Log(pvalue)", color = "", title = "KEGG Enrichment(BP)") +
    theme(plot.title = element_text(hjust = 0.5))

  CyDataPro::save_plot(p, filename = paste0(outpath, "/", gene, "_DEG_KEGG"),
    style = "g", width = plotWidth + 4, height = plotHeight + 2)
}


#' Slingshot Pseudotime Trajectory Analysis
#'
#' Performs slingshot-based pseudotime trajectory analysis on a Seurat object
#' and visualizes the lineage paths and gene dynamics over pseudotime.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param cell_start Character. Starting cell type for trajectory. Default is \code{NULL}.
#' @param cell_end Character. Ending cell type for trajectory. Default is \code{NULL}.
#' @param Reduction Character. Dimensionality reduction method: \code{"UMAP"} or \code{"TSNE"}.
#'   Default is \code{"UMAP"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 6.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves trajectory plots to disk.
#' @details
#' Output directory: \code{outdir/6_Fig/Gene_slingshot/}.
#'
#' @export
sc_slingshot <- function(scPath = NULL,
                         gene = "TP53",
                         cell_start = NULL,
                         cell_end = NULL,
                         Reduction = "UMAP",
                         outdir = "./",
                         plotWidth = 6,
                         plotHeight = 6,
                         color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/Gene_slingshot/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)
  sce <- readRDS(scPath)
  sce$celltype <- factor(sce$celltype, levels = unique(sce$celltype))
  gc()
  message("Slingshot analysis...")
  sce <- scop::RunSlingshot(
    sce,
    group.by = "celltype",
    reduction = Reduction,
    start = cell_start,
    end = cell_end
  )
  gc()
  lin <- names(sce@meta.data)[grepl("Lineage", names(sce@meta.data))]
  message("Trajectory plotting...")
  p <- scop::CellDimPlot(
    sce,
    group.by = "celltype",
    lineages = lin,
    reduction = Reduction,
    palcolor = color_dis,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ) + mytheme
  CyDataPro::save_plot(p, paste0(outpath, "/slingshot_plot"), style = "g", width = plotWidth, height = plotHeight)
  gc()
  message("Gene trajectory...")
  p1 <- scop::DynamicPlot(
    sce,
    lineages = lin,
    group.by = "celltype",
    point_palcolor = color_dis,
    features = gene, legend.position = "bottom"
  )
  CyDataPro::save_plot(p1, paste0(outpath, "/Gene_plot"), style = "g", width = plotWidth, height = plotHeight)
}


#' Dorothea Transcription Factor Activity Analysis
#'
#' Computes transcription factor activities using dorothea regulons and decoupleR
#' for a target gene-defined cell group, then visualizes top differential TFs.
#'
#' @param scPath Character. Path to the pre-processed Seurat object RDS file.
#' @param gene Character. Target gene symbol.
#' @param Tar_celltype Character. Cell type to subset and analyze. Default is \code{"Malignant"}.
#' @param topTF Numeric. Number of top TFs to display. Default is 25.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 8.
#' @param plotHeight Numeric. Plot height. Default is 4.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return None. Saves TF heatmap and result CSV to disk.
#' @details
#' Output directory: \code{outdir/6_Fig/Gene_TF/}.
#'
#' @export
sc_dorotheaTF <- function(scPath = NULL,
                          gene = "TP53",
                          Tar_celltype = "Malignant",
                          topTF = 25,
                          outdir = "./",
                          plotWidth = 8,
                          plotHeight = 4,
                          color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/6_Fig/Gene_TF/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  sce <- readRDS(scPath)
  data <- subset(sce, celltype == Tar_celltype)
  exp <- Seurat::FetchData(object = data, vars = gene)
  rm(sce)
  gc()
  id1 <- paste0(gene, "+")
  id2 <- paste0(gene, "-")
  exp$Group <- ifelse(exp[, 1] > 0, id1, id2)
  data$Group <- exp$Group
  dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A", "B", "C"))
  gc()
  mat <- Seurat::GetAssayData(data)
  gc()
  message("Running decoupleR::run_ulm...")

  acts <- decoupleR::run_ulm(mat, regulon, minsize = 5,
    .source = 'tf', .target = 'target', .mor = 'mor')
  gc()

  data[['tfsulm']] <- acts %>%
    tidyr::pivot_wider(id_cols = 'source',
      names_from = 'condition',
      values_from = 'score') %>%
    tibble::column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  Seurat::DefaultAssay(object = data) <- "tfsulm"
  data <- Seurat::ScaleData(data)
  data@assays$tfsulm@data <- data@assays$tfsulm@scale.data

  df <- t(as.matrix(data@assays$tfsulm@data)) %>%
    as.data.frame() %>%
    dplyr::mutate(cluster = data$Group) %>%
    tidyr::pivot_longer(cols = -cluster,
      names_to = "source",
      values_to = "score") %>%
    dplyr::group_by(cluster, source) %>%
    dplyr::summarise(mean = mean(score))

  tfs <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(topTF, mean) %>%
    dplyr::pull(source)

  top_acts_mat <- df %>%
    dplyr::filter(source %in% tfs) %>%
    tidyr::pivot_wider(id_cols = 'cluster',
      names_from = 'source',
      values_from = 'mean') %>%
    tibble::column_to_rownames('cluster') %>%
    as.matrix()

  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)

  p <- pheatmap::pheatmap(mat = top_acts_mat,
    color = colors.use,
    border_color = "white",
    cellwidth = 15,
    cellheight = 15,
    fontsize_row = 15,
    fontsize_col = 15,
    main = Tar_celltype,
    fontsize = 15,
    treeheight_row = 20,
    treeheight_col = 20)
  CyDataPro::save_plot(p, filename = paste0(outpath, "/Top", topTF, "_TF_heatmap"), style = "x",
    width = plotWidth, height = plotHeight)
  data.table::fwrite(df, file = paste0(outpath, "/dorotheaTF_result.csv"))
}
