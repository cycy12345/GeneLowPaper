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
