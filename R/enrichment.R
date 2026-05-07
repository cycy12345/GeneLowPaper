#' GO and KEGG Enrichment Analysis
#'
#' Performs GO and KEGG enrichment analysis on genes that are differentially
#' expressed between high and low expression groups of a target gene in TCGA data.
#'
#' @param TCGA_pro Character. TCGA project abbreviation, e.g. \code{"LUAD"}.
#' @param gene Character. Target gene symbol.
#' @param gene_type Character. Gene type: \code{"mRNA"}, \code{"lncrna"}, or \code{"mirna"}.
#'   Default is \code{"mRNA"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. (Currently unused; width is auto-calculated.)
#' @param plotHeight Numeric. Plot height. Default is 8.
#' @param logFC_cut Numeric. logFC threshold for DEG filtering. Default is 0.58.
#' @param pvalue_cut Numeric. P-value threshold for DEG filtering. Default is 0.05.
#' @param data_dir Character. Root directory of TCGA data.
#'   Default is \code{"F:/BioMed/预分析/data/TCGA"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns a list containing \code{GO_result} and \code{KEGG_result}.
#' @details
#' Results and plots are saved under \code{outdir/5_Fig/}.
#'
#' @examples
#' \dontrun{
#' Fig3_Enrich(TCGA_pro = "LUAD", gene = "IL36RN", gene_type = "mRNA")
#' }
#' @export
Fig3_Enrich <- function(TCGA_pro = "LUAD", gene = "TP53",
                        gene_type = "mRNA", outdir = "./",
                        plotWidth = NULL, plotHeight = 8,
                        logFC_cut = 0.58, pvalue_cut = 0.05,
                        data_dir = "F:/BioMed/预分析/data/TCGA",
                        color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/5_Fig/")
  tpm_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_", gene_type, "_expr_tpm.csv.gz")

  exp <- data.table::fread(tpm_path, data.table = FALSE)
  exp <- exp %>% tibble::column_to_rownames("V1")
  exp <- exp[rowMeans(exp == 0) < 0.5, ] %>% as.data.frame()
  exp_t <- CyDataPro::dectect_log(exp)

  Group <- ifelse(as.numeric(stringr::str_sub(colnames(exp_t), 14, 15)) < 10, "Tumor", "Normal")
  Tumor <- exp_t[, Group == "Tumor"]
  group <- ifelse(Tumor[gene, ] > stats::median(as.numeric(Tumor[gene, ])), "high", "low")
  group <- factor(group, levels = c("low", "high"))

  deg <- BM_DEG_analysis(
    exprset = Tumor, is_count = FALSE,
    logFC_cut = logFC_cut,
    pvalue_cut = pvalue_cut,
    deg_filter = TRUE, group = group
  )
  deg <- deg$deg_limma

  gene_df <- clusterProfiler::bitr(rownames(Tumor),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  keygene <- gene_df %>% dplyr::filter(SYMBOL %in% deg$genesymbol)

  GO <- clusterProfiler::enrichGO(
    gene = keygene$ENTREZID,
    OrgDb = org.Hs.eg.db,
    universe = gene_df$ENTREZID,
    keyType = "ENTREZID",
    ont = "ALL",
    readable = FALSE,
    minGSSize = 10
  )
  GO_sim <- clusterProfiler::simplify(GO)
  GO_sim <- clusterProfiler::setReadable(GO_sim, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  GO_result <- GO_sim@result

  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  data.table::fwrite(GO_result, file = paste0(outpath, "/GO_result.csv"))

  pd <- GO_result %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::top_n(6, -pvalue)
  wid <- round(max(stringr::str_length(pd$Description)) / 7, 0)

  p <- ggplot(pd) +
    ggforce::geom_link(aes(
      x = 0, y = Description,
      xend = -log10(pvalue), yend = Description,
      alpha = after_stat(index),
      color = ONTOLOGY,
      size = after_stat(index)
    ),
    n = 500, show.legend = FALSE
    ) +
    geom_point(aes(x = -log10(pvalue), y = Description),
      color = "black", fill = "white", size = 8, shape = 21
    ) +
    geom_text(aes(x = -log10(pvalue), y = Description), label = pd$Count, size = 5) +
    scale_color_manual(values = color_dis) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold.italic")
    ) +
    mytheme +
    facet_wrap(~ONTOLOGY, scales = "free", ncol = 1) +
    labs(x = "-log10(Pvalue)", y = "")

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/GO_plot"),
    style = "g", width = wid, height = plotHeight
  )

  options(timeout = 100000)
  KEGG <- clusterProfiler::enrichKEGG(
    gene = keygene$ENTREZID,
    organism = "hsa",
    universe = gene_df$ENTREZID,
    use_internal_data = FALSE
  )
  KEGG <- clusterProfiler::setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  KEGG_result <- KEGG@result
  KEGG_result <- KEGG_result %>% dplyr::filter(category != "Human Diseases")

  data.table::fwrite(KEGG_result, file = paste0(outpath, "/KEGG_result.csv"))
  pd <- KEGG_result %>% dplyr::top_n(20, -pvalue)
  wd <- round(max(stringr::str_length(pd$Description)) / 4.5, 0)

  p <- ggplot(pd, aes(zScore, reorder(Description, zScore))) +
    geom_point(aes(size = zScore, color = -log10(pvalue)), shape = 16) +
    scale_color_gradientn(colours = color_dis[c(2, 3, 1)]) +
    labs(x = "Zscore", y = "", title = "KEGG enrichment") +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold.italic")
    ) +
    mytheme

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/KEGG_plot"),
    style = "g", width = wd, height = plotHeight
  )

  invisible(list(GO_result = GO_result, KEGG_result = KEGG_result))
}
