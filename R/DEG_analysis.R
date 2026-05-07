#' Integrated Differential Expression Analysis
#'
#' Performs differential expression analysis using multiple methods.
#' When \code{is_count = TRUE}, runs DESeq2, limma-voom, and edgeR on raw counts.
#' When \code{is_count = FALSE}, runs limma and Wilcoxon rank-sum test on normalized data.
#'
#' @param exprset A gene expression matrix with genes as rows and samples as columns.
#' @param project Character string used as the file name prefix when \code{save = TRUE}.
#' @param is_count Logical. If \code{TRUE}, the input is treated as raw count matrix.
#'   If \code{FALSE}, the input is treated as normalized expression matrix.
#' @param logFC_cut Numeric. Log2 fold change threshold for filtering DEGs. Default is 0.
#' @param pvalue_cut Numeric. P-value threshold for filtering DEGs. Default is 1.
#' @param adjpvalue_cut Numeric. Adjusted p-value threshold for filtering DEGs. Default is 1.
#' @param deg_filter Logical. Whether to filter results by the thresholds above. Default is \code{FALSE}.
#' @param group A factor vector of group labels for each sample. Control should be the first level.
#'   If \code{NULL}, groups are inferred from TCGA barcodes (sample type 01-09 = tumor, 10-19 = normal).
#' @param save Logical. Whether to save results to disk. Default is \code{FALSE}.
#'
#' @return A named list of DEG results. For count data: \code{deg_deseq2}, \code{deg_limma},
#'   \code{deg_edger}. For normalized data: \code{deg_limma}, \code{deg_wilcoxon}.
#'
#' @details
#' The function automatically detects whether log2 transformation is needed for
#' non-count data by inspecting quantiles.
#'
#' @examples
#' \dontrun{
#' # Example with raw counts
#' res <- BM_DEG_analysis(exprset = count_matrix, project = "LUAD",
#'                        is_count = TRUE, group = group_factor, save = TRUE)
#'
#' # Example with normalized data
#' res <- BM_DEG_analysis(exprset = norm_matrix, is_count = FALSE,
#'                        logFC_cut = 1, pvalue_cut = 0.05, deg_filter = TRUE)
#' }
#'
#' @import DESeq2 limma edgeR
#' @export
BM_DEG_analysis <- function(exprset, project = NULL, is_count = TRUE,
                            logFC_cut = 0, deg_filter = FALSE,
                            pvalue_cut = 1, adjpvalue_cut = 1,
                            group = NULL, save = FALSE)
{
  if (!is.null(group)) {
    if (!is.factor(group))
      stop("The type of group should be factor.")
    metadata <- data.frame(
      sample_id = colnames(exprset),
      group = group
    )
  } else {
    group <- ifelse(as.numeric(substr(colnames(exprset), 14, 15)) < 10,
                    "tumor", "normal")
    group <- factor(group, levels = c("normal", "tumor"))
    metadata <- data.frame(
      sample_id = colnames(exprset),
      group = group
    )
  }

  res_diff <- list()

  if (is_count) {
    message("=> Running DESeq2")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = exprset,
      colData = metadata,
      design = ~group
    )
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, tidy = TRUE)
    names(res)[1] <- "genesymbol"
    deg_deseq2 <- stats::na.omit(res)
    if (deg_filter) {
      deg_deseq2 <- subset(deg_deseq2,
        abs(log2FoldChange) > logFC_cut &
          padj < adjpvalue_cut & pvalue < pvalue_cut
      )
    }
    res_diff[[1]] <- deg_deseq2
    names(res_diff)[[1]] <- "deg_deseq2"

    message("=> Running limma voom")
    y <- edgeR::DGEList(counts = exprset, group = group)
    keep <- edgeR::filterByExpr(y, group = group)
    y <- y[keep, keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    design <- stats::model.matrix(~group)
    rownames(design) <- colnames(y)
    colnames(design) <- levels(y)
    v <- limma::voom(y, design, normalize = "quantile")
    fit <- limma::lmFit(v, design)
    fit2 <- limma::eBayes(fit)
    DEG2 <- limma::topTable(fit2, coef = 2, n = Inf)
    deg_limma <- stats::na.omit(DEG2)
    deg_limma$genesymbol <- rownames(deg_limma)
    if (deg_filter) {
      deg_limma <- subset(deg_limma,
        abs(logFC) > logFC_cut &
          adj.P.Val < adjpvalue_cut & P.Value < pvalue_cut
      )
    }
    res_diff[[2]] <- deg_limma
    names(res_diff)[[2]] <- "deg_limma"

    message("=> Running edgeR")
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    DEG <- edgeR::topTags(qlf, n = Inf)
    deg_edger <- as.data.frame(DEG)
    deg_edger$genesymbol <- rownames(deg_edger)
    if (deg_filter) {
      deg_edger <- subset(deg_edger,
        abs(logFC) > logFC_cut &
          FDR < adjpvalue_cut & PValue < pvalue_cut
      )
    }
    res_diff[[3]] <- deg_edger
    names(res_diff)[[3]] <- "deg_edger"
  } else {
    exprset <- limma::normalizeBetweenArrays(exprset)
    ex <- exprset
    qx <- as.numeric(stats::quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
    LogC <- (qx[5] > 100) ||
      (qx[6] - qx[1] > 50 && qx[2] > 0) ||
      (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) {
      exprset <- log2(ex + 0.1)
      message("=> log2(x+0.1) transform finished")
    } else {
      message("=> log2 transform not needed")
    }

    message("=> Running limma")
    design <- stats::model.matrix(~group)
    fit <- limma::lmFit(exprset, design)
    fit <- limma::eBayes(fit)
    deg_limma <- limma::topTable(fit, coef = 2, number = Inf)
    deg_limma <- stats::na.omit(deg_limma)
    deg_limma$genesymbol <- rownames(deg_limma)
    if (deg_filter) {
      deg_limma <- subset(deg_limma,
        abs(logFC) > logFC_cut &
          adj.P.Val < adjpvalue_cut & P.Value < pvalue_cut
      )
    }
    res_diff[[1]] <- deg_limma
    names(res_diff)[[1]] <- "deg_limma"

    message("=> Running wilcoxon test")
    exprset <- t(exprset)
    wilcox_res <- apply(exprset, 2, function(x) {
      stats::wilcox.test(x ~ group, exact = FALSE, correct = FALSE)$p.value
    })
    deg_wilcoxon <- as.data.frame(wilcox_res)
    names(deg_wilcoxon) <- "pvalue"
    deg_wilcoxon$genesymbol <- rownames(deg_wilcoxon)
    res_diff[[2]] <- deg_wilcoxon
    names(res_diff)[[2]] <- "deg_wilcoxon"
  }

  if (save) {
    if (!dir.exists("output_diff")) {
      dir.create("output_diff")
    }
    if (is.null(project))
      stop("project should be provided!")
    save(res_diff,
      file = paste0("output_diff/", project, "_diff_results.rdata")
    )
  }
  message("=> Analysis done.")
  return(res_diff)
}
