#' Immune Infiltration Correlation Analysis
#'
#' Calculates the correlation between target gene expression and immune cell
#' infiltration scores using multiple algorithms via \pkg{IOBR}.
#'
#' @param gene Character. Target gene symbol.
#' @param TCGA_pro Character. TCGA project abbreviation, e.g. \code{"LUAD"}.
#' @param gene_type Character. Gene type: \code{"mrna"}, \code{"lncrna"}, or \code{"mirna"}. Default is \code{"mrna"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 8.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param corr_cut Numeric. Absolute correlation coefficient threshold for filtering. Default is 0.3.
#' @param Imm_method Character. Immune deconvolution method. One of
#'   \code{"mcpcounter"}, \code{"epic"}, \code{"xcell"}, \code{"cibersort"},
#'   \code{"cibersort_abs"}, \code{"quantiseq"}, or \code{"timer"}. Default is \code{"mcpcounter"}.
#' @param data_dir Character. Root directory of TCGA data. Default is \code{"F:/BioMed/预分析/data/TCGA"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the filtered correlation data frame.
#' @details
#' Results are saved under \code{outdir/5_Fig/}.
#'
#' @examples
#' \dontrun{
#' Fig5_ImmFil(gene = "IL36RN", TCGA_pro = "LUAD", Imm_method = "mcpcounter")
#' }
#' @export
Fig5_ImmFil <- function(gene,
                        TCGA_pro = "LUAD",
                        gene_type = "mrna",
                        outdir = "./",
                        plotWidth = 6,
                        plotHeight = 6,
                        corr_cut = 0.3,
                        Imm_method = c("mcpcounter", "epic", "xcell", "cibersort",
                                       "cibersort_abs", "quantiseq", "timer"),
                        data_dir = "F:/BioMed/预分析/data/TCGA",
                        color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/5_Fig/")
  Imm_method <- match.arg(Imm_method)

  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("IOBR package is not installed.")
  }

  tpm_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_", gene_type, "_expr_tpm.csv.gz")
  message("Reading TCGA expression data: ", tpm_path)
  exp <- data.table::fread(tpm_path, data.table = FALSE)
  exp <- exp %>% tibble::column_to_rownames("V1")

  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  if (any(duplicated(colnames(exp)))) {
    warning("Duplicate sample names found, adding suffixes via make.unique")
    colnames(exp) <- make.unique(colnames(exp))
  }

  exp <- exp[rowMeans(exp == 0) < 0.5, ] %>% as.data.frame()
  if (!gene %in% rownames(exp)) {
    stop("Target gene ", gene, " not found in expression matrix!")
  }
  gene_exp <- exp[gene, , drop = FALSE] %>% t() %>% as.data.frame()
  colnames(gene_exp) <- "gene_expression"

  exp_for_io <- exp
  message("Running immune infiltration analysis, method: ", Imm_method)

  tme_result <- switch(Imm_method,
    mcpcounter = IOBR::deconvo_tme(eset = exp_for_io, method = "mcpcounter"),
    epic = IOBR::deconvo_tme(eset = exp_for_io, method = "epic", arrays = FALSE),
    xcell = IOBR::deconvo_tme(eset = exp_for_io, method = "xcell", arrays = FALSE),
    cibersort = IOBR::deconvo_tme(eset = exp_for_io, method = "cibersort",
                                  arrays = FALSE, perm = 50),
    cibersort_abs = IOBR::deconvo_tme(eset = exp_for_io, method = "cibersort_abs",
                                      arrays = FALSE, perm = 20),
    quantiseq = IOBR::deconvo_tme(eset = exp_for_io, method = "quantiseq",
                                  tumor = TRUE, arrays = FALSE, scale_mrna = TRUE),
    timer = IOBR::deconvo_tme(eset = exp_for_io, method = "timer",
                              group_list = rep(tolower(TCGA_pro), ncol(exp_for_io)))
  )

  if (is.data.frame(tme_result) && ncol(tme_result) >= 2) {
    immune_mat <- tme_result
    immune_mat <- immune_mat %>% tibble::column_to_rownames("ID")
  } else {
    stop("Immune infiltration result format is abnormal.")
  }
  immune_mat <- as.data.frame(lapply(immune_mat, as.numeric))
  rownames(immune_mat) <- tme_result$ID

  remove_imm_suffix <- function(df, suffixes = c(
    "__CIBERSORT", "_MCPcounter", "_EPIC",
    "_xCell", "_quantiseq", "_CIBERSORT", "_TIMER"
  )) {
    pattern <- paste0("(", paste(suffixes, collapse = "|"), ")$")
    colnames(df) <- gsub(pattern, "", colnames(df))
    return(df)
  }
  immune_mat <- remove_imm_suffix(immune_mat)

  common_samples <- intersect(rownames(gene_exp), rownames(immune_mat))
  if (length(common_samples) == 0) stop("No common samples between gene expression and immune infiltration!")
  gene_exp_common <- gene_exp[common_samples, , drop = FALSE]
  immune_mat_common <- immune_mat[common_samples, , drop = FALSE]

  message("Calculating Spearman correlations...")
  cor_res <- data.frame(
    Cell = colnames(immune_mat_common),
    Correlation = NA,
    P_value = NA,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(ncol(immune_mat_common))) {
    test <- stats::cor.test(
      gene_exp_common$gene_expression,
      immune_mat_common[, i],
      method = "spearman",
      exact = FALSE
    )
    cor_res[i, "Correlation"] <- test$estimate
    cor_res[i, "P_value"] <- test$p.value
  }
  cor_res$Significant <- ifelse(cor_res$P_value < 0.05, "P < 0.05", "P >= 0.05")

  cor_res_filtered <- cor_res[cor_res$P_value < 0.05, ]
  if (nrow(cor_res_filtered) == 0) {
    message("No significant immune correlations found.")
    data.table::fwrite(cor_res,
      file = paste0(outdir, "/Fig5_", Imm_method, "_correlation_all.csv")
    )
    return(cor_res)
  }

  cor_res_filtered <- cor_res_filtered[order(-cor_res_filtered$Correlation), ]
  cor_res_filtered$Cell <- factor(cor_res_filtered$Cell,
    levels = rev(cor_res_filtered$Cell)
  )

  cor_res_filtered$Color <- ifelse(cor_res_filtered$Correlation > 0,
    color_dis[1], color_dis[2]
  )
  cor_res_filtered$PointSize <- abs(cor_res_filtered$Correlation) * 5

  p <- ggplot(cor_res_filtered, aes(x = Correlation, y = Cell)) +
    geom_segment(aes(
      x = 0, xend = Correlation, y = Cell, yend = Cell,
      color = Color
    ), linewidth = 1.2) +
    geom_point(aes(size = PointSize, fill = Color),
      shape = 21, stroke = 2, color = "gray"
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_size_continuous(range = c(3, 10), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
    labs(
      x = "Spearman Correlation Coefficient", y = NULL,
      title = paste0(gene, " vs Immune Infiltration (", Imm_method, ")")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15)
    ) +
    geom_text(aes(label = ifelse(P_value < 0.001, "***",
      ifelse(P_value < 0.01, "**",
        ifelse(P_value < 0.05, "*", "")
      )
    )),
    x = Inf, hjust = 1.2, size = 6, color = "black"
    ) +
    mytheme

  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/", Imm_method, "_cor_plot"),
    style = "g", width = plotWidth, height = plotHeight
  )
  data.table::fwrite(cor_res_filtered,
    file = paste0(outpath, "/", Imm_method, "_correlation_significant.csv")
  )
  message("Analysis complete! Results saved to: ", outpath)
  invisible(cor_res_filtered)
}


#' ESTIMATE and IPS Immune Score Analysis
#'
#' Computes ESTIMATE and IPS (Immunophenoscore) scores for tumor samples
#' stratified by high/low expression of the target gene.
#'
#' @param gene Character. Target gene symbol.
#' @param TCGA_pro Character. TCGA project abbreviation. Default is \code{"LUAD"}.
#' @param gene_type Character. Gene type: \code{"mrna"}, \code{"lncrna"}, or \code{"mirna"}.
#'   Default is \code{"mrna"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 6.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param data_dir Character. Root directory of TCGA data. Default is \code{"F:/BioMed/预分析/data/TCGA"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns a list with ESTIMATE and IPS results.
#' @examples
#' \dontrun{
#' Fig5_Est_IPS(gene = "IL36RN", TCGA_pro = "LUAD")
#' }
#' @export
Fig5_Est_IPS <- function(gene,
                         TCGA_pro = "LUAD",
                         gene_type = "mrna",
                         outdir = "./",
                         plotWidth = 6,
                         plotHeight = 6,
                         data_dir = "F:/BioMed/预分析/data/TCGA",
                         color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/5_Fig/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  tpm_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_", gene_type, "_expr_tpm.csv.gz")
  message("Reading TCGA expression data: ", tpm_path)
  exp <- data.table::fread(tpm_path, data.table = FALSE)
  exp <- exp %>% tibble::column_to_rownames("V1")

  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  if (any(duplicated(colnames(exp)))) {
    colnames(exp) <- make.unique(colnames(exp))
  }

  exp <- exp[rowMeans(exp == 0) < 0.5, ] %>% as.data.frame()
  if (!gene %in% rownames(exp)) {
    stop("Target gene ", gene, " not found in expression matrix!")
  }

  Group <- ifelse(as.numeric(stringr::str_sub(colnames(exp), 14, 15)) < 10, "Tumor", "Normal")
  Tumor <- exp[, Group == "Tumor", drop = FALSE]

  gene_exp_tumor <- as.numeric(Tumor[gene, ])
  group <- ifelse(gene_exp_tumor > stats::median(gene_exp_tumor, na.rm = TRUE), "high", "low")
  group <- factor(group, levels = c("low", "high"))
  message("High group: ", sum(group == "high"), "; Low group: ", sum(group == "low"))

  message("Running ESTIMATE analysis...")
  estimate_res <- IOBR::deconvo_tme(eset = Tumor, method = "estimate")
  estimate_mat <- estimate_res %>% tibble::column_to_rownames("ID")
  colnames(estimate_mat) <- gsub("_estimate$", "", colnames(estimate_mat))
  estimate_mat$Group <- group[match(rownames(estimate_mat), colnames(Tumor))]

  message("Running IPS analysis...")
  ips_res <- IOBR::deconvo_tme(eset = Tumor, method = "ips")
  ips_mat <- ips_res %>% tibble::column_to_rownames("ID")
  colnames(ips_mat) <- gsub("_IPS$", "", colnames(ips_mat))
  ips_mat$Group <- group[match(rownames(ips_mat), colnames(Tumor))]

  plot_facet_violin <- function(data, title_prefix, filename,
                                outpath, plotWidth, plotHeight, gene, color_dis) {
    pd <- data %>%
      tibble::rownames_to_column("Sample") %>%
      tidyr::pivot_longer(
        cols = -c(Sample, Group),
        names_to = "Score_Type",
        values_to = "Score"
      ) %>%
      dplyr::filter(!is.na(Group))

    stat_res <- pd %>%
      dplyr::group_by(Score_Type) %>%
      dplyr::summarise(
        p_value = stats::wilcox.test(Score ~ Group, exact = FALSE)$p.value,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        p_label = dplyr::case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )

    y_max <- pd %>%
      dplyr::group_by(Score_Type) %>%
      dplyr::summarise(y_pos = max(Score, na.rm = TRUE) * 1, .groups = "drop")
    stat_res <- dplyr::left_join(stat_res, y_max, by = "Score_Type")

    score_levels <- setdiff(colnames(data), "Group")
    pd$Score_Type <- factor(pd$Score_Type, levels = score_levels)
    stat_res$Score_Type <- factor(stat_res$Score_Type, levels = score_levels)

    p <- ggplot(pd, aes(x = Group, y = Score, fill = Group)) +
      geom_violin(trim = FALSE, alpha = 0.4, linewidth = 0.8) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, linewidth = 0.8) +
      ggbeeswarm::geom_quasirandom(aes(color = Group), size = 1.2, alpha = 0.5, width = 0.2) +
      geom_text(
        data = stat_res,
        aes(x = 1.5, y = y_pos, label = p_label),
        inherit.aes = FALSE,
        size = 5, fontface = "bold", vjust = 0
      ) +
      facet_wrap(~Score_Type, scales = "free_y", nrow = 2) +
      scale_fill_manual(values = c("low" = color_dis[2], "high" = color_dis[1])) +
      scale_color_manual(values = c("low" = color_dis[2], "high" = color_dis[1])) +
      labs(
        title = paste0(gene, " - ", title_prefix),
        x = paste0(gene, " Expression"),
        y = "Score", fill = "", color = ""
      ) +
      theme_test(base_rect_size = 1.5) +
      mytheme +
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
        strip.text = element_text(face = "bold", size = 13),
        panel.spacing = grid::unit(0.5, "cm")
      )

    CyDataPro::save_plot(p,
      filename = paste0(outpath, "/", filename),
      style = "g", width = plotWidth, height = plotHeight
    )
    list(plot = p, stats = stat_res, data = pd)
  }

  message("Plotting ESTIMATE facet...")
  estimate_out <- plot_facet_violin(
    data = estimate_mat, title_prefix = "ESTIMATE",
    filename = "ESTIMATE_facet", outpath = outpath,
    plotWidth = plotWidth, plotHeight = plotHeight,
    gene = gene, color_dis = color_dis
  )

  message("Plotting IPS facet...")
  ips_out <- plot_facet_violin(
    data = ips_mat, title_prefix = "IPS",
    filename = "IPS_facet", outpath = outpath,
    plotWidth = plotWidth, plotHeight = plotHeight,
    gene = gene, color_dis = color_dis
  )

  data.table::fwrite(estimate_res, file = paste0(outpath, "/ESTIMATE_result.csv"))
  data.table::fwrite(ips_res, file = paste0(outpath, "/IPS_result.csv"))
  data.table::fwrite(estimate_out$stats, file = paste0(outpath, "/ESTIMATE_stats.csv"))
  data.table::fwrite(ips_out$stats, file = paste0(outpath, "/IPS_stats.csv"))

  message("ESTIMATE & IPS analysis complete! Results saved to: ", outpath)
  invisible(list(
    estimate = list(results = estimate_res, stats = estimate_out$stats, plot = estimate_out$plot),
    ips = list(results = ips_res, stats = ips_out$stats, plot = ips_out$plot)
  ))
}


#' Immune Checkpoint Gene Correlation Analysis
#'
#' Computes the correlation between target gene expression and 79 immune
#' checkpoint genes (PMID: 32814346) in TCGA tumor samples.
#'
#' @param gene Character. Target gene symbol.
#' @param TCGA_pro Character. TCGA project abbreviation. Default is \code{"LUAD"}.
#' @param gene_type Character. Gene type: \code{"mrna"}, \code{"lncrna"}, or \code{"mirna"}. Default is \code{"mrna"}.
#' @param cor_method Character. Correlation method: \code{"pearson"} or \code{"spearman"}.
#'   Default is \code{"pearson"}.
#' @param p_threshold Numeric. P-value threshold. Default is 0.05.
#' @param r_threshold Numeric. Correlation coefficient threshold. Default is 0.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 6.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param data_dir Character. Root directory of TCGA data. Default is \code{"F:/BioMed/预分析/data/TCGA"}.
#' @param icc_path Character. Path to the immune checkpoint gene list CSV. Default is \code{"F:/BioMed/预分析/data/IGG79Gene.csv"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns a list with full and significant correlation tables,
#'   plus the bar plot.
#' @examples
#' \dontrun{
#' Fig5_ICC_cor(gene = "IL36RN", TCGA_pro = "LUAD")
#' }
#' @export
Fig5_ICC_cor <- function(gene,
                         TCGA_pro = "LUAD",
                         gene_type = "mrna",
                         cor_method = "pearson",
                         p_threshold = 0.05,
                         r_threshold = 0,
                         outdir = "./",
                         plotWidth = 6,
                         plotHeight = 6,
                         data_dir = "F:/BioMed/预分析/data/TCGA",
                         icc_path = "F:/BioMed/预分析/data/IGG79Gene.csv",
                         color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/5_Fig/")
  if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

  tpm_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_", gene_type, "_expr_tpm.csv.gz")
  message("Reading TCGA expression data: ", tpm_path)
  exp <- data.table::fread(tpm_path, data.table = FALSE)
  exp <- exp %>% tibble::column_to_rownames("V1")

  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  if (any(duplicated(colnames(exp)))) {
    colnames(exp) <- make.unique(colnames(exp))
  }

  exp <- exp[rowMeans(exp == 0) < 0.5, ] %>% as.data.frame()
  if (!gene %in% rownames(exp)) {
    stop("Target gene ", gene, " not found in expression matrix!")
  }
  exp <- CyDataPro::dectect_log(exp)

  Group <- ifelse(as.numeric(stringr::str_sub(colnames(exp), 14, 15)) < 10, "Tumor", "Normal")
  Tumor <- exp[, Group == "Tumor", drop = FALSE]

  icc_genes <- data.table::fread(icc_path, data.table = FALSE)
  gene_col <- grep("gene|symbol|name", colnames(icc_genes), ignore.case = TRUE, value = TRUE)
  if (length(gene_col) == 0) {
    gene_col <- colnames(icc_genes)[1]
    message("Gene column not found, using first column: ", gene_col)
  } else {
    gene_col <- gene_col[1]
    message("Using gene column: ", gene_col)
  }

  icc_list <- icc_genes[[gene_col]] %>% unique() %>% as.character()
  icc_list <- icc_list[icc_list %in% rownames(Tumor)]
  icc_list <- setdiff(icc_list, gene)

  message("Total immune checkpoint genes: ", nrow(icc_genes))
  message("Found in expression matrix: ", length(icc_list))

  if (length(icc_list) == 0) {
    stop("No available immune checkpoint genes found!")
  }

  message("Calculating correlations (method: ", cor_method, ")...")
  cor_res <- purrr::map_dfr(icc_list, function(icc_gene) {
    test <- stats::cor.test(
      as.numeric(Tumor[gene, ]),
      as.numeric(Tumor[icc_gene, ]),
      method = cor_method
    )
    tibble::tibble(
      Gene_A = gene,
      Gene_B = icc_gene,
      r = test$estimate,
      p_value = test$p.value
    )
  })

  sig_cor <- cor_res %>%
    dplyr::filter(p_value < p_threshold, abs(r) > r_threshold) %>%
    dplyr::arrange(dplyr::desc(abs(r)))

  message("Significant immune checkpoint genes: ", nrow(sig_cor))
  message("Positive: ", sum(sig_cor$r > 0), "; Negative: ", sum(sig_cor$r < 0))

  data.table::fwrite(cor_res, file = paste0(outpath, "/ICC_correlation_all.csv"))
  data.table::fwrite(sig_cor, file = paste0(outpath, "/ICC_correlation_significant.csv"))

  if (nrow(sig_cor) == 0) {
    warning("No significant immune checkpoint correlations found, skipping bar plot")
    return(invisible(list(correlation = cor_res, significant = sig_cor)))
  }

  p_bar <- sig_cor %>%
    dplyr::slice(seq_len(dplyr::n())) %>%
    dplyr::mutate(
      Gene_B = forcats::fct_reorder(Gene_B, r),
      direction = ifelse(r > 0, "positive", "negative")
    ) %>%
    ggplot2::ggplot(aes(x = Gene_B, y = r, fill = direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("positive" = color_dis[1], "negative" = color_dis[2]),
      labels = c("positive" = "Positive", "negative" = "Negative")
    ) +
    ggplot2::labs(x = "", y = "Correlation Coefficient (r)", fill = "") +
    ggplot2::coord_flip() +
    theme_test(base_rect_size = 1.5) +
    mytheme +
    ggplot2::theme(legend.position = "top")

  CyDataPro::save_plot(p_bar,
    filename = paste0(outpath, "/ICC_correlation_barplot"),
    style = "g", width = plotWidth, height = plotHeight
  )

  message("Immune checkpoint correlation analysis complete! Results saved to: ", outpath)
  invisible(list(
    correlation = cor_res,
    significant = sig_cor,
    barplot = p_bar
  ))
}
