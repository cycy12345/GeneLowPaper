#' TCGA Gene Expression and Clinical Correlation Analysis
#'
#' Analyzes the association between target gene expression and clinical features
#' (Gender, Age, T_stage, M_stage, N_stage, Stage) in TCGA data. Performs
#' univariate and multivariate Cox regression and generates forest plots.
#' If the gene is significant in multivariate Cox, a nomogram is also drawn.
#'
#' @param TCGA_pro Character. TCGA project abbreviation, e.g. \code{"LUAD"}.
#' @param gene Character. Target gene symbol.
#' @param gene_type Character. Gene type: \code{"mRNA"}, \code{"lncrna"}, or \code{"mirna"}.
#'   Default is \code{"mRNA"}.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Clinical plot width. Default is 8.
#' @param plotHeight Numeric. Clinical plot height. Default is 6.
#' @param coxWidth Numeric. Cox forest plot width. Default is 8.
#' @param coxHeight Numeric. Cox forest plot height. Default is 6.
#' @param data_dir Character. Root directory of TCGA data.
#'   Default is \code{"F:/BioMed/预分析/data/TCGA"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the multivariate Cox result data frame.
#' @details
#' Output files are saved under \code{outdir/2_Fig/Fig2D_clin/} and include:
#' alluvial plots for each clinical feature, univariate and multivariate Cox
#' CSV results and forest plots, and a nomogram (if applicable).
#'
#' @examples
#' \dontrun{
#' Fig2D_TCGA_Clin(TCGA_pro = "LUAD", gene = "IL36RN", gene_type = "mRNA")
#' }
#' @export
Fig2D_TCGA_Clin <- function(TCGA_pro = "LUAD", gene,
                            gene_type = "mRNA",
                            outdir = "./",
                            plotWidth = 8, plotHeight = 6,
                            coxWidth = 8, coxHeight = 6,
                            data_dir = "F:/BioMed/预分析/data/TCGA",
                            color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/2_Fig/Fig2D_clin/")
  tpm_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_", gene_type, "_expr_tpm.csv.gz")
  Survi_meta_path <- paste0(data_dir, "/TCGA-", TCGA_pro, "/TCGA-", TCGA_pro, "_clinical_OS.csv.gz")

  exp <- data.table::fread(tpm_path, data.table = FALSE)
  exp <- exp %>% tibble::column_to_rownames("V1")
  exp <- exp[rowMeans(exp == 0) < 0.5, ] %>% as.data.frame()

  exp_t <- CyDataPro::dectect_log(exp)
  Gene_exp <- exp_t[gene, ] %>% t() %>% as.data.frame()
  Gene_exp$sample <- stringr::str_replace_all(rownames(Gene_exp), "[.]", "-")

  clin <- data.table::fread(Survi_meta_path, data.table = FALSE)
  clin <- clin %>% dplyr::select(
    sample_submitter_id, gender, age, vital_status,
    ajcc_t, ajcc_m, ajcc_n, ajcc_stage, vitalStat, surTime
  )
  colnames(clin) <- c(
    "sample", "Gender", "Age", "Status", "T_stage",
    "M_stage", "N_stage", "Stage", "Event", "Time"
  )

  meta <- dplyr::inner_join(Gene_exp, clin)
  meta$Gender <- stringr::str_to_title(meta$Gender)
  meta$Age <- as.numeric(meta$Age)
  meta$Age <- ifelse(meta$Age >= 60, ">=60", "<60")

  meta$Stage <- stringr::str_remove_all(meta$Stage, "Stage ")
  meta <- meta %>% dplyr::filter(Stage != "" & Stage != "X")
  meta <- meta %>% dplyr::filter(!is.na(Gender))
  meta <- meta %>% dplyr::filter(!is.na(Age))
  meta <- meta %>% dplyr::filter(!is.na(Status))
  meta <- meta %>% dplyr::filter(T_stage != "" & T_stage != "TX")
  meta <- meta %>% dplyr::filter(M_stage != "" & M_stage != "MX")
  meta <- meta %>% dplyr::filter(N_stage != "" & N_stage != "NX")

  res_cat <- survminer::surv_cutpoint(meta,
    time = "Time",
    event = "Event",
    variables = gene,
    minprop = 0.3
  )
  res_cat <- survminer::surv_categorize(res_cat)
  meta$Group <- stringr::str_to_title(res_cat[, 3])

  # Helper for alluvial plot
  .plot_alluvial <- function(meta, var, gene, outpath, plotWidth, plotHeight, color_dis) {
    tmp <- meta[, c("Group", var)]
    tmp <- tmp %>% dplyr::filter(!is.na(.data[[var]]) & .data[[var]] != "" & .data[[var]] != "TX" & .data[[var]] != "MX" & .data[[var]] != "NX" & .data[[var]] != "X")
    tt <- as.data.frame(table(tmp$Group, tmp[[var]]))
    test_table <- table(tmp$Group, tmp[[var]])
    if (nrow(test_table) < 2 || ncol(test_table) < 2) {
      message(paste0(var, " table too sparse for test"))
      p <- NA
    } else {
      if (min(test_table) < 5) {
        re <- stats::fisher.test(test_table)
      } else {
        re <- stats::chisq.test(test_table)
      }
      p <- re$p.value
      tt <- tt %>%
        dplyr::group_by(Var1) %>%
        dplyr::mutate(percent = Freq / sum(Freq, na.rm = TRUE)) %>%
        dplyr::mutate(label = paste0(round(percent * 100, 2), "%"))
      tt$Var1 <- factor(tt$Var1, levels = c("Low", "High"))
      p <- ggplot(tt, aes(x = Var1, y = percent, fill = Var2, stratum = Var2, alluvium = Var2)) +
        ggalluvial::geom_stratum(width = 0.5, color = "white", size = 0.6, alpha = 0.4) +
        ggalluvial::geom_alluvium(alpha = 0.6, width = 0.5, curve_type = "linear") +
        geom_text(aes(label = label),
          position = position_fill(vjust = 0.5),
          size = 5.5
        ) +
        labs(
          x = paste0(gene, " Expression"), y = "Prob", fill = "",
          title = var, subtitle = paste0("p = ", round(p, 3))
        ) +
        ggpubr::theme_test(base_rect_size = 1.5) +
        mytheme +
        scale_fill_manual(values = color_dis) +
        theme(
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )
      if (!dir.exists(outpath)) {
        dir.create(outpath, recursive = TRUE)
      }
      CyDataPro::save_plot(p,
        filename = paste0(outpath, "/", var),
        style = "g", width = plotWidth, height = plotHeight
      )
    }
    invisible(p)
  }

  .plot_alluvial(meta, "Gender", gene, outpath, plotWidth, plotHeight, color_dis)
  .plot_alluvial(meta, "Age", gene, outpath, plotWidth, plotHeight, color_dis)
  .plot_alluvial(meta, "T_stage", gene, outpath, plotWidth, plotHeight, color_dis)
  .plot_alluvial(meta, "M_stage", gene, outpath, plotWidth, plotHeight, color_dis)
  .plot_alluvial(meta, "N_stage", gene, outpath, plotWidth, plotHeight, color_dis)
  .plot_alluvial(meta, "Stage", gene, outpath, plotWidth, plotHeight, color_dis)

  data.table::fwrite(meta, file = paste0(outpath, "Clin_meta.csv"))

  # Univariate Cox
  meta$Group <- factor(meta$Group, levels = c("Low", "High"))
  meta[, 1] <- meta$Group
  meta$Age <- factor(meta$Age, levels = c("<60", ">=60"))
  uni_cox <- ezcox::ezcox(meta,
    covariates = c(colnames(meta)[1], "Gender", "Age", "T_stage",
                   "M_stage", "N_stage", "Stage"),
    time = "Time", status = "Event"
  )

  data.table::fwrite(uni_cox, file = paste0(outdir, "/2_Fig/单变量_cox_result.csv"))

  # Forest plot helper
  .plot_forest <- function(data, title, outname, coxWidth, coxHeight) {
    data <- data %>%
      dplyr::mutate(
        HR = as.numeric(HR),
        lower_95 = as.numeric(lower_95),
        upper_95 = as.numeric(upper_95),
        p.value = as.numeric(p.value)
      ) %>%
      dplyr::filter(
        !is.na(HR), !is.na(lower_95), !is.na(upper_95),
        HR > 1e-6 & HR < 1e4,
        is.finite(lower_95), is.finite(upper_95),
        upper_95 > lower_95
      ) %>%
      dplyr::arrange(Variable, contrast_level) %>%
      dplyr::mutate(
        label = paste0(Variable, ": ", contrast_level, " vs ", ref_level),
        label = factor(label, levels = unique(label)),
        color_group = dplyr::case_when(
          HR > 1 & p.value < 0.05 ~ "Risk (HR>1, p<0.05)",
          HR < 1 & p.value < 0.05 ~ "Protective (HR<1, p<0.05)",
          TRUE ~ "Not significant (p≥0.05)"
        ),
        label_text = sprintf("%.2f (%.2f-%.2f)", HR, lower_95, upper_95)
      )

    min_lower <- min(data$lower_95, na.rm = TRUE)
    max_upper <- max(data$upper_95, na.rm = TRUE)
    x_min <- max(min_lower / 1.2, 0.01)
    x_max_expanded <- max_upper * 3.5
    text_x <- max_upper * 1.05

    breaks_raw <- scales::log_breaks(n = 6)(c(x_min, x_max_expanded))
    breaks <- sort(unique(c(breaks_raw, 1)))
    breaks <- breaks[breaks >= x_min & breaks <= x_max_expanded]

    bg_rects <- data %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(
        ymin = min(as.numeric(label)) - 0.5,
        ymax = max(as.numeric(label)) + 0.5,
        .groups = "drop"
      )

    bg_colors <- stats::setNames(
      RColorBrewer::brewer.pal(length(unique(data$Variable)), "Pastel1"),
      unique(data$Variable)
    )
    color_values <- c(
      "Risk (HR>1, p<0.05)" = "#E41A1C",
      "Protective (HR<1, p<0.05)" = "#377EB8",
      "Not significant (p≥0.05)" = "#999999"
    )

    p <- ggplot(data, aes(x = HR, y = label)) +
      geom_rect(
        data = bg_rects, inherit.aes = FALSE,
        aes(xmin = x_min, xmax = x_max_expanded, ymin = ymin, ymax = ymax, fill = Variable),
        alpha = 0.3, show.legend = FALSE
      ) +
      scale_fill_manual(values = bg_colors) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray30", linewidth = 0.6) +
      geom_errorbarh(aes(xmin = lower_95, xmax = upper_95, color = color_group),
        height = 0.2, size = 0.8
      ) +
      geom_point(aes(color = color_group), size = 3, shape = 15) +
      geom_text(aes(x = text_x, label = label_text), hjust = 0, size = 5, vjust = 0.5) +
      scale_color_manual(values = color_values, name = "") +
      ggplot2::coord_trans(x = "log10", xlim = c(x_min, x_max_expanded)) +
      scale_x_continuous(breaks = breaks, labels = sprintf("%.2f", breaks), expand = c(0, 0)) +
      labs(x = "Hazard Ratio (95% CI)", y = NULL, title = title) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        panel.grid.major.y = element_blank(),
        plot.margin = grid::unit(c(0.5, 0.1, 1, 0.5), "cm")
      ) +
      mytheme +
      guides(color = guide_legend(nrow = 1))

    CyDataPro::save_plot(p,
      filename = paste0(outdir, "/2_Fig/", outname),
      style = "g", width = coxWidth, height = coxHeight
    )
    invisible(p)
  }

  .plot_forest(uni_cox, "Univariate Cox Regression", "单变量cox_plot", coxWidth, coxHeight)

  # Multivariate Cox
  pd <- meta
  colnames(pd)[1] <- "Gene"

  multi_cox <- survival::coxph(
    survival::Surv(Time, Event) ~ Gene + Gender + Age +
      T_stage + M_stage + N_stage + Stage,
    data = pd
  )

  multi_result <- broom::tidy(multi_cox, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      Variable = dplyr::case_when(
        grepl("Gene", term) ~ "Gene",
        grepl("Gender", term) ~ "Gender",
        grepl("Age", term) ~ "Age",
        grepl("T_stage", term) ~ "T_stage",
        grepl("M_stage", term) ~ "M_stage",
        grepl("N_stage", term) ~ "N_stage",
        grepl("Stage", term) ~ "Stage"
      ),
      contrast_level = dplyr::case_when(
        term == "GeneLow" ~ "Low",
        term == "GenderMale" ~ "Male",
        term == "Age>=60" ~ ">=60",
        term == "T_stageT2" ~ "T2",
        term == "T_stageT3" ~ "T3",
        term == "T_stageT4" ~ "T4",
        term == "M_stageM1" ~ "M1",
        term == "N_stageN1" ~ "N1",
        term == "N_stageN2" ~ "N2",
        term == "N_stageN3" ~ "N3",
        term == "StageII" ~ "II",
        term == "StageIII" ~ "III",
        term == "StageIV" ~ "IV"
      ),
      ref_level = dplyr::case_when(
        Variable == "Gene" ~ "High",
        Variable == "Gender" ~ "Female",
        Variable == "Age" ~ "<60",
        Variable == "T_stage" ~ "T1",
        Variable == "M_stage" ~ "M0",
        Variable == "N_stage" ~ "N0",
        Variable == "Stage" ~ "I"
      )
    ) %>%
    dplyr::select(Variable, contrast_level, ref_level,
                  HR = estimate, lower_95 = conf.low,
                  upper_95 = conf.high, p.value)

  multi_result$Variable[1] <- gene
  data.table::fwrite(multi_result, file = paste0(outdir, "/2_Fig/多变量_cox_result.csv"))
  .plot_forest(multi_result, "Multivariable Cox Regression", "多变量cox_plot", coxWidth, coxHeight)

  # Nomogram
  multi_result$Variable[1] <- "Gene"
  p_core <- multi_result %>%
    dplyr::filter(Variable == "Gene") %>%
    dplyr::pull(p.value) %>%
    min()

  if (is.na(p_core)) {
    message(paste0("未找到 ", gene, " 变量，请检查变量名。"))
  } else if (p_core < 0.05) {
    message(sprintf(paste0(gene, " p = %.4f < 0.05，开始筛选显著变量并绘制列线图。"), p_core))
    sig_vars <- multi_result %>%
      dplyr::filter(p.value < 0.05) %>%
      dplyr::pull(Variable) %>%
      unique()

    if (length(sig_vars) == 0) {
      message("没有 p < 0.05 的变量，无法绘制列线图。")
    } else {
      message("显著变量: ", paste(sig_vars, collapse = ", "))
      var_map <- c(
        "Gene" = "Gene",
        "Gender" = "Gender",
        "Age" = "Age",
        "T_stage" = "T_stage",
        "M_stage" = "M_stage",
        "N_stage" = "N_stage",
        "Stage" = "Stage"
      )
      sig_original <- var_map[sig_vars]
      sig_original <- sig_original[!is.na(sig_original)]

      if (length(sig_original) == 0) stop("无法映射变量名")

      formula_reduced <- stats::as.formula(
        paste("survival::Surv(Time, Event) ~", paste(sig_original, collapse = " + "))
      )

      dd <- rms::datadist(pd)
      options(datadist = "dd")
      reduced_cph <- rms::cph(formula_reduced, data = pd, x = TRUE, y = TRUE, surv = TRUE)

      time_points <- c(1, 3, 5)
      pd$Time <- pd$Time / 365
      max_time <- max(pd$Time, na.rm = TRUE)
      time_points <- time_points[time_points <= max_time]

      grDevices::png(paste0(outdir, "2_Fig/Nomogram_plot.png"),
        width = coxWidth, height = coxHeight, res = 300, units = "in"
      )
      nom <- rms::nomogram(reduced_cph,
        funlabel = "3-year survival probability",
        lp = TRUE,
        conf.int = FALSE,
        maxscale = 100
      )
      plot(nom)
      grDevices::dev.off()

      grDevices::pdf(paste0(outdir, "2_Fig/Nomogram_plot.pdf"),
        width = coxWidth, height = coxHeight
      )
      plot(nom)
      grDevices::dev.off()

      message("列线图已保存为 nomogram.png 和 nomogram.pdf")
    }
  } else {
    message(sprintf(paste0(gene, " p = %.4f ≥ 0.05，不绘制列线图"), p_core))
  }

  invisible(multi_result)
}
