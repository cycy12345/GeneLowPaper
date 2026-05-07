#' Single Gene Pan-Cancer Differential Expression
#'
#' Draws a pan-cancer boxplot comparing tumor vs normal expression of a target gene
#' across all TCGA cancer types using \pkg{TCGAplot}.
#'
#' @param gene Character. Target gene symbol.
#' @param plotWidth Numeric. Plot width in inches. Default is 12.
#' @param plotHeight Numeric. Plot height in inches. Default is 6.
#' @param outdir Character. Output directory path. Default is \code{"./"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the p-value data frame for significant cancers.
#' @details
#' Saves the plot to \code{outdir/1_Fig/Fig1A_Pan_expression} and a CSV of
#' significant p-values to \code{outdir/1_Fig/Fig1A_Pan_pvalue_sig.csv}.
#'
#' @examples
#' \dontrun{
#' Fig1A_Pancer(gene = "TP53", plotWidth = 12, plotHeight = 6, outdir = "./")
#' }
#' @export
Fig1A_Pancer <- function(gene, plotWidth = 12, plotHeight = 6,
                         outdir = "./",
                         color_dis = color_dis_default) {
  p <- TCGAplot::pan_boxplot(gene)
  data <- p$data

  group_counts <- data %>%
    dplyr::group_by(Cancer) %>%
    dplyr::summarise(n_groups = dplyr::n_distinct(Group), .groups = "drop") %>%
    dplyr::filter(n_groups > 1)

  test_data <- data %>% dplyr::filter(Cancer %in% unique(group_counts$Cancer))
  p_values <- test_data %>%
    dplyr::group_by(Cancer) %>%
    dplyr::summarise(
      p = stats::wilcox.test(get(gene) ~ Group, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::filter(p < 0.05)

  group_mean <- data %>%
    dplyr::group_by(Cancer) %>%
    dplyr::summarise(mean = stats::median(get(gene))) %>%
    dplyr::arrange(dplyr::desc(mean))
  data$Cancer <- factor(data$Cancer, levels = group_mean$Cancer)

  p <- ggplot(data, aes(x = Cancer, y = log(get(gene)), fill = Group)) +
    geom_boxplot(
      width = 0.6,
      aes(color = Group),
      linewidth = 1,
      outlier.shape = NA,
      alpha = 0.4
    ) +
    scale_fill_manual(values = c("Normal" = color_dis[2], "Tumor" = color_dis[1])) +
    scale_color_manual(values = c("Normal" = color_dis[2], "Tumor" = color_dis[1])) +
    ggpubr::stat_compare_means(
      aes(group = Group), label = "p.signif",
      method = "wilcox.test",
      label.y.npc = "top",
      vjust = 1,
      hide.ns = TRUE, size = 7
    ) +
    ggpubr::theme_test(base_rect_size = 1.5) +
    mytheme +
    labs(x = "", y = "Gene Expression Level", fill = "", color = "") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "top"
    )

  if (!dir.exists(paste0(outdir, "/1_Fig"))) {
    dir.create(paste0(outdir, "/1_Fig"))
  }

  CyDataPro::save_plot(p,
    filename = paste0(outdir, "/1_Fig/Fig1A_Pan_expression"),
    style = "g", width = plotWidth, height = plotHeight
  )
  data.table::fwrite(p_values,
    file = paste0(outdir, "/1_Fig/Fig1A_Pan_pvalue_sig.csv")
  )
  invisible(p_values)
}


#' Single Gene Pan-Cancer Survival Analysis
#'
#' Performs Kaplan-Meier survival analysis for a target gene across all TCGA cancer types.
#' Uses the optimal cutpoint to stratify patients into high/low expression groups.
#'
#' @param gene Character. Target gene symbol.
#' @param plotWidth Numeric. Summary plot width in inches. Default is 6.
#' @param plotHeight Numeric. Summary plot height in inches. Default is 8.
#' @param outdir Character. Output directory path. Default is \code{"./"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the KM result table.
#' @details
#' Individual KM plots are saved per cancer type under
#' \code{outdir/1_Fig/Fig1B_kmplot/}. A summary forest-style plot and a CSV
#' of results are also saved.
#'
#' @examples
#' \dontrun{
#' Fig1B_PanKM(gene = "TP53", plotWidth = 6, plotHeight = 8, outdir = "./")
#' }
#' @export
Fig1B_PanKM <- function(gene, plotWidth = 6, plotHeight = 8,
                        outdir = "./",
                        color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/1_Fig//Fig1B_kmplot/")

  tpm <- TCGAplot::get_all_tpm()
  meta <- TCGAplot::get_all_meta()
  Cancers <- TCGAplot::get_cancers() %>%
    as.data.frame() %>%
    dplyr::pull(Var1) %>%
    unique()

  KM_table <- data.frame()
  for (i in Cancers) {
    exprSet <- subset(tpm, Group == "Tumor" & Cancer == i) %>%
      dplyr::select(dplyr::all_of(gene)) %>%
      tibble::add_column(
        ID = stringr::str_sub(rownames(.), 1, 12)
      ) %>%
      dplyr::filter(!duplicated(ID)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("ID") %>%
      dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == i))) %>%
      t() %>%
      as.matrix()

    cl <- meta[colnames(exprSet), ]
    cl$exp <- as.numeric(exprSet[1, ])
    cl$time <- cl$time / 12

    res.cut <- tryCatch({
      survminer::surv_cutpoint(cl,
        time = "time",
        event = "event",
        variables = "exp",
        minprop = 0.3
      )
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(res.cut)) next

    res_cat <- survminer::surv_categorize(res.cut)
    mysurv <- survival::Surv(res_cat$time, res_cat$event)
    group <- res_cat[, "exp"]
    group <- factor(group, levels = c("low", "high"))
    surv_dat <- data.frame(group = group)
    fit <- survminer::surv_fit(mysurv ~ group, data = surv_dat)

    data.survdiff <- survival::survdiff(mysurv ~ group)
    p.val <- 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

    HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) /
      (data.survdiff$obs[1] / data.survdiff$exp[1])
    up95 <- exp(log(HR) + stats::qnorm(0.975) *
      sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
    low95 <- exp(log(HR) - stats::qnorm(0.975) *
      sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))

    tmpd <- data.frame(Cancer = i, HR = HR, up95 = up95, low95 = low95, pvalue = p.val)
    KM_table <- rbind(KM_table, tmpd)

    p <- survminer::ggsurvplot(fit, data = surv_dat,
      conf.int = TRUE, conf.int.style = c("step"),
      censor = TRUE,
      palette = color_dis[2:1],
      legend.title = i,
      ggtheme = ggpubr::theme_test(base_rect_size = 1.5) + mytheme,
      font.x = 15, font.y = 15, font.title = 15,
      legend.labs = c(
        paste0("Low ", "(", fit$n[1], ")"),
        paste0("High ", "(", fit$n[2], ")")
      ),
      font.lenend = 20, pval = TRUE,
      xlab = "OS_time(years)",
      ylab = "Survival probablity",
      break.x.by = 3,
      break.y.by = 0.2
    )

    if (!dir.exists(outpath)) {
      dir.create(outpath, recursive = TRUE)
    }
    if (p.val < 0.05) {
      CyDataPro::save_plot(p,
        filename = paste0(outpath, "/", i, "_sig"),
        style = "x", width = 6, height = 4
      )
    } else {
      CyDataPro::save_plot(p,
        filename = paste0(outpath, "/", i, "_not"),
        style = "x", width = 6, height = 4
      )
    }
    print(i)
  }

  K1 <- (KM_table$pvalue < 0.05) & (KM_table$HR > 1)
  K2 <- (KM_table$pvalue < 0.05) & (KM_table$HR < 1)
  KM_table$change <- ifelse(K1, "Risk", ifelse(K2, "Protective", "not"))
  message("Integration plotting!!!")

  data.table::fwrite(KM_table,
    file = paste0(outdir, "/1_Fig/Fig1B_Pan_km_result.csv")
  )

  pd <- KM_table
  pd$up95 <- ifelse(pd$up95 > max(pd$HR) + 5, max(pd$HR) + 5, pd$up95)
  pd$low95 <- ifelse(pd$low95 < min(pd$HR) - 5, min(pd$HR) - 5, pd$low95)

  p <- ggplot(pd, aes(HR, reorder(Cancer, HR))) +
    geom_point(size = 3, aes(color = change)) +
    geom_errorbar(aes(xmin = low95, xmax = up95, color = change),
      width = 0.25, cex = 1
    ) +
    labs(y = "Cancer", x = "HR", color = "") +
    geom_vline(aes(xintercept = 1), size = 0.6, linetype = "dashed", colour = "gray2") +
    scale_color_manual(values = c(
      "Risk" = color_dis[1],
      "Protective" = color_dis[2],
      "not" = "gray"
    )) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    theme(legend.position = "top") +
    mytheme

  CyDataPro::save_plot(p,
    filename = paste0(outdir, "/1_Fig/Fig1B_Pan_km_plot"),
    style = "g", width = plotWidth, height = plotHeight
  )
  invisible(KM_table)
}


#' Venn Diagram of Differential Expression and Prognostic Cancers
#'
#' Reads the outputs of \code{Fig1A_Pancer} and \code{Fig1B_PanKM} and draws
#' a Venn diagram showing the overlap between cancers with significant differential
#' expression and significant prognosis.
#'
#' @param outdir Character. Output directory path. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 6.
#' @param plotHeight Numeric. Plot height. Default is 6.
#' @param show_element Logical. Whether to show the intersection table on the plot.
#'   Default is \code{TRUE}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the intersection data frame.
#' @examples
#' \dontrun{
#' Fig1D_Venn(outdir = "./", plotWidth = 6, plotHeight = 6, show_element = TRUE)
#' }
#' @export
Fig1D_Venn <- function(outdir = "./", plotWidth = 6, plotHeight = 6,
                       show_element = TRUE,
                       color_dis = color_dis_default) {
  DEG <- data.table::fread(paste0(outdir, "/1_Fig/Fig1A_Pan_pvalue_sig.csv"),
    data.table = FALSE
  )
  KM <- data.table::fread(paste0(outdir, "/1_Fig/Fig1B_Pan_km_result.csv"),
    data.table = FALSE
  ) %>%
    dplyr::filter(change != "not")

  list_data <- list(
    `Differentially Expressed` = unique(DEG$Cancer),
    `Prognostically Significant` = unique(KM$Cancer)
  )
  table_data <- data.frame(`Inner Cancers` = Reduce(intersect, list_data))

  p <- ggvenn::ggvenn(list_data,
    show_percentage = FALSE,
    text_size = 8,
    set_name_color = color_dis[1:2],
    set_name_size = 5,
    fill_color = color_dis[1:2]
  )

  if (show_element) {
    p <- p + annotation_custom(
      grob = gridExtra::tableGrob(table_data, rows = NULL),
      xmin = 0, xmax = 0, ymin = -0.7, ymax = -0.7
    )
  }

  CyDataPro::save_plot(p,
    filename = paste0(outdir, "/1_Fig/Fig1D_Vennplot"),
    style = "g", width = plotWidth, height = plotHeight
  )
  data.table::fwrite(table_data,
    file = paste0(outdir, "/1_Fig/Fig1D_VennInnerCancer.csv")
  )
  invisible(table_data)
}
