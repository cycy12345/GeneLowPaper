#' Find Validating GEO Datasets for a Gene
#'
#' Scans pre-processed GEO datasets for a given cancer type and evaluates whether
#' the target gene shows significant differential expression and prognostic value.
#'
#' @param TCGA_pro Character. TCGA project abbreviation, e.g. \code{"LUAD"}.
#' @param gene Character. Target gene symbol.
#' @param data_dir Character. Root directory containing pre-processed GEO data.
#'   Default is \code{"F:/BioMed/预分析/data/GEO"}.
#'
#' @return Invisibly returns a list with two data frames: \code{DEG_table} and
#'   \code{KM_table}.
#'
#' @details
#' The function expects a specific folder structure under \code{data_dir/TCGA_pro}:
#' \itemize{
#'   \item \code{差异/} containing \code{*_DEG.Rds} files
#'   \item \code{预后/} containing \code{*_OS.Rds} files
#' }
#' Each \code{.Rds} file should contain a list with \code{exprmatix} (expression matrix)
#' and \code{sample} (sample information data frame).
#'
#' @examples
#' \dontrun{
#' BM_Gene_Find_datasets(TCGA_pro = "LUAD", gene = "IL36RN")
#' }
#' @export
BM_Gene_Find_datasets <- function(TCGA_pro = "LUAD", gene = NULL,
                                  data_dir = "F:/BioMed/预分析/data/GEO") {
  DEG_path <- paste0(data_dir, "/", TCGA_pro, "/差异")
  Surv_path <- paste0(data_dir, "/", TCGA_pro, "/预后")

  DEG_table <- data.frame()
  data_set <- list.files(DEG_path, "DEG.Rds")
  for (i in data_set) {
    geo_id <- stringr::str_split_1(i, "_")[1]
    geo <- readRDS(paste0(DEG_path, "/", i))
    if (length(intersect(rownames(geo$exprmatix), gene)) == 0) {
      print(paste0(geo_id, "中没有该基因"))
      next
    }
    exp <- geo$exprmatix[gene, ] %>% t() %>% as.data.frame()
    exp$group <- geo$sample$Group
    wilcox_res <- stats::wilcox.test(
      exp[, gene] ~ exp[, "group"],
      exact = FALSE, correct = FALSE
    )$p.value
    if (wilcox_res > 0.05) {
      print(paste0(geo_id, "差异不显著"))
      next
    }
    group_median <- exp %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(mean_value = mean(get(gene)))
    change <- ifelse(
      group_median[group_median$group == "Tumor", 2] -
        group_median[group_median$group == "Normal", 2] > 0,
      "Up", "Down"
    )
    tmp <- data.frame(data = geo_id, deg_p = wilcox_res, deg_change = change)
    DEG_table <- rbind(DEG_table, tmp)
  }
  print(DEG_table)

  KM_table <- data.frame()
  data_set <- list.files(Surv_path, "OS.Rds")
  for (i in data_set) {
    geo_id <- stringr::str_split_1(i, "_")[1]
    geo <- readRDS(paste0(Surv_path, "/", i))
    if (length(intersect(rownames(geo$exprmatix), gene)) == 0) {
      print(paste0(geo_id, "中没有该基因"))
      next
    }

    sampleinfo <- geo$sample %>%
      dplyr::filter(Group == "Tumor") %>%
      dplyr::select("Event", "Time")
    exp <- geo$exprmatix[gene, rownames(sampleinfo)] %>% t() %>% as.data.frame()
    Surv_data <- cbind(sampleinfo, exp)
    Surv_data <- na.omit(Surv_data)
    Surv_data$Event <- as.numeric(Surv_data$Event)
    res.cut <- survminer::surv_cutpoint(Surv_data,
      time = "Time",
      event = "Event",
      variables = names(Surv_data)[3:ncol(Surv_data)],
      minprop = 0.3
    )
    res_cat <- survminer::surv_categorize(res.cut)
    mysurv <- survival::Surv(res_cat$Time, res_cat$Event)

    group <- res_cat[, gene]
    if (length(unique(group)) != 2) {
      next
    }
    group <- factor(group, levels = c("low", "high"))
    surv_dat <- data.frame(group = group)
    fit <- survminer::surv_fit(mysurv ~ group, data = surv_dat)

    data.survdiff <- survival::survdiff(mysurv ~ group)
    p.val <- 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    if (p.val > 0.05) {
      next
    }

    HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) /
      (data.survdiff$obs[1] / data.survdiff$exp[1])

    tmpd <- data.frame(data = geo_id, HR = HR, KM_pvalue = p.val)
    KM_table <- rbind(KM_table, tmpd)
  }
  print(KM_table)
  invisible(list(DEG_table = DEG_table, KM_table = KM_table))
}


#' GEO Differential Expression Validation (Box Plot)
#'
#' Validates differential expression of a target gene in a GEO dataset using
#' a violin-box plot with significance testing.
#'
#' @param GSE Character. GEO dataset ID, e.g. \code{"GSE19188"}.
#' @param gene Character. Target gene symbol.
#' @param TCGA_pro Character. TCGA project abbreviation used to locate the data.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param plotWidth Numeric. Plot width. Default is 4.
#' @param plotHeight Numeric. Plot height. Default is 4.
#' @param data_dir Character. Root directory of GEO pre-processed data.
#'   Default is \code{"F:/BioMed/预分析/data/GEO"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the ggplot object.
#' @examples
#' \dontrun{
#' Fig2A_GEO_box(GSE = "GSE19188", gene = "TP53", TCGA_pro = "LUAD")
#' }
#' @export
Fig2A_GEO_box <- function(GSE = "GSE19188", gene = "TP53",
                          TCGA_pro = "LUAD", outdir = "./",
                          plotWidth = 4, plotHeight = 4,
                          data_dir = "F:/BioMed/预分析/data/GEO",
                          color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/2_Fig/")

  geo <- readRDS(paste0(data_dir, "/", TCGA_pro, "/差异/", GSE, "_DEG.Rds"))

  exp <- geo$exprmatix[gene, ] %>% t() %>% as.data.frame()
  exp$group <- geo$sample$Group

  p <- ggplot(data = exp, aes(x = group, y = get(gene), color = group)) +
    geom_violin(size = 1.3, trim = FALSE, position = position_dodge(0.9)) +
    ggbeeswarm::geom_quasirandom(aes(color = group), size = 1.6) +
    ggsignif::geom_signif(
      comparisons = list(c("Tumor", "Normal")),
      test = "t.test", colour = "black",
      test.args = list(var.equal = TRUE, alternative = "two.sided"),
      map_signif_level = TRUE, textsize = 8,
      tip_length = c(0, 0), size = 1,
      y_position = max(exp[, 1] - 0.3),
      vjust = 0
    ) +
    stat_summary(
      aes(group = group, color = group),
      position = position_dodge(0.9),
      width = 0.3, size = 0.8, fun = mean, geom = "errorbar",
      fun.min = function(x) quantile(x, 0.25),
      fun.max = function(x) quantile(x, 0.75)
    ) +
    stat_summary(
      aes(group = group, color = group),
      position = position_dodge(0.9),
      width = 0.3, size = 0.3, show.legend = FALSE,
      fun = mean, geom = "crossbar"
    ) +
    scale_color_manual(values = c("Tumor" = color_dis[1], "Normal" = color_dis[2])) +
    labs(title = GSE, y = paste0(gene, " Expression Level"), x = "") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "black", size = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      plot.title = element_text(hjust = 0.5)
    ) +
    mytheme

  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/Fig2A_GEO_box"),
    style = "g", width = plotWidth, height = plotHeight
  )
  invisible(p)
}


#' GEO Survival Validation (KM Plot)
#'
#' Validates the prognostic value of a target gene in a GEO dataset using
#' Kaplan-Meier survival analysis.
#'
#' @param GSE Character. GEO dataset ID, e.g. \code{"GSE31210"}.
#' @param gene Character. Target gene symbol.
#' @param TCGA_pro Character. TCGA project abbreviation.
#' @param plotWidth Numeric. Plot width. Default is 4.
#' @param plotHeight Numeric. Plot height. Default is 4.
#' @param outdir Character. Output directory. Default is \code{"./"}.
#' @param data_dir Character. Root directory of GEO pre-processed data.
#'   Default is \code{"F:/BioMed/预分析/data/GEO"}.
#' @param color_dis Character vector. Color palette. Defaults to \code{color_dis_default}.
#'
#' @return Invisibly returns the ggsurvplot object.
#' @examples
#' \dontrun{
#' Fig2A_GEO_KM(GSE = "GSE31210", gene = "TP53", TCGA_pro = "LUAD")
#' }
#' @export
Fig2A_GEO_KM <- function(GSE = "GSE31210", gene = "TP53",
                         TCGA_pro = "LUAD", plotWidth = 4, plotHeight = 4,
                         outdir = "./",
                         data_dir = "F:/BioMed/预分析/data/GEO",
                         color_dis = color_dis_default) {
  outpath <- paste0(outdir, "/2_Fig/")

  geo <- readRDS(paste0(data_dir, "/", TCGA_pro, "/预后/", GSE, "_OS.Rds"))
  sampleinfo <- geo$sample %>%
    dplyr::filter(Group == "Tumor") %>%
    dplyr::select("Event", "Time")
  exp <- geo$exprmatix[gene, rownames(sampleinfo)] %>% t() %>% as.data.frame()
  Surv_data <- cbind(sampleinfo, exp)
  Surv_data <- na.omit(Surv_data)
  Surv_data$Event <- as.numeric(Surv_data$Event)

  res.cut <- survminer::surv_cutpoint(Surv_data,
    time = "Time",
    event = "Event",
    variables = gene,
    minprop = 0.3
  )
  res_cat <- survminer::surv_categorize(res.cut)
  mysurv <- survival::Surv(res_cat$Time, res_cat$Event)

  group <- res_cat[, gene]
  group <- factor(group, levels = c("low", "high"))
  surv_dat <- data.frame(group = group)
  fit <- survminer::surv_fit(mysurv ~ group, data = surv_dat)

  data.survdiff <- survival::survdiff(mysurv ~ group)
  p.val <- 1 - stats::pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)

  p <- survminer::ggsurvplot(fit, data = surv_dat,
    conf.int = TRUE, conf.int.style = c("step"),
    censor = TRUE,
    palette = color_dis[2:1],
    legend.title = GSE,
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

  CyDataPro::save_plot(p,
    filename = paste0(outpath, "/Fig2B_GEO_KM"),
    style = "x", width = plotWidth, height = plotHeight
  )
  invisible(p)
}
