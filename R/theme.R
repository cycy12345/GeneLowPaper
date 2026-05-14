#' Default ggplot2 theme for GeneLowPaper
#'
#' A customized ggplot2 theme used across all plotting functions in the package.
#' Features bold serif fonts, thicker axis lines and ticks.
#'
#' @format A \code{theme} object.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   mytheme
#' }
#' @export
mytheme <- theme(
  axis.title = element_text(size = 15, face = "bold"),
  axis.line = element_line(linewidth = 1),
  axis.ticks = element_line(linewidth = 1.5),
  axis.text = element_text(size = 15, face = "bold", colour = "black"),
  text = element_text(size = 15, colour = "black", face = "bold")
)

#' Default discrete color palette
#'
#' A vector of 11 color hex codes used as the default discrete color scale
#' throughout the package.
#'
#' @format A character vector of length 11.
#' @examples
#' color_dis_default
#' @export
color_dis_default <- c(
  "#dd6b66", "#759aa0", "#e69d87",
  "#8dc1a9", "#ea7e53", "#eedd78",
  "#73a373", "#73b9bc", "#7289ab",
  "#91ca8c", "#f49f42"
)
