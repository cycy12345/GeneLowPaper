% Generated manually
\name{mytheme}
\alias{mytheme}
\docType{data}
\title{Default ggplot2 Theme for GeneLowPaper}
\description{
A customized \code{ggplot2} theme used across all plotting functions in the package.
Features bold serif fonts, thicker axis lines and ticks.
}
\usage{mytheme}
\format{A \code{theme} object.}
\examples{
\dontrun{
library(ggplot2)
ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  mytheme
}
}
\author{GeneLowPaper package}
