#' A simulated data to illustrate the use of CoxFlex and BackSelection
#'
#' A dataframe containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format A data frame with 2307 rows and 8 variables:
#' \describe{
#'   \item{Id}{variable }
#'   \item{Event}{Binary indicator 1=event, 0=cencorsed}
#'   \item{Fup}{Length of follow up}
#'   \item{Start}{Start of the interval}
#'   \item{Stop}{Stop of the interval}
#'   \item{x1}{continuous time-varying covariate}
#'   \item{x2}{baseline continuous covariate generated from a standard normal distribution}
#'   \item{x3}{baseline continuouse covariate generated from a lognormal distribution}
#'}
#' @source This dataset is simulated using the PermAlgo package\url{http://cran.r-project.org/web/packages/PermAlgo/index.html}
"sampledata"
