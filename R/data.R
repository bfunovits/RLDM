#' Blanchard/Quah (1989) dataset
#'
#' This dataset contains a (159 x 2)-dimensional matrix with quarterly data of real GDP growth rates (first column) and the unemployment rate in the USA.
#' It starts in the second quarter of 1948 and ends in the last quarter of 1987.
#' The script for transforming the raw data (as csv-file) to the matrix BQdata is available in the *data-raw* directory.
#' This dataset was used in [Gourieroux, Monfort, Renne (2019)]( https://doi.org/10.1093/restud/rdz028)
#'
#' @format A tibble (BQdata), a ts-object (BQdata_ts), or an xts-object (BQdata_xts) with 159 rows and 2 variables (plus a timestamp, so in total 3):
#' \describe{
#'   \item{date}{Timestamp as column of type *dttm* for the tibble, and as index for the *ts* and *xts* object.}
#'   \item{rGDPgrowth_demeaned}{The real GDP growth series is demeaned with respect to two different subperiods: Till the last quarter of 1973 and from the first quarter of 1974}
#'   \item{unemp_detrended}{The unemployment rate is detrended with respect to a linear trend (and an intercept).}
#' }
#'
#' @source <https://academic.oup.com/restud/advance-article/doi/10.1093/restud/rdz028/5490841#supplementary-data>
"BQdata"

#' @rdname BQdata
"BQdata_ts"

#' @rdname BQdata
"BQdata_xts"

#' Ramey/Shapiro dataset
#'
#' This dataset contains a (248 x 7)-dimensional matrix with quarterly US data from the first quarter
#' of 1947 till the last quarter of 2008.
#' The script for transforming the raw data (as csv-file) to the matrix RSdata is available
#' in the *data-raw* directory.
#' The data set is used to identify government spending shocks and in particular to
#' investigate the response of consumption and real wages with respect to these shocks.
#' All variables are orthogonalized with respect to an intercept and a linear trend.
#' See <https://doi.org/10.1093/qje/qjq008> for more details.
#'
#' @format A matrix with 248 rows and 7 variables (plus a timestamp, so in total 8 variables).
#' Versions as tibble (RSdata), ts (RSdata_ts), or xts (RSdata_xts) object are available:
#' \enumerate{
#'   \item Date column
#'   \item the logarithm of real per capita quantities of total government spending
#'   \item real GDP
#'   \item total hours worked
#'   \item nondurable plus services consumption
#'   \item private fixed investment
#'   \item real wage (or more precisely the nominal compensation in private business divided by the deflator in private business)
#'   \item tax rate
#' }
#' @source <https://econweb.ucsd.edu/~vramey/research.html>
"RSdata"

#' @rdname RSdata
"RSdata_ts"

#' @rdname RSdata
"RSdata_xts"

