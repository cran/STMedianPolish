#' Median polish multidimensional.
#'
#' An additive model for multidimensional array is fitted, using Tukey's median polish procedure.
#' @param data object of class array, table or matrix (see details).
#' @param \dots default arguments, see \code{\link{MedianPolishM.default}}
#' @details The function MedianPolishM is generic. See the documentation for \code{\link{MedianPolishM.default}} for further details.
#' @return An object of class medpolish with the following named components in a list:
#' @return \item{residuals}{the residuals.}
#' @return \item{overall}{the fitted constant term.}
#' @return \item{effects}{the fitted every dimensions effects to array multidimensional.}
#' @return \item{iter}{number of iterations used in the range maxiter.}
#'
#' @references Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}
#' @references Hoaglin, D. C., Mosteller, F., & Tukey, J. W. (Eds.). (2011). \emph{Exploring data tables, trends, and shapes} (Vol. 101). John Wiley & Sons.\href{http://www.wiley.com/WileyCDA/WileyTitle/productCd-047004005X.html}{[link]}
#' @importFrom reshape2 melt 
#' 
#' @export 
#' 
MedianPolishM <-
function(data,...) UseMethod("MedianPolishM")
