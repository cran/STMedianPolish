#' Plot Median polish multidimensional.
#'
#' Plot the effects of an additive model for multidimensional array, using Tukey's median polish procedure.
#' @method plot MedianPolishM
#' @param x object of class MedianPolishM.
#' @param \dots ignored.
#' @details The object of class MedianPolish has a list of the contributions of every effect over data. The graphic shows for each iteration, the behavior of these components. If the median polish is applied to data of class ConstructutMPst, this method has a specific graphic for data with space - time variability.
#'
#' @references Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}
#' @references Hoaglin, D. C., Mosteller, F., & Tukey, J. W. (Eds.). (2011). \emph{Exploring data tables, trends, and shapes} (Vol. 101). John Wiley & Sons.\href{http://www.wiley.com/WileyCDA/WileyTitle/productCd-047004005X.html}{[link]}
#' @examples A<-MedianPolishM(UCBAdmissions, eps=0.1, maxiter=2, na.rm=TRUE)
#' plot(A)
#' @importFrom reshape2 melt  
#' @export 

plot.MedianPolishM <-
function(x,...)
{
a<-length(dim(x$residuals))
b<-list()
for(i in 1:a){
d<-matrix(c(1:(dim(x$residuals)[i]*x$iter)),ncol=dim(x$residuals)[i])
for(j in 1:x$iter){
d[j,]<-MedianPolishM.default(x$data,x$eps,j,na.rm=TRUE)$effects[[i]]}
b[[i]]<-d
}
b$residuals<-x$residuals
b$iter<-x$iter
b$Graphic<-x$Gr
class(b) <- "plot.MedianPolishM"
b
}
