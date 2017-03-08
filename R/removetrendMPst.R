#' Median polish trend
#'
#' Direct method to remove the trend of spatio - temporal data througth median polish.
#' @usage removetrendMPst(MPST,eps=0.01, maxiter=10L)
#'
#' @param MPST object of class \code{\link{ConstructMPst}}
#' @param eps real number greater than \code{0}, default 0.01. A tolerance for convergence of median polish.
#' @param maxiter the maximum number of iterations, default 10.
#' @return data.frame with the following fields:
#' @return \item{ET}{indicate the time of a observation.}
#' @return \item{x}{ spatial coordinate x.}
#' @return \item{y}{ spatial coordinate y.}
#' @return \item{z}{ spatial coordinate z.}
#' @return \item{Trend}{trend calculated through median polish space - time.}
#' @return \item{Value}{observed values.}
#' @return \item{Residual}{\eqn{Residual = Value-Trend}.}
#' @details Following the Berke's approach (Berke, 2001) to remove spatial trend, Martinez et al.(2017) used a structure of four dimensions to apply median polish tecnique.
#'
#' For data \eqn{\left\{Y(\mathbf{s}_{abc},t), a=1,...,U; b=1,...,V; c=1,...,W, t=1,...,n  \right\} }, a spatial and temporal process is given by:
#' \deqn{Y(\mathbf{s}_{abc},t)= \mu_{y}(\mathbf{s}_{abc},t) + \delta _{abct}}
#' where
#' \deqn{\mu _{y}(\mathbf{s}_{abc},t)= \mu +\alpha _{a} + \beta _{b} + \xi _{c} + \tau _{t}}
#' and \eqn{\delta _{abct}} is a fluctuation arising from natural variability and from the measurement process. Additionally, \eqn{\mu} is an overall mean, \eqn{\alpha_{a}} is the \eqn{a}-th row effect, \eqn{\beta_{b}} is the effect \eqn{b}-th column effect, \eqn{\xi_{c}} is the \eqn{c}-th layer effect and \eqn{\tau _{t}} is the \eqn{t}-th time effect.
#'
#' @references Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}
#' @references Berke, O. (2001). \emph{Modified median polish kriging and its application to the wolfcamp - aquifer data.} Environmetrics, 12(8):731-748.\href{http://onlinelibrary.wiley.com/doi/10.1002/env.495/abstract}{[link]}
#' @examples
#' ## Not run:
#' library(zoo)
#' data(Metadb)
#' #records of monthly precipitation from january 2007 to january 2010
#' Metadb<-Metadb[,c(1:4,89:125)]
#' x<-matrix(0,1,37)
#' for(i in 1:37){
#'  x[,i] <- 2007 + (seq(0, 36)/12)[i]
#' }
#' x<-as.Date (as.yearmon(x), frac = 1)
#' time = as.POSIXct(x, tz = "GMT")
#'
#' MPST<-ConstructMPst(Metadb[,-c(1:4)],time,pts=Metadb[,2:4],Delta=c(7,6,5))
#' residual<-removetrendMPst(MPST,eps=0.01, maxiter=2)
#' ## End(Not run)
#'
#' @importFrom reshape2 melt
#' @importFrom zoo as.yearmon
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
removetrendMPst<-function(MPST,eps=0.01, maxiter=10L){

  MpSTData<-MedianPolishM.default(MPST$Value, eps, maxiter, na.rm=TRUE)
  pts<-MPST$pts

  Delta<-MPST$Delta
  ET<-rep(MPST$time,rep(dim(pts)[1],length(MPST$time)))
  Gridx<-rep(pts[,1],length(MPST$time))
  Gridy<-rep(pts[,2],length(MPST$time))
  Gridz<-rep(pts[,3],length(MPST$time))
  Gridxyz<-data.frame(Gridx,Gridy,Gridz)
  colnames(Gridxyz)<-colnames(pts)
  Grid2<-data.frame(pts,Trend=0)

		Grid_T<-data.frame(ET,Gridxyz,Tendencia=0)

  VDelta<-list()
		Min<-c(1:length(Delta))
		Max<-c(1:length(Delta))
		for(i in 1:length(Delta)){
		Min[i]<-min(pts[,i])}
		for(i in 1:length(Delta)){
		Max[i]<-max(pts[,i])}

#Grilla minimos y maximos del cubo

	GridReal<-pts
	VDeltaGrid<-list()
        MinGrid<-c(1:length(Delta))
        MaxGrid<-c(1:length(Delta))
        for(i in 1:length(Delta)){
        MinGrid[i]<-min(GridReal[,i])}
        for(i in 1:length(Delta)){
        MaxGrid[i]<-max(GridReal[,i])}

        for(i in 1:length(Delta)){
        VDelta[1:length(Delta)][[i]]<-(seq(Min[i],Max[i],abs(Min[i]-Max[i])/(Delta[i]-1)))}
		
        VDelta2<-VDelta

        #Remplazo de Minimos
        for(j in 1:length(Delta)){
        VDelta2[[j]][1] <- MinGrid[j]-0.1}

        #Remplazo de Maximos
        for(j in 1:length(Delta)){
        VDelta2[[j]][ length(VDelta2[[j]]) ] <- MaxGrid[j]+0.1}

  pb <- txtProgressBar(min = 0, max = Delta[1]-1, style = 3)

  for(i in 1:(Delta[1]-1)){
    for(j in 1:(Delta[2]-1)){
      for(k in 1:(Delta[3]-1)){
        sub_indexi <- which(Grid2[,1] >= VDelta2[][[1]][i] & Grid2[,1] <= VDelta2[][[1]][i+1] )
        sub_indexj <- which(Grid2[,2] >= VDelta2[][[2]][j] & Grid2[,2] <= VDelta2[][[2]][j+1] )
        sub_indexk <- which(Grid2[,3] >= VDelta2[][[3]][k] & Grid2[,3] <= VDelta2[][[3]][k+1] )
        inter1<-intersect(sub_indexj,sub_indexi)
        inter2<-intersect(inter1,sub_indexk)

        if(length(inter2)==0 ){
        }else
        {
          for(l in 1:length(inter2)){
            EX<-MpSTData$effects[[1]][i] + (Grid2[,1][inter2[l]]- VDelta[][[1]][i] )/ ( VDelta[][[1]][i+1] - VDelta[][[1]][i] )*(MpSTData$effects[[1]][i+1]-MpSTData$effects[[1]][i])
            EY<-MpSTData$effects[[2]][j] + (Grid2[,2][inter2[l]]- VDelta[][[2]][j] )/ ( VDelta[][[2]][j+1] - VDelta[][[2]][j] )*(MpSTData$effects[[2]][j+1]-MpSTData$effects[[2]][j])
            EZ<-MpSTData$effects[[3]][k] + (Grid2[,3][inter2[l]]- VDelta[][[3]][k] )/ ( VDelta[][[3]][k+1] - VDelta[][[3]][k] )*(MpSTData$effects[[3]][k+1]-MpSTData$effects[[3]][k])
	    Grid2[,4][inter2[l]]<-EX + EY + EZ  + MpSTData$overall
          }
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  time_c<-data.frame(MPST$time,MpSTData$effects[[4]])

  for(tim in unique(ET) )
  {
    sub_in <- which(Grid_T[,1] == tim)
    Grid_T[sub_in,5]<- Grid2[,4] + time_c[which(time_c[,1]==tim),2]
  }

  Grid_T$Value<-as.numeric(as.matrix(MPST$valuest))
  Grid_T$Residual<-Grid_T[,6]-Grid_T[,5]

  return(Grid_T)
}
