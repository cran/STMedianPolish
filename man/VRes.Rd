\name{VRes}
\alias{VRes}
\docType{data}
\title{
Precomputed Variogram for residuals of monthly precipitation
}
\description{
Precomputed Variogram for residuals of monthly precipitation \code{\link{Metadb}}. For this, the 'variogram' [package "gstat"] function is used. 
}
\usage{data(VRes)}
\format{
  The format is: 'StVariogram' 'data.frame'
}
\examples{
#Empirical variogram
#VRes = variogram(values ~ 1, rain_residual, cutoff=90000,tlags=0:24,width=2650,
#                 assumeRegular=TRUE, na.omit=TRUE)			   		   
data(VRes)
plot(VRes)
FitPar_st = function(p, gfn, v, trace = FALSE, ...) {
  mod = gfn(v$spacelag, v$timelag,p,  ...)
  resid = v$gamma - mod
  if (trace)
    print(c(p, MSE = mean(resid^2)))
  mean(resid^2)
}
ModSpatial = function(h,p){p[2]*(1-exp(-h/p[3]))}
ModTemporal = function(u,p){p[4]*(1-exp(-u/p[5]))+ p[6]*(1-cos(pi*u/180))+
							p[7]*abs(sin(pi*u/180))}
VariogST=function(h,u,p)
		{ModTemporal(u,p)+ModSpatial(h,p)+p[8]*ModTemporal(u,p)*ModSpatial(h,p)}
#Parametros Iniciales
p1<-c(2,14.5,13900,5.9,29,1.55,3.7,-0.07)
pars.st = optim(p1, FitPar_st, method =  "BFGS",
                gfn = VariogST, v = VRes, hessian=TRUE)
fit_Variog_ST<-VRes
fit_Variog_ST$gamma<-VariogST(VRes$spacelag, VRes$timelag, pars.st$par)
plot(fit_Variog_ST)
}
\keyword{datasets}
\references{
Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}

Pebesma, E.J. (2004). \emph{Multivariable geostatistics in S: the gstat package}. Computers & Geosciences, 30: 683-691 \href{https://CRAN.R-project.org/package=gstat}{[link]}
}




