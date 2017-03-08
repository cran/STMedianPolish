#' Ordinary local Spatio - temporal Kriging
#'
#' Function for ordinary local spatio-temporal kriging
#' @usage krige0STlocalMP(data,newdata,p,model,k,stAni)
#'
#' @param data object of class 'STFDF' [package "spacetime"]. It must contain the spatio -- temporal coordinates and values.
#' @param newdata object of class 'STF' [package "spacetime"], It should contain the prediction location in space and time.
#' @param p parameters of the spatio - temporal covariance model. The first parameter must be nugget value.
#' @param model spatio -- temporal covariance model.
#' @param k defines the number of the input spatio -- temporal points that will be used to interpolate one new value.
#' @param stAni Constant of the spatio -- temporal anisotropy, assuming a metric spatio -- temporal space.
#' @return Table that contains the prediction and the prediction variance. 
#' @examples
#' library(spacetime)
#' library(sp)
#' library(gstat)
#' library(zoo)
#' library(maptools)
#' data(Metadb)
#' #records of the precipitation monthly from january 2007 to january 2010
#' Metadb<-Metadb[,c(1:4,89:125)]
#' x<-matrix(0,1,37)
#' for(i in 1:37){
#'   x[,i] <- 2007 + (seq(0, 36)/12)[i]
#' }
#' x<-as.Date (as.yearmon(x), frac = 1)
#' time = as.POSIXct(x, tz = "GMT")
#' 
#' MPST<-ConstructMPst(sqrt(0.5+Metadb[,-c(1:4)]),time,pts=Metadb[,2:4],Delta=c(7,6,5))
#' residual<-removetrendMPst(MPST,eps=0.01, maxiter=2)
#' rain.loc<-Metadb[,c("Station","East","North","Height")]
#' coordinates(rain.loc) = ~East+North+Height
#' proj4string(rain.loc) = CRS(proj4string(DemMeta))
#' rain_residual = stConstruct(data.frame(Res=residual[,7]), space = list(values = 1),
#'                             time, SpatialObj = rain.loc,interval=TRUE)
#' 
#' #NewData
#' data(HZRMeta)
#' polygon1 = polygons(HZRMeta)
#' Gridxy<- spsample(polygon1, cellsize=10000, n=1000,"regular")
#' Gridxyz<-data.frame(Gridxy,over(Gridxy,DemMeta))
#' colnames(Gridxyz)<-c("East", "North","height")
#' Grid_pred <- STF(sp=SpatialPoints(Gridxyz,CRS(proj4string(DemMeta))), time=time[c(18,19)])
#' 
#' #Product - sum covariance model generalized
#' p=c(2,12.98,13899.95,3.44,14.95,1.84,3.92,-0.07)
#' CS = function(h,p){p[2]*exp(-h/p[3])}
#' CT = function(u,p){p[4]*exp(-u/p[5])+ p[6]*cos(pi*u/180)+p[7]*(1-abs(sin(pi*u/180)))}
#' CST<-function(h,u,p){0.084*CT(u,p)+ 0.32*CS(h,p)+0.07*CT(u,p)*CS(h,p)}
#' data(VRes)
#' stAni<-estiStAni(VRes, interval=c(10, 100))
#' 
#' PredictValue<-krige0STlocalMP(data=rain_residual,newdata=Grid_pred,p,model=CST,k=10,stAni)
#' IDs = paste("ID",1)
#' mydata = data.frame(PredictValue[,5], ID=IDs)
#' wind.ST1 = STFDF(SpatialPixels(Gridxy),time[c(18,19)],mydata)
#' stplot(wind.ST1,col.regions=bpy.colors(40),par.strip.text = list(cex=0.7)
#'        ,main="Kriging ordinary residuals: Prediction surface")
#' @references Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}
#' @references Pebesma, E.J. (2004). \emph{Multivariable geostatistics in S: the gstat package}. Computers & Geosciences, 30: 683-691 \href{https://CRAN.R-project.org/package=gstat}{[link]}
#' @references Pebesma, E.J. (2012). \emph{spacetime: Spatio-Temporal Data in R.} Journal of Statistical Software, 51(7), 1-30.\href{https://CRAN.R-project.org/package=spacetime}{[link]}
#' @importFrom spacetime STFDF stplot
#' @importFrom sp spsample spDists proj4string
#' @importFrom gstat estiStAni
#' @importFrom maptools readShapePoly
#' @importFrom zoo as.yearmon
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats na.omit dist
#' @importFrom nabor knn
#' @export
#'
krige0STlocalMP <-
  function(data,newdata,p,model,k,stAni)
  {
    #Conditionals
    stopifnot(identical(proj4string(data@sp), proj4string(newdata@sp)))
    stopifnot(class(data@time) == class(newdata@time))
    stopifnot(k<=nrow(data@data))
    stopifnot(length(data@sp@coords.nrs)==2||length(data@sp@coords.nrs)==3)
    stopifnot(length(data@sp@coords.nrs)==ncol(newdata@sp@coords))
    #Evaluating number of terms of the coordinates
    ncoord<-length(data@sp@coords.nrs)
    output1 <- 0
    output2 <- 0
    C00<-model(h=0,u=0,p)+p[1]
    df= as.data.frame(data)[,c(1:ncoord,ncoord+2,ncoord+6)]
    dfcoord = as.data.frame(data)[,c(1:ncoord,ncoord+2)]
    dfcoord$time = as.numeric(dfcoord$time)*stAni
    
    query = as.data.frame(newdata)[,c(1:ncoord,ncoord+2)]
    query$time = as.numeric(query$time)*stAni
    
    nnbp <-t(apply(knn(dfcoord, query,k)[[1]],1,sort))
    query = as.data.frame(newdata)[,c(1:ncoord,ncoord+2)]
    
    pb <- txtProgressBar(min = 0, max = nrow(query), style = 3)
    for (i in 1:nrow(query)) {
      nghbrData <- df[nnbp[i,],]
      indexna<-which(is.na(nghbrData[,ncoord+2]))
      if(anyNA(nghbrData[,ncoord+2]))
      {
        nghbrData<-nghbrData[-indexna,]
      }
      #lags with regard s0,t0
      dist_s0 <- spDists(as.matrix(nghbrData[,1:ncoord]),query[i,1:ncoord])
      rez_t0<-as.matrix(abs(difftime(nghbrData[,ncoord+1],query[i,ncoord+1], units= c("days"))))
      #spatio temporal lags among neighbours
      sp_lags<- spDists(as.matrix(nghbrData[,1:ncoord]))
      temp_lags <- as.matrix(dist(as.matrix(nghbrData[,ncoord+1]))/86400)
      
      z<-nghbrData[,ncoord+2]
      v1<-matrix(rep(1,nrow(nghbrData)),nrow=nrow(nghbrData))
      R<-chol2inv(chol(model(h = sp_lags, u = temp_lags, p)))
      #======
      mgsl<-(1/(t(v1)%*%R%*%v1))*(t(v1)%*%R%*%z)
      z_star <- mgsl + (t(model(h = dist_s0, u = rez_t0, p)) %*% R)%*% matrix((z - v1%*%mgsl),ncol=1)
      #======
      ecm <-C00-t(model(h = dist_s0, u = rez_t0, p)) %*% R %*% model(h = dist_s0, u = rez_t0, p)+((1-t(v1)%*% R %*% model(h = dist_s0, u = rez_t0, p))^2)/(t(v1)%*%R%*%v1)
      output1 <- rbind(output1, z_star)
      output2 <- rbind(output2, ecm)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    outend <- data.frame(query, (output1)[-1,], output2[-1, ])
    colnames(outend)<-c(colnames(query),"predicted.values","error")
    return(outend)
  }