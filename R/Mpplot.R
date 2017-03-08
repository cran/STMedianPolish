#' Traces of the space.
#'
#' Visualization of the spatial distribution according with three perspectives. Each face has the distribution for trace x, y and z (see \code{\link{ConstructMPst}}).
#' @param MpData object of class ConstructMPst.
#' @return Graphic of the three perspectives for space data "x", "y", "z", with divisions that contain the number of points in each quadrat.
#' @examples
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
#' Mpplot(MPST)
#'
#' @importFrom graphics points
#' @importFrom sp Polygon
#' @importFrom sp SpatialPolygons
#' @importFrom sp SpatialPoints
#' @importFrom sp over
#' @importFrom sp Polygons
#' @importFrom sp plot
#' @importFrom zoo as.yearmon
#' @importFrom graphics layout text axis mtext 
#'
#' @export

Mpplot <-
function(MpData){

if (class(MpData) != "ConstructMPst" )
stop("MpData must be a class ConstructMPst")

space<-MpData$pts
Deltasp<-MpData$Delta

Coord<-list()
Coord[[1]]<-c(1,2)
Coord[[2]]<-c(1,3)
Coord[[3]]<-c(2,3)

layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
for(k in 1:3){

pts<-scale(space[,Coord[[k]]])
Delta<-Deltasp[Coord[[k]]]
VDelta<-list()
Min<-c(1:2)
Max<-c(1:2)
ADelta<-c(1:2)

for(i in 1:2){
Min[i]<-min(pts[,i]) -(max(pts[,i])-min(pts[,i]))/(Delta[i]*2)}
for(i in 1:2){
Max[i]<-max(pts[,i]) +(max(pts[,i])-min(pts[,i]))/((Delta[i])*2)}

for(i in 1:2){
VDelta[[i]]<-(seq(Min[i],Max[i],(Max[i]-Min[i])/(Delta[i])))}

for(i in 1:2){
ADelta[i]<-(Max[i]-Min[i])/Delta[i]}

Sr<-list()
m=1
for(i in 1:Delta[2]){
for(j in 1:Delta[1]){
Sr[[m]] = Polygon(cbind(c(VDelta[[1]][j],VDelta[[1]][j]+ADelta[1],VDelta[[1]][j]+ADelta[1],VDelta[[1]][j],VDelta[[1]][j]),
    c(VDelta[[2]][i],VDelta[[2]][i],VDelta[[2]][i]+ADelta[2],VDelta[[2]][i]+ADelta[2],VDelta[[2]][i])))
m=m+1
}
}

Sr1<-list()
for(i in 1:(m-1)){Sr1[[i]] = Polygons(list(Sr[[i]]),i)}
SpP = SpatialPolygons(Sr1)

#Text counts
SrCount<-matrix(0,m-1,2)
m1=0
for(i in 1:Delta[2]){
for(j in 1:Delta[1]){
SrCount[m1+1,] = cbind(VDelta[[1]][j]+ADelta[1]/2,VDelta[[2]][i]+ADelta[2]/2)
m1=m1+1
}
}

ptsp = SpatialPoints(pts)
ct<-table(over(ptsp,SpP))
Vpos<-matrix(NA,m1)
Vpos2<-matrix(NA,m1)
Vpos[as.numeric(names(ct))]<-ct[]
Vpos2[as.numeric(names(ct))]<-8

#Labels Axis X
Srlabelx<-list()
dy<-abs((Min[2]-Max[2])/(Delta[2]*2.5))
Labelx<-c(1:Delta[1])
for(j in 1:Delta[1]){
Srlabelx[[j]] = Polygon(cbind(c(VDelta[[1]][j],VDelta[[1]][j]+ADelta[1],VDelta[[1]][j]+ADelta[1],VDelta[[1]][j],VDelta[[1]][j]),
    c(VDelta[[2]][Delta[2]+1],VDelta[[2]][Delta[2]+1],VDelta[[2]][Delta[2]+1]+dy,VDelta[[2]][Delta[2]+1]+dy,VDelta[[2]][Delta[2]+1])))
Labelx[j]<-j
}
Text_x<-cbind(VDelta[[1]][1:Delta[1]]+abs(Min[1]-Max[1])/((Delta[1])*2),rep((VDelta[[2]][Delta[2]+1]+dy/2)))

#Labels Axis Y
Srlabely<-list()
dx<-abs((Min[1]-Max[1])/(Delta[1]*2))
Labely<-c(1:Delta[2])
for(i in 1:Delta[2]){
Srlabely[[i]] = Polygon(cbind(c(VDelta[[1]][1],VDelta[[1]][1]-dy,VDelta[[1]][1]-dy,VDelta[[1]][1],VDelta[[1]][1]),
   c(VDelta[[2]][i],VDelta[[2]][i],VDelta[[2]][i]+ADelta[2],VDelta[[2]][i]+ADelta[2],VDelta[[2]][i])))
Labely[i]<-i
}
Text_y<-cbind(rep(VDelta[[1]][1]-dy/2),VDelta[[2]][1:Delta[2]]+abs(Min[2]-Max[2])/((Delta[2])*2))

Srsx = Polygons(Srlabelx, "s1")
Srsy = Polygons(Srlabely, "s2")
SpPLabel = SpatialPolygons(list(Srsx,Srsy),1:2)

plot(SpPLabel,col = c(2))
plot(SpP,col = Vpos2,add=T)
points(ptsp,pch=20)
text(SrCount,labels=Vpos,col="red",cex=0.8)

axis(1,range(pts[,1]),round(range(space[,Coord[[k]][1]]),2) )
mtext(colnames(pts)[1], side=1,line = 1,cex=0.8)

axis(2,range(pts[,2]),round(range(space[,Coord[[k]][2]]),2) )
mtext(colnames(pts)[2], side=2,line = 1,cex=0.8)

text(Text_x,labels=Labelx,cex=0.8)
text(Text_y,labels=Labely,cex=0.8)
}
}
