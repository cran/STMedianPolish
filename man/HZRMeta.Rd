\name{HZRMeta}
\alias{HZRMeta}
\docType{data}
\title{
Hydrogeological zone west of Meta river.
}
\description{
Map of hydrogeological zone west of Meta river. Spatial reference system: Datum Magna Sirgas Origen Bogota.
}
\usage{data(HZRMeta)}
\format{
  The format is:
Formal class 'SpatialPolygonsDataFrame' [package "sp"] 
}
\source{
\url{http://www.arcgis.com/home/item.html?id=103b63dcc9f448acbd63f22b728b1a02}
}
\examples{
library(sp)
data(HZRMeta)
Gridxy<- spsample(HZRMeta, cellsize=2000, n=300,"regular")
plot(Gridxy)
}
\keyword{datasets}
