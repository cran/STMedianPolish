\name{DemMeta}
\alias{DemMeta}
\docType{data}
\title{
Digital Elevation Model Resolution 90 meters.
}
\description{
Digital elevation model with resolution 250 meters of Hydrogeological zone west of Meta river. Spatial reference system: Datum Magna Sirgas Origen Bogota. 
}
\usage{data(DemMeta)}
\format{
Formal class: 'SpatialGridDataFrame' [package "sp"]
}
\source{
\url{http://www.jspacesystems.or.jp/ersdac/GDEM/E/index.html}
}
\examples{
library(sp)
data(DemMeta)
Gridxy<- spsample(DemMeta, cellsize=2000, n=300,"regular")
plot(Gridxy)
}
\keyword{datasets}
