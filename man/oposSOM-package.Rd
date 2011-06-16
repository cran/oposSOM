\name{oposSOM-package}
\alias{oposSOM-package}
\docType{package}
\title{
Package for analysis and visualization of metagenes
}
\description{
This package translates microarray expression data into metagenes and provides various visualizations and second level analyses. For a given data set, a standardized analysis output is created.
}
\details{
\tabular{ll}{
Package: \tab oposSOM\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2011-06-16\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Henry Wirth
Maintainer: Henry Wirth <wirth@izbi.uni-leipzig.de>
}
\references{
Wirth, Loeffler, v.Bergen, Binder: Expression cartography of human tissues using self organizing maps. (BMC Bioinformatics 2011)
}
\keyword{ package }
\examples{
## create data matrix
m = matrix( runif(1000), 100, 10 )
## run pipeline
run.pipeline( m )
}
