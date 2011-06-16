\name{run.pipeline}
\alias{run.pipeline}
\title{Run the oposSOM pipeline.}
\author{Henry Wirth}
\description{Translates input expression data into SOM metagene structure and outputs results in standardized scheme. }
\usage{
run.pipeline( indata, dataset.name = "Unnamed", dim.som1 = 20, dim.som2 = 20 )
}
\arguments{
\item{indata}{numerical data matrix containing genes as rows and samples as columns.}
\item{dataset.name}{name of the folder which will be created to contain the results.}
\item{dim.som1}{size of the primary SOM / resolution of the expression profiles: number of tiles equals dim.som1^2.}
\item{dim.som2}{resolution of the second level SOM: number of tiles equals dim.som2^2.}
}
\details{
If \code{indata} is supplied with ENSEMBL-IDs as rownames, geneset overrepresentation analysis is carried out using predefined GO sets according to GSEA database (http://www.broadinstitute.org/gsea).
}
\value{
The pipeline will create a folder in the working directory according to \code{dataset.name} containing all results. Additionally, a HTML summary file is created which contains brief descriptions and links to all results.
}
\examples{
## create data matrix
m = matrix( runif(1000), 100, 10 )
## run pipeline
run.pipeline( m )
}
\references{
Wirth, Loeffler, v.Bergen, Binder: Expression cartography of human tissues using self organizing maps. (BMC Bioinformatics 2011)
}
