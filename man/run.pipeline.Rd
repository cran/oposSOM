\name{run.pipeline}
\alias{run.pipeline}
\title{Run the oposSOM pipeline.}
\author{Henry Wirth}
\description{This function realizes the complete pipeline functionality: single gene expression values are culstered to metagenes using a self-organizing map.
Based on these metagenes, visualizations (e.g. expression profiles) and downstreaming analysis approaches (e.g. hierarchical clustering, ICA or geneset overrepresentation analysis) are performed.
}
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
If \code{indata} is supplied with ENSEMBL-IDs as rownames, geneset overrepresentation analysis is carried out using predefined GO gene sets automatically downloaded from Ensembl database (http://www.ensembl.org). 
The results are contained in summary sheets of the integrated maps, e.g. the overexpression overview map.
}
\value{         
\item{}{The pipeline will create a folder in the working directory according to \code{dataset.name} containing all results. Graphical output is provided in respective files sole to ensure functionality in command line mode.}
\item{}{Additionally, a HTML summary file is created which contains brief descriptions and links to all results:}
\item{}{-  Expression profiles for all samples in different color modes(standard, loglog and WAD). }
\item{}{-  Supporting maps illustrating properties of the SOM as metagene population, variance and concordance.}
\item{}{-  Second level SOM providing sample-centered view on relations and similarities between the samples.}
\item{}{-  Metagene based downstreaming analysis: Hierarchical clustering, pairwise correlation matrix and spanning tree, and independent component analysis.}
\item{}{-  Report sheets for the spot summary maps containing various information of identified clusters of metagenes.}
}
\examples{
library( oposSOM )


#### Toy example ####

## create data
toy.data = matrix( runif(10000), 1000, 10 )

## run pipeline
run.pipeline( toy.data, dataset.name = "Example", dim.som1 = 20, dim.som2 = 20 )


#### Real data example - takes some time to process ####

## load data
#data( tissues )

## run pipeline
#run.pipeline( tissue.data, dataset.name="Tissue Example" )
}
\references{
Wirth, Loeffler, v.Bergen, Binder: Expression cartography of human tissues using self organizing maps. (BMC Bioinformatics 2011)

Wirth, v.Bergen, Binder: MALDI-typing of infectious algae of the genus Prototheca using SOM portraits. (Journal of microbiological methods 2012)
}
