% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_flowgraph_features.R
\name{fg_feat_edge_prop}
\alias{fg_feat_edge_prop}
\title{Generates the proportion edge feature.}
\usage{
fg_feat_edge_prop(fg, no_cores = 1, overwrite = FALSE)
}
\arguments{
\item{fg}{flowGraph object.}

\item{no_cores}{An integer indicating how many cores to parallelize on.}

\item{overwrite}{A logical variable indicating whether to
overwrite the existing proportion edge feature if it exists.}
}
\value{
flowGraph object containing the proportion edge feature.
}
\description{
Generates the proportion edge feature and returns it
 inside the flowGraph object.
}
\details{
Given a flowGraph object, \code{fg_feat_edge_prop} returns the same
 flowGraph object with an additional proportions \code{prop} \code{edge}
 feature and its meta data. The proportions feature is
 made using the node count
 feature and is the cell count of each cell population (e.g. A+B+)
 over the cell count of its parent (e.g. A+);
 each edge then corresponds with such a relationship.
 The edge feature matrix has column names <from>_<to> e.g. A+_A+B+.
}
\examples{

 no_cores <- 1
 data(fg_data_pos30)
 fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
                 prop=FALSE, specenr=FALSE,
                 no_cores=no_cores)

 fg <- fg_feat_edge_prop(fg)

}
\seealso{
\code{\link[flowGraph]{flowGraph-class}}
 \code{\link[flowGraph]{fg_feat_node_prop}}
 \code{\link[flowGraph]{fg_feat_node_specenr}}
 \code{\link[flowGraph]{fg_add_feature}}
 \code{\link[flowGraph]{fg_get_feature}}
 \code{\link[flowGraph]{fg_rm_feature}}
 \code{\link[flowGraph]{fg_get_feature_desc}}
}
