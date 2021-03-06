% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_flowgraph_accessors.R
\name{fg_get_markers}
\alias{fg_get_markers}
\title{Retrieves the markers from a given flowGraph object.}
\usage{
fg_get_markers(fg)
}
\arguments{
\item{fg}{flowGraph object.}
}
\value{
A character vector containing the markers used in a flowGraph object.
}
\description{
Retrieves the markers from a given flowGraph object.
}
\examples{

 no_cores <- 1
 data(fg_data_pos30)
 fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
                 prop=FALSE, specenr=FALSE,
                 no_cores=no_cores)
 fg_get_markers(fg)

}
\seealso{
\code{\link[flowGraph]{flowGraph-class}}
}
