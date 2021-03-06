% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_flowgraph_modifiers.R
\name{fg_replace_meta}
\alias{fg_replace_meta}
\title{Replaces sample meta.}
\usage{
fg_replace_meta(fg, meta)
}
\arguments{
\item{fg}{flowGraph object.}

\item{meta}{A data frame containing meta data; see details in
\code{\link[flowGraph]{flowGraph-class}}.}
}
\value{
A flowGraph object with an updated sample meta.
}
\description{
Replaces sample meta in a given flowGraph object.
}
\examples{

 no_cores <- 1
 data(fg_data_pos30)
 fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
                 prop=FALSE, specenr=FALSE,
                 no_cores=no_cores)
 head(fg_get_meta(fg))

 new_df <- fg_data_pos30$meta
 new_df$id[1] <- "newID"

 fg <- fg_replace_meta(fg, new_df)
 head(fg_get_meta(fg))

}
\seealso{
\code{\link[flowGraph]{flowGraph-class}}
 \code{\link[flowGraph]{fg_get_meta}}
}
