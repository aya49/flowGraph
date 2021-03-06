% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_flowgraph_modifiers.R
\name{fg_merge}
\alias{fg_merge}
\title{Merges two flowGraph objects together.}
\usage{
fg_merge(
  fg1,
  fg2,
  method_sample = c("union", "intersect", "setdiff", "none"),
  method_phenotype = c("intersect", "setdiff", "none")
)
}
\arguments{
\item{fg1}{flowGraph object.}

\item{fg2}{flowGraph object.}

\item{method_sample}{A string indicating how samples from flowGraph objects
should be merged:
\itemize{
  \item{\code{union}: keep all samples from both flowGraph objects;
   in this case \code{method_phenotype} must be \code{intersect}.}
  \item{\code{intersect}: keep only samples that exist
   in both \code{fg1} and \code{fg2}.}
  \item{\code{setdiff}: keep only samples that exist
   in \code{fg1} and not in \code{fg2}.}
  \item{\code{none}: keep all samples in \code{fg1}.}
}}

\item{method_phenotype}{A string indicating how phenotypes from
flowGraph objects should be merged:
 \itemize{
  \item{\code{intersect}: keep only phenotypes that exist in both
   \code{fg1} and \code{fg2}.}
  \item{\code{setdiff}: keep only phenotypes that exist in
   \code{fg1} and not in \code{fg2}.}
  \item{\code{none}: keep all phenotypes in \code{fg1}.}
 }}
}
\value{
flowGraph object.
}
\description{
Merges two flowGraph objects together.
}
\details{
\code{fg_merge} is a generic function that merges the samples and
 phenotypes of two flowGraph objects.
 Note that if \code{method_sample="union"}
 then \code{method_phenotype} must be set to "intersect".
}
\examples{

 no_cores <- 1
 data(fg_data_pos30)
 fg0 <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
                 prop=FALSE, specenr=FALSE,
                 no_cores=no_cores)

 fg1 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[1:5])
 fg2 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[4:7])
 fg <- fg_merge(fg1, fg2, method_sample="intersect",
                          method_phenotype="intersect")
 fg_get_feature_desc(fg)

}
\seealso{
\code{\link[flowGraph]{flowGraph-class}}
 \code{\link[flowGraph]{fg_extract_samples}}
 \code{\link[flowGraph]{fg_extract_phenotypes}}
 \code{\link[flowGraph]{fg_merge_samples}}
}
