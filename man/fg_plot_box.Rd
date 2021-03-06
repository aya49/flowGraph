% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_flowgraph_plots.R
\name{fg_plot_box}
\alias{fg_plot_box}
\title{Creates a boxplot of the values of one node/edge}
\usage{
fg_plot_box(
  fg,
  type = "node",
  index = 1,
  summary_meta = NULL,
  node_edge = 1,
  adjust_custom = "byLayer",
  p_thres = 0.05,
  filter_adjust0 = 0.5,
  filter_es = 0.5,
  filter_btwn_tpthres = 0.05,
  filter_btwn_es = 0.5,
  paired = FALSE,
  dotplot = TRUE,
  outlier = TRUE,
  all_labels = FALSE,
  show_mean = TRUE,
  main = NULL,
  path = NULL
)
}
\arguments{
\item{fg}{flowGraph object.}

\item{type}{A string indicating feature type the summary was created for
'node' or 'edge'.}

\item{index}{The user must provide \code{type} and
 additionally, one of \code{summary_meta} or \code{index}.

 \code{index} is an integer indicating the row in
 \code{fg_get_summary_desc(<flowGraph>)} of the corresponding type and
 summary the user would like to retrieve.}

\item{summary_meta}{The user must provide \code{type} and
 additionally, one of \code{summary_meta} or \code{index}.

 \code{summary_meta} is a list containing
 \code{feature} (feature name), \code{test_name} (summary statistic name),
 \code{class} (class), \code{label1}, and \code{label2} (class labels compared).
 See \code{\link[flowGraph]{fg_get_summary_desc}} for details.}

\item{node_edge}{An integer/index of or a string of the cell population (node) /
edge name (edge) the user wants to plot.}

\item{adjust_custom}{A function or a string indicating the
test adjustment method to use.
 If a string is provided, it should be one of
 \code{c("holm", "hochberg", "hommel",
 "bonferroni", "BH", "BY", "fdr", "none")} (see \code{p.adjust.methods}).
 If a function is provided, it should take as input
 a numeric vector and output the
 same vector adjusted.}

\item{p_thres}{A numeric variable indicating a p-value threshold}

\item{filter_adjust0}{A numeric variable indicating what percentage of
SpecEnr values compared (minimum) should be not close to 0.
Set to 1 to not conduct filtering.}

\item{filter_es}{A numeric variable between 0 and 1 indicating what the
Cohen's D value of the nodes/edges in question must be greater or
equal to, to be significant.}

\item{filter_btwn_tpthres}{A numeric variable between 0 and 1 indicating the
unadjusted T-test p-value threshold used to test whether the actual
and expected feature values used to calculate the specified SpecEnr
feature are significantly different for each sample class. Note this only
needs to be specified for SpecEnr features. Combined
with \code{filter_btwn_es}, we conduct three tests to understand if
there is an actual large difference between actual and expected features:
(1,2) T-test of significance between the actual and expected raw feature value
(e.g. proportion) for samples in each of the compared classes, (3) and the
T-test of significance between the differences of actual and
expected feature values of the two classes. If any two of the three tests
come out as insignificant, we set the p-value for the associated node/edge
to 1.}

\item{filter_btwn_es}{A numeric variable between 0 and 1 indicating what the
Cohen's D value of the nodes/edges in question must be greater or
equal to, to be significant -- see \code{filter_btwn_tpthres}.}

\item{paired}{A logical indicating whether the summary is paired.}

\item{dotplot}{A logical indicating whether or not to plot sample points.}

\item{outlier}{A logical indicating whether or not outliers should be plotted.}

\item{all_labels}{A logical indicating whether or not to plot samples of all
classes outside of just those used in the summary statistic test.}

\item{show_mean}{A logical indicating whether or not to label the mean.}

\item{main}{A string or the title of the plot; if left as \code{NULL},
a default title will be applied.}

\item{path}{A string indicating the path to where the function should save
the plot; leave as \code{NULL} to not save the plot. Static plots are saved
as PNG.}
}
\value{
A static boxplot.
}
\description{
Creates a boxplot comparing the
 features of samples belonging to different classes corresponding
 to an existing summary statistic using ggplot2.
}
\details{
The plot is made using the \code{ggplot2} package. The interactive
 version is the same as the static version, it is only here to support the
 shiny app.
}
\examples{

 no_cores <- 1
 data(fg_data_pos2)
 fg <- flowGraph(fg_data_pos2$count, class=fg_data_pos2$meta$class,
                 no_cores=no_cores)

 fg_plot_box(fg, type="node", summary_meta=NULL, adjust_custom="byLayer", index=1, node_edge=10)

}
\seealso{
\code{\link[flowGraph]{flowGraph-class}}
 \code{\link[flowGraph]{fg_plot}}
 \code{\link[flowGraph]{plot_gr}}
 \code{\link[flowGraph]{fg_get_feature}}
 \code{\link[flowGraph]{fg_get_summary}}
 \code{\link[flowGraph]{fg_plot_qq}}
}
