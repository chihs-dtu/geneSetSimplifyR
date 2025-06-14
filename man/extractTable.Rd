% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractTable.R
\name{extractTable}
\alias{extractTable}
\title{Extract summary table of clustering}
\usage{
extractTable(
  geneSetsList,
  resolution = NULL,
  summarize = TRUE,
  removeClusterId = TRUE
)
}
\arguments{
\item{geneSetsList}{'gsList' object with clusters and labels.}

\item{resolution}{Numeric. Resolution to plot the table for.}

\item{summarize}{Logic indicating whether to summarize table by cluster (and source). Default is TRUE.}

\item{removeClusterId}{Logical, if TRUE the cluster number will be removed from each cluster label. Default: FALSE.}
}
\value{
Table with clustering summary
}
\description{
Extract summary table of clustering
}
\examples{
extractTable(exampleGsList, 1.1)
}
