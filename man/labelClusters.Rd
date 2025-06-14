% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/labelClusters.R
\name{labelClusters}
\alias{labelClusters}
\title{Label clusters in clustered gene sets list object}
\usage{
labelClusters(
  geneSetsList,
  stopWords = geneSetSimplifyR::stopWords,
  nLabels = 10,
  verbose = TRUE
)
}
\arguments{
\item{geneSetsList}{'gsList' object that has been clustered with the 'clusterGeneSets()' function.}

\item{stopWords}{Vector with stop words. The default vector can be viewed with data(stopWords).
To add new stop words 'stopWords <- c(geneSetSimplifyR::stopWords, 'new word 1', 'new word 2', ...)}

\item{nLabels}{Number of top labels to keep. Default to 10}

\item{verbose}{a logic indicating whether to print progress statements. Defaults to TRUE.}
}
\value{
'gsList' object with labelled clusters for all clustered resolutions.
}
\description{
Label clusters in clustered gene sets list object
}
\examples{
geneSetsList <- initializeList(geneSetsList = exampleGeneSets, geneSetsDF = exampleEnrichment)
gsListObject <- clusterGeneSets(geneSetsList = geneSetsList)
gsListObject <- labelClusters(geneSetsList = gsListObject)
}
