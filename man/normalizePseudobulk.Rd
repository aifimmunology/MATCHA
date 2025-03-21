% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePseudobulkRNA.R
\name{normalizePseudobulk}
\alias{normalizePseudobulk}
\title{Normalize a ChAI scRNA object using DESEq2}
\usage{
normalizePseudobulk(rnaSE, sampleColumn = "Sample")
}
\arguments{
\item{rnaSE}{the output of makePseudobulkRNA a}

\item{sampleColumn}{The column of the Seurat object with sample information}
}
\value{
A SummarizedExperiment with normalized average expression
}
\description{
\code{normalizePseudobulk} Takes the output of makePseudobulkRNA and normalizes it.
}
\keyword{data_import}
