% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_utils.R
\name{getResiduals}
\alias{getResiduals}
\title{Get residuals from a given model object}
\usage{
getResiduals(modelSE)
}
\arguments{
\item{modelSE}{A model object output from any ChAI modeling function (whether on a single modality or the association between two modalities)}
}
\value{
A SummarizedExperiment object.
}
\description{
\code{getResiduals} Takes a model object and returns a SummarizedExperiment containing the residuals for each sample with associated metadata. Useful for modeling residuals.
}
\examples{
\dontrun{
  getResiduals(modelObj)
}

}
\keyword{model_results}
