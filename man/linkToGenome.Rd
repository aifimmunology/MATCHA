% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePseudobulkRNA.R
\name{linkToGenome}
\alias{linkToGenome}
\title{\code{linkToGenome}}
\usage{
linkToGenome(rnaSE, TxDb = NULL, OrgDb = NULL, gene_format = "SYMBOL")
}
\arguments{
\item{rnaSE}{the output of makePseudobulkRNA}

\item{TxDb}{A Transcript database to be used for identifying gene locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.}

\item{OrgDb}{An Organism database to match up gene names and locations. Must be installed already. Both TxDb and OrgDb must be provided to annote the gene locations.}

\item{gene_format}{A string, describing the gene name format of the Seurat object. This is used to annotate gene loci, and convert IDs. See documentation for AnnotationDbi::mapIds to identify formats. Default is 'SYMBOL', but 'ENSEMBL' is also common.}
}
\value{
a list of metadata by celltype and sample
}
\description{
\code{linkToGenome} Helper function for makePseudobulkRNA. Links genes to genomic loci, and tosses poorly annotated transcripts.
}
\keyword{data_import}
