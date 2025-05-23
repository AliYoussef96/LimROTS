% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Check_SummarizedExperiment.R
\name{Check_SummarizedExperiment}
\alias{Check_SummarizedExperiment}
\title{Check if SummarizedExperiment or data is correct}
\usage{
Check_SummarizedExperiment(data.exp, meta.info, group.name)
}
\arguments{
\item{data.exp}{A matrix-like object or a \code{SummarizedExperiment} containing
the data to be analyzed.}

\item{meta.info}{Data frame. Metadata associated with the samples
(columns of \code{data.exp}). If \code{data.exp} is a \code{SummarizedExperiment},}

\item{group.name}{Character. Column name in \code{meta.info} that defines the
groups or conditions for comparison.}
}
\value{
a list of \code{data} , \code{groups} and \code{meta.info}
}
\description{
Check if SummarizedExperiment or data is correct
}
