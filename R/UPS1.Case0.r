#' Spectronaut UPS1 Spiked Dataset case 0
#'
#' A \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} object containing DIA proteomics data from a UPS1-spiked E. coli protein mixture, processed using Spectronaut.
#'
#' @name UPS1.Case0
#' @docType data
#' @format An instance of the \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}} class with the following assays:
#' \describe{
#'   \item{log2}{This assay includes log2 protein intensities calculated by averaging the peptides derived from the same protein}
#' }
#' The object also contains colData and rowData:
#' \describe{
#'   \item{colData}{A DataFrame with metadata for samples.}
#'   \item{rowData}{A DataFrame with metadata for proteins.}
#' }
#'
#' The `colData` contains the following columns:
#' \describe{
#'   \item{SampleID}{Unique identifier for each sample.}
#'   \item{Conc}{Experimental condition or group for each sample , representing different conc. of UPS1-spiked proteins.}
#' }
#'
#' @source Generated with Spectronaut using a mixed mode acquisition method and FASTA mode for demonstration purposes.
#' @references
#' Gotti, C., Roux-Dalvai, F., Joly-Beauparlant, C., Mangnier, L., Leclercq, M., & Droit, A. (2022). DIA proteomics data from a UPS1-spiked E.coli protein mixture processed with six software tools. In Data in Brief (Vol. 41, p. 107829). Elsevier BV. https://doi.org/10.1016/j.dib.2022.107829
#'
#'

NULL
