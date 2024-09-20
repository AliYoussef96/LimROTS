#' BM21 Dataset case 1
#'
#' A `SummarizedExperiment` object with a dataset was obtained from the Metabolomics Workbench (ST002454). Different matrix ratios of human plasma and vegetable juice were analyzed, This dataset contains three different mixing ratios, each with three technical replicates.
#'
#' @name BM21.Case1
#' @docType data
#' @format An instance of the `SummarizedExperiment` class with the following assays:
#' \describe{
#'   \item{log2}{This assay includes log2 metabolites peaks}
#' }
#' The object also contains colData and rowData:
#' \describe{
#'   \item{colData}{A DataFrame with metadata for samples.}
#'   \item{rowData}{A DataFrame with metadata for metabolites peaks}
#' }
#'
#' The `colData` contains the following columns:
#' \describe{
#'   \item{SampleID}{Unique identifier for each sample.}
#'   \item{mixing.ratio}{Experimental condition or group for each sample , representing different different mixing ratios of human plasma and vegetable juice .}
#' }
#'
#' @source Bloody Mary (BM21)- Serial mixtures of vegetable juice/water and human plasma, Metabolomics Workbench (ST002454).
#' @references
#' This data is available at the NIH Common Fund's National Metabolomics Data Repository (NMDR) website, the Metabolomics Workbench, https://www.metabolomicsworkbench.org, where it has been assigned Project ID PR001582. The data can be accessed directly via it's Project DOI: 10.21228/M86Q7N This work is supported by NIH grant, U2C- DK119886.
#'

NULL
