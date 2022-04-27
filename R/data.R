#' Sample mRNA data
#'
#' A mRNA dataset
#' rna
#'
#' @format A data frame with 699 genes and 109 variables:
#' \describe{
#'   \item{geneSymbol}{Gene name}
#'   \item{C3L.00081}{mRNA measurement for sample C3L.00081}
#'   \item{C3L.00415}{mRNA measurement for sample C3L.00415}
#'   ...
#' }
"rna"

#' Sample Copy Number Alternation data
#'
#' A CNV dataset
#' cnv
#'
#' @format A data frame with 709 genes and 109 variables:
#' \describe{
#'   \item{geneSymbol}{Gene name}
#'   \item{C3L.00081}{CNV dosage measurement for sample C3L.00081}
#'   \item{C3L.00415}{CNV dosage measurement for sample C3L.00415}
#'   ...
#' }
"cnv"

#' Sample Mutation data
#'
#' A Mutation dataset
#' mut
#'
#' @format A data frame with 13 genes and 109 variables:
#' \describe{
#'   \item{geneSymbol}{Gene name}
#'   \item{C3L.00081}{Mutation (Yes=1; No=0) for sample C3L.00081}
#'   \item{C3L.00415}{Mutation (Yes=1; No=0) sample C3L.00415}
#'   ...
#' }
"mut"

#' Sample Protein data
#'
#' A protein dataset
#' methy
#'
#' @format A data frame with 551 genes and 109 variables:
#' \describe{
#'   \item{geneSymbol}{Gene name}
#'   \item{C3L.02665}{Protein measurement for sample C3L.02665}
#'   \item{C3L.01663}{Protein measurement for sample C3L.01663}
#'   ...
#' }
"protein"

#' Sample Phospho data
#'
#' A phospho dataset
#' phospho
#'
#' @format A data frame with 4168 phosphosites and 110 variables:
#' \describe{
#'   \item{geneSymbol}{Gene name}
#'   \item{id}{Phospho site identification}
#'   \item{C3L.02665}{Phospho measurement for C3L.02665}
#'   ...
#' }
"phospho"

#' Sample covariate data
#'

#' cov
#'
#' @format A data frame with 2 rows (covariates: age & female) and 108 variables:
#' \describe{
#'   \item{C3L.00081}{covariate measurements for sample C3L.00081}
#'   ...
#'   }
"cov"
