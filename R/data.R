#' Sample mRNA data
#'
#' A mRNA dataset
#' rna
#'
#' @format A data frame with 676 rows and 526 variables:
#' \describe{
#'   \item{Gene_ID}{Gene name}
#'   \item{TCGA-04-1331}{mRNA measurement for TCGA-04-1331}
#'   \item{TCGA-04-1332}{mRNA measurement for TCGA-04-1332}
#'   ...
#' }
"rna"

#' Sample Copy Number Alternation data
#'
#' A mRNA dataset
#' cna
#'
#' @format A data frame with 676 rows and 526 variables:
#' \describe{
#'   \item{Gene_ID}{Gene name}
#'   \item{TCGA-04-1331}{CNA measurement for TCGA-04-1331}
#'   \item{TCGA-04-1332}{CNA measurement for TCGA-04-1332}
#'   ...
#' }
"cna"

#' Sample Methylation data
#'
#' A Methylation dataset
#' methy
#'
#' @format A data frame with 1103 rows and 528 variables:
#' \describe{
#'   \item{Gene_ID}{Gene name}
#'   \item{TCGA-04-1331}{Methylation measurement for TCGA-04-1331}
#'   \item{TCGA-04-1332}{Methylation measurement for TCGA-04-1332}
#'   ...
#' }
"methy"

#' Sample Protein data
#'
#' A protein dataset
#' methy
#'
#' @format A data frame with 1103 rows and 528 variables:
#' \describe{
#'   \item{Gene_ID}{Gene name}
#'   \item{TCGA-04-1331}{Protein measurement for TCGA-04-1331}
#'   \item{TCGA-04-1332}{Protein measurement for TCGA-04-1332}
#'   ...
#' }
"protein"

#' Sample Phospho data
#'
#' A phospho dataset
#' methy
#'
#' @format A data frame with 1591 rows and 70 variables:
#' \describe{
#'   \item{Gene_ID}{Gene name}
#'   \item{TCGA-04-1331}{Phospho measurement for TCGA-04-1331}
#'   \item{TCGA-04-1332}{Phospho measurement for TCGA-04-1332}
#'   ...
#' }
"phospho"

#' First 3 principle components for Phospho data
#'
#' First 3 principle components for Phospho data
#' phospho_pc_1_3
#'
#' @format A data frame with 68 rows and 4 variables:
#' \describe{
#'   \item{subject_id}{Subject ID name}
#'   \item{phospho_pc1}{First principle component for phospho data}
#'   \item{phospho_pc2}{Second principle component for phospho data}
#'   \item{phospho_pc3}{Third principle component for phospho data}
#' }
"phospho_pc_1_3"

#' First 3 principle components for protein data
#'
#' First 3 principle components for protein data
#' protein_pc_1_3
#'
#' @format A data frame with 68 rows and 4 variables:
#' \describe{
#'   \item{subject_id}{Subject ID name}
#'   \item{protein_pc1}{First principle component for protein data}
#'   \item{protein_pc2}{Second principle component for protein data}
#'   \item{protein_pc3}{Third principle component for protein data}
#' }
"protein_pc_1_3"

#' First 3 principle components for rna data
#'
#' First 3 principle components for rna data
#' rna_pc_1_3
#'
#' @format A data frame with 68 rows and 4 variables:
#' \describe{
#'   \item{subject_id}{Subject ID name}
#'   \item{rna_pc1}{First principle component for rna data}
#'   \item{rna_pc2}{Second principle component for rna data}
#'   \item{rna_pc3}{Third principle component for rna data}
#' }
"rna_pc_1_3"
