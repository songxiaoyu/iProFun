% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iProFun.detection.R
\name{iProFun.detection}
\alias{iProFun.detection}
\title{iProFun detection using multiple criteria}
\usage{
iProFun.detection(
  reg.all,
  eFDR.all,
  FWER.all = NULL,
  filter = c(0, 0),
  NoProbButFWERIndex = NULL,
  fdr.cutoff = 0.1,
  fwer.cutoff = 0.1,
  PostPob.cutoff = 0.75,
  xType,
  yType
)
}
\arguments{
\item{reg.all}{reg.all is the regression results for all outcomes from `iProFun.reg`.}

\item{eFDR.all}{eFDR.all is eFDR results for all outcomes from `iProFun.eFDR`.}

\item{FWER.all}{FWER.all is FWER results for all outceoms from `iProFun.FWER`. Can be NULL if FWER analysis was not performed.}

\item{filter}{filter is a vector with the same length of xList, taking values of 1, -1, 0 or NULL.
"NULL" is default imposes no filtering. 1" indicates that an association is considered for
significance only if its significant associations are positive across all outcome platforms. "-1" indicates
that an association is considered
for significance only if its significant associations are negative  across all outcome platforms. "0" indicates
that an association is considered for significance only if its significant association across all outcome
platforms preserve consistent directions (either positive or negative).}

\item{NoProbButFWERIndex}{NoProbXIndex allows users to provide the index for the predictor data type(s) that are not considered
for calculating posterior probabilities of association patterns. Default=NULL. All predictors go throught posterior
probability calculation.}

\item{fdr.cutoff}{Nominal false discover rate, default at 0.1}

\item{fwer.cutoff}{family wise error rate cutoff to be considered as
significant, default at 0.1.}

\item{PostPob.cutoff}{Minimal posterior probability cutoff to be considered as
significant, default at 0.75.}

\item{xType}{A vector of string for the data types of xList, such as "mutation", "CNV" and "methylation", for output presentation.}

\item{yType}{A vector of string for the data types of xList, such as "RNA", "protein" and "phospho", for output presentation.}
}
\value{
A output table that includes gene ID, other gene info if provided, predictor data type, outcome data
type, estimate, se, p-value from Student's t-test, FWER, Posterior Association Probability, eFDR,
directional filtering, and whether it's identified in iProFun or not, in a long format.
}
\description{
iProFun detection using multiple criteria
}
