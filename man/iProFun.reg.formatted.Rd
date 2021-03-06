% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iProFun.reg.1y.R
\name{iProFun.reg.formatted}
\alias{iProFun.reg.formatted}
\title{Linear regression on all data types of outcomes (unformatted)}
\usage{
iProFun.reg.formatted(
  yList,
  xList,
  covariates,
  permutation.col = 0,
  var.ID = c("Gene_ID"),
  var.ID.additional = NULL,
  seed = NULL
)
}
\arguments{
\item{yList}{yList is a list of data matrix for outcomes.}

\item{xList}{xList is a list of data matrix for predictors.}

\item{covariates}{covariates is a list of data matrix for covariate.}

\item{permutation.col}{permutation.col provides the index of the data types that should be permuated.
permutation.col = 0 (default): No permuatation and analysis is on original data. 0 < permutate <= length of yList:
permuate the label of the corresponding data type in yList. For example, permutate =2, permute the y label of second
data matrix.}

\item{var.ID}{var.ID gives the variable name (e.g. gene/protein name) to match different data types.
If IDs are not specified, the first columns will be considered as ID variable.}

\item{var.ID.additional}{var.ID.additional allows to output additional variables in the outputs.
Often helpful if multiple rows (e.g. probes) are considered per gene to allow clear index of result.}

\item{seed}{seed allows users to externally assign seed to replicate results.}
}
\value{
list with the same length as xlist. Nested within each list, it contains
\item{betas_J:}{Coefficient estimate for predictors across J outcome data types}
\item{betas_se_J:}{Coefficent SE for predictors across J outcome data types}
\item{sigma2_J:}{Regrssion error term for predictors across J outcome data types}
\item{dfs_J:}{Regression degrees of freedom for predictors across J outcome data types}
\item{v_g_J:}{ (X^T X)^{-1} projection on predictors across J outcome data types}
\item{xName_J:}{Predictor name corresponds to each predictor-outcome pair across J outcome data types}
\item{yName_J:}{Outcome name corresponds to each predictor-outcome pair across J outcome data types}
}
\description{
Linear regression on all data types of outcomes with all types of DNA alterations (results formatted for iProFun input)
}
