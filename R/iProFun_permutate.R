#' iProFun false discovery rate assessment based on permutations
#'
#' @param ylist ylist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an outcome. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of ylist is
#' a list of mRNA, global protein and phosphoprotein.
#' @param xlist xlist is a list of data matrix. Each matrix is a data type or platform that
#' one would like to analyze as an predictor. The matrix is formated as each outcome variable
#' (such as gene) as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of xlist is
#' a list of cna and methylation.
#' @param covariates covariates is a list of data matrix. Each matrix is a data type or platform that
#' one would like to adjust as a covariate. The matrix is formated as each outcome variable
#'as a row, and each sample as a column. Multiple data matrixs need to have
#' shared sample ID and variable name for integrative iProFun analysis. Example of covariates is
#' a list of principle components.
#' @param pi prior probability
#' @param permutate_number Number of permutation, default 10
#' @param fdr False Discover Rate, default as 0.1
#' @param posterior Minimal Posterior Probabilty Cutoff, default as 0.75
#' @param filter is a vector with the same length of ylist, taking values of 1, -1, 0 or NULL. "1"
#' indicates filter in all positive results across multi-omic platforms. XXXX .... Jiayi!!!!
#' @param thresholds xxx Jiayi
#' @param seed Jiayi
#' @return list with 3 components
#' \item{Posterior Probability Cutoff:}{the cutoff values for each group based on permutation}
#' \item{iProFun Result:}{A table indicating whether a gene is identified by iProFun or not}
#' \item{iProFun Result (Negative/Positive):}{A table indicating whether a gene is identified by iProFun or not and whether the association is positive or not}
#' @export iProFun_permutate
#'
#' @examples
#' iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho), xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3), pi = rep(0.05, 3), permutate_number = 1, fdr = 0.1, posterior = 0.75, filter <- c(1,0))
iProFun_permutate = function(ylist, xlist, covariates,
                             pi = rep(0.05, 3), permutate_number = 10, fdr = 0.1, thresholds = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001))){
                             # ,filter = NULL, seed=?.Random.seed

  iprofun_result = iProFun(ylist = ylist, xlist = xlist, covariates = covariates, pi=pi, permutate = 0 )

  # need to add check function for filter function here (length of vectors equals to the length of xlist)

  for (i in 1: length(xlist)){
    if (filter[i] == 1){ # all positive
      assign(paste0("x_", i, "_filter_gene"), iprofun_result[[paste0("x_", i, "_xName")]][,1,drop = F][apply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]][,sapply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]],is.numeric)], 1, function(x) all(x > 0)),])
    }

    if (filter[i] == -1){  # all negative
      assign(paste0("x_", i, "_filter_gene"), iprofun_result[[paste0("x_", i, "_xName")]][,1,drop = F][apply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]][,sapply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]],is.numeric)], 1, function(x) all(x < 0)),])
    }

    if (filter[i] == 0){ # all positive or all negative
      assign(paste0("x_", i, "_filter_gene"), c(iprofun_result[[paste0("x_", i, "_xName")]][,1,drop = F][apply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]][,sapply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]],is.numeric)], 1, function(x) all(x > 0)),], iprofun_result[[paste0("x_", i, "_xName")]][,1,drop = F][apply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]][,sapply(iprofun_result[[paste0("x_", i, "_iProFun_gene_beta")]],is.numeric)], 1, function(x) all(x < 0)),]))
    }
    if (is.null(filter[i]) ){ # no requirement
      assign(paste0("x_", i, "_filter_gene"), iprofun_result[[paste0("x_", i, "_xName")]][,1,drop = F])
    }
  }
  # define the vectors to put the results
  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      for (p in 1:length(permutate_number)) {
        assign(
          paste0("x_", k, "_iprofun_perm_", j, "_fdr_", p, "_permutate"),
          vector("numeric", length(thresholds))
        )
      }
    }
  }

  Q = makeQ(1:length(ylist))

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_iprofun_perm_", j, "_fdr"),
             NULL)

    }
  }
  # The following codes run the the whole permutation "permutate_number" of times (each permutaion consist of "length(ylist)" times of sub-permutation)
  for (i in 1 : permutate_number) {
    # set.seed(seed)
    # In each permutation we run "length(ylist)" times of sub-permutation
    set.seed(i)
    for (j in 1:length(ylist)){
      assign(paste0("iprofun_perm_", j), iProFun(ylist = ylist, xlist = xlist, covariates = covariates, permutate = j))
    }

    # define the vectors to put the results
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
          assign(paste0("x_", k, "_iprofun_perm_", j, "_count"),
                 NULL)
          assign(paste0("x_", k, "_iprofun_perm_", j, "_original"),
                 NULL)
          assign(paste0("x_", k, "_iprofun_perm_", j, "_fdr"),
                 NULL)
        }
    }
    # define the vectors to put the results at the threshold level
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        for (p in 1:length(thresholds)) {
          assign(paste0("x_", k, "_iprofun_perm_", j, "_count_", p, "_thresholds"),
                 vector("numeric", length(thresholds)))
          assign(paste0("x_", k, "_iprofun_perm_", j, "_original_", p, "_thresholds"),
                 vector("numeric", length(thresholds)))
        }
      }
    }
    # calculate the number of genes significant at each threshold level
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        for (p in 1:length(thresholds)) {
          assign(paste0("x_", k, "_iprofun_perm_", j, "_count_", p, "_thresholds"),
                 sum(iprofun_result[[paste0("x_", k, "_xName")]][apply(eval(parse(text = paste0(
                   "iprofun_perm_", j
                 )))[[paste0("x_", k, "_iProFun_gene_posterior_probability_only")]][, Q[, j] ==
                                                                                      1], 1, sum) > thresholds[p], , drop = F][,1] %in% eval(parse(text = paste0("x_",k, "_filter_gene")))))
          assign(paste0("x_", k, "_iprofun_perm_", j, "_original_", p, "_thresholds"),
                 sum(iprofun_result[[paste0("x_", k, "_xName")]][apply(iprofun_result[[paste0("x_", k, "_iProFun_gene_posterior_probability_only")]][, Q[, j] ==
                                                                                                                                                       1], 1, sum) > thresholds[p], , drop = F][,1] %in% eval(parse(text = paste0("x_",k, "_filter_gene")))))
        }
      }
    }
    # paste the numbers at threshold level into a vector
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        for (p in 1:length(thresholds)) {
          assign(paste0("x_", k, "_iprofun_perm_", j, "_count"),
                 c(eval(parse(
                   text = paste0("x_", k, "_iprofun_perm_", j, "_count")
                 )), eval(parse(
                   text = paste0("x_", k, "_iprofun_perm_", j, "_count_", p, "_thresholds")
                 ))))

          assign(paste0("x_", k, "_iprofun_perm_", j, "_original"),
                 c(eval(parse(
                   text = paste0("x_", k, "_iprofun_perm_", j, "_original")
                 )), eval(parse(
                   text = paste0("x_", k, "_iprofun_perm_", j, "_original_", p, "_thresholds")
                 ))))
        }
      }
    }
    # calculte the empirical FDR at each threshold (vectorized division)
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        assign(paste0("x_", k, "_iprofun_perm_", j, "_fdr_", i, "_permutate"),
               eval(parse(text = paste0("x_", k, "_iprofun_perm_", j, "_count")))/eval(parse(text = paste0("x_", k, "_iprofun_perm_", j, "_original"))))
      }
    }
    print(paste("Finish Seed", i))
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      for (p in 1:permutate_number) {
        assign(paste0("x_", k, "_iprofun_perm_", j, "_fdr"),
               cbind(eval(parse(
                 text = paste0("x_", k, "_iprofun_perm_", j, "_fdr")
               )), eval(parse(
                 text = paste0("x_", k, "_iprofun_perm_", j, "_fdr_", p, "_permutate")
               ))))
      }
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_iprofun_perm_", j, "_fdr_overall"),
             apply(eval(parse(
               text = paste0("x_", k, "_iprofun_perm_", j, "_fdr")
             )), 1, mean, na.rm = T))
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_iprofun_", j, "_cutoff"),
             thresholds[min(which(eval(parse(
               text = paste0("x_", k, "_iprofun_perm_", j, "_fdr_overall")
             )) < fdr))])
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_iprofun_", j, "_cutoff"),
             thresholds[min(which(eval(parse(
               text = paste0("x_", k, "_iprofun_perm_", j, "_fdr_overall")
             )) < fdr))])
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("gene_x_", k, "_iprofun_", j),
             iprofun_result[[paste0("x_", k, "_xName")]][apply(eval(parse(text = paste0(
               "iprofun_perm_", j
             )))[[paste0("x_", k, "_iProFun_gene_posterior_probability_only")]][, Q[, j] == 1], 1, sum) > eval(parse(text = paste0("x_", k, "_iprofun_", j, "_cutoff"))), 1])
    }
  }

  for (i in 1: dim(Q)[1]){
    group = c(group, paste(Q[i,], collapse = ""))
  }

  cutoff_names <- NULL
  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
  cutoff_names <- c(cutoff_names, paste0("x_", j, "_y_", k, "_cutoff"))
    }
  }
  cutoff <- NULL
  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      cutoff <- c(cutoff, eval(parse(
        text = paste0("x_", k, "_iprofun_", j, "_cutoff"))))
    }
  }

  cutoff <- data.frame(cutoff_names, cutoff)

  result_permutation <- NULL

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      result_permutation <- append(result_permutation, list(eval(parse(text = paste0("gene_x_", k, "_iprofun_", j)))))
    }
  }

  names_result_permutation <- NULL
  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      names_result_permutation <- c(names_result_permutation, paste0("x_", j, "_y_", k))
    }
  }
  result_permutation <- append(result_permutation, list(cutoff))
  names(result_permutation) <- c(names_result_permutation, "cutoff")

  return(result_permutation)
}
