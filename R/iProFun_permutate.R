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
#' @param fdr false discover rate, default as 0.1
#' @param PostCut PostCut specifies minimal posterior probabilty cutoff to be considered as
#' significant, default as 0.75
#' @param filter filter is a vector with the same length of xlist, taking values of 1, -1, 0 or NULL.
#' "NULL" is default imposes no filtering. 1" indicates that an association is considered for
#' significance only if in positive directions across all outcome platforms. "-1" indicates
#' that an association is considered
#' for significance only if in negative directions across all outcome platforms. "0" indicates
#' that an association is considered for significance only if association across all outcome
#' platforms preserve consistent directions (either positive or negative).

#' @param grids grids specify the searching grids to find significant associations
#' @param seed seed allows users to externally assign seed to replicate results.
#' @return list with 3 components
#' \item{Posterior Probability Cutoff:}{the cutoff values for each group based on permutation}
#' \item{iProFun Result:}{A table indicating whether a gene is identified by iProFun or not}
#' \item{iProFun Result (Negative/Positive):}{A table indicating whether a gene is identified
#' by iProFun or not and whether the association is positive or not}
#' @export iProFun_permutate
#'
#' @examples
#' iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho),
#' xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
#' pi = rep(0.05, 3), permutate_number = 1, fdr = 0.1, PostCut = 0.75, filter <- c(1,0)),
#' grids = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001)),
#' seed=123
iProFun_permutate = function(ylist, xlist, covariates,
                             pi = rep(0.05, 3), permutate_number = 10, fdr = 0.1, PostCut=0.75,
                             grids = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001)),
                             filter = NULL, seed=NULL,
                             ID=NULL, x.ID=NULL, y.ID=NULL,  sub.ID.common="TCGA",
                             colum.to.keep=c("phospho_ID", "Hybridization", "chr"),
                             missing.rate.filter=NULL){
  xlength=length(xlist)
  ylength=length(ylist)

  # Obtain posterior probabiliy from original data
  iprofun_result = iProFun(ylist = ylist, xlist = xlist, covariates = covariates,
                           pi=pi, permutate = 0, ID=ID, x.ID=x.ID, y.ID=y.ID, sub.ID.common = sub.ID.common,
                           colum.to.keep=colum.to.keep, missing.rate.filter=missing.rate.filter)
  post_original=vector("list", xlength)
  for (k in 1:xlength) {post_original[[k]]=iprofun_result[[k]]$PostProb}

  # filter genes

  x_filter_gene=vector("list", xlength)
  for (i in 1: xlength){
    x_filter_gene[[i]]= iprofun_result[[i]]$xName_J[,1,drop = F][apply(iprofun_result[[i]]$betas_J, 1, function(x) all(x > 0)),]
    if (filter[i] == 1){ # all positive
      x_filter_gene[[i]]= iprofun_result[[i]]$xName_J[,1,drop = F][apply(iprofun_result[[i]]$betas_J, 1, function(x) all(x > 0)),]
    }

    if (filter[i] == -1){  # all negative
      x_filter_gene[[i]]= iprofun_result[[i]]$xName_J[,1,drop = F][apply(iprofun_result[[i]]$betas_J, 1, function(x) all(x < 0)),]
    }

    if (filter[i] == 0){ # all positive or all negative
      x_filter_gene[[i]]= iprofun_result[[i]]$xName_J[,1,drop = F][apply(iprofun_result[[i]]$betas_J, 1, function(x) (all(x > 0) | all(x<0))==T ),]
    }
    if (is.null(filter[i]) ){ # no requirement
      x_filter_gene[[i]]=iprofun_result[[i]]$xName_J[,1,drop = F]
    }
  }

  # define the vectors to put the results
      count_perm_grid=vector("list", xlength);
     for (k in 1:xlength) {
         count_perm_grid[[k]]=vector("list", ylength)
     }


#   for (j in 1:ylength) {
#     for (k in 1:xlength) {
#       for (p in 1:permutate_number) {
#         assign(
#           paste0("x_", k, "_perm_outcome_", j, "_num_", p, "_permutate"),
#           vector("numeric", length(grids))
#         )
#       }
#     }
#   }
#
   Q = makeQ(1:ylength)
#
#   for (j in 1:ylength) {
#     for (k in 1:xlength) {
#       assign(paste0("x_", k, "_perm_outcome_", j, "_fdr"),
#              NULL)
#     }
#   }
  # The following codes run the the whole permutation "permutate_number" of times (each permutaion consist of "length(ylist)" times of sub-permutation)
  for (i in 1 : permutate_number) {

    set.seed(seed+i)
    post_perm=vector("list", xlength); for (k in 1:xlength) {post_perm[[k]]=vector("list", ylength)}
    for (j in 1:ylength){
      iprofun_perm = iProFun(ylist = ylist, xlist = xlist,
                                                 covariates = covariates, permutate = j,
                                                 pi=pi, ID=ID, x.ID=x.ID, y.ID=y.ID,
                                                 sub.ID.common = sub.ID.common,
                                                 colum.to.keep=colum.to.keep,
                                                 missing.rate.filter=missing.rate.filter)
      for (k in 1:xlength) {post_perm[[k]][[j]]=iprofun_perm[[k]]$PostProb}
    }


    # calculate the number of genes significant at each threshold level


    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        temp=NULL
        for (p in 1:length(grids)) {
        temp=c(temp, sum(apply(post_perm[[k]][[j]][, Q[, j] ==1],1,sum)>grids[p]))
        count_perm_grid[[k]][[j]]=rbind(count_perm_grid[[k]][[j]], temp)
        }
      }
    }


        for (p in 1:length(grids)) {

                 sum(iprofun_result[[k]]$xName_J
                     [apply(eval(parse(text = paste0("perm_outcome_", j )))
                            [[k]]$PostProb[, Q[, j] ==1], 1, sum) > grids[p],1] %in%
                   x_filter_gene[[k]] ) )

          assign(paste0("x_", k, "_perm_outcome_", j, "_original_", p, "_grids"),
                 sum(iprofun_result[[k]]$xName_J
                     [apply(iprofun_result[[k]]$PostProb
                            [, Q[, j] ==1], 1, sum) > grids[p],1] %in%
                       x_filter_gene[[k]] ))
        }
      }
    }
    # paste the numbers at threshold level into a vector
    for (j in 1:length(ylist)) {
      for (k in 1:length(xlist)) {
        for (p in 1:length(grids)) {
          assign(paste0("x_", k, "_perm_outcome_", j, "_count"),
                 c(eval(parse(
                   text = paste0("x_", k, "_perm_outcome_", j, "_count")
                 )), eval(parse(
                   text = paste0("x_", k, "_perm_outcome_", j, "_count_", p, "_grids")
                 ))))

          assign(paste0("x_", k, "_perm_outcome_", j, "_original"),
                 c(eval(parse(
                   text = paste0("x_", k, "_perm_outcome_", j, "_original")
                 )), eval(parse(
                   text = paste0("x_", k, "_perm_outcome_", j, "_original_", p, "_grids")
                 ))))
        }
      }
    }
    # calculte the empirical FDR at each threshold (vectorized division)
    for (j in 1:ylength) {
      for (k in 1:xlength) {
        assign(paste0("x_", k, "_perm_outcome_", j, "_fdr_", i, "_permutate"),
               eval(parse(text = paste0("x_", k, "_perm_outcome_", j, "_count")))
               /eval(parse(text = paste0("x_", k, "_perm_outcome_", j, "_original"))))
      }
    }
    print(paste("Finish permutation", i))
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
      assign(paste0("x_", k, "_perm_outcome_", j, "_fdr_overall"),
             apply(eval(parse(
               text = paste0("x_", k, "_perm_outcome_", j, "_fdr")
             )), 1, mean, na.rm = T))
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_outcome_", j, "_cutoff"),
             grids[min(which(eval(parse(
               text = paste0("x_", k, "_perm_outcome_", j, "_fdr_overall")
             )) < fdr))])
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("x_", k, "_outcome_", j, "_cutoff"),
             grids[min(which(eval(parse(
               text = paste0("x_", k, "_perm_outcome_", j, "_fdr_overall")
             )) < fdr))])
    }
  }

  for (j in 1:length(ylist)) {
    for (k in 1:length(xlist)) {
      assign(paste0("gene_x_", k, "_outcome_", j),
             iprofun_result[[paste0("x_", k, "_xName")]][apply(eval(parse(text = paste0(
               "iprofun_perm_", j
             )))[[paste0("x_", k, "_iProFun_gene_posterior_probability_only")]][, Q[, j] == 1], 1, sum) > eval(parse(text = paste0("x_", k, "_iprofun_", j, "_cutoff"))), 1])
    }
  }
  group <- NULL
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
