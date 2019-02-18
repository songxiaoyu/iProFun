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
#' @return list with 4 components
#' \item{fdr_grid:}{FDR by posterior probability cutoffs considered in grids}
#' \item{fdr_cutPob:}{the cutoff values for prespecified FDR rate and posterior probability for each X and Y pairs based on permutation}
#' \item{NoIdentified:}{the number of identified variables for each X and Y pair}
#' \item{Gene_fdr:}{A table indicating whether a gene is significantly identified by iProFun or not. "1" indicates the gene
#' is significantly identified; "0" indicates the gene is not.}
#' @export iProFun_permutate
#'
#' @examples
#' iprofun_permutate_result <- iProFun_permutate(ylist = list(rna, protein, phospho),
#' xlist = list(cna, methy), covariates = list(rna_pc_1_3, protein_pc_1_3, phospho_pc_1_3),
#' pi = rep(0.05, 3), permutate_number = 1, fdr = 0.1, PostCut = 0.75, filter <- c(1,0),
#' grids = c(seq(0.75, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001)),
#' seed=123)
iProFun_permutate = function(ylist, xlist, covariates, pi,
                             permutate_number = 10, fdr = 0.1, PostCut=0.75,
                             grids = seq(0.01, 0.99, by=0.01),
                             filter = NULL, seed=.Random.seed[1],
                             ID=NULL, x.ID=NULL, y.ID=NULL,  sub.ID.common="TCGA",
                             colum.to.keep=c("phospho_ID", "Hybridization", "chr"),
                             missing.rate.filter=NULL){
  xlength=length(xlist)
  ylength=length(ylist)

  # Obtain posterior probabiliy from original data
  iprofun_result = iProFun(ylist = ylist, xlist = xlist, covariates = covariates,
                           pi=pi, permutate = 0, ID=ID, x.ID=x.ID, y.ID=y.ID, sub.ID.common = sub.ID.common,
                           colum.to.keep=colum.to.keep, missing.rate.filter=missing.rate.filter, verbose=F)
  PostProb_original=vector("list", xlength)
  for (k in 1:xlength) {PostProb_original[[k]]=iprofun_result[[k]]$PostProb}

  # filter genes
  x_filter_gene=vector("list", xlength)
  for (k in 1: xlength){
    if (is.null(filter[k]) ){ # no requirement
      x_filter_gene[[k]]=seq(1, nrow(iprofun_result[[k]]$betas_J))
    } else if (filter[k] == 1){ # all positive
      x_filter_gene[[k]]= which(apply(iprofun_result[[k]]$betas_J, 1, function(x) all(x > 0)) )
    } else if (filter[k] == -1){  # all negative
      x_filter_gene[[k]]= which(apply(iprofun_result[[k]]$betas_J, 1, function(x) all(x < 0)))
    } else if (filter[k] == 0){ # all positive or all negative
      x_filter_gene[[k]]= which(apply(iprofun_result[[k]]$betas_J, 1, function(x) (all(x > 0) | all(x<0))==T ))
    }

  }


  # No of genes that passes threshold grids in original data
  Q = makeQ(1:ylength)
  count_original_grid=vector("list", xlength); for (k in 1:xlength) {count_original_grid[[k]]=vector("list", ylength)}

  for (k in 1:xlength) {
    for (j in 1:ylength) {
      Post_filter=PostProb_original[[k]][x_filter_gene[[k]],]
      temp=NULL
      for (p in 1:length(grids)) {
        temp=c(temp, sum(apply(Post_filter[, Q[, j] ==1],1,sum)>grids[p]))
      }
      count_original_grid[[k]][[j]]=rbind(count_original_grid[[k]][[j]], temp)
    }
  }

  # No of genes that passes threshold grids in permutation data
    count_perm_grid=vector("list", xlength);
    for (k in 1:xlength) {
      count_perm_grid[[k]]=vector("list", ylength)
    }


  # The following codes run the the whole permutation "permutate_number" of times (each permutaion consist of "length(ylist)" times of sub-permutation)
  for (i in 1 : permutate_number) {
    set.seed(seed+i)
    PostProb_perm=vector("list", xlength); for (k in 1:xlength) {PostProb_perm[[k]]=vector("list", ylength)}
    for (j in 1:ylength){
      iprofun_perm = iProFun(ylist = ylist, xlist = xlist,
                                                 covariates = covariates, permutate = j,
                                                 pi=pi, ID=ID, x.ID=x.ID, y.ID=y.ID,
                                                 sub.ID.common = sub.ID.common,
                                                 colum.to.keep=colum.to.keep,
                                                 missing.rate.filter=missing.rate.filter, verbose=F)
      for (k in 1:xlength) {PostProb_perm[[k]][[j]]=iprofun_perm[[k]]$PostProb}
    }


    # calculate the number of genes significant at each threshold level
    for (k in 1:xlength) {
      for (j in 1:ylength) {
        Post_filter=PostProb_perm[[k]][[j]][x_filter_gene[[k]],]
        temp=NULL
        for (p in 1:length(grids)) {
        temp=c(temp, sum(apply(Post_filter[, Q[, j] ==1],1,sum)>grids[p]))
        }
        count_perm_grid[[k]][[j]]=rbind(count_perm_grid[[k]][[j]], temp)
      }
    }

    print(paste("Finish permutation", i))
  }


 # calculate FDR based cutoff
    fdr_grid=vector("list", xlength);

    for (k in 1:xlength) {
      for (j in 1:ylength) {
        fdr_grid[[k]]=rbind(fdr_grid[[k]], apply( count_perm_grid[[k]][[j]], 2, mean, na.rm = T)/ count_original_grid[[k]][[j]])
        }
      colnames(fdr_grid[[k]])=paste0("Prob", grids)
      rownames(fdr_grid[[k]])=paste0("Y", 1:ylength)
    }

    AboveCut=which(grids>=PostCut)
    grid_PostCut= AboveCut[which.min(abs(AboveCut- PostCut))]

    fdr_cut=lapply(fdr_grid, function(f)
      apply(f, 1, function(f2) max(min(which(f2<fdr)), grid_PostCut) ))
    fdr_cutPob=lapply(1:xlength, function(f)
      sapply(1:ylength, function(f2) grids[fdr_cut[[f]][f2] ] ))

    NoIdentified=lapply(1:xlength, function(f)
      sapply(1:ylength, function(f2) count_original_grid[[f]][[f2]][fdr_cut[[f]][f2]]))


    Gene_fdr=vector("list", xlength);
    for (k in 1:length(xlist)) {
      for (j in 1:ylength) {
        temp=rep(0, nrow(iprofun_result[[k]]$xName_J))
        sig=intersect(x_filter_gene[[k]], which(apply(PostProb_original[[k]][, Q[, j] ==1],1,sum)>fdr_cutPob[[k]][j]))
        temp[sig]=1
        Gene_fdr[[k]]=cbind(Gene_fdr[[k]], temp)
      }
      colnames(Gene_fdr[[k]])=paste0("Y", 1:ylength)
      Gene_fdr[[k]]=cbind(iprofun_result[[k]]$xName_J,  Gene_fdr[[k]])
    }

    permutation_result=list(fdr_grid=fdr_grid, fdr_cutPob=fdr_cutPob, NoIdentified=NoIdentified, Gene_fdr=Gene_fdr )

  return(permutation_result)
}
