# `iProFun`

## An *i*ntegrative analysis tool to screen for *Pro*teogenomic *Fun*ctional traits perturbed by DNA copy number alterations (CNA) and DNA methylation

The goal is to characterize functional consequences of DNA copy number and methylation alterations in tumors and to facilitate screening for cancer drivers contributing to tumor initiation and progression, since CNAs and DNA methylations that preserve functional consequences are more likely to be cancer drivers. We applied iProFun to the ovarian high-grade serous carcinoma tumor data from The Cancer Genome Atlas and Clinical Proteomic Tumor Analysis Consortium, and identified a large number of CNAs and methylations affecting the RNA expression and/or protein/phosphoprotein abundances. The power gain is notable. For example, iProFun identified 130 genes whose CNAs were associated with phosphoprotein abundances by leveraging mRNA expression level and global protein abundance information, while analyses based on phosphoprotein data alone identified none. Genes merged from iProFun results could serve as potential drug targets for ovarian cancer.

This paper describes the functionality in the package.

---------------
  
### Installation
  
To install the latest patched version directly from Github, please use `devtools::install_github("iProFun")` for `iProFun` 

To install the developer version with experimental features directly from Github, please use `devtools::install_github("iProFun", ref="devel")`.
