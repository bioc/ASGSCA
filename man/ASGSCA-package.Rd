\name{ASGSCA-package}
\alias{ASGSCA-package}
\alias{ASGSCA}
\docType{package}
\title{Association Studies for multiple SNPs and multiple traits using Generalized Structured Component Analysis}
\description{The package provides tools to model and test the association between multiple genotypes and multiple traits, taking into account the prior biological knowledge. Functional genomic regions, e.g., genes, and clinical pathways are incorporated in the model as latent variables that are not directly observed.  The method is based on Generalized Structured Component Analysis (GSCA).
}
\details{
\tabular{ll}{
Package: \tab ASGSCA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-30\cr
License: \tab GPL-3\cr
}

}
\author{
Hela Romdhani, Stepan Grinek, Heungsun Hwang and Aurelie Labbe.

Maintainer: Hela Romdhani <hela.romdhani@mcgill.ca>
}
\references{
Romdhani, H., Hwang, H., Paradis, G., Roy-Gagnon, M.-H. and Labbe, A. (2014). Pathway-based Association Study of Multiple Candidate Genes and Multiple Traits Using Structural Equation Models. Submitted.
}
\examples{
data(GenPhen)
W0 <- matrix(c(rep(1,2),rep(0,8),rep(1,2),rep(0,8),rep(1,3),rep(0,7),rep(1,2)),nrow=8,ncol=4)
B0 <- matrix(c(rep(0,8),rep(1,2),rep(0,3),1,rep(0,2)),nrow=4,ncol=4)

#Estimation only
GSCA(GenPhen,W0, B0,estim=TRUE,path.test=FALSE)
#Estimation and test for all the path coefficients in the model
GSCA(GenPhen,W0, B0,estim=TRUE,path.test=TRUE)
#Test only
GSCA(GenPhen,W0, B0,estim=FALSE,path.test=TRUE)
#Give names to the latent variables
GSCA(GenPhen,W0, B0,latent.names=c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2"),
estim=TRUE,path.test=TRUE)
#Testing only a subset of path coefficients
GSCA(GenPhen,W0, B0,estim=FALSE,path.test=TRUE,path=matrix(c(1,2,3,4),ncol=2))
}
