\name{GSCA}
\alias{GSCA}
\title{Association test for multiple genotypes and multiple traits using Generalized Structured Component Analysis (GSCA)}

\description{For a specified structural equation model with latent variables relating the traits and the genotypes, the function GSCA gives estimates of the parameters of the model and performs permutation tests for the association between multiple genotypes and multiple traits (see Romdhani et al., 2014).}

\usage{GSCA(data,W0, B0,latent.names=NULL,estim=TRUE,path.test=TRUE,path=NULL,nperm=1000)}
\arguments{
\item{data}{data frame containing the observed variables (genotypes and traits).}
\item{W0}{matrix with 0's and 1's indicating connections between the observed variables (genotypes and traits) and the latent variables (genes and clinical pathways). The rows correspond to the observed variables in the same order as in data; the columns to the latent variables. A value of 1 indicates an arrow from the observed variable in the row to the latent variable in the column.}
\item{B0}{square matrix with 0's and 1's indicating connections among latent variables (genes and clinical pathways). Both rows and columns correspond to the latent variables. A value of 1 indicates an arrow directed from the latent variable in the row to the latent variable in the column.}
\item{latent.names}{optional vector of characters containing names for the latent variables that will be displayed in the results. If NULL c("Latent1","Latent2",\ldots) will be used. Default is NULL.} 
\item{estim}{logical. If TRUE the estimates of the weight and path coefficient are returned. Default is TRUE.}
\item{path.test}{logical. If TRUE tests for path coefficients are performed. Default is TRUE.}
\item{path}{an optional matrix with 2 columns indicating particular connections to be tested. Each row contains the indices of the two latent variables (the gene and the clinical pathway) corresponding to the connection to be tested.  If NULL, the test is performed for all gene-clinical pathway connections specified in the model. Default is NULL.}
\item{nperm}{number of permutations. Default is 1000.}
}


\value{If estim is TRUE, returns a list with 2 items:
\item{Weight}{Matrix of the same dimension as W0 with 1's replaced by weight coefficients estimates.}
\item{Path}{Matrix of the same dimension as B0 with 1's replaced by path coefficients estimates.}
If, path.test is TRUE and path is NULL, the function returns a matrix of the same dimensions as B0 with 1's replaced by the corresponding p-values and 0's replaced by NA's. If path.test is TRUE and path is not NULL, only p-values for the specified path coefficients are returned.
}

%\section*{Note}
\author{Hela Romdhani, Stepan Grinek, Heungsun Hwang and Aurelie Labbe.}

\references{
Romdhani, H., Hwang, H., Paradis, G., Roy-Gagnon, M.-H. and Labbe, A. (2014). Pathway-based Association Study of Multiple Candidate Genes and Multiple Traits Using Structural Equation Models. Submitted.
}



\examples{
#Scenario (g) in Romdhani et al. (2014): 4 SNPs mapped to 2 genes and 4 
#traits involved in 2 clinical pathways 
#In total: 8 observed variables and 4 latent variables.
#One of the traits is involved in both clinical pathways.
#One gene is connected to one of the clinical pathways and
#the other to both of them.

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
