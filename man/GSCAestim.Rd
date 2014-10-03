\name{GSCAestim}
\alias{GSCAestim}
\title{Structural Equation Models for multiple genotypes and multiple traits using Generalized Structured Component Analysis.}

\description{GSCAestim fits the Generalized Structured Component Analysis (GSCA) model to data on multiple genetic variants and multiple traits (see Romdhani et al., 2014). An Alternating Least-Squares algorithm (ALS) (de Leeuw, Young and Takane, 1976) is used to minimize a global least squares criterion. The ALS algorithm alternates between two main steps until convergence. In the first step, the weight coefficients are fixed, and the path coefficients are updated in the least-squares sense. In the second step, the weights are updated in the least-squares sense for fixed path coefficients. }


\usage{GSCAestim(data,W0,B0)}
\arguments{
\item{data}{Data frame containing the observed variables.}
\item{W0}{Matrix with 0's and 1's indicating connections between the observed variables (genotypes and traits) and the latent variables (genes and clinical pathways). The rows correspond to the observed variables in the same order as in data; the columns to the latent variables. A value of 1 indicates an arrow from the observed variable in the row to the latent variable in the column.}
\item{B0}{square matrix with 0's and 1's indicating connections among latent variables (genes and clinical pathways). Both rows and columns correspond to the latent variables. A value of 1 indicates an arrow directed from the latent variable in the row to the latent variable in the column.}
}

\value{Returns a list with 2 items.
\item{Weight}{Matrix of the same dimension as W0 with 1's replaced by weight coefficients estimates.}
\item{Path}{Matrix of the same dimension as B0 with 1's replaced by path coefficients estimates.}
}

%\section*{Note}
\author{Hela Romdhani, Stepan Grinek, Heungsun Hwang and Aurelie Labbe}

\references{
de Leeuw, J., Young, F. W., and Takane, Y. (1976). Additive structure in qualitative data: An alternating least squares method with optimal scaling features. \emph{ Psychometrika}, 41, 471-503. 

Romdhani, H., Hwang, H., Paradis, G., Roy-Gagnon, M.-H. and Labbe, A. (2014). Pathway-based Association Study of Multiple Candidate Genes and Multiple Traits Using Structural Equation Models. Submitted.
}

\examples{
#Scenario (g) in Romdhani et al. (2014): 4 SNPs mapped to 2 genes 
#and 4 traits involved in 2 clinical pathways.
#In total: 8 observed variables and 4 latent variables.
#One of the traits is involved in both clinical pathways.
#One gene is connected to one of the clinical pathways and the other
#to both of them.
data(GenPhen)
W0 <- matrix(c(rep(1,2),rep(0,8),rep(1,2),rep(0,8),rep(1,3),rep(0,7),rep(1,2)),nrow=8,ncol=4)
B0 <- matrix(c(rep(0,8),rep(1,2),rep(0,3),1,rep(0,2)),nrow=4,ncol=4)
res<-GSCAestim(data=GenPhen,W0,B0)
}
