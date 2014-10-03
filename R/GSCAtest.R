GSCA<-function(data,W0, B0,latent.names=NULL,estim=TRUE,path.test=TRUE,path=NULL,nperm=1000){
  
 if (!is.logical(estim) | !is.logical(path.test)) stop('The arguments estim and path.test must be logical') else {if (!estim & !path.test) {stop('At least one of the two arguments estim and path.test must be equal to TRUE') }}

 if (!is.data.frame(data)){stop('The argument data must be a data frame')}
 
 Z0=as.matrix(data)

 if(ncol(Z0)<2){stop('The data must contain at least two observed variables')}
 
 N <- nrow(Z0)
 W0 <- as.matrix(W0)
 B0 <- as.matrix(B0)
 if(nrow(W0)!=ncol(Z0)) stop('The number of rows of W0 is not equal to the number of observed variables')
 if (nrow(B0)!=ncol(B0)) stop('B0 must be a square matrix')

 #path is a vector of 2 elements: gene index and Clinical pathway index
 if (is.null(latent.names)){
    latent.names=0
    for(i in 1:nrow(B0)){latent.names[i]=sprintf('Latent%d',i)}
    }  else {if (!is.vector(latent.names)) {stop('The argument latent.names must be a vector of characters')
 }else {if(length(latent.names)!=nrow(B0)){stop('The length of the argument latent.names does not match the dimension of B0')
 }else latent.names=as.character(latent.names)}}

 #Estimate the original parameters
 GSCAestim_res <- GSCAestim(Z0, W0, B0)
 Weight=GSCAestim_res$Weight
 colnames(Weight)=latent.names
 rownames(Weight)=colnames(data)
 
 Path= GSCAestim_res$Path
 colnames(Path)=latent.names
 rownames(Path)=latent.names
 

 if (path.test){

 BIND <- which(B0 == 1, arr.ind = T)

 if (is.null(path)){path=BIND
 } else if (!is.matrix(path) ){stop('The argument path must be NULL or a matrix with 2 columns')
 }  else if (ncol(path)!=2)     {stop('The argument path must be NULL or a matrix with 2 columns')
 } else{ 
  is.path=0
  for (i in 1:nrow(path)) {is.path[i]=sum(apply(BIND, 1, function(matrow, vec) isTRUE(all.equal(as.vector(matrow), vec)), as.vector(path[i,])) ) }
  if ( sum(is.path)!=nrow(path))  {stop('One or more rows of the argument path do not correspond to a  connection specified in the model')}
  }

 
  indGenes=unique(path[,1])
  Ngenes=length(indGenes)
  Genes=vector("list",Ngenes)
  indpath_per_gene=vector("list",Ngenes)

  for (g in indGenes){
    Genes[[g]]= which(W0[,g]==1)
    ind_path_tested=which(path[,1]==g)
    indpath_per_gene[[g]]=path[ind_path_tested,2]
 }

 Pval=matrix(nrow=nrow(B0),ncol=ncol(B0),dimnames=list(latent.names,latent.names))
   for(g in indGenes){
    gene=Genes[[g]]
    n <- length(gene)

    indpath_gene_g=indpath_per_gene[[g]]

      #Permutations
    Bperm=matrix(nrow=nperm,ncol=length(indpath_gene_g))

      for (perm in 1:nperm){
        Z0perm <- Z0
        indices <-  sample(1:N, N, replace = FALSE)

        Z0perm[ , gene] <- Z0[indices, gene]

        GSCAestim_res_perm <- GSCAestim(Z0perm, W0, B0)
        Bperm[perm,] <- GSCAestim_res_perm$Path[g,indpath_gene_g]
      }

      for(p in 1:length(indpath_gene_g)){
       Pval[g,indpath_gene_g[p]] <- mean(abs(GSCAestim_res$Path[g,indpath_gene_g[p]]) <= abs(Bperm[,p]))
      }
      rownames(Pval)=latent.names
      colnames(Pval)=latent.names
     }
 }

 if (estim & path.test) return (list("Weight" = Weight,"Path" = Path,"pvalues" =Pval))
 else if (estim & !path.test)  return(list("Weight" = Weight,"Path" = Path))
 else if (!estim & path.test)  return(Pval)
 }

