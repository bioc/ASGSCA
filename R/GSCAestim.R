GSCAestim<-function(data, W0, B0){

  
  Z0 <- as.matrix(data); W0 <- as.matrix(W0); B0 <- as.matrix(B0);
  
  N <- nrow(Z0)
  
  J <- nrow(W0); P <- ncol(W0)
  
  JP <- J + P
  
  W <- W0
  B <- B0
  C <- matrix(0, P, J)
  
  for (p in c(1:P)){
    WINDp <- which(W0[ , p] != 0)  #find the non null elements
    #print(WINDp)
    if (length(WINDp) > 1){    
      W0[WINDp, p] <- 99}
    }
  
  B0[B0 == 1] <- 99
  
  WIND <- which(W0 == 99, arr.ind = T)
  
  BIND <- which(B0 == 99, arr.ind = T)
  
  W[WIND] <- runif(nrow(WIND))         
  B[BIND] <- runif(nrow(BIND))
  
  V <- cbind(diag(J), W)
  
  Z <- scale(Z0)*sqrt(N)/sqrt(N-1)
  ZZ<-Z
  
  if (rankMatrix(t(Z)%*%Z) == J){
    Z <- chol(t(Z)%*%Z)}
  sizeZ <- dim(Z)[1]
  Gamma <- Z%*%W
  Psi <- Z%*%V
  
  f0 <- 100000000
  imp <- 100000000
  
  it <- 0             
  
    while(abs(imp) > 0.000001){
      it <- it+1
      #Step 1: Update B
      tr_b <- 0
        for (p in 1:P){ 
          ee <- matrix(0, 1, P)
          ee[p] <- 1
          LL <- diag(P)
          LL[p, p] <- 0
          b0 <- B0[, p]
          bindex_p <- which(b0 == 99)
          YY <- Gamma - Gamma%*%B%*%LL
          if (length(bindex_p)!=0){
             B[bindex_p, p] <- solve(t(Gamma[ , bindex_p])%*%Gamma[ , bindex_p])%*%((t(Gamma[ , bindex_p])%*%YY)%*%t(ee))
            tr_b <- tr_b + t(B[bindex_p, p])%*%B[bindex_p, p]
          }
        }
      A <- cbind(C, B)
      
      # Step 2: Update W
      tr_w <- 0
      for(p in 1:P){
        t <- J + p
        windex_p <- which(W0[, p] == 99)
        m <- matrix(0, 1, JP)
        m[t] <- 1
        a <- A[p, ]
        beta <- m - a
        H1 <- diag(P)
        H2 <- diag(JP)
        H1[p,p] <- 0
        H2[t,t] <- 0
        Delta <- W%*%H1%*%A - V%*%H2 
        Zp <- Z[ , windex_p]
        if (length(windex_p)!=0){        
          #browser()
          
          #theta <- solve(as.numeric(beta%*%t(beta))*(t(Zp)%*%Zp))%*%(t(Zp)%*%(Z%*%Delta)%*%t(beta))
           theta <- ginv(as.numeric(beta%*%t(beta))*(t(Zp)%*%Zp))%*%(t(Zp)%*%(Z%*%Delta)%*%t(beta))
         
          zw <- Zp%*%theta        
          theta <- sqrt(N)/sqrt(sum(zw^2))*theta
          W[windex_p, p] <- theta
          V[windex_p, t] <- theta
          tr_w <- tr_w + t(theta)%*%theta
        }
      }
      Gamma <- Z%*%W
      Psi <- Z%*%V
      dif <- Psi-Gamma%*%A                    
      f <- sum(diag(t(dif)%*%dif))                   
      imp <- f0-f
      f0 <- f
    }
  #NPAR <- nrow(WIND) + nrow(BIND)
    
  Westim <- W
  Bestim <- B
  
  #vecWestim <- W[WIND]
  #vecBestim <- B[BIND]
  return(list("Weight" = Westim,"Path" = Bestim))
}


