
######## Choosing penalty parameters using rGCV ##########

GCV <- function(Y, X, t, s, nubasis, ntbasis, nsbasis, lam){
  
  if(nubasis == 0){
    gcv <- rep(NA, length(lam))
    for (i in 1 : length(lam)){
      gcv[i] <- MEst(Y, X, t, s, nubasis, ntbasis, nsbasis, 0, lam[i])$GCV
    }
    optimal <- which(gcv == min(gcv))
    opgcv <- gcv[optimal]
    oplam <- lam[optimal]
    gcvlam <- list(optlam = oplam, Optgcv = opgcv, gcv = gcv)
  }
  else{
    gcv <- matrix(NA, length(lam), length(lam))
    for (i in 1 : length(lam)){
      for (j in 1 : length(lam)){
        gcv[i, j] <- MEst(Y, X, t, s, nubasis, ntbasis, nsbasis, lam[i], lam[j])$GCV
      }
    }
    optimal <- which(gcv == min(gcv), arr.ind=T)
    alphalam <- lam[optimal[1]]
    betalam <- lam[optimal[2]]
    rgcv <- gcv[optimal[1], optimal[2]] 
    gcvlam <- list(alphalam = alphalam, betalam = betalam, gcv = rgcv)
  }
  return(gcvlam)
}
