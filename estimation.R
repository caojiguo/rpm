
######### computing the estimates using IRLS  #############

Huberfun <- function(u, tao = 1.345, escale){
  c <- tao * escale
  ifelse (abs(u) <= c, 0.5 * u^2, c * abs(u) - 0.5 * c^2)
}

irls <- function(Y, X, Pmat, iniPar, tao = 1.345, tol = 1e-04){
  weighfun <- function(u, tao, escale){
    c <- tao * escale
    ifelse (abs(u) <= c, 1, c / abs(u))
  }
  Par <- iniPar
  res <- as.vector(Y - X %*% Par)
  s <- mad(res)
  if (s == 0) 
    stop("cannot estimate scale: MAD is zero for this sample")
  repeat {
    W <- diag(weighfun(res, tao, s))
    Par1 <- solve(t(X) %*% W %*% X + Pmat, t(X) %*% W %*% Y)
    res1 <- as.vector(Y - X %*% Par1)
    s1 <- mad(res1)
    if (max(abs(Par1 - Par)) < tol) 
      break
    Par<- Par1
    res <- res1
    s <- s1
  }
  list(Par = Par, res = res, escale = s, weight = W)
} 

MEst <- function(Y, X, t, s, nubasis, ntbasis, nsbasis, lam0, lam1){
  n <- nrow(Y)
  nygrid <- length(t)
  nxgrid <- length(s)
  if (ncol(Y) != nygrid){
    stop("Number of observations in Y do not match with nygrid")
  }
  if (ncol(X) != nxgrid){
    stop("Number of observations in X do not match with nxgrid")
  }
  y <- as.vector(t(Y))
  tRange <- range(t)
  sRange <- range(s)
  
  if(nubasis == 0){
    betabasis.s <- create.bspline.basis(sRange, nsbasis)
    betabasis.t <- create.bspline.basis(tRange, ntbasis)
    betabasisMa.s <- eval.basis(s, betabasis.s)  
    betabasisMa.t <- eval.basis(t, betabasis.t)  
    betaX <- X %*% betabasisMa.s * (s[2] - s[1])  
    Xdesign <- kronecker(betaX, betabasisMa.t)   
    
    betaPenaMat0.t <- eval.penalty(betabasis.t, 0) 
    betaPenaMat2.t <- eval.penalty(betabasis.t, 2) 
    betaPenaMat0.s <- eval.penalty(betabasis.s, 0) 
    betaPenaMat2.s <- eval.penalty(betabasis.s, 2)
    betaPenaMat.t <- kronecker(betaPenaMat0.s, betaPenaMat2.t)
    betaPenaMat.s <- kronecker(betaPenaMat2.s, betaPenaMat0.t)
    lamPmat <- as.matrix(lam1 * (betaPenaMat.t + betaPenaMat.s))
    
    olsPar <- solve(t(Xdesign) %*% Xdesign + lamPmat, t(Xdesign)) %*% y
    HM <- irls(y, Xdesign, lamPmat, iniPar = olsPar)
    HMPar <- HM$Par
    HMbetaPar <- HMPar
    betaParMat <- matrix(HMbetaPar, nrow = ntbasis, ncol = nsbasis)
    HMbeta <- betabasisMa.t %*% betaParMat %*% t(betabasisMa.s)
    HMyfit <- as.vector(Xdesign %*% HMPar)
    HMres <- HM$res
    escale <- HM$escale
    Weigh <- HM$weight
    Hmat <- Xdesign %*% solve(t(Xdesign) %*% Weigh %*% Xdesign + lamPmat, t(Xdesign) %*% Weigh)
    rGCV <- (t(HMres) %*% Weigh %*% HMres) / (length(y) * (1-mean(diag(Hmat)))^2)
    
    HMest <- list(beta = HMbeta, yfit = HMyfit, residuals = HMres, GCV = rGCV)
  }
  else{
    alpbasis <- create.bspline.basis(tRange, nubasis)
    betabasis.s <- create.bspline.basis(sRange, nsbasis)
    betabasis.t <- create.bspline.basis(tRange, ntbasis)
    alpbasisMa <- eval.basis(t, alpbasis)  
    betabasisMa.s <- eval.basis(s, betabasis.s)  
    betabasisMa.t <- eval.basis(t, betabasis.t)  
    betaX <- X %*% betabasisMa.s * (s[2] - s[1])  
    X1design <- kronecker(betaX, betabasisMa.t)   
    X0design <- kronecker(rep(1, n), alpbasisMa)
    Xdesign <- cbind(X0design, X1design)
    
    alpPenaMat <- eval.penalty(alpbasis, 2)
    betaPenaMat0.t <- eval.penalty(betabasis.t, 0) 
    betaPenaMat2.t <- eval.penalty(betabasis.t, 2) 
    betaPenaMat0.s <- eval.penalty(betabasis.s, 0) 
    betaPenaMat2.s <- eval.penalty(betabasis.s, 2)
    
    betaPenaMat.t <- kronecker(betaPenaMat0.s, betaPenaMat2.t)
    betaPenaMat.s <- kronecker(betaPenaMat2.s, betaPenaMat0.t)
    lamPmat <- as.matrix(bdiag(lam0 * alpPenaMat, 
                               lam1 * (betaPenaMat.t + betaPenaMat.s)))
    
    olsPar <- solve(t(Xdesign) %*% Xdesign + lamPmat, t(Xdesign)) %*% y
    HM <- irls(y, Xdesign, lamPmat, iniPar = olsPar)
    HMPar <- HM$Par
    HMalphaPar <- HMPar[1 : nubasis]
    HMalpha <- alpbasisMa %*% HMalphaPar
    HMbetaPar <- HMPar[-c(1 : nubasis)]
    betaParMat <- matrix(HMbetaPar, nrow = ntbasis, ncol = nsbasis)
    HMbeta <- betabasisMa.t %*% betaParMat %*% t(betabasisMa.s)
    HMyfit <- as.vector(Xdesign %*% HMPar)
    HMres <- HM$res
    escale <- HM$escale
    Weigh <- HM$weight
    Hmat <- Xdesign %*% solve(t(Xdesign) %*% Weigh %*% Xdesign + lamPmat, t(Xdesign) %*% Weigh)
    rGCV <- (t(HMres) %*% Weigh %*% HMres) / (length(y) * (1-mean(diag(Hmat)))^2)
    
    HMest <- list(alpha = HMalpha, beta = HMbeta, yfit = HMyfit, residuals = HMres, GCV = rGCV)
  }
  
  return(HMest)
}  

