
require(fda)
require(refund)
require(MASS)

source("generateSample.R")
source("estimation.R")
source("rGCV.R")

alpha <- function(t) 2 * exp(-(t - 1)^2)
beta <- function(t, s) {
  4 * cos(2 * pi * t) * sin(pi * s)
}


simulationNorm <- function(n){
  nxgrid <- 50
  nygrid <- 60
  s <- (1 : nxgrid) / nxgrid
  t <- (1 : nygrid) / nygrid
  alpha.t <- alpha(t)
  beta.ts <- outer(t, s, beta)
  
  sd.error = 0.01
  data <- NormSample(n, NumberOfphi = 10, ygrid = t, xgrid = s, sd.error)
  Y <- data$Y
  X <- data$X
  Ytrue <- data$Ytrue
  
  nubasis = 15
  ntbasis = 10
  nsbasis = 10
  lam <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 10)
  GCV <- GCV(Y, X, t, s, nubasis, ntbasis, nsbasis, lam)
  lam0 <- GCV$alphalam
  lam1 <- GCV$betalam
  
  HMfit <- MEst(Y, X, t, s, nubasis = 15, ntbasis = 10, nsbasis = 10, lam0, lam1)
  HMalpha <- HMfit$alpha
  HMbeta <- HMfit$beta
  HMISE_alpha <- mean((alpha.t - HMalpha)^2)
  HMISE_beta <- mean((beta.ts - HMbeta)^2)
  
  ntest = 100
  dataTest <- NormSample(ntest, NumberOfphi = 10, ygrid = t, xgrid = s, sd.error)
  Ytest <- dataTest$Ytrue
  Xtest <- dataTest$X
  HMalpha.Ma <- matrix(HMalpha, nrow = ntest, ncol = nygrid, byrow = TRUE)
  HMYpred <- HMalpha.Ma + Xtest %*% t(HMbeta) * (s[2] - s[1])
  HMMSPE <- mean((HMYpred - Ytest)^2)
  
  ols <- pffr(Y ~ ff(X, xind = s), yind = t, data = data, 
              bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
              bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
  olsfit <- plot(ols, select = 0)
  olsalpha <- olsfit[[1]]$fit + ols$coefficients[1]
  olsalpgrid <- olsfit[[1]]$x
  olsbetafit <- olsfit[[2]]$fit
  olsgrids <- olsfit[[2]]$x
  olsgridt <- olsfit[[2]]$y
  olsbeta <- t(matrix(olsbetafit, length(olsgrids), length(olsgridt)))
  olsalpha.t <- approx(olsalpgrid, olsalpha, t)$y 
  olsbeta.ts <- outer(olsgridt, olsgrids, beta)
  olsISE_alpha <- mean((olsalpha.t - alpha.t)^2)
  olsISE_beta <- mean((olsbeta.ts - olsbeta)^2)
  
  olsYpred <- predict(ols, newdata = list(X = Xtest))
  olsMSPE <- mean((olsYpred - Ytest)^2)

  return(list(Ralpha = HMalpha, Rbeta = HMbeta, RISE_alpha = HMISE_alpha, 
              RISE_beta = HMISE_beta, RMSPE = HMMSPE, olsalpha = olsalpha.t, 
              olsbeta = olsbeta, olsxgrid = olsgrids, olsygrid = olsgridt, 
              olsISE_alpha = olsISE_alpha, olsISE_beta = olsISE_beta, olsMSPE = olsMSPE))
}

simulationBump <- function(n, ContamProb){
  nxgrid <- 50
  nygrid <- 60
  s <- (1 : nxgrid) / nxgrid
  t <- (1 : nygrid) / nygrid
  alpha.t <- alpha(t)
  beta.ts <- outer(t, s, beta)
  
  sd.error = 0.01
  data <- BumpOutlierSample(n, NumberOfphi = 10, ygrid = t, xgrid = s, 
                            sd.error, SamContamProb = ContamProb, sl = 5, su = 8, l = 0.4)
  Y <- data$Y
  X <- data$X
  Ytrue <- data$Ytrue
  
  nubasis = 15
  ntbasis = 10
  nsbasis = 10
  lam <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 10)
  GCV <- GCV(Y, X, t, s, nubasis, ntbasis, nsbasis, lam)
  lam0 <- GCV$alphalam
  lam1 <- GCV$betalam
  
  HMfit <- MEst(Y, X, t, s, nubasis = 15, ntbasis = 10, nsbasis = 10, lam0, lam1)
  HMalpha <- HMfit$alpha
  HMbeta <- HMfit$beta
  HMISE_alpha <- mean((alpha.t - HMalpha)^2)
  HMISE_beta <- mean((beta.ts - HMbeta)^2)
  
  ntest = 100
  dataTest <- BumpOutlierSample(ntest, NumberOfphi = 10, ygrid = t, xgrid = s, 
                                sd.error, SamContamProb = ContamProb, sl = 5, su = 8, l = 0.4)
  Ytest <- dataTest$Ytrue
  Xtest <- dataTest$X
  HMalpha.Ma <- matrix(HMalpha, nrow = ntest, ncol = nygrid, byrow = TRUE)
  HMYpred <- HMalpha.Ma + Xtest %*% t(HMbeta) * (s[2] - s[1])
  HMMSPE <- mean((HMYpred - Ytest)^2)
  
  ols <- pffr(Y ~ ff(X, xind = s), yind = t, data = data, 
              bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
              bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
  olsfit <- plot(ols, select = 0)
  olsalpha <- olsfit[[1]]$fit + ols$coefficients[1]
  olsalpgrid <- olsfit[[1]]$x
  olsbetafit <- olsfit[[2]]$fit
  olsgrids <- olsfit[[2]]$x
  olsgridt <- olsfit[[2]]$y
  olsbeta <- t(matrix(olsbetafit, length(olsgrids), length(olsgridt)))
  olsalpha.t <- approx(olsalpgrid, olsalpha, t)$y 
  olsbeta.ts <- outer(olsgridt, olsgrids, beta)
  olsISE_alpha <- mean((olsalpha.t - alpha.t)^2)
  olsISE_beta <- mean((olsbeta.ts - olsbeta)^2)
  
  olsYpred <- predict(ols, newdata = list(X = Xtest))
  olsMSPE <- mean((olsYpred - Ytest)^2)
  
  return(list(Ralpha = HMalpha, Rbeta = HMbeta, RISE_alpha = HMISE_alpha, 
              RISE_beta = HMISE_beta, RMSPE = HMMSPE, olsalpha = olsalpha.t, 
              olsbeta = olsbeta, olsxgrid = olsgrids, olsygrid = olsgridt, 
              olsISE_alpha = olsISE_alpha, olsISE_beta = olsISE_beta, olsMSPE = olsMSPE))
}
  
simulationPart <- function(n, ContamProb){
  nxgrid <- 50
  nygrid <- 60
  s <- (1 : nxgrid) / nxgrid
  t <- (1 : nygrid) / nygrid
  alpha.t <- alpha(t)
  beta.ts <- outer(t, s, beta)
  
  sd.error = 0.01
  data <- PartOutlierSample(n, NumberOfphi = 10, ygrid = t, xgrid = s, 
                            sd.error, SamContamProb = ContamProb, sl = 5, su = 8)
  Y <- data$Y
  X <- data$X
  Ytrue <- data$Ytrue
  
  nubasis = 15
  ntbasis = 10
  nsbasis = 10
  lam <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 10)
  GCV <- GCV(Y, X, t, s, nubasis, ntbasis, nsbasis, lam)
  lam0 <- GCV$alphalam
  lam1 <- GCV$betalam
  
  HMfit <- MEst(Y, X, t, s, nubasis = 15, ntbasis = 10, nsbasis = 10, lam0, lam1)
  HMalpha <- HMfit$alpha
  HMbeta <- HMfit$beta
  HMISE_alpha <- mean((alpha.t - HMalpha)^2)
  HMISE_beta <- mean((beta.ts - HMbeta)^2)
  
  ntest = 100
  dataTest <- PartOutlierSample(ntest, NumberOfphi = 10, ygrid = t, xgrid = s, 
                                sd.error, SamContamProb = ContamProb, sl = 5, su = 8)
  Ytest <- dataTest$Ytrue
  Xtest <- dataTest$X
  HMalpha.Ma <- matrix(HMalpha, nrow = ntest, ncol = nygrid, byrow = TRUE)
  HMYpred <- HMalpha.Ma + Xtest %*% t(HMbeta) * (s[2] - s[1])
  HMMSPE <- mean((HMYpred - Ytest)^2)
  
  ols <- pffr(Y ~ ff(X, xind = s), yind = t, data = data, 
              bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
              bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
  olsfit <- plot(ols, select = 0)
  olsalpha <- olsfit[[1]]$fit + ols$coefficients[1]
  olsalpgrid <- olsfit[[1]]$x
  olsbetafit <- olsfit[[2]]$fit
  olsgrids <- olsfit[[2]]$x
  olsgridt <- olsfit[[2]]$y
  olsbeta <- t(matrix(olsbetafit, length(olsgrids), length(olsgridt)))
  olsalpha.t <- approx(olsalpgrid, olsalpha, t)$y 
  olsbeta.ts <- outer(olsgridt, olsgrids, beta)
  olsISE_alpha <- mean((olsalpha.t - alpha.t)^2)
  olsISE_beta <- mean((olsbeta.ts - olsbeta)^2)
  
  olsYpred <- predict(ols, newdata = list(X = Xtest))
  olsMSPE <- mean((olsYpred - Ytest)^2)
  
  return(list(Ralpha = HMalpha, Rbeta = HMbeta, RISE_alpha = HMISE_alpha, 
              RISE_beta = HMISE_beta, RMSPE = HMMSPE, olsalpha = olsalpha.t, 
              olsbeta = olsbeta, olsxgrid = olsgrids, olsygrid = olsgridt, 
              olsISE_alpha = olsISE_alpha, olsISE_beta = olsISE_beta, olsMSPE = olsMSPE))
}
  

simulationShiOut <- function(n, ContamProb){
  nxgrid <- 50
  nygrid <- 60
  s <- (1 : nxgrid) / nxgrid
  t <- (1 : nygrid) / nygrid
  alpha.t <- alpha(t)
  beta.ts <- outer(t, s, beta)
  
  sd.error = 0.01
  data <- ShiOutlierSample(n, NumberOfphi = 10, ygrid = t, xgrid = s, sd.error,
                           SamContamProb = ContamProb, sl = 2, su = 4)
  Y <- data$Y
  X <- data$X
  Ytrue <- data$Ytrue
  
  nubasis = 15
  ntbasis = 10
  nsbasis = 10
  lam <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 10)
  GCV <- GCV(Y, X, t, s, nubasis, ntbasis, nsbasis, lam)
  lam0 <- GCV$alphalam
  lam1 <- GCV$betalam
  
  HMfit <- MEst(Y, X, t, s, nubasis = 15, ntbasis = 10, nsbasis = 10, lam0, lam1)
  HMalpha <- HMfit$alpha
  HMbeta <- HMfit$beta
  HMISE_alpha <- mean((alpha.t - HMalpha)^2)
  HMISE_beta <- mean((beta.ts - HMbeta)^2)
  
  ntest = 100
  dataTest <- ShiOutlierSample(ntest, NumberOfphi = 10, ygrid = t, xgrid = s, sd.error,
                               SamContamProb = ContamProb, sl = 2, su = 4)
  Ytest <- dataTest$Ytrue
  Xtest <- dataTest$X
  HMalpha.Ma <- matrix(HMalpha, nrow = ntest, ncol = nygrid, byrow = TRUE)
  HMYpred <- HMalpha.Ma + Xtest %*% t(HMbeta) * (s[2] - s[1])
  HMMSPE <- mean((HMYpred - Ytest)^2)
  
  ols <- pffr(Y ~ ff(X, xind = s), yind = t, data = data, 
              bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
              bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
  olsfit <- plot(ols, select = 0)
  olsalpha <- olsfit[[1]]$fit + ols$coefficients[1]
  olsalpgrid <- olsfit[[1]]$x
  olsbetafit <- olsfit[[2]]$fit
  olsgrids <- olsfit[[2]]$x
  olsgridt <- olsfit[[2]]$y
  olsbeta <- t(matrix(olsbetafit, length(olsgrids), length(olsgridt)))
  olsalpha.t <- approx(olsalpgrid, olsalpha, t)$y 
  olsbeta.ts <- outer(olsgridt, olsgrids, beta)
  olsISE_alpha <- mean((olsalpha.t - alpha.t)^2)
  olsISE_beta <- mean((olsbeta.ts - olsbeta)^2)
  
  olsYpred <- predict(ols, newdata = list(X = Xtest))
  olsMSPE <- mean((olsYpred - Ytest)^2)
  
  return(list(Ralpha = HMalpha, Rbeta = HMbeta, RISE_alpha = HMISE_alpha, 
              RISE_beta = HMISE_beta, RMSPE = HMMSPE, olsalpha = olsalpha.t, 
              olsbeta = olsbeta, olsxgrid = olsgrids, olsygrid = olsgridt, 
              olsISE_alpha = olsISE_alpha, olsISE_beta = olsISE_beta, olsMSPE = olsMSPE))
}

simulationCau <- function(n){
  nxgrid <- 50
  nygrid <- 60
  s <- (1 : nxgrid) / nxgrid
  t <- (1 : nygrid) / nygrid
  alpha.t <- alpha(t)
  beta.ts <- outer(t, s, beta)
  
  sd.error <- 0.01
  cauScale <- 0.5
  data <- NormCauchySample(n, NumberOfphi = 10, ygrid = t, xgrid = s, 
                           sd.error, cauScale, MitureProb = 0.2)
  Y <- data$Y
  X <- data$X
  Ytrue <- data$Ytrue
  
  nubasis = 15
  ntbasis = 10
  nsbasis = 10
  lam <- c(10^{-3}, 10^{-2}, 10^{-1}, 1, 10)
  GCV <- GCV(Y, X, t, s, nubasis, ntbasis, nsbasis, lam)
  lam0 <- GCV$alphalam
  lam1 <- GCV$betalam
  
  HMfit <- MEst(Y, X, t, s, nubasis = 15, ntbasis = 10, nsbasis = 10, lam0, lam1)
  HMalpha <- HMfit$alpha
  HMbeta <- HMfit$beta
  HMISE_alpha <- mean((alpha.t - HMalpha)^2)
  HMISE_beta <- mean((beta.ts - HMbeta)^2)
  
  ntest = 100
  dataTest <- NormCauchySample(ntest, NumberOfphi = 10, ygrid = t, xgrid = s, 
                               sd.error, cauScale, MitureProb = 0.2)
  Ytest <- dataTest$Ytrue
  Xtest <- dataTest$X
  HMalpha.Ma <- matrix(HMalpha, nrow = ntest, ncol = nygrid, byrow = TRUE)
  HMYpred <- HMalpha.Ma + Xtest %*% t(HMbeta) * (s[2] - s[1])
  HMMSPE <- mean((HMYpred - Ytest)^2)
  
  ols <- pffr(Y ~ ff(X, xind = s), yind = t, data = data, 
              bs.int = list(bs = "ps", k = 15, m = c(2, 2)),
              bs.yindex = list(bs = "ps", k = 10, m = c(2, 2)))
  olsfit <- plot(ols, select = 0)
  olsalpha <- olsfit[[1]]$fit + ols$coefficients[1]
  olsalpgrid <- olsfit[[1]]$x
  olsbetafit <- olsfit[[2]]$fit
  olsgrids <- olsfit[[2]]$x
  olsgridt <- olsfit[[2]]$y
  olsbeta <- t(matrix(olsbetafit, length(olsgrids), length(olsgridt)))
  olsalpha.t <- approx(olsalpgrid, olsalpha, t)$y 
  olsbeta.ts <- outer(olsgridt, olsgrids, beta)
  olsISE_alpha <- mean((olsalpha.t - alpha.t)^2)
  olsISE_beta <- mean((olsbeta.ts - olsbeta)^2)
  
  olsYpred <- predict(ols, newdata = list(X = Xtest))
  olsMSPE <- mean((olsYpred - Ytest)^2)
  
  return(list(Ralpha = HMalpha, Rbeta = HMbeta, RISE_alpha = HMISE_alpha, 
              RISE_beta = HMISE_beta, RMSPE = HMMSPE, olsalpha = olsalpha.t, 
              olsbeta = olsbeta, olsxgrid = olsgrids, olsygrid = olsgridt, 
              olsISE_alpha = olsISE_alpha, olsISE_beta = olsISE_beta, olsMSPE = olsMSPE))
}

Norm50 <- replicate(100, simulationNorm(n = 50), simplify = "array")
IntO50 <- replicate(100, simulationBump(n = 50, ContamProb = 0.1), simplify = "array")
StepO50 <- replicate(100, simulationPart(n = 50, ContamProb = 0.1), simplify = "array")
Shiout50 <- replicate(100, simulationShiOut(n = 50, ContamProb = 0.1), simplify = "array")
Cau50 <- replicate(100, simulationCau(n = 50), simplify = "array")


############## plots  of beta ###############

IntO <- simulationBump(n = 50, ContamProb = 0.1)
StepO <- simulationPart(n = 50, ContamProb = 0.1)
Shiout <- simulationShiOut(n = 50, ContamProb = 0.1)
Cau <- simulationCau(n = 50)

HMbeta1 <- IntO$Rbeta
HMbeta2 <- StepO$Rbeta
HMbeta3 <- Shiout$Rbeta
HMbeta4 <- Cau$Rbeta
olsbeta1 <- IntO$olsbeta
olsbeta2 <- StepO$olsbeta
olsbeta3 <- Shiout$olsbeta
olsbeta4 <- Cau$olsbeta
olsgridt <- IntO$olsygrid
olsgrids <- IntO$olsxgrid

beta <- function(t, s) {
  4 * cos(2 * pi * t) * sin(pi * s)
}
nxgrid <- 50
nygrid <- 60
s <- (1 : nxgrid) / nxgrid
t <- (1 : nygrid) / nygrid
beta.ts <- outer(t, s, beta)

library(fields)

pdf(file = "beta.pdf")
drape.plot(t, s, beta.ts,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(beta.ts), zlab = "beta(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "HMbetaBump.pdf")
drape.plot(t, s, HMbeta1,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(HMbeta1), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "HMbetaPart.pdf")
drape.plot(t, s, HMbeta2,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(HMbeta1), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "HMbetaShit.pdf")
drape.plot(t, s, HMbeta3,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(HMbeta1), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "HMbetaCau.pdf")
drape.plot(t, s, HMbeta4,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(HMbeta4), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "olsbetaBump.pdf")
drape.plot(olsgridt, olsgrids, olsbeta1,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(olsbeta1), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "olsbetaPart.pdf")
drape.plot(olsgridt, olsgrids, olsbeta2,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(olsbeta2), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "olsbetaShit.pdf")
drape.plot(olsgridt, olsgrids, olsbeta3,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(olsbeta3), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

pdf(file = "olsbetaCau.pdf")
drape.plot(olsgridt, olsgrids, olsbeta4,
           xlab = "t", ylab = "s", ticktype = "detailed", phi = 30, theta = 30, 
           zlim = range(olsbeta4), zlab = "hat(beta)(t, s)", add.legend = FALSE, 
           col = tim.colors(64))
dev.off()

