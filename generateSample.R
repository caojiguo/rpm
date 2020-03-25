
X.phi1 <- function(NumberOfphi, xgrid) {
  phi <- matrix(0, NumberOfphi, length(xgrid))
  for (j in 1 : NumberOfphi)
    phi[j, ] <- (j^{-2}) * sqrt(2) * sin(j * pi * xgrid)
  return(phi)
}

X.phi2 <- function(NumberOfphi, xgrid) {
  phi <- matrix(0, NumberOfphi, length(xgrid))
  for (j in 1 : NumberOfphi)
    phi[j, ] <- (j^{-2}) * sqrt(2) * cos(j * pi * xgrid)
  return(phi)
}

rX.s <- function(NumberOfphi, xgrid){
  xsi <- rnorm(2 * NumberOfphi, 0, 1)
  X <- xsi * rbind(X.phi1(NumberOfphi, xgrid), X.phi2(NumberOfphi, xgrid)) 
  Xs <- colSums(X)
  return(Xs)
}

alpha <- function(t) 2 * exp(-(t - 1)^2)
beta <- function(t, s) {
  4 * cos(2 * pi * t) * sin(pi * s)
}

NormSample <- function(n, NumberOfphi = 10, ygrid, xgrid, sd.error){
  data <- list()
  coeffun <- list()
  s <- xgrid
  t <- ygrid
  nxgrid <- length(xgrid)
  nygrid <- length(ygrid)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = nygrid, byrow = TRUE)
  #  alpha.t <- matrix(1 + dbeta(t, 2, 7), nrow = n, ncol = nygrid, byrow = TRUE)
  beta.ts <- outer(t, s, beta)   ## nygrid * nxgrid
  data$X <- I(t(replicate(n, rX.s(NumberOfphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])    ## 100*60
  data$Xbeta <- I(Xbeta.t)
  eps <- matrix(rnorm(n * nygrid, 0, sd.error), n, nygrid)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  data$Y <- I(alpha.t + Xbeta.t + eps)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

BumpOutlierSample <- function(n, NumberOfphi = 10, ygrid, xgrid, 
                              sd.error, SamContamProb, sl, su, l){
  data <- list()
  coeffun <- list()
  s <- xgrid
  t <- ygrid
  nxgrid <- length(xgrid)
  nygrid <- length(ygrid)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = nygrid, byrow = TRUE)
  beta.ts <- outer(t, s, beta)   ## nygrid * nxgrid
  data$X <- I(t(replicate(n, rX.s(NumberOfphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])    ## 100*60
  eps <- matrix(rnorm(n * nygrid, 0, sd.error), n, nygrid)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  Y <- alpha.t + Xbeta.t + eps
  n.out <- floor(SamContamProb * n)
  n.real <- n - n.out
  SamContaVec <- rep(0, n)
  SamRamLoca <- sample(1 : n, n.out, replace = FALSE)
  SamContaVec[SamRamLoca] <- 1
  Ranl <- floor(runif(n, 0, nygrid - l * nygrid))
  Ranu <- Ranl + l * nygrid
  Outlier1 <- runif(n * nygrid, sl, su)    
  Outlier2 <- runif(n * nygrid, -su, -sl) 
  Prob <- rbinom(n * nygrid, size = 1, prob = 0.5)
  Outliers <- Prob * Outlier1 + (1 - Prob) * Outlier2
  OutMa <- matrix(Outliers, n, nygrid)
  for (i in 1 : n){
    OutMa[i, -c(Ranl[i] : Ranu[i])] <- 0
  }
  StepOuty <- Y + SamContaVec * OutMa
  data$Y <- I(StepOuty)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

PartOutlierSample <- function(n, NumberOfphi = 10, ygrid, xgrid, 
                              sd.error, SamContamProb, sl, su){
  data <- list()
  coeffun <- list()
  s <- xgrid
  t <- ygrid
  nxgrid <- length(xgrid)
  nygrid <- length(ygrid)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = nygrid, byrow = TRUE)
  beta.ts <- outer(t, s, beta)   ## nygrid * nxgrid
  data$X <- I(t(replicate(n, rX.s(NumberOfphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])    ## 100*60
  eps <- matrix(rnorm(n * nygrid, 0, sd.error), n, nygrid)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  Y <- alpha.t + Xbeta.t + eps
  n.out <- floor(SamContamProb * n)
  n.real <- n - n.out
  SamContaVec <- rep(0, n)
  SamRamLoca <- sample(1 : n, n.out, replace = FALSE)
  SamContaVec[SamRamLoca] <- 1
  Ranl <- floor(runif(n, 0, nygrid))
  Ranu <- nygrid
  StOutlier1 <- runif(n * nygrid, sl, su)    
  StOutlier2 <- runif(n * nygrid, -su, -sl) 
  StProb <- rbinom(n * nygrid, size = 1, prob = 0.5)
  StOutliers <- StProb * StOutlier1 + (1 - StProb) * StOutlier2
  StOutMa <- matrix(StOutliers, n, nygrid)
  for (i in 1 : n){
    StOutMa[i, -c(Ranl[i] : Ranu)] <- 0
  }
  StepOuty <- Y + SamContaVec * StOutMa
  data$Y <- I(StepOuty)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

ShiOutlierSample <- function(n, NumberOfphi = 10, ygrid, xgrid, sd.error,
                             SamContamProb, sl, su){
  data <- list()
  coeffun <- list()
  s <- xgrid
  t <- ygrid
  nxgrid <- length(xgrid)
  nygrid <- length(ygrid)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = nygrid, byrow = TRUE)
  beta.ts <- outer(t, s, beta)   ## nygrid * nxgrid
  data$X <- I(t(replicate(n, rX.s(NumberOfphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])    ## 100*60
  eps <- matrix(rnorm(n * nygrid, 0, sd.error), n, nygrid)
  data$Ytrue <- I(alpha.t + Xbeta.t)
  Y <- alpha.t + Xbeta.t + eps
  n.out <- floor(SamContamProb * n)
  n.real <- n - n.out
  SamContaVec <- rep(0, n)
  SamRamLoca <- sample(1 : n, n.out, replace = FALSE)
  SamContaVec[SamRamLoca] <- 1
  outlier1 <- runif(n, sl, su)    
  outlier2 <- runif(n, -su, -sl) 
  Prob <- rbinom(n, size = 1, prob = 0.5)
  outliers <- Prob * outlier1 + (1 - Prob) * outlier2
  outy <- Y + SamContaVec * outliers
  data$Y <- I(outy)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

NormCauchySample <- function(n, NumberOfphi = 10, ygrid, xgrid, 
                             sd.error, cauScale, MitureProb){
  data <- list()
  coeffun <- list()
  s <- xgrid
  t <- ygrid
  nxgrid <- length(xgrid)
  nygrid <- length(ygrid)
  alpha.t <- matrix(alpha(t), nrow = n, ncol = nygrid, byrow = TRUE)
  beta.ts <- outer(t, s, beta)   ## nygrid * nxgrid
  data$X <- I(t(replicate(n, rX.s(NumberOfphi, s))))
  Xbeta.t <- (data$X) %*% t(beta.ts) * (s[2] - s[1])    ## 100*60
  MixVec <- rbinom(n, size = 1, MitureProb)
  Mixeps <- (1 - MixVec) * matrix(rnorm(n * nygrid, 0, sd.error), n, nygrid) + 
             MixVec * matrix(rcauchy(n * nygrid, 0, cauScale), n, nygrid)
  Cauy <- alpha.t + Xbeta.t + Mixeps
  data$Ytrue <- I(alpha.t + Xbeta.t)
  data$Y <- I(Cauy)
  coeffun$alpha <- alpha.t
  coeffun$beta <- beta.ts
  return(structure(as.data.frame(data, row.names = 1 : n), xindex = s, yindex = t,
                   truth = coeffun))
}

