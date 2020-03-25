
rm(list=ls())
setwd("C:/Users/cx/Desktop/Submitted rPFFRcodes")


library(fda); library(refund); library(MASS); library(fda.usc)

source("estimation.R")
source("rGCV.R")

y <- read.csv("oxygen.csv", header = T)
x <- read.csv("temper.csv", header = T)
y <- as.matrix(y); x <- as.matrix(x)
n <- nrow(y); m <- ncol(y)
rownames(y) <- 1 : n; rownames(x) <- 1 : n

yy <- scale(y, scale = FALSE)
xx <- scale(x, scale = FALSE)
s = t = seq(0, 1, length = 100)
nubasis <- 0; ntbasis = nsbasis = 20 

######### outlier detect in oxygen ###########

yfdat <- fdata(y, argvals = seq(2, 200, by = 2))
xfdat <- fdata(x, argvals = seq(2, 200, by = 2))

plot(yfdat, lty = 1, ylab = "Oxygen", xlab = "Depth below the sea surface", 
     main = "Oxygen", mgp = c(2, 0.5, 0))

plot(xfdat, lty = 1, ylab = "Temperature", xlab = "Depth below the sea surface", 
     main = "Temperature", mgp = c(2, 0.5, 0))

nb <- 100; smo <- 0.05; trim <- 0.05
out.mode <- outliers.depth.trim(yfdat, dfunc=depth.mode, nb=nb, smo=smo, trim=trim)
out <- out.mode

plot(yfdat, col="gray65", lty = 1, ylab = "Oxygen", xlab = "Depth below the sea surface",
     main = "Outlier detection", mgp = c(2, 0.5, 0))
lines(yfdat[out[[1]]], col = 2, lty = 2, lwd = 2)

############ coefficient surface estimation ############

lamPar <- c(0.01, 0.1, 1)
mgcv <- GCV(yy, xx, t, s, nubasis, ntbasis, nsbasis, lamPar)
lam <- mgcv$optlam
HM <- MEst(yy, xx, t, s, nubasis, ntbasis, nsbasis, lam0 = 0, lam1 = lam)
HMbeta <- HM$beta

image.plot(200 * t, 200 * s, HMbeta, main = expression(hat(beta)(t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))


############### prediction over 100 repeats #####################
predfun <- function(yy, xx){
  n = nrow(yy); ntrain = floor(0.7 * n); ntest = n - ntrain; m_s = m_t = ncol(yy)
  ind.n = sample(c(1 : n), ntrain, replace = FALSE)
  ytrain <- yy[ind.n, ]; xtrain <- xx[ind.n, ]
  ytest <- yy[-ind.n, ]; xtest <- xx[-ind.n, ]
  
  lamPar <- c(0.001, 0.01, 0.1, 1)
  mgcv <- GCV(ytrain, xtrain, t, s, nubasis, ntbasis, nsbasis, lamPar)
  lam <- mgcv$optlam
  HM <- MEst(ytrain, xtrain, t, s, nubasis, ntbasis, nsbasis, lam0 = 0, lam1 = lam)
  HMbeta <- HM$beta
  HMYpred <- (xtest %*% t(HMbeta) * (s[2]-s[1]))
  HMMSPE <- mean((HMYpred - ytest)^2)
  
  ols <- pffr(ytrain ~ ff(xtrain, xind = s) - 1, yind = t,
              bs.yindex = list(bs = "ps", k = 20, m = c(2, 2)))
  olsYpred <- predict(ols, newdata = list(xtrain = xtest))
  olsMSPE <- mean((olsYpred - ytest)^2)
  pred <- c(HMMSPE, olsMSPE)
  names(pred) <- c("HMMSPE", "olsMSPE")
  return(pred)
}

pred <- replicate(100, predfun(yy, xx))
PredMean <- apply(pred, 1, mean)
Predsd <- apply(pred, 1, sd)
