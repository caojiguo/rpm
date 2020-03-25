
source("estimation.R")
source("rGCV.R")

library(fda); library(refund); library(MASS)

y <- as.matrix(na.omit(yfull))
x <- as.matrix(na.omit(xfull))
rownames(y) <- 1 : nrow(y); rownames(x) <- 1 : nrow(y)
Y <- log(y + 1)

yy <- scale(Y, scale = FALSE)
xx <- scale(x, scale = FALSE)
s = t = seq(0, 1, length = 24)
nubasis <- 0; ntbasis = nsbasis = 24 

######### outlier detect in Y ###########

library(fda.usc)
yfdat <- fdata(Y)
xfdat <- fdata(x)

nb <- 100; smo <- 0.05; trim <- 0.05
out.mode <- outliers.depth.trim(yfdat, dfunc=depth.mode, nb=nb, smo=smo, trim=trim)
out <- out.mode

plot(yfdat, col="gray65", lty = 1, ylab = "log(1+Count)", xlab = "Hour", 
     main = "Outlier detection", mgp = c(2, 0.5, 0))
lines(yfdat[out[[1]]], col = 2, lty = 2, lwd = 2)

ols <- pffr(yy ~ ff(xx, xind = s) - 1, yind = t,
            bs.yindex = list(bs = "ps", k = 24, m = c(2, 2)))
olsfit <- predict(ols)
olsPE <- apply((yy - olsfit)^2, 1, mean)

plot(density(olsPE), xlab = "Residuals", main = "Density of residuals", mgp = c(2, 0.5, 0))

############ coefficient surface estimation ############

lamPar <- c(0.5, 1, 5)
mgcv <- GCV(yy, xx, t, s, nubasis, ntbasis, nsbasis, lamPar)
lam <- mgcv$optlam
HM <- MEst(yy, xx, t, s, nubasis, ntbasis, nsbasis, lam0 = 0, lam1 = lam)
HMbeta <- HM$beta

image.plot(24 * t, 24 * s, HMbeta, main = expression(hat(beta)(t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))


############### prediction over 100 repeats #####################
predfun <- function(yy, xx){
  ntrue = 102; n = 70; ntest = 32; m_s = m_t = 24
  ind.n = sample(c(1 : ntrue), n, replace = FALSE)
  ytrain <- yy[ind.n, ]; xtrain <- xx[ind.n, ]
  ytest <- yy[-ind.n, ]; xtest <- xx[-ind.n, ]
  
  lamPar <- c(0.001, 0.01, 0.1, 0.5, 1, 5)
  mgcv <- GCV(ytrain, xtrain, t, s, nubasis, ntbasis, nsbasis, lamPar)
  lam <- mgcv$optlam
  HM <- MEst(ytrain, xtrain, t, s, nubasis, ntbasis, nsbasis, lam0 = 0, lam1 = lam)
  HMbeta <- HM$beta
  HMYpred <- (xtest %*% t(HMbeta) * (s[2]-s[1]))
  HMMSPE <- mean((HMYpred - ytest)^2)
  
  ols <- pffr(ytrain ~ ff(xtrain, xind = s) - 1, yind = t,
              bs.yindex = list(bs = "ps", k = 24, m = c(2, 2)))
  olsYpred <- predict(ols, newdata = list(xtrain = xtest))
  olsMSPE <- mean((olsYpred - ytest)^2)
  pred <- c(HMMSPE, olsMSPE)
  names(pred) <- c("HMMSPE", "olsMSPE")
  return(pred)
}

pred <- replicate(100, predfun(yy, xx))
PredMean <- apply(pred, 1, mean)
Predsd <- apply(pred, 1, sd)

