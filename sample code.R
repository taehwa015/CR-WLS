##############################################################################
## Weighted least-squares regression with competing risks data
## By Sangbum Choi, Taehwa Choi, Hyunsoon Cho and Dipankar Bandyopadhyay
## Update: July 13, 2021
##############################################################################
## Required libraries
library(crrSC)
library(ranger)
library(MASS)
library(cmprsk)

## User defined functions
weight = function(d, type) {
  if (type == "km") {
    km = survfit(Surv(Y, cause == 0) ~ 1, d)
    suv = approx(x = km$time,
                 y = km$surv,
                 xout = d$Y)$y
  } else if (type == "cox") {
    cox = survfit(coxph(Surv(Y, cause == 0) ~ ., d))
    suv = approx(x = cox$time,
                 y = cox$surv,
                 xout = d$Y)$y
  } else if (type == "rf") {
    rf = ranger(Surv(Y, cause == 0) ~ ., d)
    msurv = apply(predictions(rf), 2, mean)
    suv = approx(x = ranger::timepoints(rf),
                 y = msurv,
                 xout = d$Y)$y
  }
  suv[suv < 0.0001] = 1
  I(d$cause == 1) / suv
}

lmcrr.cor = function(d, type = c("km", "cox", "rsf"), alpha = 1) {
  Qfunc = function(t, Y, cause, Z, w, beta, m) {
    c(t(Z) %*% c((Y>t) * w/m * (log(Y)-Z%*%beta)))/sum(Y>=t)
  }
  m = as.numeric(table(d$id))
  m = unlist(sapply(m, function(a) rep(a, a)))^alpha
  Y = d$Y
  cause = d$cause
  id = d$id
  Z = cbind(1, z1 = d$z1, z2 = d$z2)
  n = length(unique(id))
  p = ncol(Z)
  nn = nrow(Z)
  wdt = data.frame(Y, cause, Z)
  w = d$w = weight(wdt, type = type)
  beta = lm(log(Y) ~ z1 + z2, weights = w/m, d)$coef
  
  A = t(Z) %*% (Z*w/m) 
  Q = matrix(0, nn, p)
  for (i in 1:nn) Q[i,] = Qfunc(Y[i], Y, cause, Z, w, beta, m = m)
  
  Bmat = Z * c(w/m*(log(Y)-Z%*%beta)) - (cause==0)*Q
  for (i in 1:nn) {
    for (j in 1:nn) {
      if (cause[j] == 0)
        Bmat[i,] = Bmat[i,] - (Y[i]>=Y[j]) * Q[j,]/sum(Y>=Y[j])
    }
  }
  B = t(Bmat) %*% Bmat
  vv = ginv(A) %*% B %*% ginv(A)
  se = sqrt(diag(vv))
  data.frame(beta = beta[-1], se = se[-1])
}


resfun = function(fit) {
  est = fit$beta
  round(cbind(est, fit$se), 2)
}

## Data analysis
## Cause 1
data(center)
d = with(center,
         data.frame(
           Y = ftime,
           cause = fstatus,
           z1 = fm,
           z2 = cells,
           id = id
         ))
d = d[complete.cases(d), ]
d = d[d$id %in% names(which(table(d$id) != 1)), ]

# Cause 1, cluster adjusted AFT
fit1 = lmcrr.cor(d = d, type = "km")
fit2 = lmcrr.cor(d = d, type = "cox")
fit3 = lmcrr.cor(d = d, type = "rf")

# Cause 1, cluster unadjusted AFT
fit4 = lmcrr.cor(d = d, type = "km", alpha = 0)
fit5 = lmcrr.cor(d = d, type = "cox", alpha = 0)
fit6 = lmcrr.cor(d = d, type = "rf", alpha = 0)

# Cause 1, cluster adjusted SPH
cov.test = cbind(center$fm, center$cells)
fit7 = crrc(ftime = center[, 1], fstatus = center[, 2], 
            cov1 = cov.test, cluster = center$id)
fit7 = list(beta = fit7$coef, se = sqrt(diag(fit7$var)))

resfun(fit1)
resfun(fit2)
resfun(fit3)
resfun(fit4)
resfun(fit5)
resfun(fit6)
resfun(fit7)


## Cause 2
d = with(center,
         data.frame(
           Y = ftime,
           cause = fstatus,
           z1 = fm,
           z2 = cells,
           id = id
         ))
d$cause[d$cause == 1] = 2
d$cause[d$cause == 2] = 1
d = d[complete.cases(d), ]

# Cause 2, cluster adjusted AFT
fit1 = lmcrr.cor(d = d, type = "km")
fit2 = lmcrr.cor(d = d, type = "cox")
fit3 = lmcrr.cor(d = d, type = "rf")

# Cause 2, cluster unadjusted AFT
fit4 = lmcrr.cor(d = d, type = "km", alpha = 0)
fit5 = lmcrr.cor(d = d, type = "cox", alpha = 0)
fit6 = lmcrr.cor(d = d, type = "rf", alpha = 0)

# Cause 2, cluster adjusted SPH
cov.test = cbind(d$z1, d$z2)
fit7 = crrc(ftime = d[, 1], fstatus = d[, 2], 
            cov1 = cov.test, cluster = d$id)
fit7 = list(beta = fit7$coef, se = sqrt(diag(fit7$var)))

resfun(fit1)
resfun(fit2)
resfun(fit3)
resfun(fit4)
resfun(fit5)
resfun(fit6)
resfun(fit7)