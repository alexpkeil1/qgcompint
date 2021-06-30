library(qgcomp)
library(qgcompint)
set.seed(40)
# linear model
dat <- data.frame(y=runif(50),
                  x1=runif(50),
                  x2=runif(50),
                  z=rbinom(50,1,0.5),
                  r=rbinom(50,1,0.5))
(qfit <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
(qfitemmboot <- qgcomp.emm.boot(f=y ~ z + x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
(qfitemmb <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, bayes=TRUE, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
# check
getstratweights(qfitemm, emmval=0.0)
getstrateffects(qfitemm, emmval=0.0)
getstratweights(qfitemm, emmval=1.0)
getstrateffects(qfitemm, emmval=1.0)
getstratweights(qfitemmb, emmval=1.0)
getstrateffects(qfitemmb, emmval=1.0)

dat$z=as.factor(dat$z)
(qfitemmf <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))

checkres <- function(qfitemm){
  q0 <- q1 <- qfitemm$fit$data
  #q0$z <- q1$z <- 1
  q0$x1 <- q0$x2 <-  0
  q1$x1 <- q1$x2 <-  1
  #qfitemm$fit$coefficients
  qfitemm
  mean(predict(qfitemm, newdata = q1[q1$z==0,]))-
    mean(predict(qfitemm, newdata = q0[q0$z==0,]))
  mean(predict(qfitemm, newdata = q1[q1$z==1,]))-
    mean(predict(qfitemm, newdata = q0[q0$z==1,]))
}



# categorical modifier
dat <- data.frame(y=runif(50),
                  x1=runif(50),
                  x2=runif(50),
                  z=as.factor(sample(1:3, 50,replace=TRUE)),
                  r=rbinom(50,1,0.5))
(qfit <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
(qfitemm <- qgcomp.emm.noboot(f=y ~ x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian()))
(qfitemmboot <- qgcomp.emm.boot(f=y ~ x1 + x2, degree=1, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian()))
(qfitemmb <- qgcomp.emm.noboot(f=y ~ x1 + x2, bayes=TRUE, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
# check
getstratweights(qfitemm, emmval=0.0)
getstrateffects(qfitemm, emmval=0.0)
getstratweights(qfitemm, emmval=1.0)
getstrateffects(qfitemm, emmval=1.0)
getstratweights(qfitemmb, emmval=1.0)



##### continuous emm, assessing type 1 error only
n = 200
dat <- data.frame(y=runif(n),
                  x1=runif(n),
                  x2=runif(n),
                  x3=runif(n),
                  x4=runif(n),
                  x5=runif(n),
                  z=rnorm(n,1,0.5),
                  r=rbinom(n,1,0.5))
(qfit <- qgcomp.noboot(f=y ~ ., expnms = paste0("x",1:5), data=dat, q=4, family=gaussian()))
(qfitemm <- qgcomp.emm.noboot(f=y ~x1+x2+x3+x4+x5, emmvar="z", expnms = paste0("x",1:5), data=dat, q=2, family=gaussian()))


qfitemm$fit
getstratweights(qfitemm, emmval=0.0)
getstratweights(qfitemm, emmval=1.0)
getstratweights(qfitemm, emmval=2.0)
getstrateffects(qfitemm, emmval=0.0)
getstrateffects(qfitemmb, emmval=1.0)
getstratweights(qfitemm, emmval=-2.0)
getstrateffects(qfitemmb, emmval=-2.0)


# -----------------------------------
# long run simulations (none of these are actually run in testing)
# -----------------------------------

## type 1 error true linear model, binary modifier ##
rft <- function(i,n=50){
  dat <- data.frame(y=runif(n),
                    x1=runif(n),
                    x2=runif(n),
                    z=rbinom(n,1,0.5),
                    r=rbinom(n,1,0.5))
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + r, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  (qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + r, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  c(qfit$pval[2], qfitemm$pval[c(2,4)], qfit$psi,qfitemm$psi,qfitemm$psiint, qfit$covmat.psi, qfitemm$covmat.psi, qfitemm$covmat.psiint)
}

runsims <- function(){
  resl = lapply(1:1000, rft)
  res = do.call(rbind, resl)

  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  apply(res[,4:6], 2, mean)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
}
#runsims()

## type 1 error for logistic model, binary modifier ##
rftb <- function(i,n=50){
  dat <- data.frame(y=rbinom(n,1,0.5),
                    x1=runif(n),
                    x2=runif(n),
                    z=rbinom(n,1,0.5),
                    r=rbinom(n,1,0.5))
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + r, expnms = c('x1', 'x2'), data=dat, q=2, family=binomial()))
  (qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + r, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=binomial()))
  c(qfit$pval[2], qfitemm$pval[c(2,4)], qfit$psi,qfitemm$psi,qfitemm$psiint, qfit$covmat.psi, qfitemm$covmat.psi, qfitemm$covmat.psiint)
}

runsimsb <- function(){
  resl = lapply(1:1000, rftb, n=100)
  res = do.call(rbind, resl)

  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  apply(res[,4:6], 2, mean)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
}


## type 1 error for true linear model, continuous modifier ##
rft2 <- function(i,n=50){
  dat <- data.frame(y=runif(n),
                    x1=runif(n),
                    x2=runif(n),
                    z=rnorm(n,1,0.5),
                    r=rbinom(n,1,0.5))
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + r, expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + r, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  c(qfit$pval[2], qfitemm$pval[c(2,4)], qfit$psi,qfitemm$psi,qfitemm$psiint, qfit$covmat.psi, qfitemm$covmat.psi, qfitemm$covmat.psiint)
}

#not run
runsims2 <- function(){
  resl = lapply(1:10000, rft2, n=100)
  res = do.call(rbind, resl)

  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  apply(res[,4:6], 2, mean)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
}
#runsims2()

## bias, coverage, power for true linear model, binary modifier ##
rft3 <- function(i, n=100){
  dat <- qgcompint:::.dgm_quantized_linear_emm(
    N = n,
    b0=0,
    mainterms=  c(0.3,0.0,0.3,0.0),
    prodterms = c(0.3,0.1,0.1,0.0),
    ztype = "binary",
    ncor=0,
    corr=0.75
  )
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + x3 + x4, emmvar="z", expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[2],
    pz = qfitemm$pval[c(2,4)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint,
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint,
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}

runsims3 <- function(){
  resl = lapply(1:1000, rft3, n=100)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  (bias <- est - truth)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}
#runsims3()


## bias, coverage, power for true linear model, continuous modifier ##
rft4 <- function(i, n=100){
  dat <- qgcompint:::.dgm_quantized_linear_emm(
    N = n,
    b0=0,
    mainterms=  c(0.3,0.0,0.3,0.0),
    prodterms = c(0.3,0.1,0.1,0.0),
    ztype = "cont",
    ncor=0,
    corr=0.75
  )
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + x3 + x4, emmvar="z", expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[2],
    pz = qfitemm$pval[c(2,4)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint,
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint,
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}

runsims4 <- function(){
  resl = lapply(1:1000, rft4, n=100)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  (bias <- est - truth)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}

#runsims4()



## bias, coverage, power for true logistic model, binary modifier ##
rft5 <- function(i, n=100){
  dat <- qgcompint:::.dgm_quantized_logistic_emm(
    N = n,
    b0=0,
    mainterms=  c(0.3,0.0,0.3,0.0),
    prodterms = c(0.3,0.1,0.1,0.0),
    ztype = "binary",
    ncor=0,
    corr=0.75
  )
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = paste0("x",1:4), data=dat, q=NULL, family=binomial()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + x3 + x4, emmvar="z", expnms = paste0("x",1:4), data=dat, q=NULL, family=binomial()))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[2],
    pz = qfitemm$pval[c(2,4)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint,
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint,
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}

runsims5 <- function(){
  resl = lapply(1:1000, rft5, n=200)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  apply(res[,4:6], 2, median)
  (bias <- est - truth)
  plot(density(res[,4]))
  plot(density(res[,5]))
  plot(density(res[,6]))
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}
#runsims5()


## bias, coverage, power for true logistic model, binary modifier, balancing modification ##
rft6 <- function(i, n=100){
  dat <- qgcompint:::.dgm_quantized_logistic_emm(
    N = n,
    b0=0,
    mainterms=  c(0.3,0.0,0.3,0.0),
    prodterms = c(-0.3,0.15,0.15,0.0),
    ztype = "binary",
    ncor=0,
    corr=0.75
  )
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = paste0("x",1:4), data=dat, q=NULL, family=binomial()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + x3 + x4, emmvar="z", expnms = paste0("x",1:4), data=dat, q=NULL, family=binomial()))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[2],
    pz = qfitemm$pval[c(2,4)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint,
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint,
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}

runsims6 <- function(){
  resl = lapply(1:1000, rft6, n=200)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  apply(res[,4:6], 2, median)
  (bias <- est - truth)
  #plot(density(res[,4]))
  #plot(density(res[,5]))
  #plot(density(res[,6]))
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}
#runsims6()



## bias, coverage, power for true Cox model, binary modifier, balancing modification ##
rft7 <- function(i, n=1000){
  dat = qgcompint:::.dgm_quantized_survival_emm(
    N = n,
    b0=0,
    mainterms=c(1,0,0,0),
    prodterms = c(1,0,0,0),
    ztype = "bin",
    ncor=0,
    corr=0.75
  )
  f = survival::Surv(time, d)~x1 + x2 + x3 + x4+z
  expnms = c("x1", "x2", "x3", "x4")
  #(fit1 <- survival::coxph(f, data = dat))
  (qfit <- qgcomp.cox.noboot(f, expnms = expnms, data = dat, q=4))
  suppressMessages(qfitemm <- qgcomp.emm.cox.noboot(f, expnms = expnms, emmvar="z", data = dat, q=4))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[1],
    pz = qfitemm$pval[c(1,3)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint,
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint,
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}
runsims7 <- function(){
  resl = lapply(1:1000, rft7, n=200)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  apply(res[,4:6], 2, median)
  (bias <- est - truth)
  #plot(density(res[,4]))
  #plot(density(res[,5]))
  #plot(density(res[,6]))
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}
#runsims7()




## bias, coverage, power for true linear model, categorical modifier ##
rft8 <- function(i, n=100){
  dat <- qgcompint:::.dgm_quantized_linear_emm(
    N = n,
    b0=0,
    mainterms=  c(0.3,0.0,0.3,0.0),
    prodterms = c(0.3,0.1,0.1,0.0),
    ztype = "cat",
    ncor=0,
    corr=0.75
  )
  dat$z = as.factor(dat$z)
  (qfit <- qgcomp.noboot(f=y ~ x1 + x2 + x3 + x4, expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  suppressMessages(qfitemm <- qgcomp.emm.noboot(f=y ~ z + x1 + x2 + x3 + x4, emmvar="z", expnms = paste0("x",1:4), data=dat, q=NULL, family=gaussian()))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[2],
    pz = qfitemm$pval[c(2,4)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint[1],
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint[1],
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}

runsims8 <- function(){
  resl = lapply(1:1000, rft8, n=200)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  (bias <- est - truth)
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}

#runsims8()


## bias, coverage, power for true Cox model, categorical modifier, balancing modification ##
rft9 <- function(i, n=1000){
  dat = qgcompint:::.dgm_quantized_survival_emm(
    N = n,
    b0=0,
    mainterms=c(1,0,0,0),
    prodterms = c(1,0,0,0),
    ztype = "cat",
    ncor=0,
    corr=0.75
  )
  dat$z = as.factor(dat$z)
  f = survival::Surv(time, d)~x1 + x2 + x3 + x4
  expnms = c("x1", "x2", "x3", "x4")
  #(fit1 <- survival::coxph(f, data = dat))
  (qfit <- qgcomp.cox.noboot(f, expnms = expnms, data = dat, q=4))
  suppressMessages(qfitemm <- qgcomp.emm.cox.noboot(f, expnms = expnms, emmvar="z", data = dat, q=4))
  truth = attr(dat, "truecoef")
  c(pmarg = qfit$pval[1],
    pz = qfitemm$pval[c(1,3)],
    psixmarg = qfit$psi,
    psix0 = qfitemm$psi,
    psix1 = qfitemm$psiint[1],
    psivarxmarg = qfit$covmat.psi,
    psivarx0 = qfitemm$covmat.psi,
    psivarx0 = qfitemm$covmat.psiint[1],
    truth0 = sum(truth$mainterms),
    truthint = sum(truth$prodterms)
  )
}
runsims9 <- function(){
  resl = lapply(1:1000, rft9, n=200)
  res = as.data.frame(do.call(rbind, resl))
  (truth <- c(NA, mean(res$truth0), mean(res$truthint)))
  pow = (res[,1:3]<0.05)*1.0
  apply(pow, 2, mean)
  (est <- apply(res[,4:6], 2, mean))
  apply(res[,4:6], 2, median)
  (bias <- est - truth)
  #plot(density(res[,4]))
  #plot(density(res[,5]))
  #plot(density(res[,6]))
  apply(res[,4:6], 2, var)
  apply(res[,7:9], 2, mean)
  #
  varest = res[,7:9]
  ll = res[,4:6] + sqrt(varest)*qnorm(.05/2)
  ul = res[,4:6] + sqrt(varest)*qnorm(1-.05/2)
  cover = (t(ll) < truth) & (t(ul) > truth)
  apply(cover, 1, mean)
}
#runsims9()
