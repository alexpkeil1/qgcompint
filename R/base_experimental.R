predict.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  if(is.null(newdata)){
    pred <- predict(object$fit, type=type, ...)
  }
  if(!is.null(newdata)){
    if(is.null(expnms[1])) expnms = object$expnms # testing
    oldZ = object$fit$data[,object$call$emmvar]
    oldl = length(oldZ)
    zdata = zproc(c(oldZ,newdata[,object$call$emmvar]))
    #emmvars = names(zdata)
    newdata = cbind(newdata, zdata[-seq_len(oldl)])
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...)
  }
  return(pred)
}

predictmsm <- function (object, ...)
  UseMethod("predictmsm")


predictmsm.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  if(is.null(newdata)){
    pred <- predict(object$msmfit, type=type, ...)
  }
  if(!is.null(newdata)){
    if(is.null(expnms[1])) expnms = object$expnms # testing
    oldZ = object$fit$data[,object$call$emmvar]
    oldl = length(oldZ)
    zdata = zproc(c(oldZ,newdata[,object$call$emmvar]))
    #emmvars = names(zdata)
    newdata = cbind(newdata, zdata[-seq_len(oldl)])
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...)
  }
  return(pred)
}



qgcomp.survcurve.boot <- function(x, ...){
  stop("Not yet implemented")
  namespaceImport("survival")
  rootdat <- as.data.frame(x$fit$x)
  psidat <- data.frame(psi=0)
  rootfun <- function(idx, df){
    df[,x$expnms] <- idx
    df
  }
  rootfun2 <- function(idx, df){
    df[,"psi"] <- idx
    df[,"psi1"] <- idx
    df[,"psi2"] <- idx^2
    df[,"psi3"] <- idx^3
    df[,"psi4"] <- idx^4
    df
  }
  newmarg = lapply(0:(x$q-1), rootfun2, df=psidat)
  margdf = data.frame(do.call("rbind", newmarg))
  newcond = lapply(0:(x$q-1), rootfun, df=rootdat)
  conddf = data.frame(do.call("rbind", newcond))
  msmobj = survival::survfit(x$msmfit, newdata=margdf)
  gcompobj = survival::survfit(x$fit, newdata=conddf)
  #
  mdfl = lapply(seq_len(x$q), function(zz) with(survival::survfit(x$msmfit, newdata=newmarg[[zz]]), data.frame(time=time, surv=surv, q=zz)))
  cdfl = lapply(seq_len(x$q), function(zz) with(survival::survfit(x$fit, newdata=newcond[[zz]][1,]), data.frame(time=time, surv=surv, q=zz)))
  mdfq = do.call(rbind, mdfl)
  cdfq = do.call(rbind, cdfl)
  mdf = with(msmobj, data.frame(time=time, surv=apply(surv, 1, mean)))
  cdf = with(gcompobj, data.frame(time=time, surv=apply(surv, 1, mean)))
  list(
    mdfpop = mdf, #
    cdfpop = cdf,
    mdfq = mdfq,
    cdfq = cdfq
  )
}

getjointeffects <- function(x, emmval=1.0, ...){
  #' @title Calculate joint effect of mixture effect and modifier vs. common referent
  #'
  #' @description A standard qgcomp fit with effect measure modification
  #' only estimates effects at the referent (0) level of the modifier (psi1).
  #' This function can be used to estimate a "common referent" parameter that
  #' estimates the effect of being in a non-referent category of the modifier and
  #' increasing exposure by one quantile, relative to no change in exposure in the
  #' referent category of the modifier. This is generally useful for binary exposures
  #' (for a mixture with a set of binary exposures,
  #' this would be the "effect" of being exposed and at the index level of the mediator,
  #' relative to being unexposed in the referent level of the mediator), but it may also
  #' be of interest with more general exposures.
  #'
  #'
  #' @param x "qgcompemmfit" object from qgcomp.emm.noboot
  #' function
  #' @param emmval numerical: value of effect measure modifier at which weights are generated
  #' @param ... unused
  #' @seealso \code{\link[qgcompint]{qgcomp.emm.noboot}} \code{\link[qgcompint]{getstrateffects}}
  #' @concept variance mixtures
  #' @return
  #' An object of class "qgcompemmeffects", which inherits from "qgcompemmfit" and "list"
  #'
  #' This class contains the `emmval`-stratum specific effect estimates of the mixture. By default, this prints a coefficient table, similar to objects of type "qgcompemmfit" which displays the stratum specific joint effects from a "qgcompemmfit" model.
  #'
  #' @export
  #' @examples
  #' library(qgcompint)
  #' n = 500
  #' dat <- data.frame(y=rbinom(n,1,0.5), cd=runif(n), pb=runif(n),
  #'                   raceth=factor(sample(c("WNH", "BNH", "AMIND"), n, replace=TRUE),
  #'                           levels = c("BNH", "WNH", "AMIND")))
  #' (qfit <- qgcomp.emm.noboot(f=y ~cd + pb, emmvar="raceth",
  #'                            expnms = c('cd', 'pb'), data=dat, q=4,
  #'                            family=binomial()))
  #'
  #'
  #' # first level of the stratifying variable should be the referent category,
  #' #  which you can set with the "levels" argument to "factor" when
  #' #  cleaning/generating data
  #' levels(dat$raceth)
  #'
  #' # stratum specific mixture log-odds ratios
  #' # this one comes straight from the model (psi 1)
  #' getjointeffects(qfit, emmval = "BNH")
  #' # this will coincide with joint effects, since it is in the referent category
  #' getstrateffects(qfit, emmval = "BNH")
  #'
  #' # the stratum specific effect for a non-referent category of the EMM
  #' #  will not coincide with the joint effect
  #' getjointeffects(qfit, emmval = "AMIND")
  #' getstrateffects(qfit, emmval = "AMIND")
  #'

  #expnms = x$expnms
  #addedintsord =  x$intterms  zvar = x$fit$data[,x$call$emmvar]
  #if(x$bootstrap) stop("This method does not work for bootstrapped fits. If using a linear parameterization, then stratified effects can be estimated using non-bootstrapped methods.")
  if(x$degree>1) stop("not implemented for non-linear fits")
  zvar = x$fit$data[,x$call$emmvar]
  res = .calcjointffects(x,emmval=emmval, zvar=zvar)
  class(res) <- "qgcompemmeffects"
  res
}



.calcjointffects <- function(x, emmval=1.0, zvar){
  #x$call$emmvar
  whichintterms = x$intterms
  if(is.factor(zvar)){
    whichlevels = zproc(zvar[which(zvar==emmval)][1], znm = x$call$emmvar)
    whichvar = names(whichlevels)[which(whichlevels==1)]
    whichmainterms = whichvar
    whichintterms = NULL
    if(length(whichvar)>0) whichintterms = grep(whichvar, x$intterms, value = TRUE)
  }

  #lnx = length(x$expnms)
  #lnxz = length(whichintterms)
  mod = summary(x$fit)
  if( x$fit$family$family=="cox" ){
    covmat = as.matrix(x$fit$var)
    colnames(covmat) <- rownames(covmat) <- names(x$fit$coefficients)
  } else{
    covmat = as.matrix(mod$cov.scaled)
  }
  #stopifnot(lnx == lnxz)
  if(is.factor(zvar)){
    indeffects =
      x$fit$coefficients[x$expnms]
    if(!is.null(whichintterms)){
      indeffects =
        indeffects +
        x$fit$coefficients[whichintterms]
    }
  } else{
    indeffects =
      x$fit$coefficients[x$expnms] +
      x$fit$coefficients[x$intterms]*emmval
  }
  if(length(whichmainterms)>1)
    stop("getjointeffects: length(whichmainterms)>1, which generally means something is wrong in code")
  maineffects = 0
  if(length(whichmainterms)==1)
    maineffects = x$fit$coefficients[whichmainterms] # this
  effectatZ <- sum(indeffects)  + maineffects
  expidx <- which(colnames(covmat) %in% x$expnms)
  mainidx <- which(colnames(covmat) %in% whichmainterms)
  intidx <- which(colnames(covmat) %in% whichintterms)
  effgrad = 0*x$fit$coefficients
  effgrad[expidx] <- 1
  effgrad[mainidx] <- 1
  if(is.factor(zvar)){
    effgrad[intidx] <- 1.0
  } else effgrad[intidx] <- emmval
  seatZ <-  se_comb2(c(x$expnms, whichmainterms, x$intterms),
                     covmat = covmat,
                     grad = effgrad
  )
  ciatZ <- cbind(
    effectatZ + seatZ * qnorm(x$alpha / 2),
    effectatZ + seatZ * qnorm(1 - x$alpha / 2)
  )
  res <- list(
    effectmat = rbind(
      terms.emm = x$fit$coefficients[whichmainterms],
      terms.main = x$fit$coefficients[c(x$expnms)],
      terms.prod = x$fit$coefficients[x$intterms],
      indeffects = indeffects
    )
    , # main effect + product term
    eff = effectatZ,
    se = seatZ,
    ci = ciatZ,
    emmvar = x$call$emmvar,
    emmlev = x$emmlev,
    emmval = emmval
  )
  res
}
