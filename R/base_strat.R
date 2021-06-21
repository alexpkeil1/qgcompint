vc_comb <- function (aname = "(Intercept)", expnms, covmat, grad = NULL) {
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  if (!is.null(grad))
    grad = NULL
  if (is.null(grad))
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  outcov = matrix(NA, nrow = 2, ncol = 2)
  acol = which(colnames(as.matrix(covmat)) %in% aname)
  bcol = which(colnames(as.matrix(covmat)) %in% expnms)
  outcov[1, 1] <- covmat[acol, acol]
  outcov[1, 2] <- outcov[2, 1] <- sum(covmat[acol, bcol])
  outcov[2, 2] <- weightvec %*% covmat %*% weightvec
  outcov
}

se_comb <- function (expnms, covmat, grad = NULL) {
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  if (is.null(grad))
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  if (!is.null(grad))
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  var <- weightvec %*% covmat %*% weightvec
  sqrt(var)[1, , drop = TRUE]
}

.intmaker <- function(
  f,
  expnms,
  emmvar
){
  rightside = as.character(f)[3]
  trms = strsplit(gsub(" ", "", rightside), "+",fixed=TRUE)[[1]]
  expidx <- which(trms %in% expnms)
  newtrms = paste(paste(trms[expidx], "*", emmvar), collapse = "+")
  newrightside = paste(rightside, "+", newtrms)
  newf <- as.formula(paste0(as.character(f)[2],as.character(f)[1], newrightside))
  newf
}

.intchecknames <- function(terms,emmvar){
  #nonlin <- ifelse(sum(grep("\\(|\\:|\\^", terms)) > 0, TRUE,
  #                 FALSE)
  nonlin <- any(attr(terms,"order")>1)
  if (nonlin) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

qgcomp.emm.noboot <- function(f,
                          data,
                          expnms=NULL,
                          emmvar=NULL,
                          q=4,
                          breaks=NULL,
                          id=NULL,
                          weights,
                          alpha=0.05,
                          bayes=FALSE,
                          ...){
  #' @title Quantile g-computation for continuous, binary, and count outcomes under linearity/additivity
  #'
  #' @description This function estimates a linear dose-response parameter representing a one quantile
  #' increase in a set of exposures of interest. This function is limited to linear and additive
  #' effects of individual components of the exposure. This model estimates the parameters of a marginal
  #' structural model (MSM) based on g-computation with quantized exposures. Note: this function is
  #' valid only under linear and additive effects of individual components of the exposure, but when
  #' these hold the model can be fit with very little computational burden.
  #'
  #' @details For continuous outcomes, under a linear model with no
  #' interaction terms, this is equivalent to g-computation of the effect of
  #' increasing every exposure by 1 quantile. For binary/count outcomes
  #' outcomes, this yields a conditional log odds/rate ratio(s) representing the
  #' change in the expected conditional odds/rate (conditional on covariates)
  #' from increasing every exposure by 1 quantile. In general, the latter
  #' quantity is not equivalent to g-computation estimates. Hypothesis test
  #' statistics and confidence intervals are based on using the delta
  #' estimate variance of a linear combination of random variables.
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param expnms character vector of exposures of interest
  #' @param emmvar character vector of effect measure modifier
  #' @param q NULL or number of quantiles used to create quantile indicator variables
  #' representing the exposure variables. If NULL, then gcomp proceeds with un-transformed
  #' version of exposures in the input datasets (useful if data are already transformed,
  #' or for performing standard g-computation)
  #' @param breaks (optional) NULL, or a list of (equal length) numeric vectors that
  #' characterize the minimum value of each category for which to
  #' break up the variables named in expnms. This is an alternative to using 'q'
  #' to define cutpoints.
  #' @param id (optional) NULL, or variable name indexing individual units of
  #' observation (only needed if analyzing data with multiple observations per
  #' id/cluster). Note that qgcomp.noboot will not produce cluster-appropriate
  #' standard errors (this parameter is essentially ignored in qgcomp.noboot).
  #' Qgcomp.boot can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights "case weights" - passed to the "weight" argument of
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param alpha alpha level for confidence limit calculation
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated
  #' exposures. Note: this does not lead to fully Bayesian inference in general,
  #' so results should be interpreted as frequentist.
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the
  #'  weights/standardized coefficients in the positive (pos.weights) and
  #'  negative (neg.weights) directions.
  #' @concept variance mixtures
  #' @import stats arm
  #' @export
  #' @examples
  #' set.seed(50)
  #' # linear model
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
  #' qfit <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian())
  #' # logistic model
  #' dat2 <- data.frame(y=rbinom(50, 1,0.5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat2, q=2, family=binomial())
  #' # poisson model
  #' dat3 <- data.frame(y=rpois(50, .5), x1=runif(50), x2=runif(50), z=runif(50))
  #' qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat3, q=2, family=poisson())
  #' # weighted model
  #' N=5000
  #' dat4 <- data.frame(y=runif(N), x1=runif(N), x2=runif(N), z=runif(N))
  #' dat4$w=runif(N)*2
  #' qdata = quantize(dat4, expnms = c("x1", "x2"))$data
  #' (qgcfit <- qgcomp.noboot(f=y ~ z + x1 + x2, expnms = c('x1', 'x2'), data=dat4, q=4,
  #'                          family=gaussian(), weights=w))
  #' qgcfit$fit
  #' glm(y ~ z + x1 + x2, data = qdata, weights=w)
  require("qgcomp")
  requireNamespace("qgcomp")
  if (is.null(expnms)) {
    stop("'expnms' must be specified explicitly\n")
  }
  if (is.null(emmvar)) {
    stop("'emmvar' must be specified explicitly\n")
  }
  emmlev <- length(unique(data[,emmvar]))
  newform_oldform <- terms(f, data = data)
  f = .intmaker(f,expnms,emmvar)
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(newform_oldform, "term.labels"))
  addedmain <- setdiff(addedterms, grep(":",addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  if (length(addedmain)>0) {
    message(paste0("Adding main term for ",emmvar," to the model\n"))
  }
  oord <- order(expnms)
  # order interaction terms in same order as main terms
  s0 <- gsub(paste0("^", emmvar,":"), "",
        gsub(paste0(":", emmvar,"$"), "", addedints))
  intord = order(s0)
  equalord = all.equal(oord, intord)
  if( equalord ) addedintsord = addedints
  if( !equalord ){
    neword = match(s0, expnms)
    addedintsord = addedints[neword]
  }
  nobs = nrow(data)
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  #m <- match(c("f", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE
  #
  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)
  #
  lin = .intchecknames(expnms)
  if(!lin) stop("Model appears to be non-linear: use qgcomp.boot instead")
  if (!is.null(q) | !is.null(breaks)){
    ql <- quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  } else{
    qdata <- data
    br <- breaks
  }
  if(is.null(id)) {
    # not yet implemented
    id = "id__"
    qdata$id__ = seq_len(dim(qdata)[1])
  }

  if(!bayes) fit <- glm(newform, data = qdata,
                        weights=weights,
                        ...)
  if(bayes){
    requireNamespace("arm")
    fit <- bayesglm(newform, data = qdata,
                    weights=weights,
                    ...)
  }
  mod <- summary(fit)
  if(length(setdiff(expnms, rownames(mod$coefficients)))>0){
    stop("Model aliasing occurred, likely due to perfectly correlated quantized exposures.
           Try one of the following:
             1) set 'bayes' to TRUE in the qgcomp function (recommended)
             2) set 'q' to a higher value in the qgcomp function (recommended)
             3) check correlation matrix of exposures, and drop all but one variable in each highly correlated set  (not recommended)
           ")
  }
  # terms for reference level of emmvar
  estb <- c(fit$coefficients[1], sum(mod$coefficients[expnms,1, drop=TRUE]))
  # terms for index level of emmvar
  estb1 <- c(fit$coefficients[emmvar], sum(mod$coefficients[addedints,1, drop=TRUE]))
  seb <- c(sqrt(mod$cov.scaled[1,1]), se_comb(expnms, covmat = mod$cov.scaled))
  seb1 <- c(sqrt(mod$cov.scaled[emmvar,emmvar]),
           se_comb(addedints, covmat = mod$cov.scaled))
  tstat <- estb / seb
  tstat1 <- estb1 / seb1
  df <- mod$df.null - length(expnms) - length(addedints) - 1
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pval1 <- 2 - 2 * pt(abs(tstat1), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  pvalz1 <- 2 - 2 * pnorm(abs(tstat1))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  ci1 <- cbind(estb1 + seb1 * qnorm(alpha / 2), estb1 + seb1 * qnorm(1 - alpha / 2))
  #
  effect1 <- sum(
    fit$coefficients[expnms] +
    fit$coefficients[addedintsord] +
      0 #fit$coefficients[emmvar]
  )
  seeff1 <-  se_comb(c(expnms,emmvar,addedints), covmat = mod$cov.scaled)
  cieff1 <- cbind(effect1 + seeff1 * qnorm(alpha / 2), effect1 + seeff1 * qnorm(1 - alpha / 2))
  # 'weights' at referent level
  wcoef <- fit$coefficients[expnms]
  # weights calculated based on effect size at emmvar=1, rather than coefficients
  wcoef1 <-
    fit$coefficients[expnms] +
    fit$coefficients[addedintsord] +
    0 #fit$coefficients[emmvar]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  names(wcoef1) <- gsub("_q", "", names(wcoef1))
  pos.coef <- which(wcoef > 0)
  neg.coef <- which(wcoef <= 0)
  pos.coef1 <- which(wcoef1 > 0)
  neg.coef1 <- which(wcoef1 <= 0)
  pos.weights <- abs(wcoef[pos.coef]) / sum(abs(wcoef[pos.coef]))
  neg.weights <- abs(wcoef[neg.coef]) / sum(abs(wcoef[neg.coef]))
  pos.weights1 <- abs(wcoef1[pos.coef1]) / sum(abs(wcoef1[pos.coef1]))
  neg.weights1 <- abs(wcoef1[neg.coef1]) / sum(abs(wcoef1[neg.coef1]))
  # 'post-hoc' positive and negative estimators
  # similar to constrained gWQS
  pos.psi <- sum(wcoef[pos.coef])
  neg.psi <- sum(wcoef[neg.coef])
  pos.psi1 <- sum(wcoef1[pos.coef1])
  neg.psi1 <- sum(wcoef1[neg.coef1])
  #nmpos <- names(pos.weights)
  #nmneg <- names(neg.weights)
  #se.pos.psi <- se_comb(nmpos, covmat = mod$cov.scaled)
  #se.neg.psi <- se_comb(nmneg, covmat = mod$cov.scaled)
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  names(estb) <- c('(Intercept)', "psi1")
  names(estb1) <- c(emmvar, paste0(emmvar, ":mixture"))
  res <- list(
    qx = qx,
    fit = fit,
    psi = estb[-1],
    psiint = estb1[-1],
    effect1 = effect1,
    var.effect1 = seeff1 ^ 2,
    var.psi = seb[-1] ^ 2,
    var.psiint = seb1[-1] ^ 2,
    covmat.psi=c('psi1' = seb[-1]^2),
    covmat.psiint=c('psiint' = seb1[-1]^2),
    ci = ci[-1,],
    ciint = ci1[-1,],
    cieff = cieff1,
    coef = c(estb, estb1),
    var.coef = c(seb ^ 2, seb1 ^ 2),
    #covmat.coef=c('(Intercept)' = seb[1]^2, 'psi1' = seb[2]^2),
    covmat.coef=vc_comb(aname="(Intercept)", expnms=expnms, covmat = mod$cov.scaled),
    ci.coef = rbind(ci, ci1),
    #ci.coefint = ci1,
    expnms=expnms,
    q=q,
    breaks=br,
    degree=1,
    pos.psi = pos.psi,
    neg.psi = neg.psi,
    pos.psi1 = pos.psi1,
    neg.psi1 = neg.psi1,
    pos.weights = sort(pos.weights, decreasing = TRUE),
    neg.weights = sort(neg.weights, decreasing = TRUE),
    pos.weights1 = sort(pos.weights1, decreasing = TRUE),
    neg.weights1 = sort(neg.weights1, decreasing = TRUE),
    pos.size = sum(abs(wcoef[pos.coef])),
    neg.size = sum(abs(wcoef[neg.coef])),
    pos.size1 = sum(abs(wcoef1[pos.coef1])),
    neg.size1 = sum(abs(wcoef1[neg.coef1])),
    bootstrap=FALSE,
    cov.yhat=NULL,
    alpha=alpha,
    call=origcall,
    emmlev = emmlev
  )
  if(fit$family$family=='gaussian'){
    res$tstat <- c(tstat, tstat1)
    res$df <- df
    res$pval <- c(pval,pval1)
  }
  if(fit$family$family %in% c('binomial', 'poisson')){
    res$zstat <- c(tstat, tstat1)
    res$pval <- c(pvalz, pvalz1)
  }
  attr(res, "class") <- c("qgcompemmfit", "qgcompfit")

  res
}


##### generics ######
vcov.qgcompemmfit <- function(object, ...){
  #' @importFrom stats vcov
  #' @export
  stop("not implemented")
  object$covmat.coef
}


confint.qgcompemmfit <- function(object, ...){
  #' @importFrom stats anova
  #' @export
  message("not yet implemented")
  anova(object$fit)
}


print.qgcompemmfit <- function(x, showweights=TRUE, ...){
  #' @title Default printing method for a qgcompemmfit object
  #'
  #' @description Gives variable output depending on whether `qgcomp.noboot` or `qgcomp.boot`
  #' is called. For `qgcomp.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients). For `qgcomp.boot`,
  #' the marginal effect is given, but no weights are reported since this approach
  #' generally incorporates non-linear models with interaction terms among exposures,
  #' which preclude weights with any useful interpretation.
  #'
  #' @param x "qgcompemmfit" object from `qgcomp`, `qgcomp.noboot` or `qgcomp.boot`
  #' function
  #' @param showweights logical: should weights be printed, if estimated?
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.emm.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
  #' @concept variance mixtures
  #' @export
  emmvar <- x$call$emmvar
  rnm = c("(Intercept)", 'psi1', emmvar, paste0(emmvar,"*mix"))
  fam <- x$fit$family$family
  if(!is.null(x$pos.size) & showweights) {
    cat(paste0("Scaled effect size (positive direction, sum of positive coefficients = ", signif(x$pos.size, 3) , ")\n"))
    if (length(x$pos.weights) > 0) {
      print(x$pos.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$neg.size) & showweights) {
    cat(paste0("Scaled effect size (negative direction, sum of negative coefficients = ", signif(-x$neg.size, 3) , ")\n"))
    if (length(x$neg.weights) > 0) {
      print(x$neg.weights, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$pos.size1) & showweights) {
    cat(paste0("Scaled effect size (positive direction, sum of positive effects = ", signif(x$pos.size1, 3) , ")\n"))
    if (length(x$pos.weights1) > 0) {
      print(x$pos.weights1, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if(!is.null(x$neg.size1) & showweights) {
    cat(paste0("Scaled effect size (negative direction, sum of negative effects = ", signif(-x$neg.size1, 3) , ")\n"))
    if (length(x$neg.weights1) > 0) {
      print(x$neg.weights1, digits = 3)
    } else cat("None\n")
    cat("\n")
  }
  if (fam == "binomial"){
    estimand <- 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (fam == "poisson"){
    #message("Poisson family still experimental: use with caution")
    estimand <- 'RR'
    cat(paste0("Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (fam == "gaussian"){
    cat(paste0("Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "t"
    x$zstat = x$tstat
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(x$ci.coef))){
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "pval"=x$pval)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[,1], "Upper CI"=x$ci.coef[,2], "test"=x$zstat, "pval"=x$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  .printeffects(x, digits=5)
  invisible(x)
}

.printeffects <- function(x, digits=2){
  emmv = x$call$emmvar
  x$alpha
  if( x$emmlev == 2){
    eff <- signif(x$effect1, digits=digits)
    ci <- signif(x$cieff, digits=digits)
    l1 <- paste0("Estimate (CI), ", emmv, "=1: \n")
    l2 <- paste0(eff, " (", ci[1], ", ", ci[2], ")")
    cat("\n");cat(l1);cat(l2)
  } else if( x$emmlev > 2 ){
    TRUE
  }
}



