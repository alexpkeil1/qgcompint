
.intmaker <- function(
  f,
  expnms,
  emmvars
){
  rightside = as.character(f)[3]
  trms = strsplit(gsub(" ", "", rightside), "+",fixed=TRUE)[[1]]
  expidx <- which(trms %in% expnms)
  newtrmsl = lapply(emmvars, function(x) paste(paste0(trms[expidx], "*", x), collapse = "+"))
  newtrms = paste0(newtrmsl, collapse="+")
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

qgcomp.emm.noboot <- function(
  f,
  data,
  expnms=NULL,
  emmvar=NULL,
  q=4,
  breaks=NULL,
  id=NULL,
  weights,
  alpha=0.05,
  bayes=FALSE,
  errcheck = TRUE,
  ...){
  #' @title EMM for Quantile g-computation for continuous, binary, and count outcomes under linearity/additivity
  #'
  #' @description This function fits a quantile g-computation model, allowing
  #' effect measure modification by a binary or continuous covariate. This allows
  #' testing of statistical interaction as well as estimation of stratum specific effects.
  #'
  #'
  #' @param f R style formula
  #' @param data data frame
  #' @param expnms character vector of exposures of interest
  #' @param emmvar (character) name of effect measure modifier in dataset (if categorical, must be coded as a factor variable)
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
  #' @param errcheck (logical, default=TRUE) include some basic error checking. Slightly
  #' faster if set to false (but be sure you understand risks)
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.noboot}}
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
  #' (qfit <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z", expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  #' # logistic model
  #' dat2 <- data.frame(y=rbinom(50, 1,0.5), x1=runif(50), x2=runif(50), z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
  #' (qfit2 <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",, expnms = c('x1', 'x2'), data=dat2, q=2, family=binomial()))
  #' # categorical modifier
  #' dat3 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50), z=as.factor(sample(0:2, 50,replace=TRUE)), r=rbinom(50,1,0.5))
  #' (qfit3 <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",, expnms = c('x1', 'x2'), data=dat3, q=2, family=gaussian()))
  require("qgcomp")
  requireNamespace("qgcomp")
  if(errcheck){
    # basic argument checks
    if (is.null(expnms)) {
      stop("'expnms' must be specified explicitly\n")
    }
    if (is.null(emmvar)) {
      stop("'emmvar' must be specified explicitly\n")
    }
    #if (is.factor(data[,emmvar])) {
      #stop("'emmvar' must be numeric\n")
    #}
  }
  # housekeeping
  allemmvals<- unique(data[,emmvar])
  emmlev <- length(allemmvals)
  ## process to expand factors if needed
  zdata = zproc(data[,emmvar])
  emmvars = names(zdata)
  data = cbind(data, zdata)
  ### end new
  if(errcheck){
    #if( emmlev == 2 && !all.equal(range(allemmvals), c(0,1))){
    #  stop(paste0("Variable ", emmvar, " should only take on 0/1 values"))
    #}
    #if(min(data[,emmvar])>0 || max(data[,emmvar])<0){
    #  message(paste0("Note: default weights reported are for exposure effects at ", emmvar, " = 0"))
    #}
  }
  # keep track of added terms by remembering old model
  newform_oldform <- terms(f, data = data)
  #f = .intmaker(f,expnms,emmvar) # create necessary interaction terms with exposure
  (f <- qgcompint:::.intmaker(f,expnms,emmvars)) # create necessary interaction terms with exposure
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(newform_oldform, "term.labels"))
  addedmain <- setdiff(addedterms, grep(":",addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, value = TRUE))

  #if (length(addedmain)>0) {
  #  message(paste0("Adding main term for ",emmvar," to the model\n"))
  #}
  #oord <- order(expnms)
  ## order interaction terms in same order as main terms
  #s0 <- gsub(paste0("^", emmvar,":"), "",
  #      gsub(paste0(":", emmvar,"$"), "", addedints))
  #intord = order(s0)
  #equalord = all.equal(oord, intord)
  addedintsord = addedints
  #if( equalord ) addedintsord = addedints
  #if( !equalord ){
  #  neword = match(s0, expnms)
  #  addedintsord = addedints[neword]
  #}
  nobs = nrow(data)
  # a convoluted way to handle arguments that correspond to variable names in the data frame
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
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
  if(!lin) stop("Model appears to be non-linear: this is not yet implemented")
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
  #
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
  # intercept, main effect term
  estb <- c(
    fit$coefficients[1],
    sum(mod$coefficients[expnms,1, drop=TRUE])
  )
  names(estb) <- c('(Intercept)', "psi1")
  seb <- c(
    sqrt(mod$cov.scaled[1,1]),
    se_comb2(expnms, covmat = mod$cov.scaled)
  )
  tstat <- estb / seb
  df <- mod$df.null - length(expnms) - length(addedints) - 1
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(
    estb + seb * qnorm(alpha / 2),
    estb + seb * qnorm(1 - alpha / 2)
  )

  # modifier main term, product term
  estb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
    fit$coefficients[emmvars[x]],
    sum(mod$coefficients[addedintsl[[x]],1, drop=TRUE])
  )))
  names(estb.prod) <- do.call(c, lapply(1:length(emmvars), function(x) c(emmvars[x], paste0(emmvars[x], ":mixture"))))
  seb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
    sqrt(mod$cov.scaled[emmvars[x],emmvars[x]]),
    qgcompint:::se_comb2(addedintsl[[x]], covmat = mod$cov.scaled)
  )))
  #estb.prod <- c(
  #  fit$coefficients[emmvar],
  #  sum(mod$coefficients[addedints,1, drop=TRUE])
  #)
  #names(estb.prod) <- c(emmvar, paste0(emmvar, ":mixture"))
  #seb.prod <- c(
  #  sqrt(mod$cov.scaled[emmvar,emmvar]),
  #  se_comb2(addedints, covmat = mod$cov.scaled)
  #  )
  tstat.prod <- estb.prod / seb.prod
  pval.prod <- 2 - 2 * pt(abs(tstat.prod), df = df)
  pvalz.prod <- 2 - 2 * pnorm(abs(tstat.prod))
  ci.prod <- cbind(
    estb.prod + seb.prod * qnorm(alpha / 2),
    estb.prod + seb.prod * qnorm(1 - alpha / 2)
  )
  # 'weights' at referent level
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  pos.coef <- which(wcoef > 0)
  neg.coef <- which(wcoef <= 0)
  pos.weights <- abs(wcoef[pos.coef]) / sum(abs(wcoef[pos.coef]))
  neg.weights <- abs(wcoef[neg.coef]) / sum(abs(wcoef[neg.coef]))
  # 'post-hoc' positive and negative estimators
  # similar to constrained gWQS
  pos.psi <- sum(wcoef[pos.coef])
  neg.psi <- sum(wcoef[neg.coef])
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- list(
    qx = qx,
    fit = fit,
    psi = estb[2],
    psiint = estb.prod[2*(1:length(emmvars))],
    var.psi = seb[2] ^ 2,
    var.psiint = seb.prod[2*(1:length(emmvars))] ^ 2,
    covmat.psi=c('psi1' = seb[2]^2),
    covmat.psiint=c('psiint' = seb.prod[2*(1:length(emmvars))]^2),
    ci = ci[2,],
    ciint = ci.prod[2*(1:length(emmvars)),],
    coef = c(estb, estb.prod),
    var.coef = c(seb ^ 2, seb.prod ^ 2),
    #covmat.coef=c('(Intercept)' = seb[1]^2, 'psi1' = seb[2]^2),
    covmat.coef=vc_comb(aname="(Intercept)", expnms=expnms, covmat = mod$cov.scaled), # todo: fix this
    ci.coef = rbind(ci, ci.prod),
    #ci.coefint = ci1,
    expnms=expnms,
    intterms = addedintsord,
    q=q,
    breaks=br,
    degree=1,
    pos.psi = pos.psi,
    neg.psi = neg.psi,
    pos.weights = sort(pos.weights, decreasing = TRUE),
    neg.weights = sort(neg.weights, decreasing = TRUE),
    pos.size = sum(abs(wcoef[pos.coef])),
    neg.size = sum(abs(wcoef[neg.coef])),
    bootstrap=FALSE,
    cov.yhat=NULL,
    alpha=alpha,
    call=origcall,
    emmlev = emmlev
  )
  # include some extra things by default for binary modifier (convenience only)
  if(emmlev==2){
    ww = getweightsemm(res, emmval = 1)
    ff = .calcstrateffects(res, emmval = 1)
    res = c(res,
            list(
              pos.weights1 = ww$pos.weights,
              neg.weights1 = ww$neg.weights,
              pos.psi1 = ww$pos.psi,
              neg.psi1 = ww$neg.psi,
              pos.size1 = ww$pos.size,
              neg.size1 = ww$neg.size
            ),
            list(
              effect = ff$eff,
              vareffect = ff$se^2,
              cieffect = ff$ci
            )
    )
  }
  if(fit$family$family=='gaussian'){
    res$tstat <- c(tstat, tstat.prod)
    res$df <- df
    res$pval <- c(pval,pval.prod)
  }
  if(fit$family$family %in% c('binomial', 'poisson')){
    res$zstat <- c(tstat, tstat.prod)
    res$pval <- c(pvalz, pvalz.prod)
  }
  attr(res, "class") <- c("qgcompemmfit", "qgcompfit")
  res
}






