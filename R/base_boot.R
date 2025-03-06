msm.emm.fit <- function(f,
                    qdata,
                    intvals,
                    expnms,
                    emmvar, # original z variable
                    emmvars, #  z variable possibly split into indicator variables
                    rr=TRUE,
                    main=TRUE,
                    degree=1,
                    id=NULL,
                    weights,
                    bayes=FALSE,
                    MCsize=nrow(qdata), hasintercept=TRUE, ...){

  newform <- terms(f, data = qdata)
  nobs = nrow(qdata)
  thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("qdata", "data", names(thecall))
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE

  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    qdata$weights <- as.vector(model.weights(thecalle))
  } else qdata$weights = rep(1, nobs)

  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  # conditional outcome regression fit
  nidx = which(!(names(qdata) %in% id))
  if(!bayes) fit <- glm(newform, data = qdata,
                        weights=weights,
                        ...)
  if(bayes){
    requireNamespace("arm")
    fit <- bayesglm(f, data = qdata[,nidx,drop=FALSE],
                    weights=weights,
                    ...)
  }
  if(fit$family$family %in% c("gaussian", "poisson")) rr=FALSE
  ###
  # get predictions (set exposure to 0,1,...,q-1)
  if(is.null(intvals)){
    intvals <- (seq_len(length(table(qdata[expnms[1]])))) - 1
  }
  predit <- function(idx, newdata){
    #newdata <- qdata
    newdata[,expnms] <- idx
    suppressWarnings(predict(fit, newdata=newdata, type='response'))
  }
  if(MCsize==nrow(qdata)){
    newdata <- qdata
  }else{
    newids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), MCsize,
                                          #probs=weights, #bootstrap sampling with weights works with fixed weights, but not time-varying weights
                                          replace = TRUE
    )))
    names(newids) <- id
    newdata <- merge(qdata,newids, by=id, all.x=FALSE, all.y=TRUE)[seq_len(MCsize),]
  }
  predmat = lapply(intvals, predit, newdata=newdata)
  # fit MSM using g-computation estimates of expected outcomes under joint
  #  intervention
  #nobs <- dim(qdata)[1]
  msmdat <- data.frame(
    cbind(
      #Ya = unlist(predmat),
      Ya = do.call(c,predmat),
      psi = rep(intvals, each=MCsize),
      weights = rep(newdata$weights, times=length(intvals))
      #times=length(table(qdata[expnms[1]])))
    )
  )
  msmdat[,emmvars] <- newdata[,emmvars]
  polydat <- as.data.frame(poly(msmdat$psi, degree=degree, raw=TRUE))
  newexpnms <- paste0("psi",1:degree)
  names(polydat) <- newexpnms
  msmdat <- cbind(msmdat, polydat)
  msmf <- paste0("Ya ~ ",
                 ifelse(hasintercept, "1 +", "-1 +"),
                 paste0(c(newexpnms, emmvars), collapse = "+"))
  # TODO: categorical Z
  msmform <- .intmaker(as.formula(msmf), expnms = newexpnms, emmvars = emmvars, emmvar)
  class(msmform) <- "formula"
  newterms <- terms(msmform)
  nterms = length(attr(newterms, "term.labels"))
  nterms = nterms + attr(newterms, "intercept")

  # to do: allow functional form variations for the MSM via specifying the model formula
  if(bayes){
    if(!rr) suppressWarnings(msmfit <- bayesglm(msmform, data=msmdat,
                                                weights=weights, x=TRUE,
                                                ...))
    if(rr)  suppressWarnings(msmfit <- bayesglm(msmform, data=msmdat,
                                                family=binomial(link='log'), start=rep(-0.0001, nterms),
                                                weights=weights, x=TRUE))
  }
  if(!bayes){
    if(!rr) suppressWarnings(msmfit <- glm(msmform, data=msmdat,
                                           weights=weights, x=TRUE,
                                           ...))
    if(rr)  suppressWarnings(msmfit <- glm(msmform, data=msmdat,
                                           family=binomial(link='log'), start=rep(-0.0001, nterms),
                                           weights=weights, x=TRUE))
  }
  res <- list(fit=fit, msmfit=msmfit)
  if(main) {
    res$Ya <- msmdat$Ya   # expected outcome under joint exposure, by gcomp
    res$Yamsm <- as.numeric(predict(msmfit, type='response'))
    res$Yamsml <- as.numeric(predict(msmfit, type="link"))
    res$A <- msmdat$psi # joint exposure (0 = all exposures set category with
    res[[emmvar]] <- do.call(c, lapply(intvals, function(x) newdata[,emmvar,drop=TRUE]))
    # upper cut-point as first quantile)
  }
  #newterms <- terms(msmform)
  #prodterms <- do.call(c, lapply(1:length(emmvars), function(x) c(emmvars[x], paste0(emmvars[x], ":mixture"))))
  newtermlabels <- attr(newterms, "term.labels")
  #newtermlabels[(degree+1):length(newtermlabels)] <- prodterms
  for(emmv in emmvars){
    newtermlabels <- gsub(paste0("psi([0-9]):", emmv), paste0(emmv,":","mixture", "^\\1"), newtermlabels)
  }
  newtermlabels <- gsub("\\^1", "", newtermlabels)
  attr(res, "term.labels") <- newtermlabels
  res
}



qgcomp.emm.glm.boot <- function(
  f,
  data,
  expnms=NULL,
  emmvar="",
  q=4,
  breaks=NULL,
  id=NULL,
  weights,
  alpha=0.05,
  B=200,
  rr=TRUE,
  degree=1,
  seed=NULL,
  bayes=FALSE,
  MCsize=nrow(data),
  parallel=FALSE,
  parplan = FALSE,
  errcheck=FALSE,
  ...){
  #' @title EMM for Quantile g-computation for continuous, binary, and count outcomes under linearity/additivity
  #'
  #' @description This function fits a quantile g-computation model, allowing
  #' effect measure modification by a binary or continuous covariate. This allows
  #' testing of statistical interaction as well as estimation of stratum specific effects.
  #' This particular implementation formally fits a marginal structural model using
  #' a Monte Carlo-based g-computation method, utilizing bootstrapping for variance
  #' estimates. Because this approach allows for non-linear/non-additive effects of
  #' exposures, it does not report weights nor EMM stratum specific effects.
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
  #' id/cluster). Note that qgcomp.emm.glm.noboot will not produce cluster-appropriate
  #' standard errors (this parameter is essentially ignored in qgcomp.emm.glm.noboot).
  #' Qgcomp.emm.boot can be used for this, which will use bootstrap
  #' sampling of clusters/individuals to estimate cluster-appropriate standard
  #' errors via bootstrapping.
  #' @param weights "case weights" - passed to the "weight" argument of
  #' \code{\link[stats]{glm}} or \code{\link[arm]{bayesglm}}
  #' @param alpha alpha level for confidence limit calculation
  #' @param B integer: number of bootstrap iterations (this should typically be >=200,
  #'  though it is set lower in examples to improve run-time).
  #' @param rr logical: if using binary outcome and rr=TRUE, qgcomp.boot will
  #'   estimate risk ratio rather than odds ratio
  #' @param degree polynomial bases for marginal model (e.g. degree = 2
  #'  allows that the relationship between the whole exposure mixture and the outcome
  #'  is quadratic (default = 1).
  #' @param seed integer or NULL: random number seed for replicable bootstrap results
  #' @param bayes use underlying Bayesian model (`arm` package defaults). Results
  #' in penalized parameter estimation that can help with very highly correlated
  #' exposures. Note: this does not lead to fully Bayesian inference in general,
  #' so results should be interpreted as frequentist.
  #' @param MCsize integer: sample size for simulation to approximate marginal
  #'  zero inflated model parameters. This can be left small for testing, but should be as large
  #'  as needed to reduce simulation error to an acceptable magnitude (can compare psi coefficients for
  #'  linear fits with qgcomp.noboot to gain some intuition for the level of expected simulation
  #'  error at a given value of MCsize). This likely won't matter much in linear models, but may
  #'  be important with binary or count outcomes.
  #' @param parallel use (safe) parallel processing from the future and future.apply packages
  #' @param parplan (logical, default=FALSE) automatically set future::plan to plan(multisession) (and set to existing plan after bootstrapping)
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
  #' @importFrom  qgcomp quantize
  #' @importFrom future plan
  #' @importFrom future.apply future_lapply
  #' @export
  #' @examples
  #' set.seed(50)
  #' # linear model, binary modifier
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
  #'   z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
  #' (qfit <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
  #'   expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian()))
  #' # set B larger for real examples
  #' (qfit2 <- qgcomp.emm.boot(f=y ~ z + x1 + x2, emmvar="z",
  #'   degree = 1,
  #'   expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian(), B=10))
  #' # categorical modifier
  #' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
  #'   z=sample(0:2, 50,replace=TRUE), r=rbinom(50,1,0.5))
  #' dat2$z = as.factor(dat2$z)
  #' (qfit3 <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
  #'   expnms = c('x1', 'x2'), data=dat2, q=4, family=gaussian()))
  #' # set B larger for real examples
  #' (qfit4 <- qgcomp.emm.boot(f=y ~ z + x1 + x2, emmvar="z",
  #'   degree = 1,
  #'   expnms = c('x1', 'x2'), data=dat2, q=4, family=gaussian(), B=10))
  oldq = NULL
  if(is.null(seed)) seed = round(runif(1, min=0, max=1e8))

  if(errcheck){
    # basic argument checks
    if (is.null(expnms)) {
      stop("'expnms' must be specified explicitly\n")
    }
    if (is.null(emmvar)) {
      stop("'emmvar' must be specified explicitly\n")
    }
  }
  # housekeeping
  allemmvals<- unique(data[,emmvar,drop=TRUE])
  emmlev <- length(allemmvals)
  zdata = zproc(data[,emmvar], znm = emmvar)
  emmvars = names(zdata)
  data = cbind(data, zdata)
  data = data[,unique(names(data)),drop=FALSE]
  if(errcheck){
    # placeholder
  }
  # keep track of added terms by remembering old model
  originalform <- terms(f, data = data)
  hasintercept = as.logical(attr(originalform, "intercept"))

  #f = .intmaker(f,expnms,emmvar) # create necessary interaction terms with exposure
  (f <- .intmaker(f,expnms,emmvars, emmvar)) # create necessary interaction terms with exposure
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(originalform, "term.labels"))
  addedmain <- setdiff(addedterms, grep(":",addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, value = TRUE))
  addedintsord = addedints

  class(newform) <- "formula"

  nobs = nrow(data)
  ##### #######
  origcall <- thecall <- match.call(expand.dots = FALSE)
  names(thecall) <- gsub("f", "formula", names(thecall))
  m <- match(c("formula", "data", "weights", "offset"), names(thecall), 0L)
  hasweights = ifelse("weights" %in% names(thecall), TRUE, FALSE)
  thecall <- thecall[c(1L, m)]
  thecall$drop.unused.levels <- TRUE

  thecall[[1L]] <- quote(stats::model.frame)
  thecalle <- eval(thecall, parent.frame())
  if(hasweights){
    data$weights <- as.vector(model.weights(thecalle))
  } else data$weights = rep(1, nobs)


  if (is.null(expnms)) {

    #expnms <- attr(terms(f, data = data), "term.labels")
    expnms <- attr(newform, "term.labels")

    message("Including all model terms as exposures of interest\n")
  }
  lin = .intchecknames(expnms)
  if(!lin) stop("Model appears to be non-linear and I'm having trouble parsing it:
                  please use `expnms` parameter to define the variables making up the exposure")
  if (!is.null(q) & !is.null(breaks)){
    # if user specifies breaks, prioritize those
    oldq = q
    q <- NULL
  }
  if (!is.null(q) | !is.null(breaks)){
    ql <- qgcomp::quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
    if(is.null(q)){
      # rare scenario with user specified breaks and q is left at NULL
      nvals <- length(br[[1]])-1
    } else{
      nvals <- q
    }
    intvals <- (seq_len(nvals))-1
  } else {
    # if( is.null(breaks) & is.null(q)) # also includes NA
    qdata <- data[,unique(names(data)),drop=FALSE]
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # : allow user specification of this
    nvals = length(table(unlist(data[,expnms])))
    if(nvals < 10){
      message("\nNote: using all possible values of exposure as the
              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[,expnms]))))

      br <- lapply(seq_len(p), function(x) c(-1e16, intvals[2:nvals]-1e-16, 1e16))
    }else{
      message("\nNote: using quantiles of all exposures combined in order to set
          proposed intervention values for overall effect (25th, 50th, 75th %ile)
        You can ensure this is valid by scaling all variables in expnms to have similar ranges.")
      intvals = as.numeric(quantile(unlist(data[,expnms]), c(.25, .5, .75)))
      br <- NULL
    }
  }
  if(is.null(id)) {
    id <- "id__"
    qdata$id__ <- seq_len(dim(qdata)[1])
  }
  ###
  msmfit <- msm.emm.fit(newform, qdata, intvals, emmvar=emmvar, emmvars=emmvars, expnms=expnms, rr, main=TRUE,degree=degree, id=id,
                    weights,
                    bayes,
                    MCsize=MCsize,
                    ...)
  # main estimate
  msmcoefnames <- attr(msmfit, "term.labels")
  estb <- as.numeric(msmfit$msmfit$coefficients)
  #bootstrap to get std. error
  nobs <- dim(qdata)[1]
  nids <- length(unique(qdata[,id, drop=TRUE]))
  starttime = Sys.time()
  psi.emm.only <- function(i=1, f=f, qdata=qdata, intvals=intvals, emmvar=emmvar, emmvars=emmvars, expnms=expnms, rr=rr, degree=degree,
                       nids=nids, id=id,
                       weights,MCsize=MCsize,
                       ...){
    if(i==2 & !parallel){
      timeiter = as.numeric(Sys.time() - starttime)
      if((timeiter*B/60)>0.5) message(paste0("Expected time to finish: ", round(B*timeiter/60, 2), " minutes \n"))
    }
    bootids <- data.frame(temp=sort(sample(unique(qdata[,id, drop=TRUE]), nids,
                                           replace = TRUE
    )))
    names(bootids) <- id
    qdata_ <- merge(qdata,bootids, by=id, all.x=FALSE, all.y=TRUE)
    ft = msm.emm.fit(f, qdata_, intvals=intvals, expnms=expnms, emmvar=emmvar, emmvars=emmvars, rr, main=FALSE, degree, id, weights=weights, bayes, MCsize=MCsize,
                 ...)
    yhatty = data.frame(yhat=predict(ft$msmfit), psi=ft$msmfit$data[,"psi"])
    as.numeric(
      # the yhat estimates will be identical across individuals due to this being a marginal model
      c(
        with(yhatty, tapply(yhat, psi, mean)), # yhats
        ft$msmfit$coefficients                 # coefficients
        )
    )
  }
  set.seed(seed)
  if(parallel){
    #Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
    if (parplan) {
        oplan <- future::plan(strategy = future::multisession)
        on.exit(future::plan(oplan), add = TRUE)
    }

    bootsamps <- future.apply::future_lapply(X=seq_len(B), FUN=psi.emm.only,f=f, qdata=qdata, intvals=intvals,
                                             emmvar=emmvar, emmvars=emmvars, expnms=expnms, rr=rr, degree=degree, nids=nids, id=id,
                                             weights=qdata$weights,MCsize=MCsize,
                                             future.seed=TRUE,
                                             ...)


  }else{
    bootsamps <- lapply(X=seq_len(B), FUN=psi.emm.only,f=f, qdata=qdata, intvals=intvals,
                        emmvar=emmvar, emmvars=emmvars, expnms=expnms, rr=rr, degree=degree, nids=nids, id=id,
                        weights=weights, MCsize=MCsize,
                        ...)

  }
  bootsamps = do.call("cbind", bootsamps)
  # these are the linear predictors (predictions at each discrete value of exposure)
  hatidx = seq_len(length(intvals))
  hats = t(bootsamps[hatidx,])
  # covariance of the linear predictors
  cov.yhat = cov(hats)
  # coefficients
  bootsamps = bootsamps[-hatidx,]
  #rownames(bootsamps) = msmcoefnames
  seb <- apply(bootsamps, 1, sd)
  covmat.coef <- cov(t(bootsamps))
  #print(msmcoefnames)
  #print(names(estb))
  #print(rownames(bootsamps))
  colnames(covmat.coef) <- rownames(covmat.coef) <- names(estb) <- c("(Intercept)", msmcoefnames)
  #colnames(covmat) <- rownames(covmat) <- names(estb) <- c("(intercept)", paste0("psi", seq_len(nrow(bootsamps)-1)))

  tstat <- estb / seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(estb + seb * qnorm(alpha / 2), estb + seb * qnorm(1 - alpha / 2))
  # outcome 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity,
  #   then weights will vary with level of exposure)
  if (!is.null(oldq)){
    q = oldq
  }
  psidx = 1:(hasintercept+1)
  qx <- qdata[, expnms]
  res <- .qgcompemm_object(
    qx = qx, fit = msmfit$fit, msmfit = msmfit$msmfit,
    psi = estb[-1],
    var.psi = seb[-1] ^ 2,
    covmat.psi=covmat.coef["psi1", "psi1"],
    covmat.psiint=covmat.coef[grep("mixture", colnames(covmat.coef)), grep("mixture", colnames(covmat.coef))], # to fix
    ci = ci[-1,],
    coef = estb, var.coef = seb ^ 2, covmat.coef=covmat.coef, ci.coef = ci,
    expnms=expnms,
    intterms = addedintsord,
    q=q, breaks=br, degree=degree,
    pos.psi = NULL, neg.psi = NULL,
    pos.weights = NULL, neg.weights = NULL, pos.size = NULL,neg.size = NULL, bootstrap=TRUE,
    y.expected=msmfit$Ya, y.expectedmsm=msmfit$Yamsm, index=msmfit$A,
    emmvar.msm = msmfit[[emmvar]],
    bootsamps = bootsamps,
    cov.yhat=cov.yhat,
    alpha=alpha,
    call=origcall,
    emmlev = emmlev
  )
  if(msmfit$fit$family$family=='gaussian'){
    res$tstat <- tstat
    res$df <- df
    res$pval <- pval
  }
  if(msmfit$fit$family$family %in% c('binomial', 'poisson')){
    res$zstat <- tstat
    res$pval <- pvalz
  }
  res
}

#' @rdname qgcomp.emm.glm.boot
#' @export
qgcomp.emm.boot <- qgcomp.emm.glm.boot


