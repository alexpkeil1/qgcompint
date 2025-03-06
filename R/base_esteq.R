.esteq_qgc <- qgcomp:::.esteq_qgc

.esteq_qgcemmdf <- function(f, data, theta, family, intvals, expnms, emmvar, emmvars, hasintercept, weights, degree=1,rr=FALSE,offset=0, delta=-Inf, ...){
  # family specific estimating equation based on input dataframes
  # used for "B" part of sandwich variance: V = solve(A) %*% B %&% t(solve(A)) because it facilitates
  # individual level calculations.
  fam = family$family
  lnk = family$link
  X = model.matrix(f, data) # model.frame
  #modframe = model.frame(f, data=data)
  Y = model.response(data)
  Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {data[,expnms] = x; model.matrix(f,data)})))
  Xmsm = poly(Xint[, expnms[1]], degree=degree, raw=TRUE) # intercept and constant exposure
  if(hasintercept){
    Xmsm = cbind(Xint[,colnames(Xint)[1]], Xmsm)
  }
  Zmsm = do.call(rbind, lapply(intvals, function(x) X[,emmvars,drop=FALSE]))
  msmdf = cbind(data.frame(Xmsm), Zmsm, data.frame(fakey=1))
  msmf = fakey~X1
  (msmf <- .intmaker(msmf,"X1",emmvars, emmvar)) # create necessary interaction terms with exposure
  newmsmform <- terms(msmf, data = msmdf)
  hasintercept = as.logical(attr(newmsmform, "intercept"))
  class(newmsmform) <- "formula"
  msmvars = get_all_vars(newmsmform, data=msmdf) # variables before being processed by model.frame
  msmmodframe = model.frame(msmvars, data=msmdf)
  Xmsm = model.matrix(newmsmform, msmmodframe)

  binlink = ifelse(rr, "logitlog", "logit")
  FUN <- switch(fam,
                gaussian = qgcomp:::.esteq_qgclin,
                tobit = qgcomp:::.esteq_qgctobit,
                poisson = qgcomp:::.esteq_qgcpoisson,
                #binomial = .esteq_qgclogit(theta, Y, X, Xint, Xmsm, weights),
                binomial = switch(binlink,
                                  logit = qgcomp:::.esteq_qgclogit,
                                  logitlog = qgcomp:::.esteq_qgclogitlog
                ),
                qgcomp:::.esteq_error
  )
  FUN(theta=theta, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=weights, offset=offset, delta=-Inf, ...)
}

#' @title EMM for Quantile g-computation for continuous, binary, and count outcomes under linearity/additivity
#'
#' @description This function estimates a dose-response parameter representing a one quantile
#' increase in a set of exposures of interest at levels of a binary, factor, or continuous covariate.
#' This allows
#' testing of statistical interaction as well as estimation of stratum specific effects.
#' This model estimates the parameters of a marginal
#' structural model (MSM) based on g-computation with quantized exposures.
#' Note: this function allows clustering of data and/or sampling weights and yields cluster-robust standard
#' errors for all estimates
#'
#' Estimating equation methodology is used as the underlying estimation scheme. This allows that observations can be correlated, and is fundementally identical to some implementations of "generalized estimating equations" or GEE. Thus, it allows for a more general set of longitudinal data structures than does the qgcomp.glm.noboot function, because it allows that the outcome can vary over time within an individual. Interpretation of parameters is similar to that of a GEE: this function yields population average estimates of the effect of exposure over time.
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
#' @param rr (logical, default=TRUE) estimate a risk ratio from the MSM when using an underlying logistic model
#' @param degree polynomial degree for non-linearity of MSM (values other than 1 are not supported, currently)
#' @param includeX (logical, default=TRUE) include data in the output?
#' @param verbose (logical, default=TRUE) include informative messages
#' @param errcheck (logical, default=TRUE) include some basic error checking. Slightly
#' faster if set to false (but be sure you understand risks)
#' @param ... arguments to glm (e.g. family)
#' @seealso \code{\link[qgcomp]{qgcomp.noboot}}
#' @return a eeqgcompfit/qgcompemmfit/qgcompfit/list object, which contains information about the effect
#'  measure of interest (psi) and associated variance (var.psi), as well
#'  as information on the conditional/underlying model fit (fit) and the marginal structural model fit (msmfit).
#' @concept variance mixtures
#' @importFrom rootSolve multiroot
#' @importFrom numDeriv jacobian
#' @export
#' @examples
#' set.seed(50)
#' # linear model, binary modifier
#' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#'   z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
#' (qfit <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#'   expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
#' (qfitee <- qgcomp.emm.glm.ee(f=y ~ z + x1 + x2, emmvar="z",
#'   expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
#' # logistic model, continuous modifier
#' dat2 <- data.frame(y=rbinom(50, 1,0.5), x1=runif(50), x2=runif(50),
#'   z=runif(50), r=rbinom(50,1,0.5))
#' (qfit2 <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#'   expnms = c('x1', 'x2'), data=dat2, q=2, family=binomial()))
#' # get weights and stratum specific effects at specific value of Z
#' #  (note that when Z=0, the effect is equal to psi1)
qgcomp.emm.glm.ee <- function(
    f,
    data,
    expnms=NULL,
    emmvar=NULL,
    q=4,
    breaks=NULL,
    id=NULL,
    weights,
    offset=NULL,
    alpha=0.05,
    rr=TRUE,
    degree=1,
    includeX=TRUE,
    verbose=TRUE,
    errcheck = TRUE,
    ...
){
  if(degree>1){
    message("Degree > 1 not supported yet in this function, setting to 1")
    degree = 1
  }
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
  allemmvals<- unique(data[,emmvar,drop=TRUE])
  emmlev <- length(allemmvals)
  ## process to expand factors if needed
  zdata = zproc(data[,emmvar,drop=TRUE], znm = emmvar)
  emmvars = names(zdata)
  data = cbind(data, zdata)

  # keep track of added terms by remembering old model
  originalform <- terms.formula(f, data = data)
  # need to explicitly expand factors because there is no automatic processing by lm or glm
 mm = model.matrix(f, data=data[1,,drop=FALSE])
 expandedterms = setdiff(colnames(mm), "(Intercept)")
 expandedfright = paste0(expandedterms, collapse="+")
 updatedf = as.formula(paste0(originalform[[2]], "~", expandedfright))
 #originalform = terms.formula(updatedf, data = data)

  hasintercept = as.logical(attr(originalform, "intercept"))
  #f = .intmaker(f,expnms,emmvar) # create necessary interaction terms with exposure
  (f <- .intmaker(updatedf,expnms,emmvars, emmvar)) # create necessary interaction terms with exposure
  #(f <- .intmaker(f,expnms,emmvars, emmvar)) # create necessary interaction terms with exposure
  newform <- terms(f, data = data)
  addedterms <- setdiff(attr(newform, "term.labels"), attr(originalform, "term.labels"))
  addedmain <- setdiff(addedterms, grep(":",addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, value = TRUE))
  addedintsord = addedints

  oldq = NULL

  newform <- terms(f, data = data)
  hasintercept = as.logical(attr(newform, "intercept"))
  class(newform) <- "formula"

  cc = match.call(expand.dots = TRUE)

  if(is.null(cc$delta)){
    delta=-Inf
  } else{
    delta = eval(cc$delta)
  }
  if(!is.null(cc$family) && cc$family != "tobit"){
    testfit <- glm(y~., data=data.frame(y=c(0,1)), eval(cc$family))
    family = testfit$family
  }
  if(!is.null(cc$family) && cc$family == "tobit"){
    family  = tobit()
  }
  if(is.null(cc$family)){
    family=gaussian()
  }
  #family = testfit$family

  #  famlist = c("binomial", "gaussian", "poisson")
  #  if(!(family$family %in% famlist))
  #    stop(paste0("Distribution (family) `", family$family, "` not (yet?) supported."))


  nobs = nrow(data)
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
    if(verbose) {
      message("Including all model terms as exposures of interest\n")
    }
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
    qdata <- data
    # if no transformation is made (no quantiles, no breaks given)
    # then draw distribution values from quantiles of all the exposures
    # pooled together
    # : allow user specification of this
    nvals = length(table(unlist(data[,expnms])))
    if(nvals < 10){
      if(verbose) message("\nNote: using all possible values of exposure as the
              intervention values\n")
      p = length(expnms)
      intvals <- as.numeric(names(table(unlist(data[,expnms]))))

      br <- lapply(seq_len(p), function(x) c(-1e16, intvals[2:nvals]-1e-16, 1e16))
    }else{
      if(verbose) message("\nNote: using quantiles of all exposures combined in order to set
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
  if(is.null(offset)) {
    offset <- "offset__"
    qdata$offset__ <- rep(0,dim(qdata)[1])
  }



  # simultaneous estimation of conditional and
  basevars = get_all_vars(newform, data=qdata) # variables before being processed by model.frame
  #modframe = model.frame(newform, data=basevars)
  modframe = model.frame(newform, data=qdata)
  X = model.matrix(newform, modframe)
  Y = model.response(modframe)
  # nesting model.frame within the model.matrix function seems to be necessary to get interaction terms to propogate after setting exposures
  #  Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {modframe[,expnms] = x; model.matrix(newform,model.frame(newform, modframe))})))
  #Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {modframe[,expnms] = x; model.matrix(newform,modframe)})))
  #Xint = as.matrix(do.call(rbind,lapply(intvals, function(x) {mf2 = modframe; mf2[,expnms] = x; model.matrix(newform,model.frame(newform, data=mf2))}))) # works in non-linear setting
  Xint = as.matrix(model.matrix(newform, do.call(rbind,lapply(intvals, function(x) {mf2 = basevars; mf2[,expnms] = x; model.frame(newform, data=mf2)})))) # works in non-linear setting

  # TODO: make this work for an interaction term

  Xmsm = poly(Xint[, expnms[1]], degree=degree, raw=TRUE) # intercept and constant exposure
  colnames(Xmsm) = paste0("mixture")
  if(hasintercept){
    Xmsm = cbind(Xint[,colnames(Xint)[1]], Xmsm)
  }
  Zmsm = do.call(rbind, lapply(intvals, function(x) modframe[,emmvars,drop=FALSE]))
  msmdf = cbind(data.frame(Xmsm), Zmsm, data.frame(fakey=1))
  msmf = fakey~mixture



  (msmf <- .intmaker(msmf,"mixture",emmvars,emmvar)) # create necessary interaction terms with exposure
  newmsmform <- terms(msmf, data = msmdf)
  hasintercept = as.logical(attr(newmsmform, "intercept"))
  class(newmsmform) <- "formula"


  msmvars = get_all_vars(newmsmform, data=msmdf) # variables before being processed by model.frame
  msmmodframe = model.frame(msmvars, data=msmdf)
  Xmsm = model.matrix(newmsmform, msmmodframe)
  #Y = model.response(modframe)


  # point estimates
  #.esteq_qgc(parminits, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights, rr=FALSE)
  np = ncol(X)
  npmsm = ncol(Xmsm)

  startvals <- function(family,Y,X,np,npmsm,offset=0){
    fam = family$family
    res = switch(fam,
                 binomial = c(log(mean(Y))-(mean(offset)), rep(0, np-1), log(mean(Y))-(mean(offset)), rep(0, npmsm-1)),
                 poisson = c(log(mean(Y))-(mean(offset)), rep(0, np-1), log(mean(Y))-(mean(offset)), rep(0, npmsm-1)),
                 tobit = c(rep(0, np+npmsm),1),
                 rep(0, np+npmsm)
    )
    res
  }
  parminits = startvals(family,Y,X,np,npmsm)
  #.esteq_qgc(parminits, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights, rr=FALSE)
  eqfit <- rootSolve::multiroot(.esteq_qgc, start=parminits, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights, rr=rr, offset=qdata$offset__, delta=delta)
  # "bread" of the sandwich covariance
  A = numDeriv::jacobian(func=.esteq_qgc, x=eqfit$root, family=family, Y=Y, X=X, Xint=Xint, Xmsm=Xmsm, weights=qdata$weights,method="Richardson", rr=rr,offset=qdata$offset__, delta=delta)
  #

  # "meat" of the sandwich covariance
  uid =   unique(qdata[,id,drop=TRUE])
  psii = lapply(uid, function(x){
    selidx = which(qdata[,id,drop=TRUE] == x)
    .esteq_qgcemmdf(newform, data=modframe[selidx,,drop=FALSE], theta=eqfit$root, family, intvals, expnms, emmvar, emmvars, hasintercept, weights=qdata$weights[selidx], degree=degree, rr=rr,offset=qdata$offset__[selidx], delta=delta)
  } )
  Bi = lapply(psii, function(x) x%*%t(x))
  n = length(Y)
  B = Bi[[1]]
  for(i in 2:length(Bi)){
    B = B + Bi[[i]]
  }

  # sandwich covariance
  ibread = solve(A)
  (fullcovmat = ibread %*% B %*% t(ibread))

  condidx = 1:np
  msmidx = (np+1):(np+npmsm)
  estb <- as.numeric(eqfit$root[msmidx])
  nobs <- dim(qdata)[1]
  covmat <- fullcovmat[msmidx,msmidx,drop=FALSE]
  seb <- sqrt(diag(covmat))
  cnms = colnames(Xmsm)#c(paste0("psi", seq_len(length(msmidx)-1)))
  cnms[which(cnms == "mixture")] = "psi1"
  allest = eqfit$root
  names(allest) <- c(colnames(X), cnms)
  rownames(fullcovmat) <- colnames(fullcovmat) <- names(allest)


  #if(hasintercept)
  #  cnms = c("(Intercept)", cnms)
  colnames(covmat) <- rownames(covmat) <- names(estb) <- cnms

  tstat <- estb / seb
  df <- nobs - length(attr(terms(f, data = data), "term.labels")) - 1 - degree # df based on obs - gcomp terms - msm terms
  pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- cbind(
    estb + seb * qnorm(alpha / 2),
    estb + seb * qnorm(1 - alpha / 2)
  )

  # modifier main term, product term
  condcovmat = fullcovmat[condidx,condidx]
  #estb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
  #  coef(fit)[emmvars[x]],
  #  #sum(mod$coefficients[addedintsl[[x]],1, drop=TRUE])
  #  sum(coef(fit)[addedintsl[[x]]])
  #)))
  #names(estb.prod) <- do.call(c, lapply(1:length(emmvars), function(x) c(emmvars[x], paste0(emmvars[x], ":mixture"))))
  #seb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
  #  sqrt(condcovmat[emmvars[x],emmvars[x]]),
  #  se_comb2(addedintsl[[x]], covmat = condcovmat)
  #)))
  #tstat.prod <- estb.prod / seb.prod
  #pval.prod <- 2 - 2 * pt(abs(tstat.prod), df = df)
  #pvalz.prod <- 2 - 2 * pnorm(abs(tstat.prod))
  #ci.prod <- cbind(
  #  estb.prod + seb.prod * qnorm(alpha / 2),
  #  estb.prod + seb.prod * qnorm(1 - alpha / 2)
  #)
  # qgcomp 'weights' not applicable in this setting, generally (i.e. if using this function for non-linearity,
  #   then weights will vary with level of exposure)
  if (!is.null(oldq)){
    q = oldq
  }

  msmfamily = family
  if(msmfamily$family == "binomial" & rr == TRUE){
    msmfamily = binomial(link="log")
  }

  fit = list(formula = newform, est=allest[condidx], vcov=fullcovmat[condidx,condidx], family=family, type="conditional", data = qdata)
  msmfit = list(est=allest[msmidx], vcov=fullcovmat[msmidx,msmidx], family=msmfamily, type="msm")
  if(includeX){
    fit$X = X
    msmfit$X = Xmsm
  }

  attr(fit, "class") <- c("eefit", attr(fit, "class"))
  attr(msmfit, "class") <- c("eefit", attr(msmfit, "class"))

  psiidx = seq_len(degree)+ifelse(hasintercept,1,0)
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  res <- .qgcompemm_object(
    qx = qx, fit=fit, msmfit = msmfit,
    psi = estb[psiidx], var.psi = seb[psiidx] ^ 2,
    covmat.psi=covmat[psiidx,psiidx, drop=FALSE],
    ci = ci[psiidx,],
    coef = estb,
    var.coef = seb ^ 2,
    covmat.coef=covmat,
    ci.coef = ci,
    expnms=expnms,
    intterms = addedintsord,
    q=q,
    breaks=br,
    degree=degree,
    bootstrap=FALSE,
    y.expected = fit$family$linkinv(Xint %*% coef(fit)), index=Xint[,expnms[1]], # predictions from conditional fit at intervention data
    y.expectedmsm=predict(msmfit),
    bread = A, meat = B,
    covmat.all_robust = fullcovmat,
    alpha=alpha,
    call=origcall,
    emmlev = emmlev,
    hasintercept=hasintercept
  )
  # include some extra things by default for binary modifier (convenience only)

  if(emmlev==2){
    #ww = getstratweights(res, emmval = 1)
    ff = getstrateffects(res, emmval = 1)
    cl = class(res)
    res = c(res,
            list(
              effect = ff$eff,
              vareffect = ff$se^2,
              cieffect = ff$ci
            )
    )
    class(res) = cl
  }
  if(fit$family$family=='gaussian'){
    res$tstat <- c(tstat)
    res$df <- df
    res$pval <- c(pval)
  }
  if(fit$family$family %in% c('binomial', 'poisson')){
    res$zstat <- c(tstat)
    res$pval <- c(pvalz)
  }
  res
  attr(res, "class") <- c("eeqgcompfit", attr(res, "class"))
  res
}


#' @rdname qgcomp.emm.glm.ee
#' @export
qgcomp.emm.ee <- qgcomp.emm.glm.ee
