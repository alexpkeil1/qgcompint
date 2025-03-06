qgcomp.emm.cox.noboot <- function (
  f,
  data,
  expnms = NULL,
  emmvar=NULL,
  q = 4,
  breaks = NULL,
  id=NULL,
  weights,
  cluster=NULL,
  alpha=0.05,
  errcheck=TRUE,
  ...
) {
  #' @title EMM for Quantile g-computation with survival outcomes under linearity/additivity
  #'
  #' @description This function performs quantile g-computation in a survival
  #' setting, , allowing
  #' effect measure modification by a binary, categorical or continuous covariate. This allows
  #' testing of statistical interaction as well as estimation of stratum specific effects.
  #'
  #' @param f R style survival formula, which includes \code{\link[survival]{Surv}}
  #'   in the outcome definition. E.g. \code{Surv(time,event) ~ exposure}. Offset
  #'   terms can be included via \code{Surv(time,event) ~ exposure + offset(z)}
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
  #' id/cluster)
  #' @param weights "case weights" - passed to the "weight" argument of
  #' \code{\link[survival]{coxph}}
  #' @param cluster not yet implemented
  #' @param alpha alpha level for confidence limit calculation
  #' @param errcheck (logical, default=TRUE) include some basic error checking. Slightly
  #' faster if set to false (but be sure you understand risks)
  #' @param ... arguments to glm (e.g. family)
  #' @seealso \code{\link[qgcomp]{qgcomp.cox.noboot}}
  #' @return a qgcompfit object, which contains information about the effect
  #'  measure of interest (psi) and associated variance (var.psi), as well
  #'  as information on the model fit (fit) and information on the
  #'  weights/standardized coefficients in the positive (pos.weights) and
  #'  negative (neg.weights) directions.
  #' @concept variance mixtures
  #' @import survival
  #' @importFrom  qgcomp quantize
  #' @export
  #' @examples
  #' set.seed(5)
  #' N=200
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))),
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N), z=runif(N))
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2+z
  #' (fit1 <- survival::coxph(f, data = dat))
  #' (obj <- qgcomp.emm.cox.noboot(f, expnms = expnms, emmvar="z", data = dat))
  #'
  #' #categorical emm
  #' dat <- data.frame(time=(tmg <- pmin(.1,rweibull(N, 10, 0.1))),
  #'                 d=1.0*(tmg<0.1), x1=runif(N), x2=runif(N),
  #'                 z=sample(0:2, N, replace=TRUE))
  #'  dat$z = as.factor(dat$z)
  #' expnms=paste0("x", 1:2)
  #' f = survival::Surv(time, d)~x1 + x2+z
  #' (obj2 <- qgcomp.emm.cox.noboot(f, expnms = expnms, emmvar="z", data = dat))
  #'
  if(errcheck){
    # basic argument checks
    if (is.null(expnms)) {
      stop("'expnms' must be specified explicitly\n")
    }
    if (is.null(emmvar)) {
      stop("'emmvar' must be specified explicitly\n")
    }
  }
  allemmvals<- unique(data[,emmvar,drop=TRUE])
  emmlev <- length(allemmvals)
  ## process to expand factors if needed
  zdata = zproc(data[,emmvar], znm = emmvar)
  emmvars = names(zdata)
  data = cbind(data, zdata)
  ### end new
  # housekeeping
  #of <- f
  # keep track of added terms by remembering old model
  originalform <- terms(f, data = data)
  #f = .intmaker(f,expnms,emmvar) # create necessary interaction terms with exposure
  (f <- .intmaker(f,expnms,emmvars, emmvar)) # create necessary interaction terms with exposure
  newform <- terms(f, data = data)
  class(newform) <- "formula"
  addedterms <- setdiff(attr(newform, "term.labels"), attr(originalform, "term.labels"))
  addedmain <- setdiff(addedterms, grep(":",addedterms, value = TRUE))
  addedints <- setdiff(addedterms, addedmain)
  addedintsl <- lapply(emmvars, function(x) grep(x, addedints, value = TRUE))
  if (length(addedmain)>0) {
    message(paste0("Adding main term for ",emmvar," to the model\n"))
  }
  #oord <- order(expnms)
  # order interaction terms in same order as main terms
  #s0 <- gsub(paste0("^", emmvar,":"), "",
  #           gsub(paste0(":", emmvar,"$"), "", addedints))
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
  if (is.null(expnms)) {
    message("Including all model terms as exposures of interest")
    expnms <- attr(newform, "term.labels")
  }
  lin = .intchecknames(expnms)
  if(!lin) stop("Model appears to be non-linear: this is not yet implemented")
  if (!is.null(q) | !is.null(breaks)) {
    ql <- qgcomp::quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  }
  else {
    qdata <- data
    br <- breaks
  }
  # original fit
  environment(newform) <- list2env(list(qdata=qdata))
  #newform = Surv(time, d) ~ x1 + x2 + x1 * z + x2 * z
  fit <- coxph(newform, data = qdata,
               weights=weights, x=FALSE,y=FALSE,
               #cluster=cluster,
               ...)
  fit$data = data.frame(qdata)
  coxfam = list(family='cox', link='log', linkfun=log)
  class(coxfam) = "family"
  fit[['family']] = coxfam # kludge for print function
  mod <- summary(fit)
  #
  estb <- sum(
    mod$coefficients[expnms, 1]
    )
  names(estb) = "psi1"
  covMat = fit$var
  colnames(covMat) <- rownames(covMat) <- names(mod$coefficients[, 1])
  seb <- se_comb2(
    expnms, covmat = covMat
    )
  tstat <- estb/seb
  #df <- mod$df.null - length(expnms)
  #pval <- 2 - 2 * pt(abs(tstat), df = df)
  pvalz <- 2 - 2 * pnorm(abs(tstat))
  ci <- c(
    estb + seb * qnorm(alpha/2),
    estb + seb * qnorm(1 -alpha/2)
  )
  # modifier main term, product term
  estb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
    fit$coefficients[emmvars[x]],
    sum(mod$coefficients[addedintsl[[x]],1, drop=TRUE])
  )))
  names(estb.prod) <- do.call(c, lapply(1:length(emmvars), function(x) c(emmvars[x], paste0(emmvars[x], ":mixture"))))
  seb.prod <- do.call(c, lapply(1:length(emmvars), function(x) c(
    sqrt(covMat[emmvars[x],emmvars[x]]),
    se_comb2(addedintsl[[x]], covmat = covMat)
  )))
  #estb.prod <- c(
  #  mod$coefficients[emmvar, 1],
  #  sum(mod$coefficients[addedints,1, drop=TRUE])
  #)
  #names(estb.prod) <- c(emmvar, paste0(emmvar, ":mixture"))
  #seb.prod <- c(
  #  sqrt(covMat[emmvar,emmvar]),
  #  se_comb2(addedints, covmat = covMat)
  #)
  tstat.prod <- estb.prod / seb.prod
  #pval.prod <- 2 - 2 * pt(abs(tstat.prod), df = df)
  pvalz.prod <- 2 - 2 * pnorm(abs(tstat.prod))
  ci.prod <- cbind(
    estb.prod + seb.prod * qnorm(alpha / 2),
    estb.prod + seb.prod * qnorm(1 - alpha / 2)
  )
  # 'weights' at referent level
  wcoef <- fit$coefficients[expnms]
  names(wcoef) <- gsub("_q", "", names(wcoef))
  poscoef <- which(wcoef > 0)
  negcoef <- which(wcoef <= 0)
  pos.weights <- abs(wcoef[poscoef])/sum(abs(wcoef[poscoef]))
  neg.weights <- abs(wcoef[negcoef])/sum(abs(wcoef[negcoef]))
  # 'post-hoc' positive and negative estimators
  # similar to constrained gWQS
  pos.psi <- sum(wcoef[poscoef])
  neg.psi <- sum(wcoef[negcoef])
  qx <- qdata[, expnms]
  names(qx) <- paste0(names(qx), "_q")
  covmat.coef = vc_multiscomb(inames = NULL, emmvars=emmvars,
                              expnms=expnms,addedintsl=addedintsl, covmat=covMat, grad = NULL
  )
  colnames(covmat.coef) <- rownames(covmat.coef) <- names(c(estb, estb.prod))

  #
  res <- .qgcompemm_object(
    qx = qx, fit = fit,
    psi = estb,
    psiint = estb.prod[2*(1:length(emmvars))],
    var.psi = seb^2,
    var.psiint = seb.prod[2*(1:length(emmvars))] ^ 2,
    covmat.psi=covmat.coef["psi1", "psi1"],
    covmat.psiint=covmat.coef[grep("mixture", colnames(covmat.coef)), grep("mixture", colnames(covmat.coef))], # to fix
    ci = ci,
    ciint = ci.prod[2*(1:length(emmvars)),],
    coef = c(estb, estb.prod),
    var.coef = c(seb ^ 2, seb.prod ^ 2),
    covmat.coef = covmat.coef,
    ci.coef = rbind(ci, ci.prod),
    expnms = expnms,
    intterms = addedintsord,
    q = q,
    breaks = br,
    degree = 1,
    pos.psi = pos.psi,
    neg.psi = neg.psi,
    pos.weights = sort(pos.weights, decreasing = TRUE),
    neg.weights = sort(neg.weights, decreasing = TRUE),
    pos.size = sum(abs(wcoef[poscoef])),
    neg.size = sum(abs(wcoef[negcoef])),
    bootstrap = FALSE,
    zstat = c(tstat,tstat.prod),
    pval = c(pvalz,pvalz.prod),
    alpha=alpha,
    call=origcall,
    emmlev = emmlev
  )
  if(emmlev==2){
    ww = getstratweights(res, emmval = 1)
    ff = getstrateffects(res, emmval = 1)
    cl = class(res)
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
    class(res) = cl
  }
  res$hasintercept = FALSE
  attr(res, "class") <- c( "survqgcompemmfit", "survqgcompfit", attr(res, "class"))
  res
}
