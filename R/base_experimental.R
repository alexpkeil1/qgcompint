#' @export
predict.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  message("Experimental feature, not validated")
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

#' @export
predictmsm.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...){
  message("Experimental feature, not validated")
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

#' @exportS3Method stats::anova
anova.eeqgcompfit = function (object, ..., dispersion = NULL, test = NULL)
{
  # based on geepack::anova.geeglm
  stop("Not yet implemented")
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("The following arguments to anova.glm(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.eefit <- unlist(lapply(dotargs, function(x) inherits(x,
                                                        "eeqgcompfit")))
  dotargs <- dotargs[is.eefit]
  if (length(dotargs) > 0)
    return(anova.eeqgcompfit(c(list(object), dotargs), dispersion = dispersion,
                            test = test))
  #varlist <- attr(object$terms, "variables")
  varlist = names(coef(object))
  x <- if (n <- match("x", names(object), 0)) {
    object[[n]]
  }else {
    #model.matrix(object)
    object$msmfit$X
  }
  varseq <- attr(x, "assign")
  nvars <- max(0, varseq)
  betaList <- vbetaList <- NULL
  if (nvars > 1) {
    method <- object$method
    if (!is.function(method))
      method <- get(method, mode = "function", envir = parent.frame())
    for (i in 1:(nvars - 1)) {
      eprint("calling fit....")
      # TODO: START HERE: fit model with decreasing number of variables (should just be underlying fit)
      fit <- method(x = x[, varseq <= i, drop = FALSE],
                    y = object$y, weights = object$prior.weights,
                    corstr = object$corstr, start = object$start,
                    offset = object$offset, id = object$id, family = object$family,
                    control = object$control)
      betaList <- c(betaList, list(fit$beta))
      vbetaList <- c(vbetaList, list(fit$vbeta))
    }
  }
  betaList <- c(betaList, list(object$geese$beta))
  vbetaList <- c(vbetaList, list(object$geese$vbeta))
  hasIntercept <- (length(grep("(Intercept)", names(betaList[[1]]))) !=
                     0)
  dimVec <- unlist(lapply(betaList, length))
  if (hasIntercept) {
    dfVec <- dimVec[1] - 1
  }
  else {
    dfVec <- dimVec[1]
  }
  if (length(dimVec) > 1) {
    for (i in 2:length(dimVec)) dfVec <- c(dfVec, dimVec[i] -
                                             dimVec[i - 1])
  }
  X2Vec <- NULL
  for (i in 1:length(dfVec)) {
    beta <- betaList[[i]]
    vbeta <- vbetaList[[i]]
    beta0 <- rep(1, length(beta))
    beta0[1:dfVec[i]] <- 0
    beta0 <- rev(beta0)
    zeroidx <- beta0 == 0
    X2 <- t(beta[zeroidx]) %*% solve(vbeta[zeroidx, zeroidx,
                                           drop = FALSE]) %*% beta[zeroidx]
    X2Vec <- c(X2Vec, X2)
  }
  resdf <- dfVec
  resdev <- X2Vec
  tab <- data.frame(resdf, resdev, 1 - pchisq(resdev, resdf))
  colnames(tab) <- c("Df", "X2", "P(>|Chi|)")
  tl <- attr(object$terms, "term.labels")
  if (length(tl) == 0)
    tab <- tab[1, , drop = FALSE]
  if (length(tl))
    rownames(tab) <- c(tl)
  title <- paste("Analysis of 'Wald statistic' Table", "\nModel: ",
                 object$family$family, ", link: ", object$family$link,
                 "\nResponse: ", as.character(varlist[-1])[1], "\nTerms added sequentially (first to last)\n",
                 sep = "")
  structure(tab, heading = title, class = c("anova", "data.frame"))
}


