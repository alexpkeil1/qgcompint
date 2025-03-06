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

anova.eeqgcompfit = function (object, ..., dispersion = NULL, test = NULL)
  #' @exportS3Method stats::anova
{
  # based on geepack:::anova.geeglm
  stop("Not yet implemented")
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("The following arguments to anova.glm(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.eefit <- unlist(lapply(dotargs, function(x) inherits(x, "eeqgcompfit")))
  dotargs <- dotargs[is.eefit]
  if (length(dotargs) > 0)
    return(anova.eeqgcompfitlist(c(list(object), dotargs), dispersion = dispersion,
                            test = test))
  else{
    stop("Two model fits are needed")
  }
}


anova.eeqgcompfitlist <- function (object, ..., dispersion = NULL, test = NULL)
  #' @exportS3Method stats::anova
{
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2]])
  }))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]),
            " removed because response differs from ", "model 1")
  }
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1]))
    stop("models were not all fitted to the same size of dataset")
  objects <- list(object, ...)
  m1 <- objects[[1]][[1]]
  if (length(objects[[1]]) > 1)
    m2 <- objects[[1]][[2]]
  else m2 <- NULL
  value <- anovageePrim2(m1, m2)
  return(value)
}
