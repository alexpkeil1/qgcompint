#' @export
predict.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...) {
  message("Experimental feature, not validated")
  if (is.null(newdata)) {
    pred <- predict(object$fit, type=type, ...)
  }
  if (!is.null(newdata)) {
    if (is.null(expnms[1])) expnms = object$expnms # testing
    oldZ = object$fit$data[, object$call$emmvar]
    oldl = length(oldZ)
    zdata = zproc(c(oldZ, newdata[, object$call$emmvar]))
    #emmvars = names(zdata)
    newdata = cbind(newdata, zdata[-seq_len(oldl)])
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...)
  }
  return(pred)
}

predictmsm <- function(object, ...)
  UseMethod("predictmsm")

#' @export
predictmsm.qgcompemmfit <- function(object, expnms=NULL, newdata=NULL, type="response", ...) {
  message("Experimental feature, not validated")
  if (is.null(newdata)) {
    pred <- predict(object$msmfit, type=type, ...)
  }
  if (!is.null(newdata)) {
    if (is.null(expnms[1])) expnms = object$expnms # testing
    oldZ = object$fit$data[, object$call$emmvar]
    oldl = length(oldZ)
    zdata = zproc(c(oldZ, newdata[, object$call$emmvar]))
    #emmvars = names(zdata)
    newdata = cbind(newdata, zdata[-seq_len(oldl)])
    newqdata <- quantize(newdata, expnms, q=NULL, object$breaks)$data
    pred <- predict(object$fit, newdata=newqdata, type=type, ...)
  }
  return(pred)
}



qgcomp.survcurve.boot <- function(x, ...) {
  stop("Not yet implemented")
  namespaceImport("survival")
  rootdat <- as.data.frame(x$fit$x)
  psidat <- data.frame(psi=0)
  rootfun <- function(idx, df) {
    df[, x$expnms] <- idx
    df
  }
  rootfun2 <- function(idx, df) {
    df[, "psi"] <- idx
    df[, "psi1"] <- idx
    df[, "psi2"] <- idx^2
    df[, "psi3"] <- idx^3
    df[, "psi4"] <- idx^4
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
  cdfl = lapply(seq_len(x$q), function(zz) with(survival::survfit(x$fit, newdata=newcond[[zz]][1, ]), data.frame(time=time, surv=surv, q=zz)))
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

# TODO: move to qgcomp
anova.eeqgcompfit = function(object, ..., dispersion = NULL, test = NULL)
  #' @exportS3Method stats::anova
{
  # based on geepack:::anova.geeglm
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("The following arguments to anova(..) are invalid and dropped: ",
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

# TODO: move to qgcomp
anova.eeqgcompfitlist <- function(object, ..., dispersion = NULL, test = NULL)
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
  ns <- sapply(object, function(x) length(x$y.expected))
  if (any(ns != ns[1]))
    stop("models were not all fitted to the same size of dataset")
  objects <- list(object, ...)
  m1 <- objects[[1]][[1]]
  if (length(objects[[1]]) > 1)
    m2 <- objects[[1]][[2]]
  else m2 <- NULL
  value <- anovaqgcompgee(m1, m2)
  return(value)
}

# TODO: move to qgcomp
anovaqgcompgee <- function(m1, m2, ...) {
  mm1 <- model.matrix(m1)
  mm2 <- model.matrix(m2)
  P1 <- mm1 %*% solve(t(mm1) %*% mm1) %*% t(mm1)
  P2 <- mm2 %*% solve(t(mm2) %*% mm2) %*% t(mm2)
  e2 <- mm2 - P1 %*% mm2
  e1 <- mm1 - P2 %*% mm1
  m2inm1 <- all(apply(e2, 2, var) < 1e-10)
  m1inm2 <- all(apply(e1, 2, var) < 1e-10)
  if (!any(c(m2inm1, m1inm2)))
    cat("Models not nested\n")
  else if (all(c(m2inm1, m1inm2)))
    cat("Models are identical\n")
  else {
    if (m1inm2) {
      tmp <- m1
      m1 <- m2
      m2 <- tmp
    }
    mm1 <- model.matrix(m1)
    mm2 <- model.matrix(m2)
    m1emm = m1$call$emmvar
    m2emm = m2$call$emmvar
    mf1 <- paste(paste(formula(m1))[c(2, 1, 3)], collapse = " ")
    mf2 <- paste(paste(formula(m2))[c(2, 1, 3)], collapse = " ")
    if (!any(is.null(m1$expnms))) mf1 = paste0(mf1, ", expnms: ", paste(m1$expnms, collapse=", "))
    if (!any(is.null(m2$expnms))) mf2 = paste0(mf2, ", expnms: ", paste(m2$expnms, collapse=", "))
    if (!any(is.null(m1emm))) mf1 = paste0(mf1, ", EMM: ", m1emm)
    if (!any(is.null(m2emm))) mf2 = paste0(mf2, ", EMM: ", m2emm)


    mm <- cbind(mm2, mm1)
    qmm <- qr(mm)
    qmmq <- qr.Q(qmm)
    nymm1 <- as.data.frame(qmmq[, 1:qmm$rank])
    colnames(nymm1) <- paste("parm", 1:ncol(nymm1), sep = ".")
    nymm2 <- nymm1[, 1:ncol(mm2), drop = FALSE]
    formula1 <- formula(paste(formula(m1)[[2]], formula(m1)[[1]],
                              paste(c("-1", colnames(nymm1)), collapse = "+"),
                              collapse = ""))
    beta = coef(m1$fit)
    vbeta = vcov(m1$fit)
    df <- dim(mm1)[2] - dim(mm2)[2]
    rbeta <- rep(1, length(beta))
    rbeta[1:df] <- 0
    beta0 <- rev(rbeta)
    zeroidx <- beta0 == 0
    V0 <- vbeta[zeroidx, zeroidx, drop = FALSE]
    b0 <- beta[zeroidx]
    #X2 <- as.numeric(t(b0) %*% ginv(V0) %*% b0) # MASS::ginv is not in the dependencies, reverting to solve
    X2 <- as.numeric(t(b0) %*% solve(V0) %*% b0)
    ev <- eigen(V0, only.values = TRUE)$values
    df.real <- sum(ev > 1e-12)
    topnote <- paste("Model 1", mf1, "\nModel 2", mf2)
    title <- "Analysis of 'Wald statistic' Table\n"
    table <- data.frame(Df = df.real, X2 = X2, p = 1 - pchisq(X2,
                                                              df.real))
    dimnames(table) <- list("1", c("Df", "X2", "P(>|Chi|)"))
    val <- structure(table, heading = c(title, topnote),
                     class = c("anova", "data.frame"))
    return(val)
  }
}


# TODO: move to qgcomp
model.matrix.eeqgcompfit <- function(object, ...) {
  #' @exportS3Method stats::model.matrix
  object$fit$X
}

# TODO: move to qgcomp
formula.eeqgcompfit <- function(x, ...) {
  #' @exportS3Method stats::formula
  x$call$f
}

