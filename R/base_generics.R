
print.qgcompemmweights <- function(x, ...) {
  #' @exportS3Method base::print
  cat(paste0("## Qgcomp weights/partial effects at ", x$emmvar, " = ", as.character(x$emmval), "\n"))
  #cat(paste0("## ", x$emmvar, " = ", x$emmval, "\n"))
  cat(paste0("Scaled effect size (positive direction, sum of positive effects = ", signif (x$pos.psi, 3) , ")\n"))
  if (length(x$pos.weights) > 0) {
    print(x$pos.weights, digits = 3)
  } else cat("None\n")
  cat(paste0("Scaled effect size (negative direction, sum of negative effects = ", signif (x$neg.psi, 3) , ")\n"))
  if (length(x$neg.weights) > 0) {
    print(x$neg.weights, digits = 3)
  } else cat("None\n")
  cat("\n")
}

print.qgcompemmeffects <- function(x, ..., digits = 2) {
#' @exportS3Method base::print
  cat(paste0("Joint effect at ", x$emmvar, "=", x$emmval, "\n"))
  zz = x$eff/x$se
  pval <- 2 - 2 * pnorm(abs(zz))
  pdat <- cbind(Estimate=x$eff, "Std. Error"=x$se, "Lower CI"=x$ci[1], "Upper CI"=x$ci[2], "z value"=zz, "Pr(>|z|)"=pval)
  rownames(pdat) <- "Mixture"
  printCoefmat(pdat, has.Pvalue=TRUE, tst.ind=5L, signif.stars=FALSE, cs.ind=1L:2)
}


.printeffects <- function(x, ..., digits=2) {
  emmv = x$call$emmvar
  x$alpha
  if ( x$emmlev == 2 ) {
    eff <- signif (x$effect, digits=digits)
    ci <- signif (x$cieffect, digits=digits)
    l1 <- paste0("Estimate (CI), ", emmv, "=1: \n")
    l2 <- paste0(eff, " (", ci[1], ", ", ci[2], ")")
    cat("\n");cat(l1);cat(l2);cat("\n")
  } else if ( x$emmlev > 2 ) {
    TRUE
  }
}

print.getstrateffects <- function(x, ..., digits=2) {
  #' @exportS3Method base::print
  cat("Independent effects\n")
  print(x$effectmat)

  cat("Joint effects\n")
  eff <- signif (x$eff, digits=digits)
  ci <- signif (x$ci, digits=digits)
  l1 <- paste0("Estimate (CI) Std. Error, ", x$emmv, "=", x$emmlev, ": \n")
  l2 <- paste0(eff, " (", ci[1], ", ", ci[2], "), x$se")
  cat("\n");cat(l1);cat(l2);cat("\n")
}


print.qgcompemmfit <- function(x, showweights=TRUE, ...) {
  #' @title Default printing method for a qgcompemmfit object
  #'
  #' @description Prints output depending the model used. `qgcomp.emm.glm.noboot` and `qgcomp.emm.cox.noboot` will output final estimate of joint exposure
  #' effect, as well
  #' as estimates of the 'weights' (scaled coefficients). `qgcomp.emm.glm.boot` and `qgcomp.emm.glm.ee` methods will only output final effect estimates.
  #'
  #' @param x "qgcompemmfit" object from `qgcomp.emm.glm.noboot`
  #' function
  #' @param showweights logical: should weights be printed, if estimated?
  #' @param ... unused
  #' @return Invisibly returns x. Called primarily for side effects.
  #' @seealso \code{\link[qgcompint]{qgcomp.emm.glm.noboot}}, \code{\link[qgcompint]{getstratweights}}
  #' @concept variance mixtures
  #' @exportS3Method base::print
  emmvar <- x$call$emmvar
  isboot <- x$bootstrap
  isee <- inherits(x, "eeqgcompfit")
  isbinemm <- x$emmlev == 2
  #rnm = c("(Intercept)", 'psi1', emmvar, paste0(emmvar, ":mixture"))

  cilabel = ifelse(isboot, " (bootstrap CI)", ifelse(isee, " (robust CI)", " (delta method CI)"))

  rnm = names(x$coef)
  fam <- x$fit$family$family
  if (showweights && !isboot && !isee) {
    if (inherits(x$emmvals[1], "factor")) {
      ww = getstratweights(x, emmval=levels(x$emmvals)[1])
    } else{
      ww = getstratweights(x, emmval=0.0)
    }
    print(ww)
    cat("\n")
  }
  if (!is.null(x$pos.size1) && showweights && isbinemm && !isboot && !isee) {
    ww = getstratweights(x, emmval=levels(x$emmvals)[1])
    print(ww)
    cat("\n")
  }
  if (fam == "binomial") {
    estimand <- 'OR'
    if (x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("## Mixture log(", estimand, ")", cilabel, ":\n\n"))
    testtype = "Z"
  }
  if (fam == "poisson") {
    estimand <- 'RR'
    cat(paste0("## Mixture log(", estimand, ")", cilabel, ":\n\n"))
    testtype = "Z"
  }
  if (fam == "gaussian") {
    cat(paste0("## Mixture slope parameters", cilabel, ":\n\n"))
    testtype = "t"
    x$zstat = x$tstat
  }
  if (fam == "cox") {
    cat(paste0("Mixture log(hazard ratio)", cilabel, ":\n\n"))
    testtype = "Z"
    rnm = rnm#[-1]
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if (is.null(dim(x$ci.coef))) {
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "pval"=x$pval)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[, 1], "Upper CI"=x$ci.coef[, 2], "test"=x$zstat, "pval"=x$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  printCoefmat(pdat, has.Pvalue=TRUE, tst.ind=5L, signif.stars=FALSE, cs.ind=1L:2)
  if ( isbinemm  && !isboot) .printeffects(x, digits=5)
  invisible(x)
}


##### generics ######

vcov.qgcompemmfit <- function(object, ...) {
  #' @exportS3Method stats::vcov
  object$covmat.coef
}

confint.qgcompemmfit <- function(object, ...) {
  #' @exportS3Method stats::confint
  message("not yet implemented")
  anova(object$fit)
}

coef.qgcompemmfit <- function(object, ...) {
  #' @exportS3Method stats::coef
  object$coef
}

formula.qgcompemmfit <- function(x, ...) {
  #' @exportS3Method stats::formula
  x$call$f
}

model.matrix.qgcompemmfit <- function(object, ...) {
  #' @exportS3Method stats::model.matrix
  object$fit$X
}


##########################################
# section below will be moved to qgcomp
##########################################

anova.eeqgcompfit = function(object, ..., dispersion = NULL, test = NULL)
{
  #' @exportS3Method stats::anova
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


anova.eeqgcompfitlist <- function(object, ..., dispersion = NULL, test = NULL)
{
  #' @exportS3Method stats::anova
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

model.matrix.eeqgcompfit <- function(object, ...) {
  #' @exportS3Method stats::model.matrix
  object$fit$X
}

formula.eeqgcompfit <- function(x, ...) {
  #' @exportS3Method stats::formula
  x$call$f
}
