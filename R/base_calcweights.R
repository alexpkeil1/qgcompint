.calcstrateffects <- function(x, emmval=1.0){
  lnx = length(x$expnms)
  lnxz = length(x$intterms)
  mod = summary(x$fit)
  if( x$fit$family$family=="cox" ){
    covmat = as.matrix(x$fit$var)
    colnames(covmat) <- rownames(covmat) <- names(x$fit$coefficients)
  } else{
    covmat = as.matrix(mod$cov.scaled)
  }
  stopifnot(lnx == lnxz)
  indeffects =
    x$fit$coefficients[x$expnms] +
    x$fit$coefficients[x$intterms]*emmval
  effectatZ <- sum(indeffects)
  expidx <- which(colnames(covmat) %in% x$expnms)
  intidx <- which(colnames(covmat) %in% x$intterms)
  effgrad = 0*x$fit$coefficients
  effgrad[expidx] <- 1
  effgrad[intidx] <- emmval
  seatZ <-  se_comb2(c(x$expnms,x$intterms),
                    covmat = covmat,
                    grad = effgrad
                    )
  ciatZ <- cbind(
    effectatZ + seatZ * qnorm(x$alpha / 2),
    effectatZ + seatZ * qnorm(1 - x$alpha / 2)
    )
  res <- list(
    effectmat = rbind(
      terms.main = x$fit$coefficients[x$expnms],
      terms.prod = x$fit$coefficients[x$intterms],
      indeffects = indeffects
    )
  , # main effect + product term
    eff = effectatZ,
    se = seatZ,
    ci = ciatZ,
    emmvar = x$call$emmvar,
    emmlev = x$emmlev
  )
  res
}

#.calcstrateffects(lst)


getweightsemm <- function(x, emmval=1.0){
  #' @title Calculate weights at a set value of effect measure modifier
  #'
  #' @description A standard qgcomp fit with effect measure modification
  #' only estimates weights at the referent (0) level of the modifier.
  #' This function can be used to estimate weights at arbitrary levels of the modifier
  #'
  #'
  #' @param x "qgcompemmfit" object from qgcomp.emm.noboot
  #' function
  #' @param emmval numerical: value of effect measure modifier at which weights are generated
  #' @param ... unused
  #' @seealso \code{\link[qgcompintf]{qgcomp.emm.noboot}}
  #' @concept variance mixtures
  #' @export
  #expnms = x$expnms
  #addedintsord =  x$intterms

  #fit <- x$fit
  wcoef1 <-
    x$fit$coefficients[x$expnms] +
    x$fit$coefficients[x$intterms]*emmval

  pos.coef1 <- which(wcoef1 > 0)
  neg.coef1 <- which(wcoef1 <= 0)
  pos.size1 = sum(abs(wcoef1[pos.coef1]))
  neg.size1 = sum(abs(wcoef1[neg.coef1]))
  pos.weights1 <- abs(wcoef1[pos.coef1]) / sum(abs(wcoef1[pos.coef1]))
  neg.weights1 <- abs(wcoef1[neg.coef1]) / sum(abs(wcoef1[neg.coef1]))
  pos.psi1 <- sum(wcoef1[pos.coef1])
  neg.psi1 <- sum(wcoef1[neg.coef1])
  res <- list(
    pos.weights = sort(pos.weights1, decreasing = TRUE) ,
    neg.weights = sort(neg.weights1, decreasing = TRUE),
    pos.psi = pos.psi1,
    neg.psi = neg.psi1,
    pos.size = pos.size1,
    neg.size = neg.size1,
    emmvar = x$call$emmvar,
    emmval = emmval
  )
  class(res) <- "qgcompemmweights"
  res
}
