getstrateffects <- function(x, emmval=1.0, ...){
#' @title Calculate mixture effect at a set value of effect measure modifier
#'
#' @description A standard qgcomp fit with effect measure modification
#' only estimates effects at the referent (0) level of the modifier (psi1).
#' This function can be used to estimate effects at arbitrary levels of the modifier
#'
#'
#' @param x "qgcompemmfit" object from qgcomp.emm.noboot
#' function
#' @param emmval numerical: value of effect measure modifier at which weights are generated
#' @param ... unused
#' @seealso \code{\link[qgcompint]{qgcomp.emm.noboot}} \code{\link[qgcompint]{getstratweights}}
#' @concept variance mixtures
#' @export
#' @examples
#' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#'   z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
#' (qfit <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#'   expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
#' getstrateffects(qfit, emmval = 0)
#' getstrateffects(qfit, emmval = 1)

#expnms = x$expnms
#addedintsord =  x$intterms  zvar = x$fit$data[,x$call$emmvar]
  zvar = x$fit$data[,x$call$emmvar]
  res = .calcstrateffects(x,emmval=emmval, zvar=zvar)
  class(res) <- "qgcompemmeffects"
  res
}


.calcstrateffects <- function(x, emmval=1.0, zvar){
  #x$call$emmvar
  whichintterms = x$intterms
  if(is.factor(zvar)){
    whichlevels = zproc(zvar[which(zvar==emmval)][1])
    whichvar = names(whichlevels)[which(whichlevels==1)]
    whichintterms = NULL
    if(length(whichvar)>0) whichintterms = grep(whichvar, x$intterms, value = TRUE)
  }

  #lnx = length(x$expnms)
  #lnxz = length(whichintterms)
  mod = summary(x$fit)
  if( x$fit$family$family=="cox" ){
    covmat = as.matrix(x$fit$var)
    colnames(covmat) <- rownames(covmat) <- names(x$fit$coefficients)
  } else{
    covmat = as.matrix(mod$cov.scaled)
  }
  #stopifnot(lnx == lnxz)
  if(is.factor(zvar)){
    indeffects =
      x$fit$coefficients[x$expnms]
    if(!is.null(whichintterms)){
      indeffects =
        indeffects +
        x$fit$coefficients[whichintterms]
    }
  } else{
    indeffects =
      x$fit$coefficients[x$expnms] +
      x$fit$coefficients[x$intterms]*emmval
  }
  effectatZ <- sum(indeffects)
  expidx <- which(colnames(covmat) %in% x$expnms)
  intidx <- which(colnames(covmat) %in% whichintterms)
  effgrad = 0*x$fit$coefficients
  effgrad[expidx] <- 1
  if(is.factor(zvar)){
    effgrad[intidx] <- 1.0
  } else effgrad[intidx] <- emmval
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


getstratweights <- function(x, emmval=1.0, ...){
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
  #' @seealso \code{\link[qgcompint]{qgcomp.emm.noboot}}
  #' @concept variance mixtures
  #' @export
  #' @examples
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
  #'   z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
  #' (qfit <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",
  #'   expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  #' getstratweights(qfit, emmval = 0)
  #' getstratweights(qfit, emmval = 1)

  #fit <- x$fit
  zvar = x$fit$data[,x$call$emmvar]
  if(!is.factor(zvar)){
    wcoef1 <-
      x$fit$coefficients[x$expnms] +
      x$fit$coefficients[x$intterms]*emmval
  }
  if(is.factor(zvar)){
    whichlevels = zproc(zvar[which(zvar==emmval)][1])
    whichvar = names(whichlevels)[which(whichlevels==1)]
    whichintterms = NULL
    if(length(whichvar)>0) whichintterms = grep(whichvar, x$intterms, value = TRUE)
    # select only the involved indicator terms
    wcoef1 <-
      x$fit$coefficients[x$expnms]
    if(!is.null(whichintterms))
      wcoef1 <-
      wcoef1 +
      x$fit$coefficients[whichintterms]
  }

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
