#   (f <- .intmaker(f,expnms,emmvars)) # create necessary interaction terms with exposure
# y ~ z + x1 + x2 + x3 + x4
.intmaker <- function(
    f,
    expnms,
    emmvars,
    emmvar
){
  rightside = as.character(f)[3]
  trms = strsplit(gsub(" ", "", rightside), "+",fixed=TRUE)[[1]]
  # drop emmvar so it isn't duplicated when it's a factor
  trms = setdiff(trms, emmvar)
  rightside = paste0(trms, collapse="+")
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


.logit <- function(p) log(p) - log(1-p)
.expit <- function(mu) 1/(1+exp(-mu))
.safelog <- function(x,eps=1e-10) ifelse(x==0, log(x+eps), log(x))

zproc <- function(z, znm="z"){
  znames = ifelse(is.null(names(z)), znm, names(z))
  zres = data.frame(model.matrix(~z)[,-1,drop=FALSE])
  names(zres) <- gsub("z", znames[1], names(zres))
  zres
}


vc_comb <- function (aname = "(Intercept)", expnms, covmat, grad = NULL) {
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  if (!is.null(grad))
    grad = NULL
  if (is.null(grad))
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  outcov = matrix(NA, nrow = 2, ncol = 2)
  acol = which(colnames(as.matrix(covmat)) %in% aname)
  bcol = which(colnames(as.matrix(covmat)) %in% expnms)
  outcov[1, 1] <- covmat[acol, acol]
  outcov[1, 2] <- outcov[2, 1] <- sum(covmat[acol, bcol])
  outcov[2, 2] <- weightvec %*% covmat %*% weightvec
  outcov
}

vc_multiscomb <- function (
  inames = c("(Intercept)"),
  emmvars,
  expnms,
  addedintsl,
  covmat,
  grad = NULL
) {
  #  construct new covariance matrix as linear combination (e.g. with all grad=1 we have)
  #       x1 x2  z 1z 2z
  # x1  | 11 12 13 14 15 |        x1+x2         z    z*(x1+x2)
  # x2  | 21 22 23 24 25 |     | 11+12+21+22  13+23 14+15+24+25 |
  #  z  | 31 32 33 34 35 |  -> | 31+32        33    34+35       |
  # 1z  | 41 42 43 44 45 |     | 41+42+51+52  43+53 44+45+54+55 |
  # 2z  | 51 52 53 54 55 |

  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  expidx = ifelse(is.null(inames), 1, 2)
  # eventual dimension
  dimnew <- expidx + length(emmvars) + length(addedintsl)
  dimold <- dim(covmat)[1]
  # initialize "weight" vector
  if (!is.null(grad[1])) # will fix later to allow non-null gradients
    grad = NULL
  if (is.null(grad[1]))
    grad <- 1
  # order of variables
  nms = list(expnms)
  if(!is.null(inames))
    nms = c(inames, nms)
  for(i in seq_len(length(emmvars))){
    nms = c(nms, emmvars[i])
    nms = c(nms, addedintsl[i])
  }
  weightvec <- list()
  for(j in seq_len(dimnew)){
    weightvec[[j]] = rep(0, dimold)
    vars = nms[[j]]
    if(j == expidx){
      weightvec[[j]][which(colnames(covmat) %in% vars)] <- grad
    } else{
      weightvec[[j]][which(colnames(covmat) %in% vars)] <- 1
    }
  }
  outcov = matrix(NA, nrow = dimnew, ncol = dimnew)
  for(jj in seq_len(dimnew)){
    for(ii in jj:dimnew){
      outcov[jj,ii] <- outcov[ii,jj] <- weightvec[[jj]] %*% covmat %*% weightvec[[ii]]
    }
  }
  outcov
}


se_comb2 <- function (expnms, covmat, grad = NULL) {
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  if (is.null(grad)){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- 1
  } else if (!is.null(grad) && length(grad)==1){
    weightvec[which(colnames(as.matrix(covmat)) %in% expnms)] <- grad
  } else if (!is.null(grad) && length(grad)==length(weightvec)){
    weightvec <- grad
  }
  var <- weightvec %*% covmat %*% weightvec
  sqrt(var)[1, , drop = TRUE]
}





.qgcompemm_object <- function(...){
  res = list(...)
  nms = names(res)
  if(is.na(match("hasintercept", nms))) res$hasintercept = TRUE
  if(is.na(match("bootstrap", nms))) res$bootstrap=FALSE
  if(is.na(match("cov.yhat", nms))) res$cov.yhat=NULL
  if(is.na(match("degree", nms))) res$degree=1
  if(is.na(match("pos.psi", nms))) res$pos.psi = NULL
  if(is.na(match("neg.psi", nms))) res$neg.psi = NULL
  if(is.na(match("pos.weights", nms))) res$pos.weights = NULL
  if(is.na(match("neg.weights", nms))) res$neg.weights = NULL
  if(is.na(match("pos.size", nms))) res$pos.size = NULL
  if(is.na(match("neg.size", nms))) res$neg.size = NULL
  if(is.na(match("df", nms))) res$df = NULL
  if(is.na(match("covmat.all_robust", nms))) res$covmat.all_robust = NULL
  attr(res, "class") <- c("qgcompemmfit", "qgcompfit", "list")
  res
}

.qgc.require <- function (package, message = paste("loading required package (",
                                                   package, ") failed", sep = "")){
  if (!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

.devinstall <- function (...)
{
  .qgc.require("devtools")
  devtools::install_github("alexpkeil1/qgcompint", ...)
}


getjointeffects <- function(x, emmval=1.0, ...){
  #' @title Calculate joint effect of mixture effect and modifier vs. common referent
  #'
  #' @description A standard qgcomp fit with effect measure modification
  #' only estimates effects at the referent (0) level of the modifier (psi1).
  #' This function can be used to estimate a "common referent" parameter that
  #' estimates the effect of being in a non-referent category of the modifier and
  #' increasing exposure by one quantile, relative to no change in exposure in the
  #' referent category of the modifier. This is generally useful for binary exposures
  #' (for a mixture with a set of binary exposures,
  #' this would be the "effect" of being exposed and at the index level of the mediator,
  #' relative to being unexposed in the referent level of the mediator), but it may also
  #' be of interest with more general exposures.
  #'
  #'
  #' @param x "qgcompemmfit" object from qgcomp.emm.glm.noboot
  #' function
  #' @param emmval numerical: value of effect measure modifier at which weights are generated
  #' @param ... unused
  #' @seealso \code{\link[qgcompint]{qgcomp.emm.glm.noboot}} \code{\link[qgcompint]{getstrateffects}}
  #' @concept variance mixtures
  #' @return
  #' An object of class "qgcompemmeffects", which inherits from "qgcompemmfit" and "list"
  #'
  #' This class contains the `emmval`-stratum specific effect estimates of the mixture. By default, this prints a coefficient table, similar to objects of type "qgcompemmfit" which displays the stratum specific joint effects from a "qgcompemmfit" model.
  #'
  #' @export
  #' @examples
  #' library(qgcompint)
  #' n = 500
  #' dat <- data.frame(y=rbinom(n,1,0.5), cd=runif(n), pb=runif(n),
  #'                   raceth=factor(sample(c("WNH", "BNH", "AMIND"), n, replace=TRUE),
  #'                           levels = c("BNH", "WNH", "AMIND")))
  #' (qfit <- qgcomp.emm.glm.noboot(f=y ~cd + pb, emmvar="raceth",
  #'                            expnms = c('cd', 'pb'), data=dat, q=4,
  #'                            family=binomial()))
  #'
  #'
  #' # first level of the stratifying variable should be the referent category,
  #' #  which you can set with the "levels" argument to "factor" when
  #' #  cleaning/generating data
  #' levels(dat$raceth)
  #'
  #' # stratum specific mixture log-odds ratios
  #' # this one comes straight from the model (psi 1)
  #' getjointeffects(qfit, emmval = "BNH")
  #' # this will coincide with joint effects, since it is in the referent category
  #' getstrateffects(qfit, emmval = "BNH")
  #'
  #' # the stratum specific effect for a non-referent category of the EMM
  #' #  will not coincide with the joint effect
  #' getjointeffects(qfit, emmval = "AMIND")
  #' getstrateffects(qfit, emmval = "AMIND")
  #'

  #expnms = x$expnms
  #addedintsord =  x$intterms  zvar = x$fit$data[,x$call$emmvar]
  #if(x$bootstrap) stop("This method does not work for bootstrapped fits. If using a linear parameterization, then stratified effects can be estimated using non-bootstrapped methods.")
  if(x$degree>1) stop("not implemented for non-linear fits")
  zvar = x$fit$data[,x$call$emmvar]
  res = .calcjointffects(x,emmval=emmval, zvar=zvar)
  class(res) <- "qgcompemmeffects"
  res
}



.calcjointffects <- function(x, emmval=1.0, zvar){
  #x$call$emmvar
  whichintterms = x$intterms
  if(is.factor(zvar)){
    whichlevels = zproc(zvar[which(zvar==emmval)][1], znm = x$call$emmvar)
    whichvar = names(whichlevels)[which(whichlevels==1)]
    whichmainterms = whichvar
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
  if(length(whichmainterms)>1)
    stop("getjointeffects: length(whichmainterms)>1, which generally means something is wrong in code")
  maineffects = 0
  if(length(whichmainterms)==1)
    maineffects = x$fit$coefficients[whichmainterms] # this
  effectatZ <- sum(indeffects)  + maineffects
  expidx <- which(colnames(covmat) %in% x$expnms)
  mainidx <- which(colnames(covmat) %in% whichmainterms)
  intidx <- which(colnames(covmat) %in% whichintterms)
  effgrad = 0*x$fit$coefficients
  effgrad[expidx] <- 1
  effgrad[mainidx] <- 1
  if(is.factor(zvar)){
    effgrad[intidx] <- 1.0
  } else effgrad[intidx] <- emmval
  seatZ <-  se_comb2(c(x$expnms, whichmainterms, x$intterms),
                     covmat = covmat,
                     grad = effgrad
  )
  ciatZ <- cbind(
    effectatZ + seatZ * qnorm(x$alpha / 2),
    effectatZ + seatZ * qnorm(1 - x$alpha / 2)
  )
  res <- list(
    effectmat = rbind(
      terms.emm = x$fit$coefficients[whichmainterms],
      terms.main = x$fit$coefficients[c(x$expnms)],
      terms.prod = x$fit$coefficients[x$intterms],
      indeffects = indeffects
    )
    , # main effect + product term
    eff = effectatZ,
    se = seatZ,
    ci = ciatZ,
    emmvar = x$call$emmvar,
    emmlev = x$emmlev,
    emmval = emmval
  )
  res
}


getstrateffects <- function(x, emmval=1.0, ...){
  #' @title Calculate mixture effect at a set value of effect measure modifier
  #'
  #' @description A standard qgcomp fit with effect measure modification
  #' only estimates effects at the referent (0) level of the modifier (psi1).
  #' This function can be used to estimate effects at arbitrary levels of the modifier
  #'
  #'
  #' @param x "qgcompemmfit" object from qgcomp.emm.glm.noboot
  #' function
  #' @param emmval numerical: value of effect measure modifier at which weights are generated
  #' @param ... unused
  #' @seealso \code{\link[qgcompint]{qgcomp.emm.glm.noboot}} \code{\link[qgcompint]{getstratweights}}
  #' @concept variance mixtures
  #' @return
  #' An object of class "qgcompemmeffects", which inherits from "qgcompemmfit" and "list"
  #'
  #' This class contains the `emmval`-stratum specific effect estimates of the mixture. By default, this prints a coefficient table, similar to objects of type "qgcompemmfit" which displays the stratum specific joint effects from a "qgcompemmfit" model.
  #'
  #' @export
  #' @examples
  #' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
  #'   z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
  #' (qfit <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
  #'   expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
  #' getstrateffects(qfit, emmval = 0)
  #' strateffects = getstrateffects(qfit, emmval = 1)
  #'

  #expnms = x$expnms
  #addedintsord =  x$intterms  zvar = x$fit$data[,x$call$emmvar]
  #if(x$bootstrap) stop("This method does not work for bootstrapped fits. If using a linear parameterization, then stratified effects can be estimated using non-bootstrapped methods.")
  if(x$degree>1) stop("not implemented for non-linear fits")
  zvar = x$fit$data[,x$call$emmvar]
  res = .calcstrateffects(x,emmval=emmval, zvar=zvar)
  class(res) <- "qgcompemmeffects"
  res
}


.calcstrateffects <- function(x, emmval=1.0, zvar){
  #x$call$emmvar
  whichintterms = x$intterms
  if(is.factor(zvar)){
    whichlevels = zproc(zvar[which(zvar==emmval)][1], znm = x$call$emmvar)
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
  } else if(any(class(x$fit)=="eefit")){
    covmat = x$fit$vcov
  }
  else{
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
      coef(x$fit)[x$expnms] +
      coef(x$fit)[x$intterms]*emmval
  }
  effectatZ <- sum(indeffects)
  expidx <- which(colnames(covmat) %in% x$expnms)
  intidx <- which(colnames(covmat) %in% whichintterms)
  effgrad = 0*coef(x$fit)
  effgrad[expidx] <- 1
  if(is.factor(zvar)){
    effgrad[intidx] <- 1.0
  } else effgrad[intidx] <- emmval
  seatZ <- se_comb2(c(x$expnms,x$intterms),
                    covmat = covmat,
                    grad = effgrad
  )
  ciatZ <- cbind(
    effectatZ + seatZ * qnorm(x$alpha / 2),
    effectatZ + seatZ * qnorm(1 - x$alpha / 2)
  )
  res <- list(
    effectmat = rbind(
      terms.main = coef(x$fit)[x$expnms],
      terms.prod = coef(x$fit)[x$intterms],
      indeffects = indeffects
    )
    , # main effect + product term
    eff = effectatZ,
    se = seatZ,
    ci = ciatZ,
    emmvar = x$call$emmvar,
    emmlev = x$emmlev,
    emmval = emmval
  )
  res
}

#.calcstrateffects(lst)

