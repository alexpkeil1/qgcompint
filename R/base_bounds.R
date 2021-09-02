###############################################################################
#
# Pointwise bounds/ confidence intervals for pointwise comparisons
#
###############################################################################

.pointwise.lin <- function(q, py, se.diff, alpha, pwr){
  # mean, mean differences
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             mean.diff = py - py[pwr],
             se.diff = se.diff, # standard error on link scale
             ll.diff =  py - py[pwr] + qnorm(alpha/2) * se.diff,
             ul.diff =  py - py[pwr] + qnorm(1-alpha/2) * se.diff
  )
}



.pointwise.log <- function(q, py, se.diff, alpha, pwr){
  # risk, risk ratios / prevalence ratios
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             rr = exp(py - py[pwr]),
             se.lnrr = se.diff, # standard error on link scale
             ll.rr = exp(py - py[pwr] + qnorm(alpha/2) * se.diff),
             ul.rr = exp(py - py[pwr] + qnorm(1-alpha/2) * se.diff)
  )
}



.pointwise.logit <- function(q, py, se.diff, alpha, pwr){
  # odds, odds ratios
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log odds
             or = exp(py - py[pwr]),
             se.lnor = se.diff, # standard error on link scale
             ll.or = exp(py - py[pwr] + qnorm(alpha/2) * se.diff),
             ul.or = exp(py - py[pwr] + qnorm(1-alpha/2) * se.diff)
  )
}


.pointwise.lin.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.lin(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - pw$mean.diff + pw$ll.diff
  pw$ul.linpred = pw$hx - pw$mean.diff + pw$ul.diff
  pw
}

.pointwise.log.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.log(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - log(pw$rr) + log(pw$ll.rr)
  pw$ul.linpred = pw$hx - log(pw$rr) + log(pw$ul.rr)
  pw
}

.pointwise.logit.boot <- function(q, py, se.diff, alpha, pwr){
  pw = .pointwise.logit(q, py, se.diff, alpha, pwr)
  pw$ll.linpred = pw$hx - log(pw$or) + log(pw$ll.or)
  pw$ul.linpred = pw$hx - log(pw$or) + log(pw$ul.or)
  pw
}

.makenewdesign <- function(x, qvals, emmval=0.0, degree=1,...){
  expnms = x$expnms
  emmvar = x$call$emmvar
  zvar = x$fit$data[,emmvar]
  if(is.factor(zvar))
    zdat = zproc(zvar[zvar==emmval], znm = emmvar)
  if(!is.factor(zvar))
    zdat = zproc(zvar*0 + emmval, znm = emmvar)

  zinmodel = zdat[1,,drop=FALSE][rep(1,length(qvals)),, drop=FALSE]
  coefnm = names(coef(x))
  cfi = coefnm[which(tolower(coefnm) != "(intercept)")]
  cfi <- gsub("psi([0-9])", "I(q^\\1)", cfi)
  cfi <- gsub("mixture([\\^0-9]{0,2})", "I(q\\1)", cfi)
  cfi <- gsub("mixture", "q", cfi)
  cfi <- gsub("I\\(q\\)|I\\(q\\^1\\)", "q", cfi)
  cfi <- c(ifelse(x$hasintercept, "~1", "~-1"), cfi)
  ff <- as.formula(paste(cfi, collapse="+"))
  df <- model.frame(formula=ff,data = data.frame(q=qvals, zinmodel))
  res <- model.matrix.default(ff, data = df)
  row.names(res) <- NULL
  res
}


#' Estimating pointwise comparisons for qgcompint fits
#' @description
#' Calculates: expected outcome (on the link scale), mean difference (link scale) and the standard error of the mean difference (link scale) for pointwise comparisons
#' @details The comparison of interest following a qgcomp fit is often comparisons of model predictions at various values of the joint-exposures (e.g. expected outcome at all exposures at the 1st quartile vs. the 3rd quartile). The expected outcome at a given joint exposure and at a given level of non-exposure covariates (W) is given as E(Y|S,W=w), where S takes on integer values 0 to q-1. Thus, comparisons are of the type E(Y|S=s,W=w) - E(Y|S=s2,W=w) where s and s2 are two different values of the joint exposures (e.g. 0 and 2). This function yields E(Y|S,W=w) as well as E(Y|S=s,W=w) - E(Y|S=p,W=w) where s is any value of S and p is the value chosen via "pointwise ref" - e.g. for binomial variables this will equal the risk/ prevalence difference at all values of S, with the referent category S=p-1. For the non-boostrapped version of quantile g-computation (under a linear model)
#' Note that function only works with standard "qgcompfit" objects from qgcomp.noboot (so it doesn't work with zero inflated, hurdle, or Cox models)
#' Variance for the overall effect estimate is given by: \eqn{transpose(G) Cov(\beta) G}
#' Where the "gradient vector" G is given by
#' \deqn{G = [\partial(f(\beta))/\partial\beta_1 = 1, ..., \partial(f(β))/\partial\beta_3k= 1]}
#' \eqn{f(\beta) = \sum_i^p \beta_i, and \partial y/ \partial x} denotes the partial derivative/gradient. The vector G takes on values that equal the difference in quantiles of S for each pointwise comparison (e.g. for a comparison of the 3rd vs the 5th category, G is a vector of 2s)
#' This variance is used to create pointwise confidence intervals via a normal approximation: (e.g. upper 95% CI = psi + variance*1.96)
#'
#' @param x qgcompemmfit object from qgcomp.emm.noboot
#'
#' @param alpha alpha level for confidence intervals
#' @param pointwiseref referent quantile (e.g. 1 uses the lowest joint-exposure category as the referent category for calculating all mean differences/standard deviations)
#' @param emmval fixed value for effect measure modifier at which pointwise comparisons are calculated
#' @param ... not used
#' @return A data frame containing
#'  \describe{
#'  \item{hx: }{The "partial" linear predictor \eqn{\beta_0 + \psi\sum_j X_j^q w_j}, or the effect of the mixture + intercept after
#'  conditioning out any confounders. This is similar to the h(x) function in bkmr. This is not a full prediction of the outcome, but
#'  only the partial outcome due to the intercept and the confounders}
#'  \item{rr/or/mean.diff: }{The canonical effect measure (risk ratio/odds ratio/mean difference) for the marginal structural model link}
#'  \item{se....: }{the stndard error of the effect measure}
#'  \item{ul..../ll....: }{Confidence bounds for the effect measure}
#' }
#' @seealso \code{\link[qgcompint]{qgcomp.emm.noboot}}, \code{\link[qgcomp]{pointwisebound.noboot}}
#' @export
#' @examples
#' set.seed(50)
#' # linear model, binary modifier
#' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#' z=rbinom(50,1,0.5), r=rbinom(50,1,0.5))
#' (qfit <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#'   expnms = c('x1', 'x2'), data=dat, q=4, family=gaussian()))
#' pointwisebound(qfit, pointwiseref = 2, emmval = 0.1)
#' # linear model, categorical modifier
#' dat3 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#' z=as.factor(sample(0:2, 50,replace=TRUE)), r=rbinom(50,1,0.5))
#' (qfit3 <- qgcomp.emm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#' expnms = c('x1', 'x2'), data=dat3, q=5, family=gaussian()))
#' pointwisebound(qfit3, pointwiseref = 2, emmval = 0)
#' pointwisebound(qfit3, pointwiseref = 2, emmval = 1)
#' pointwisebound(qfit3, pointwiseref = 2, emmval = 2)
pointwisebound <- function(x, alpha = 0.05, pointwiseref = 1, emmval=0.0, ...)
  UseMethod("pointwisebound")



#' @export
pointwisebound.qgcompemmfit <- function (x, alpha = 0.05, pointwiseref = 1, emmval=0.0, ...)
{
  #if (x$bootstrap || inherits(x, "ziqgcompfit") || inherits(x,"survqgcompfit")) {
    #stop("This function is only for qgcomp.emm.noboot objects")
  #}
  link = x$fit$family$link
  qvals = c(1:x$q)-1
  designdf = as.data.frame(.makenewdesign(x, qvals, emmval=emmval))

  vc = vcov(x)
  coefnm = colnames(vc)
  cfi <- gsub("psi1|mixture", "q", coefnm)
  rownames(vc) <- colnames(vc) <- cfi
  vcord = 1:dim(vc)[1]
  designnm = names(designdf)
  newnames = designnm
  designord = numeric(length(vcord))
  for(vci in vcord){
    designord[vci] = match(cfi[vci], designnm)
    if(is.na(designord[vci])){
      subnms = gsub("([^:]+):([^:]+)", "\\2:\\1", designnm)
      designord[vci] = match(cfi[vci], subnms)
      newnames[vci] = subnms[vci]
    }
  }
  designdf = designdf[,designord]
  names(designdf) = newnames



  refrow = designdf[pointwiseref,]
  nrows = nrow(designdf)
  se.diff = numeric(nrows)
  for(nr in seq_len(nrows)){
    grad = as.numeric(designdf[nr,] - refrow)
    #whichgrad = (grad!=0)
    se.diff[nr] = se_comb2(names(grad), covmat = vc, grad)
  }


  #x$fit$data
  #####
  py = as.matrix(designdf) %*% coef(x)
  if(!x$bootstrap){
    res = switch(link,
                 identity = .pointwise.lin(x$q, py, se.diff,alpha, pointwiseref),
                 log = .pointwise.log(x$q, py, se.diff,alpha, pointwiseref),
                 logit = .pointwise.logit(x$q, py, se.diff, alpha, pointwiseref))
  }
  if(x$bootstrap){
    if(x$degree>1) stop("not implemented for non-linear fits")
    res = switch(link,
                 identity = .pointwise.lin.boot(x$q, py, se.diff,alpha, pointwiseref),
                 log = .pointwise.log.boot(x$q, py, se.diff,alpha, pointwiseref),
                 logit = .pointwise.logit.boot(x$q, py, se.diff, alpha, pointwiseref))
  }
  res$emm_level = emmval
  attr(res, "link") = link
  res
}

###############################################################################
#
# Model bounds/ confidence bounds on regression lines # ----
#
###############################################################################


.modelwise.lin <- function(q, py, se.diff, alpha, ll, ul){
  # mean, mean differences
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py,
             m = py,
             se.pw = se.diff, # standard error on link scale
             ll.pw = py + qnorm(alpha/2) * se.diff,
             ul.pw = py + qnorm(1-alpha/2) * se.diff,
             ll.simul= ifelse(is.na(ll), NA, ll),
             ul.simul=ifelse(is.na(ul), NA, ul)
  )
}

.modelwise.log <- function(q, py, se.diff, alpha, ll, ul){
  # risk
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log risk
             r = exp(py), # risk
             se.pw = se.diff, # standard error on link scale
             ll.pw = exp(py + qnorm(alpha/2) * se.diff),# risk
             ul.pw = exp(py + qnorm(1-alpha/2) * se.diff), # risk
             ll.simul= ifelse(is.na(ll), NA, exp(ll)),
             ul.simul=ifelse(is.na(ul), NA, exp(ul))
  )
}

.modelwise.logit <- function(q, py, se.diff, alpha, ll, ul){
  # odds
  data.frame(quantile= (seq_len(q)) - 1,
             quantile.midpoint=((seq_len(q)) - 1 + 0.5)/(q),
             hx = py, # log odds
             o = exp(py), # odds
             se.pw = se.diff, # standard error on link scale
             ll.pw = exp(py + qnorm(alpha/2) * se.diff),
             ul.pw = exp(py + qnorm(1-alpha/2) * se.diff), # odds
             ll.simul= ifelse(is.na(ll), NA, exp(ll)),
             ul.simul=ifelse(is.na(ul), NA, exp(ul))
  )
}


#' @title Estimating qgcomp regression line confidence bounds
#'
#' @description Calculates: expected outcome (on the link scale), and upper and lower
#'  confidence intervals (both pointwise and simultaneous)
#'
#' @details This method leverages the bootstrap distribution of qgcomp model coefficients
#' to estimate pointwise regression line confidence bounds. These are defined as the bounds
#' that, for each value of the independent variable X (here, X is the joint exposure quantiles)
#' the 95% bounds (for example) for the model estimate of the regression line E(Y|X) are expected to include the
#' true value of E(Y|X) in 95% of studies. The "simultaneous" bounds are also calculated, and the 95%
#' simultaneous bounds contain the true value of E(Y|X) for all values of X in 95% of studies. The
#' latter are more conservative and account for the multiple testing implied by the former. Pointwise
#' bounds are calculated via the standard error for the estimates of E(Y|X), while the simultaneous
#' bounds are estimated using the bootstrap method of Cheng (reference below). All bounds are large
#' sample bounds that assume normality and thus will be underconservative in small samples. These
#' bounds may also inclue illogical values (e.g. values less than 0 for a dichotomous outcome) and
#' should be interpreted cautiously in small samples.
#'
#'
#' Reference:
#'
#' Cheng, Russell CH. "Bootstrapping simultaneous confidence bands."
#' Proceedings of the Winter Simulation Conference, 2005.. IEEE, 2005.
#'
#' @param x "qgcompemmfit" object from `qgcomp.emm.boot`,
#' @param emmval fixed value for effect measure modifier at which pointwise comparisons are calculated
#' @param alpha alpha level for confidence intervals
#' @param pwonly logical: return only pointwise estimates (suppress simultaneous estimates)
#' @param ... not used
#' @return A data frame containing
#'  \describe{
#'  \item{linpred: }{The linear predictor from the marginal structural model}
#'  \item{r/o/m: }{The canonical measure (risk/odds/mean) for the marginal structural model link}
#'  \item{se....: }{the stndard error of linpred}
#'  \item{ul..../ll....: }{Confidence bounds for the effect measure, and bounds centered at the canonical measure (for plotting purposes)}
#' }
#' The confidence bounds are either  "pointwise" (pw) and "simultaneous" (simul) confidence
#' intervals at each each quantized value of all exposures.
#' @seealso \code{\link[qgcompint]{qgcomp.emm.boot}}
#' @export
#' @examples
#' set.seed(12)
#' \dontrun{
#' }
modelbound <- function(x, emmval=0.0, alpha=0.05, pwonly=FALSE, ...){
  UseMethod("modelbound")
}


#' @export
modelbound.qgcompemmfit <- function(x, emmval=0.0, alpha=0.05, pwonly=FALSE, ...){
  if(!x$bootstrap || inherits(x, "survqgcompfit")){
    stop("This function does not work with this type of qgcomp fit")
  }
  if(x$degree>1) stop("not implemented for non-linear fits")
  link = x$fit$family$link
    qvals = c(1:x$q)-1
    designmat = .makenewdesign(x, qvals, emmval=emmval)
    coef.boot = x$bootsamps
    #as.matrix(designmat) %*% coef(x$msmfit) # expected given new design matrix
    ypred <- t(designmat %*% coef.boot) # linearpredictor at each bootstrap rep, given new design matrix
    #
    #py = tapply(x$y.expectedmsm, x$index, mean)
    py <- as.numeric(designmat %*% coef(x$msmfit))
    ycovmat = cov(ypred) # covariance at specific value of emmval
    #  ycovmat = x$cov.yhat # bootstrap covariance matrix of E(y|x) from MSM
    pw.vars = diag(ycovmat)

  if(!pwonly){
    boot.err = t(coef.boot - x$coef)
    iV = solve(qr(x$covmat.coef, tol=1e-20))
    chi.boots = vapply(seq_len(nrow(boot.err)), function(i) boot.err[i,] %*% iV %*% boot.err[i,], 0.0)
    chicrit = qchisq(1-alpha, length(x$coef))
    C.set = t(coef.boot[,which(chi.boots<chicrit)])
    fx = function(coef){
      reps = designmat %*% coef
    }
    fullset = apply(C.set, 1, fx)
    ll = apply(fullset, 1, min)
    ul = apply(fullset, 1, max)
  } else{
    ll=NA
    ul=NA
  }
  res = switch(link,
               identity = .modelwise.lin(x$q, py, sqrt(pw.vars), alpha, ll, ul),
               log = .modelwise.log(x$q, log(py), sqrt(pw.vars), alpha, ll, ul),
               logit = .modelwise.logit(x$q, .logit(py), sqrt(pw.vars), alpha, ll, ul)
               #,
               #zi = .modelwise.zi(x$q, py, NULL, alpha, ll, ul, bootY)
  )
  fix = which(names(res)=="hx")
  names(res)[fix] = "linpred"
  res$emm_level = emmval
  attr(res, "link") = link
  res
}
