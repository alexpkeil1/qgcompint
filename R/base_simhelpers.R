.expit <- function(mu) 1/(1+exp(-mu))

#' @title Simulate quantized exposures for testing methods
#'
#' @description Simulate quantized exposures for testing methods
#'
#'
#' @details Simulate continuous (normally distributed errors), binary (logistic function), or event-time outcomes as a linear function
#'
#' @param outcometype Character variable that is one of c("continuous", "logistic", "survival"). Selects what type of outcome should be simulated (or how). continuous = normal, continous outcome, logistic= binary outcome from logistic model, survival = right censored survival outcome from Weibull model.
#' @param n Sample size
#' @param corr NULL, or vector of correlations between the first exposure and subsequent exposures (if length(corr) < (length(coef)-1), then this will be backfilled with zeros)
#' @param b0  (continuous, binary outcomes) model intercept
#' @param mainterms beta coefficients for X in the outcome model at referent (0) level of interacting variable
#' @param prodterms product term coefficients for interacting variable
#' @param ztype type of interacting variable: "continuous", "binary", "categorical"
#' @param q Number of levels or "quanta" of each exposure
#' @param yscale (continuous outcomes) error scale (residual error) for normally distributed outcomes
#' @param shape0 (survival outcomes) baseline shape of weibull distribution \link[stats]{rweibull}
#' @param scale0 (survival outcomes) baseline scale of weibull distribution \link[stats]{rweibull}
#' @param censtime (survival outcomes) administrative censoring time
#' @param ncheck (logical, default=TRUE) adjust sample size if needed so that exposures are exactly evenly distributed (so that qgcomp::quantize(exposure) = exposure)
#' @param ... unused
#'
#' @seealso \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp.noboot}}
#' @return a data frame
#' @export
#' @examples
#' set.seed(50)
#' qdat = simdata_quantized_emm(
#'   outcomtype="continuous",
#'   n=10000, corr=c(.9,.3,0,0), mainterms=c(1,1,0,0), prodterms=c(1,1,0,0),
#'   q = 8
#' )
#' cor(qdat)
#' qdat = simdata_quantized_emm(
#'   outcomtype="continuous",
#'   n=10000, corr=c(-.9,.3,0,0), mainterms=c(1,2,0,0), prodterms=c(1,1,0,0),
#'   q = 4
#' )
#' cor(qdat)
#' table(qdat$x1)
#' qgcomp.emm.glm.noboot(y~x1+x2+x3+x4,expnms = c("x1", "x2", "x3", "x4"), emmvar = "z", data=qdat)
simdata_quantized_emm <- function(
  outcometype=c("continuous", "logistic", "survival"),
  n = 100,
  corr=NULL,
  b0=0,
  mainterms=c(1,0,0,0),     #
  prodterms = c(1,0,0,0),   #
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  q = 4,
  yscale = 1.0,
  shape0 = 3.0,
  scale0 = 5.0,
  censtime = 4.0,
  ncheck=TRUE,
  ...
){
  stopifnot(length(prodterms)==length(mainterms))
  rem = n %% q
  if(rem != 0 & ncheck){
    n = n+rem
    message("Adjusting sample size to have evenly distributed exposure values")
  }
  otype=substr(outcometype[1], 1, 5)
  dat <- switch(otype,
                conti = .dgm_quantized_linear_emm(N = n,b0=b0,mainterms=mainterms, prodterms=prodterms,corr=corr,yscale=yscale,q=q,ztype=ztype),
                logis = .dgm_quantized_logistic_emm(N = n,b0=b0,mainterms=mainterms, prodterms=prodterms,corr=corr,q=q,ztype=ztype),
                survi = .dgm_quantized_survival_emm(N = n,b0=0, mainterms=mainterms, prodterms=prodterms,corr=corr,shape0=shape0,scale0=shape0,censtime=censtime,q=q,ztype=ztype)
  )
  dat
}


# design matrix maker for a quantized version of exposures and a single modifier
.quantized_design_emm <- function(
  N = 100,                  # sample size
  b0=0,                     # baseline expected outcome (model intercept)
  mainterms=c(1,0,0,0),     # beta coefficients for X in the outcome model at referent (0) level of Z
  prodterms = c(1,0,0,0),   # product term coefficients for effect measure modifier (Z)
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  corr=0.75,                 # Pearson/spearman (the same here) correlation
  q=4
){
  ztype = substr(ztype, 1,2)
  Z = switch(ztype,
             co = rnorm(N, 0, 1),
             ca = sample(0:2, size=N, replace=TRUE),
             bi = rbinom(N, 1, 0.5)
  )
  p = length(mainterms)
  #if(ncor >= p) ncor = p-1
  corv = rep(0, length(coef))
  corv[1] = 1.0
  if(is.null(corr)) corr = corv[-1]
  if(length(corr)>1) corv[1+(1:length(corr))] = corr
  X = matrix(nrow=N, ncol=p)
  xtemplate = sample(rep(0:(q-1), length.out=N), N, replace=FALSE)
  for(k in 1:p){
    newx = numeric(N)
    #c1 = as.logical(rbinom(N, 1, sqrt(corv[k])))
    if(corv[k]>=0){
      c1 = as.logical(rbinom(N, 1, corv[k]))
      newx[which(c1)] = xtemplate[which(c1)]
      newx[which(!c1)] = sample(xtemplate[which(!c1)])
    }
    if(corv[k]<0){
      c1 = as.logical(rbinom(N, 1, -corv[k]))
      newx[which(c1)] = q-1-xtemplate[which(c1)]
      newx[which(!c1)] = sample(q-1-xtemplate[which(!c1)])
    }
    X[,k] = newx
  }
  colnames(X) <- paste0("x", 1:p)
  Design = cbind(X,X * Z)
  mu <- b0 + Design %*% c(mainterms,prodterms)
  list(mu=mu,X=X,Z=Z)
}

.dgm_quantized_linear_emm <- function(
  N = 100,                  # sample size
  b0=0,                     # baseline expected outcome (model intercept)
  mainterms=c(1,0,0,0),     # beta coefficients for X in the outcome model at referent (0) level of Z
  prodterms = c(1,0,0,0),   # product term coefficients for effect measure modifier (Z)
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  corr=0.75,                # Pearson/spearman (the same here) correlation
  yscale = 1.0,              # standard deviation of error term in outcome
  q=4
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N=N,b0=b0,mainterms=mainterms,prodterms=prodterms,ztype=ztype,corr=corr,q=q)
  y = rnorm(N,0,yscale) + lst$mu
  res = data.frame(z=lst$Z,lst$X,y)
  attr(res, "truecoefs") = list(intercept=b0,mainterms=mainterms,prodterms=prodterms)
  res
}

.dgm_quantized_logistic_emm <- function(
  N = 100,                  # sample size
  b0=0,                     # baseline expected outcome (model intercept)
  mainterms=c(1,0,0,0),     # beta coefficients for X in the outcome model at referent (0) level of Z
  prodterms = c(1,0,0,0),   # product term coefficients for effect measure modifier (Z)
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  corr=0.75,                 # Pearson/spearman (the same here) correlation
  q=4
  ){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N=N,b0=b0,mainterms=mainterms,prodterms=prodterms,ztype=ztype,corr=corr,q=q)
  py <- .expit(lst$mu)
  pextreme = mean(py>.995) + mean(py<0.005)
  if(pextreme > .10) warning("model implies > 10% of observations with very high/low (<0.5%) outcome probability, which may distort estimates")
  y = rbinom(N, 1, py)
  res = data.frame(z=lst$Z,lst$X,y)
  attr(res, "truecoefs") = list(intercept=b0,mainterms=mainterms,prodterms=prodterms)
  res
}


.dgm_quantized_survival_emm <- function(
  N = 100,                  # sample size
  b0=0,                     # baseline expected outcome (model intercept)
  mainterms=c(1,0,0,0),     # beta coefficients for X in the outcome model at referent (0) level of Z
  prodterms = c(1,0,0,0),   # product term coefficients for effect measure modifier (Z)
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  corr=0.75,                 # Pearson/spearman (the same here) correlation
  shape0=shape0,
  scale0=scale0,
  censtime=censtime,
  q=4
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N=N,b0=b0,mainterms=mainterms,prodterms=prodterms,ztype=ztype,corr=corr,q=q)
  #py <- .expit(lst$mu)
  t0 <- rweibull(N, shape = shape0, scale = scale0)
  #t1 <- exp(log(t0) - log(HR)/(shape0))
  tmg <- pmin(censtime, exp(log(t0) -lst$mu/(shape0)))

  #tmg <- pmin(.1,rweibull(N, 10, 0.1))
  d=1.0*(tmg<censtime)
  if(mean(d) < .10) warning("model implies > 10% of observations experience the event of interest")
  res = data.frame(z=lst$Z,lst$X,time=tmg,d=d)
  attr(res, "truecoefs") = list(intercept=b0,mainterms=mainterms,prodterms=prodterms)
  res
}
