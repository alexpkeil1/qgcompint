.expit <- function(mu) 1/(1+exp(-mu))

# design matrix maker for a quantized version of exposures and a single modifier
.quantized_design_emm <- function(
  N = 100,                  # sample size
  b0=0,                     # baseline expected outcome (model intercept)
  mainterms=c(1,0,0,0),     # beta coefficients for X in the outcome model at referent (0) level of Z
  prodterms = c(1,0,0,0),   # product term coefficients for effect measure modifier (Z)
  ztype = "binary",         # continuous, binary, categorical (or 2+ letter abbreviations)
  ncor=0,                   # Number of correlated exposures
  corr=0.75                 # Pearson/spearman (the same here) correlation
){
  p = length(mainterms)
  if(ncor >= p) ncor = p-1
  X = matrix(nrow=N, ncol=p)
  ztype = substr(ztype, 1,2)
  Z = switch(ztype,
             co = rnorm(N, 0, 1),
             ca = sample(1:4, size=N, replace=TRUE),
             bi = rbinom(N, 1, 0.5)
  )
  xmaster = sample(rep(0:3, length.out=N), N, replace=FALSE)
  for(k in 1:p){
    newx = numeric(N)
    c1 = as.logical(rbinom(N, 1, sqrt(corr)))
    newx[which(c1)] = xmaster[which(c1)]
    newx[which(!c1)] = sample(xmaster[which(!c1)])
    if(k<=(ncor+1)){
      X[,k] = newx
    } else X[,k] = sample(xmaster)
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
  ncor=0,                   # Number of correlated exposures
  corr=0.75,                # Pearson/spearman (the same here) correlation
  yscale = 1.0              # standard deviation of error term in outcome
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N,b0,mainterms,prodterms,ztype,ncor,corr)
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
  ncor=0,                   # Number of correlated exposures
  corr=0.75                 # Pearson/spearman (the same here) correlation
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N,b0,mainterms,prodterms,ztype,ncor,corr)
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
  ncor=0,                   # Number of correlated exposures
  corr=0.75                 # Pearson/spearman (the same here) correlation
){
  #'
  # simulate under data structure where WQS/qgcomp is the truth:
  #  e.g. a multivariate exposure with multinomial distribution
  #  and an outcome that is a linear function of exposure scores
  lst = .quantized_design_emm(N,b0,mainterms,prodterms,ztype,ncor,corr)
  #py <- .expit(lst$mu)
  shape0 = 3
  scale0 = 5
  censtime = 4.0
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
