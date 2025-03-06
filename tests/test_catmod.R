library(qgcomp)
library(qgcompint)
set.seed(23)
dat2 <- simdata_quantized_emm(
  outcometype="logistic",
  # sample size
  n = 100,
  # correlation between x1 and x2,x3,...
  corr=c(0.8,0.6,0.3,-0.3,-0.3,-0.3),
  # model intercept
  b0=-2,
  # linear model coefficients for x1,x2,... at referent level of interacting variable
  mainterms=c(0.3,-0.1,0.1,0.0,0.3,0.1,0.1),
  # linear model coefficients for product terms between x1,x2,... and interacting variable
  prodterms = c(1.0,0.0,0.0,0.0,0.2,0.2,0.2),
  # type of interacting variable
  ztype = "categorical",
  # number of levels of exposure
  q = 4,
  # residual variance of y
  yscale = 2.0
)

head(dat2)
table(dat2$z)
table(dat2$y)
dat2$z = as.factor(dat2$z)
qfit2 <- qgcomp.emm.glm.noboot(y~x1,
                           data = dat2,
                           expnms = paste0("x",1:1),
                           emmvar = "z",
                           q = 4)
qfit2
qfit2$fit
