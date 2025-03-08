# testing bugs
library(qgcomp)
library(qgcompint)
set.seed(23)
dat <- simdata_quantized_emm(
  outcometype="logistic",
  # sample size
  n = 100,
  # correlation between x1 and x2, x3, ...
  corr=c(0.8, 0.6, 0.3, -0.3, -0.3, -0.3),
  # model intercept
  b0=-2,
  # linear model coefficients for x1, x2, ... at referent level of interacting variable
  mainterms=c(0.3, -0.1, 0.1, 0.0, 0.3, 0.1, 0.1),
  # linear model coefficients for product terms between x1, x2, ... and interacting variable
  prodterms = c(1.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1),
  # type of interacting variable
  ztype = "categorical",
  # number of levels of exposure
  q = 4,
  # residual variance of y
  yscale = 2.0
)

dat$zspace <- paste("cat", dat$z)
dat$zspace <- ifelse(dat$z==1, "cat1", dat$zspace)

dat$zspace2 = as.factor(dat$zspace)
# testing with character
qfit1 <- qgcomp.emm.glm.noboot(y~x1+x2,
                           data = dat,
                           expnms = paste0("x", 1:2),
                           emmvar = "zspace",
                           q = 4)
# testing with factor
qfit2 <- qgcomp.emm.glm.noboot(y~x1+x2,
                           data = dat,
                           expnms = paste0("x", 1:2),
                           emmvar = "zspace2",
                           q = 4)
qfit1
qfit2

# testing tibble (not included in package dependencies, so commenting out)
#library(tibble)

#dat2 = as_tibble(dat)

#qgcomp.emm.glm.noboot(y~x1+x2,
#                  data = dat2,
#                  expnms = paste0("x", 1:2),
#                  emmvar = "zspace",
#                  q = 4)

