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
  attr(res, "class") <- c("qgcompemmfit", "qgcompfit", "list")
  res
}

