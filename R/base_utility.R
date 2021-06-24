zproc <- function(z){
  data.frame(model.matrix(~z)[,-1,drop=FALSE])
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




