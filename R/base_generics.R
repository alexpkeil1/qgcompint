
print.qgcompemmweights <- function(x, ...){
  #' @export
  cat(paste0("## Qgcomp weights/partial effects at ", x$emmvar, " = ", as.character(x$emmval), "\n"))
  #cat(paste0("## ", x$emmvar, " = ",x$emmval, "\n"))
  cat(paste0("Scaled effect size (positive direction, sum of positive effects = ", signif(x$pos.psi, 3) , ")\n"))
  if (length(x$pos.weights) > 0) {
    print(x$pos.weights, digits = 3)
  } else cat("None\n")
  cat(paste0("Scaled effect size (negative direction, sum of negative effects = ", signif(x$neg.psi, 3) , ")\n"))
  if (length(x$neg.weights) > 0) {
    print(x$neg.weights, digits = 3)
  } else cat("None\n")
  cat("\n")
}



.printeffects <- function(x, digits=2){
  emmv = x$call$emmvar
  x$alpha
  if( x$emmlev == 2 ){
    eff <- signif(x$effect, digits=digits)
    ci <- signif(x$cieffect, digits=digits)
    l1 <- paste0("Estimate (CI), ", emmv, "=1: \n")
    l2 <- paste0(eff, " (", ci[1], ", ", ci[2], ")")
    cat("\n");cat(l1);cat(l2);cat("\n")
  } else if( x$emmlev > 2 ){
    TRUE
  }
}

print.qgcompemmfit <- function(x, showweights=TRUE, ...){
  #' @title Default printing method for a qgcompemmfit object
  #'
  #' @description Prints output depending for `qgcomp.emm.noboot` will output final estimate of joint exposure
  #' effect (similar to the 'index' effect in weighted quantile sums), as well
  #' as estimates of the 'weights' (standardized coefficients).
  #'
  #' @param x "qgcompemmfit" object from `qgcomp.emm.noboot`
  #' function
  #' @param showweights logical: should weights be printed, if estimated?
  #' @param ... unused
  #' @seealso \code{\link[qgcomp]{qgcomp.emm.noboot}}, \code{\link[qgcomp]{getweightsemm}}
  #' @concept variance mixtures
  #' @export
  emmvar <- x$call$emmvar
  isbinemm <- x$emmlev == 2
  rnm = c("(Intercept)", 'psi1', emmvar, paste0(emmvar,":mixture"))
  fam <- x$fit$family$family
  if(showweights) {
    ww = getweightsemm(x, emmval=0.0)
    print(ww)
    cat("\n")
  }
  if(!is.null(x$pos.size1) & showweights & isbinemm) {
    ww = getweightsemm(x, emmval=1.0)
    print(ww)
    cat("\n")
  }
  if (fam == "binomial"){
    estimand <- 'OR'
    if(x$bootstrap && x$msmfit$family$link=='log') estimand = 'RR'
    cat(paste0("## Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (fam == "poisson"){
    estimand <- 'RR'
    cat(paste0("## Mixture log(",estimand,")", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
  }
  if (fam == "gaussian"){
    cat(paste0("## Mixture slope parameters", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "t"
    x$zstat = x$tstat
  }
  if (fam == "cox"){
    cat(paste0("Mixture log(hazard ratio)", ifelse(x$bootstrap, " (bootstrap CI)", " (Delta method CI)"), ":\n\n"))
    testtype = "Z"
    rnm = rnm[-1]
  }
  plab = ifelse(testtype=="Z", "Pr(>|z|)", "Pr(>|t|)")
  if(is.null(dim(x$ci.coef))){
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[1], "Upper CI"=x$ci.coef[2], "test"=x$zstat, "pval"=x$pval)
  } else{
    pdat <- cbind(Estimate=coef(x), "Std. Error"=sqrt(x$var.coef), "Lower CI"=x$ci.coef[,1], "Upper CI"=x$ci.coef[,2], "test"=x$zstat, "pval"=x$pval)
  }
  colnames(pdat)[which(colnames(pdat)=="test")] = eval(paste(testtype, "value"))
  colnames(pdat)[which(colnames(pdat)=="pval")] = eval(paste(plab))
  rownames(pdat) <- rnm
  printCoefmat(pdat,has.Pvalue=TRUE,tst.ind=5L,signif.stars=FALSE, cs.ind=1L:2)
  if( isbinemm ) .printeffects(x, digits=5)
  invisible(x)
}


##### generics ######
vcov.qgcompemmfit <- function(object, ...){
  #' @importFrom stats vcov
  #' @export
  stop("not implemented")
  object$covmat.coef
}


confint.qgcompemmfit <- function(object, ...){
  #' @importFrom stats anova
  #' @export
  message("not yet implemented")
  anova(object$fit)
}
