######## underlying functions ##########



.plot.md.mod.bounds <- function(x,alpha, emmval=0){
  # working
  ymin <- ymax <- v <- w <- y <- NULL
  modbounds = modelbound(x, emmval=emmval, pwonly = TRUE, alpha = alpha)
  geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax,
                  fill="Model confidence band"),
              data=data.frame(ymin=modbounds$ll.pw, ymax=modbounds$ul.pw,
                              x=modbounds$quantile.midpoint))
}


.plot.rr.mod.bounds <- .plot.md.mod.bounds
.plot.or.mod.bounds <- .plot.md.mod.bounds

.plot.linear.smooth.line <- function(x, emmval = 0){
  ymin <- ymax <- v <- w <- y <- NULL
  emmvar = x$emmvar.msm
  # this will work for binary
  emmidx = which(emmvar==emmval)
  yexp = x$y.expected[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if(length(emmdidx)<5){
    message("Too few observed values at emmval, suppressing smooth fit")
    ret = theme()
  } else{
    ret = geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
                method = "gam", formula = y ~ s(x, k=length(table(ind))-1, bs = "cs"),
                data=data.frame(y=yexp, x=(ind+0.5)/max(ind+1)))
  }
  ret
}

.plot.rr.smooth.line <- .plot.linear.smooth.line

.plot.or.smooth.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  emmvar = x$emmvar.msm
  # this will work for binary
  emmidx = which(emmvar==emmval)
  yexp = x$y.expected[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if(length(emmdidx)<5){
    message("Too few observed values at emmval, suppressing smooth fit")
    ret = theme()
  } else{
    ret = geom_smooth(aes(x=x,y=y, color="Smooth conditional fit"),se = FALSE,
              method = "gam", formula = y ~ s(x, k=length(table(yexp))-1, bs = "cs"),
              data=data.frame(y=(yexp)/(1-yexp), x=(yexp+0.5)/max(yexp+1)))
  }
  ret
}



.plot.linear.line <- function(x, emmval){
  ymin <- ymax <- v <- w <- y <- NULL
  emmvar = x$emmvar.msm
  # this will work for binary
  emmidx = which(emmvar==emmval)
  yexp = x$y.expectedmsm[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if(length(emmdidx)<5){
    message("Too few observed values at emmval, suppressing linear fit")
    ret = theme()
  } else{
    ret = geom_line(aes(x=x,y=y, color="MSM fit"),
            data=data.frame(y=yexp, x=(ind+0.5)/max(ind+1)))
  }
  ret
}

.plot.loglin.line <- .plot.linear.line

.plot.logitlin.line <- function(x){
  ymin <- ymax <- v <- w <- y <- NULL
  emmvar = x$emmvar.msm
  # this will work for binary
  emmidx = which(emmvar==emmval)
  yexp = x$y.expectedmsm[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if(length(emmdidx)<5){
    message("Too few observed values at emmval, suppressing linear fit")
    ret = theme()
  } else{
    ret = geom_line(aes(x=x,y=y, color="MSM fit"),
            data=data.frame(y=(yexp/(1-yexp)), x=(ind+0.5)/max(ind+1)))
  }
  ret
}



.plot.md.pw.boot <- function(x, alpha, pointwiseref, emmval=0){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound(x, alpha=alpha, pointwiseref=pointwiseref, emmval=emmval)
  py = pwbdat$hx
  ll = pwbdat$ll.linpred
  ul = pwbdat$ul.linpred
  list(
    geom_point(aes(x=x,y=y,
                   color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
               data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax,
                      color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                  data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.rr.pw.boot <- function(x, alpha, pointwiseref, emmval=0){
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound(x, alpha=alpha, pointwiseref=pointwiseref, emmval=emmval)
  py = exp(pwbdat$hx)
  ll = exp(pwbdat$ll.linpred)
  ul = exp(pwbdat$ul.linpred)
  list(
    geom_point(aes(x=x,y=y,
                   color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")),
               data=data.frame(y=py, x=pwbdat$quantile.midpoint)) ,
    geom_errorbar(aes(x=x,ymin=ymin,ymax=ymax,
                      color=paste0("Pointwise ",as.character(100*(1-alpha)),"% CI")), width = 0.03,
                  data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint))
  )
}

.plot.or.pw.boot <- .plot.rr.pw.boot

.plfun <- function(plt){
  grid::grid.newpage()
  grid::grid.draw(plt)
}

.plot.boot.gaussian <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0){
  if(!(x$msmfit$family$link == "identity")) stop("Plotting not implemented for this link function")
  p <- p + labs(x = "Joint exposure quantile", y = "Y") + lims(x=c(0,1))
  #
  if(modelband)     p <- p + .plot.md.mod.bounds(x,alpha=alpha, emmval=emmval) # : add alpha to main function
  if(flexfit)       p <- p + .plot.linear.smooth.line(x, emmval=emmval)
  if(modelfitline)  p <- p + .plot.linear.line(x, emmval=emmval)
  if(pointwisebars) p <- p + .plot.md.pw.boot(x,alpha,pointwiseref, emmval=emmval)
  p
}


.plot.boot.binomial <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0){
  if(!(x$msmfit$family$link %in% c("log", "logit"))) stop("Plotting not implemented for this link function")
  #
  p <- p + scale_y_log10()
  if(x$msmfit$family$link == "logit"){
    p <- p + labs(x = "Joint exposure quantile", y = "Odds(Y=1)") + lims(x=c(0,1))
    if(modelband) p <- p + .plot.or.mod.bounds(x,alpha, emmval=emmval)
    if(flexfit)   p <- p + .plot.or.smooth.line(x, emmval=emmval)
    if(modelfitline) p <- p + .plot.logitlin.line(x, emmval=emmval)
    if(pointwisebars) p <- p + .plot.or.pw.boot(x,alpha,pointwiseref, emmval=emmval)
  } else if(x$msmfit$family$link=="log"){
    p <- p + labs(x = "Joint exposure quantile", y = "Pr(Y=1)") + lims(x=c(0,1))
    if(modelband) p <- p + .plot.rr.mod.bounds(x,alpha, emmval=emmval)
    if(flexfit)   p <- p + .plot.rr.smooth.line(x, emmval=emmval)
    if(modelfitline) p <- p + .plot.loglin.line(x, emmval=emmval)
    if(pointwisebars) p <- p + .plot.rr.pw.boot(x,alpha,pointwiseref, emmval=emmval)
  }
  p
}


.plot.boot.poisson <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0){
  if(!(x$msmfit$family$link == "log")) stop("Plotting not implemented for this link function")
  p <- p + scale_y_log10()
  if(x$msmfit$family$link == "log"){
    p <- p + labs(x = "Joint exposure quantile", y = "E(Y)") + lims(x=c(0,1))
    if(modelband) p <- p + .plot.rr.mod.bounds(x,alpha, emmval=emmval)
    if(flexfit)   p <- p + .plot.rr.smooth.line(x, emmval=emmval)
    if(modelfitline) p <- p + .plot.loglin.line(x, emmval=emmval)
    if(pointwisebars) p <- p + .plot.rr.pw.boot(x,alpha,pointwiseref, emmval=emmval)
  }
  p
}

.plot.boot.cox <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=emmval){
  # : make the plot more configurable
  surv <- NULL
  stop("Not yet implemented for survival models")
  scl = qgcomp.survcurve.boot(x, emmval=emmval)
  cdf0 = scl$cdfq[scl$cdfq$q==1,]
  cdfmax = scl$cdfq[scl$cdfq$q==x$q,]
  mdf0 = scl$mdfq[scl$mdfq$q==1,]
  mdfmax = scl$mdfq[scl$mdfq$q==x$q,]
  p <- p +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Average (all quantiles)"), data=scl$mdfpop)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Average (all quantiles)"), data=scl$cdfpop) +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Lowest quantile"), data=mdf0)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Lowest quantile"), data=cdf0) +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Highest quantile"), data=mdfmax)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Highest quantile"), data=cdfmax) +
    scale_y_continuous(name="Survival", limits=c(0,1)) +
    scale_x_continuous(name="Time") +
    scale_linetype_discrete(name="")+
    theme(legend.position = c(0.01, 0.01), legend.justification = c(0,0))
  return(p)
}



#' @title Default plotting method for a qgcompfit object
#'
#' @description Plot a quantile g-computation object. For qgcomp.noboot, this function will
#' create a butterfly plot of weights. For qgcomp.boot, this function will create
#' a box plot with smoothed line overlaying that represents a non-parametric
#' fit of a model to the expected outcomes in the population at each quantile
#' of the joint exposures (e.g. '1' represents 'at the first quantile for
#' every exposure')
#'
#' @param x "qgcompfit" object from `qgcomp.noboot`,  `qgcomp.boot`,
#'   `qgcomp.cox.noboot`,  `qgcomp.cox.boot`, `qgcomp.zi.noboot` or `qgcomp.zi.boot` functions
#' @param emmval fixed value for effect measure modifier at which pointwise comparisons are calculated
#' @param suppressprint If TRUE, suppresses the plot, rather than printing it
#'   by default (it can be saved as a ggplot2 object (or list of ggplot2 objects if x is from a zero-
#'   inflated model) and used programmatically)
#'   (default = FALSE)
#' @param ... unused
#' @seealso \code{\link[qgcomp]{qgcomp.noboot}}, \code{\link[qgcomp]{qgcomp.boot}}, and \code{\link[qgcomp]{qgcomp}}
#' @import ggplot2 grid gridExtra qgcomp
#' @importFrom grDevices gray
#' @export
#' @examples
#' set.seed(12)
#' \dontrun{
#' }
plot.qgcompemmfit <- function(x,emmval=0.0, suppressprint=FALSE,...){
  if(!x$bootstrap){
    zwts <- qgcompint::getstratweights(x,emmval=emmval)
    x$pos.weights <- zwts$pos.weights
    x$neg.weights <- zwts$neg.weights
    x$pos.psi <- zwts$pos.psi
    x$neg.psi <- zwts$neg.psi
    class(x) <- setdiff(class(x), "qgcompemmfit")
    # this should plot from the qgcomp package
    p <- plot(x, suppressprint=suppressprint, ...)
  }

  if(x$bootstrap){
    #stop("not yet implemented")
    # TODO: implement some qgcomp defaults
    # variance based on delta method and knowledge that non-linear
    #functions will always be polynomials in qgcomp
    # default plot for bootstrap results (no weights obtained)
    p <- ggplot()
    if(is.null(x$msmfit$family)) {
      # ZI model
      stop("Not implemented for this model")
      #p <- .plot.boot.zi(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref, alpha=0.05)
    } else{
      # Cox model or standard GLM
      if(x$msmfit$family$family=='cox') p <-
          .plot.boot.cox(p, x, ..., alpha=x$alpha, emmval=emmval)
      if(x$msmfit$family$family=='gaussian') p <-
          .plot.boot.gaussian(p, x, ..., alpha=x$alpha, emmval=emmval)
      if(x$msmfit$family$family=='binomial') p <-
          .plot.boot.binomial(p, x, ..., alpha=x$alpha, emmval=emmval)
      if(x$msmfit$family$family=='poisson') p <-
          .plot.boot.poisson(p, x, ..., alpha=x$alpha, emmval=emmval)
    }

    p <- p + scale_fill_grey(name="", start=.9) +
      scale_colour_grey(name="", start=0.0, end=0.6) +
      theme_classic()
    if(!suppressprint) print(p)
  }
  if(suppressprint) return(p)
}
#
#qgcomp::
#
#
#plot(qfit2, emmval=0)
#plot(qfit2, emmval=1)
#plot(qfit2, emmval=2)
#
