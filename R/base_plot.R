# note: this was an attempt to make the plots functions work on proto objects by using
# lists rather than ggplot_add(); I couldn't afford the testing time, but it's worth pursuing.




######## underlying functions ##########



.plot.md.mod.bounds <- function(x, alpha, emmval=0) {
  # working
  ymin <- ymax <- v <- w <- y <- NULL
  modbounds = modelbound(x, emmval=emmval, pwonly = TRUE, alpha = alpha)
  ret = geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax,
                  fill=col),
              data=data.frame(ymin=modbounds$ll.pw, ymax=modbounds$ul.pw,
                              x=modbounds$quantile.midpoint,
                              col = "Model confidence band"))
  list(ribbon=ret)
}


.plot.rr.mod.bounds <- .plot.md.mod.bounds
.plot.or.mod.bounds <- .plot.md.mod.bounds

.plot.linear.smooth.line <- function(x, emmval = 0) {
  ymin <- ymax <- v <- w <- y <- NULL
    emmvar = x$emmvar.msm
  if (is.factor(x$emmvals)){
    zexp = zproc(factor(emmval, levels(x$emmvals)), x$call$emmvar)
    emmidx = which(sapply(1:nrow(emmvar), function(x) all(emmvar[x,] == zexp)))
  } else{
    # this will work for binary
    emmidx = which(emmvar==emmval) # need to modify with proc
  }
  yexp = x$y.expected[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if (length(emmidx)<5) {
    message("Too few observed values at emmval, suppressing smooth fit")
    ret = theme()
  } else{
    ret = geom_smooth(aes(x=x, y=y, color=col), se = FALSE,
                method = "gam", formula = y ~ s(x, k=length(table(ind))-1, bs = "cs"),
                data=data.frame(y=yexp, x=(ind+0.5)/max(ind+1),
                                col="Smooth conditional fit"))
  }
  list(smooth=ret)
}

.plot.rr.smooth.line <- .plot.linear.smooth.line

.plot.or.smooth.line <- function(x, emmval = 0) {
  ymin <- ymax <- v <- w <- y <- NULL
    emmvar = x$emmvar.msm
  if (is.factor(x$emmvals)){
    zexp = zproc(factor(emmval, levels(x$emmvals)), x$call$emmvar)
    emmidx = which(sapply(1:nrow(emmvar), function(x) all(emmvar[x,] == zexp)))
  } else{
    # this will work for binary
    emmidx = which(emmvar==emmval) # need to modify with proc
  }
  yexp = x$y.expected[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if (length(emmidx)<5) {
    message("Too few observed values at emmval, suppressing smooth fit")
    ret = theme()
  } else{
    ret = geom_smooth(aes(x=x, y=y, color=col), se = FALSE,
              method = "gam", formula = y ~ s(x, k=length(table(yexp))-1, bs = "cs"),
              data=data.frame(y=(yexp)/(1-yexp), x=(yexp+0.5)/max(yexp+1),
                              col="Smooth conditional fit"))
  }
  list(smooth=ret)
}



.plot.linear.line <- function(x, emmval = 0) {
  ymin <- ymax <- v <- w <- y <- NULL
    emmvar = x$emmvar.msm
  if (is.factor(x$emmvals)){
    zexp = zproc(factor(emmval, levels(x$emmvals)), x$call$emmvar)
    emmidx = which(sapply(1:nrow(emmvar), function(x) all(emmvar[x,] == zexp)))
  } else{
    # this will work for binary
    emmidx = which(emmvar==emmval) # need to modify with proc
  }
  yexp = x$y.expectedmsm[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if (length(emmidx)<5) {
    message("Too few observed values at emmval, suppressing linear fit")
    ret = theme()
  } else{
    ret = geom_line(aes(x=x, y=y, color=col),
            data=data.frame(y=yexp, x=(ind+0.5)/max(ind+1), col="MSM fit"))
  }
  list(line=ret)
}

.plot.loglin.line <- .plot.linear.line

.plot.logitlin.line <- function(x, emmval = 0) {
  ymin <- ymax <- v <- w <- y <- NULL
    emmvar = x$emmvar.msm
  if (is.factor(x$emmvals)){
    zexp = zproc(factor(emmval, levels(x$emmvals)), x$call$emmvar)
    emmidx = which(sapply(1:nrow(emmvar), function(x) all(emmvar[x,] == zexp)))
  } else{
    # this will work for binary
    emmidx = which(emmvar==emmval) # need to modify with proc
  }
  yexp = x$y.expectedmsm[emmidx]
  ind = x$index[emmidx]
  #predict(x$msmfit, newdata = data.frame(M=0))
  #
  if (length(emmidx)<5) {
    message("Too few observed values at emmval, suppressing linear fit")
    ret = theme()
  } else{
    ret = geom_line(aes(x=x, y=y, color=col),
            data=data.frame(y=(yexp/(1-yexp)), x=(ind+0.5)/max(ind+1),
                            col="MSM fit"))
  }
  list(line=ret)
}



.plot.md.pw.boot <- function(x, alpha, pointwiseref, emmval = NULL) {
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound(x, alpha=alpha, pointwiseref=pointwiseref, emmval=emmval)
  py = pwbdat$hx
  md = pwbdat$mean.diff
  ll = py-md+pwbdat$ll.diff
  ul = py-md+pwbdat$ul.diff
  list(
    point = geom_point(aes(x=x, y=y, color=col),
               data=data.frame(y=py, x=pwbdat$quantile.midpoint,
                               col=paste0("Pointwise ", 100*(1-alpha), "% CI"))) ,
    errorbar = geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax, color=col),
                  width = 0.03,
                  data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint,
                                  col=paste0("Pointwise ", 100*(1-alpha), "% CI")))
  )
}

.plot.rr.pw.boot <- function(x, alpha, pointwiseref, emmval=0) {
  ymin <- ymax <- v <- w <- y <- NULL
  # plot actual risk or odds bounds
  pwbdat = pointwisebound(x, alpha=alpha, pointwiseref=pointwiseref, emmval=emmval)
  py = exp(pwbdat$hx)
  ll = exp(pwbdat$ll.linpred)
  ul = exp(pwbdat$ul.linpred)
  list(
    point = geom_point(aes(x=x, y=y, color=col),
               data=data.frame(y=py, x=pwbdat$quantile.midpoint,
                               col=paste0("Pointwise ", 100*(1-alpha), "% CI"))) ,
    errorbar = geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax, color=col),
                  width = 0.03,
                  data=data.frame(ymin=ll, ymax=ul, x=pwbdat$quantile.midpoint,
                                  col=paste0("Pointwise ", 100*(1-alpha), "% CI")))
  )
}

.plot.or.pw.boot <- .plot.rr.pw.boot

.plfun <- function(plt) {
  grid::grid.newpage()
  grid::grid.draw(plt)
}

.plot.boot.gaussian <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0) {
  if (!(x$msmfit$family$link == "identity")) stop("Plotting not implemented for this link function")
  p <- c(p, list(labs(x = "Joint exposure quantile", y = "Y"),  lims(x=c(0, 1))))
  #
  if (modelband)     p <- c(p, .plot.md.mod.bounds(x, alpha=alpha, emmval=emmval) )# : add alpha to main function
  if (flexfit)       p <- c(p, .plot.linear.smooth.line(x, emmval=emmval))
  if (modelfitline)  p <- c(p, .plot.linear.line(x, emmval=emmval))
  if (pointwisebars) p <- c(p, .plot.md.pw.boot(x, alpha, pointwiseref, emmval=emmval))
  p
}


.plot.boot.binomial <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0) {
  if (!(x$msmfit$family$link %in% c("log", "logit"))) stop("Plotting not implemented for this link function")
  #
  p <- c(p, scale_y_log10())
  if (x$msmfit$family$link == "logit") {
    p <- c(p, list(labs(x = "Joint exposure quantile", y = "Odds(Y=1)"),  lims(x=c(0, 1))))
    if (modelband) p <- c(p, .plot.or.mod.bounds(x, alpha, emmval=emmval))
    if (flexfit)   p <- c(p, .plot.or.smooth.line(x, emmval=emmval))
    if (modelfitline) p <- c(p, .plot.logitlin.line(x, emmval=emmval))
    if (pointwisebars) p <- c(p, .plot.or.pw.boot(x, alpha, pointwiseref, emmval=emmval))
  } else if (x$msmfit$family$link=="log") {
    p <- c(p, list(labs(x = "Joint exposure quantile", y = "Pr(Y=1)"),  lims(x=c(0, 1))))
    if (modelband) p <- c(p, .plot.rr.mod.bounds(x, alpha, emmval=emmval))
    if (flexfit)   p <- c(p, .plot.rr.smooth.line(x, emmval=emmval))
    if (modelfitline) p <- c(p, .plot.loglin.line(x, emmval=emmval))
    if (pointwisebars) p <- c(p, .plot.rr.pw.boot(x, alpha, pointwiseref, emmval=emmval))
  }
  p
}


.plot.boot.poisson <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=0) {
  if (!(x$msmfit$family$link == "log")) stop("Plotting not implemented for this link function")
  p <- c(p, scale_y_log10())
  if (x$msmfit$family$link == "log") {
    p <- c(p, list(labs(x = "Joint exposure quantile", y = "E(Y)"),  lims(x=c(0, 1))))
    if (modelband) p <- c(p, .plot.rr.mod.bounds(x, alpha, emmval=emmval))
    if (flexfit)   p <- c(p, .plot.rr.smooth.line(x, emmval=emmval))
    if (modelfitline) p <- c(p, .plot.loglin.line(x, emmval=emmval))
    if (pointwisebars) p <- c(p, .plot.rr.pw.boot(x, alpha, pointwiseref, emmval=emmval))
  }
  p
}

.plot.boot.cox <- function(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, emmval=emmval) {
  # : make the plot more configurable
  surv <- NULL
  stop("Not yet implemented for survival models")
  scl = qgcomp.survcurve.boot(x, emmval=emmval)
  cdf0 = scl$cdfq[scl$cdfq$q==1, ]
  cdfmax = scl$cdfq[scl$cdfq$q==x$q, ]
  mdf0 = scl$mdfq[scl$mdfq$q==1, ]
  mdfmax = scl$mdfq[scl$mdfq$q==x$q, ]
  p <- p +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Average (all quantiles)"), data=scl$mdfpop)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Average (all quantiles)"), data=scl$cdfpop) +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Lowest quantile"), data=mdf0)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Lowest quantile"), data=cdf0) +
    geom_step(aes(x=time, y=surv, color="MSM", linetype="Highest quantile"), data=mdfmax)+
    geom_step(aes(x=time, y=surv, color="Conditional", linetype="Highest quantile"), data=cdfmax) +
    scale_y_continuous(name="Survival", limits=c(0, 1)) +
    scale_x_continuous(name="Time") +
    scale_linetype_discrete(name="")+
    theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0))
  return(p)
}



#' @title Default plotting method for a qgcompemmfit object
#'
#' @description Plot a quantile g-computation object from qgcompint.
#' For qgcomp.emm.glm.noboot, this function will create a butterfly plot of
#' weights. For qgcomp.emm.glm.boot and qgcomp.emm.glm.ee, this function will create
#' a box plot with smoothed line overlaying that represents a non-parametric
#' fit of a model to the expected outcomes in the population at each quantile
#' of the joint exposures (e.g. '0' represents 'at the lowest quantile for
#' every exposure')
#'
#' @param x "qgcompemmfit" object from `qgcomp.emm.glm.noboot`, `qgcomp.emm.glm.boot`,
#'   `qgcomp.emm.glm.ee`, or `qgcomp.emm.cox.noboot` functions
#' @param emmval fixed value for effect measure modifier at which pointwise comparisons are calculated
#' @param suppressprint If TRUE, suppresses the plot, rather than printing it
#'   by default (it can be saved as a ggplot2 object (or list of ggplot2 objects if x is from a zero-
#'   inflated model) and used programmatically)
#' @param geom_only If TRUE, returns only the geometry (i.e. does not contain the entire plot object). Used for overlays. Only used for `.ee` and `.boot` methods.
#'   (default = FALSE)
#' @param pointwisebars (boot/ee only) If TRUE (TRUE=default), adds 95% error bars for pointwise comparisons of E(Y|joint exposure) to the smooth regression line plot
#' @param modelfitline (boot/ee only) If TRUE (FALSE=default), adds fitted (MSM) regression line of E(Y|joint exposure) to the smooth regression line plot
#' @param modelband (boot/ee only) If TRUE (FALSE=default), adds 95% prediction bands for E(Y|joint exposure) (the MSM fit)
#' @param flexfit (boot/ee only) if TRUE (FALSE=default), adds flexible interpolation of predictions from underlying (conditional) model
#' @param pointwiseref (boot/ee only) integer (0=default): which category of exposure (from 1 to q) should serve as the referent category for pointwise comparisons? (default=1)
#' @param alpha alpha level for all confidence intervals
#' @param ... Arguments (listed under ``details'' ) to underlying functions:
#'
#'
#' @details
#' The `...` argument calls underlying plot functions, with arguments given here (similar to the plot function in the `qgcomp` package)
#' @return
#'
#' If suppressprint=FALSE, then this function prints a plot specific to a "qgcompemmfit" object.
#' If no bootstrapping is used, it will print a butterfly plot of the weights at the specified value of the modifier (set via `emmval` parameter)
#' If bootstrapping is used, it will print a joint regression line for all exposures at the specified value of the modifier (set via `emmval` parameter)
#'
#'  If suppressprint=TRUE, then this function returns a "gg" (regression line) or "gtable" (butterfly plot) object (from ggplot2 package or gtable/grid packages), which can be used to print a ggplot figure and modify either of the above figures (see example below)
#'
#'
#'
#' @seealso \code{\link[qgcompint]{qgcomp.emm.glm.noboot}}, \code{\link[qgcompint]{qgcomp.emm.glm.boot}},  \code{\link[qgcompint]{qgcomp.emm.glm.ee}}, \code{\link[qgcompint]{qgcomp.emm.cox.noboot}}
#' @import ggplot2 grid gridExtra qgcomp
#' @importFrom grDevices gray
#' @exportS3Method base::plot
#' @examples
#' set.seed(50)
#' # linear model, binary modifier
#' dat <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#' z=rbinom(50, 1, 0.5), r=rbinom(50, 1, 0.5))
#' (qfit <- qgcomp.emm.glm.noboot(f=y ~ z + x1 + x2, emmvar="z",
#' expnms = c('x1', 'x2'), data=dat, q=2, family=gaussian()))
#' plot(qfit, emmval = 1)
#' #
#' library(ggplot2)
#'
#' # example with estimating equation approach
#' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#' z=sample(0:2, 50, replace=TRUE), r=rbinom(50, 1, 0.5))
#' dat2$z = as.factor(dat2$z)
#' (qfit4ee <- qgcomp.emm.glm.ee(f=y ~ z + x1 + x2, emmvar="z",
#'                           degree = 1,
#'                          expnms = c('x1', 'x2'), data=dat2, q=4, family=gaussian()))
#' pp0ee = plot(qfit4ee, emmval=0, suppressprint=TRUE)
#' pp1ee = plot(qfit4ee, emmval=1, suppressprint=TRUE)
#' pp2ee = plot(qfit4ee, emmval=2, suppressprint=TRUE)
#' pp1ee + theme_linedraw() # can use with other ggplot functions
#' # example with bootstrapping
#' dat2 <- data.frame(y=runif(50), x1=runif(50), x2=runif(50),
#' z=sample(0:2, 50, replace=TRUE), r=rbinom(50, 1, 0.5))
#' dat2$z = as.factor(dat2$z)
#' (qfit4 <- qgcomp.emm.glm.boot(f=y ~ z + x1 + x2, emmvar="z",
#'                           degree = 1, B = 20,
#'                          expnms = c('x1', 'x2'), data=dat2, q=4, family=gaussian()))
#' pp0 = plot(qfit4, emmval=0, suppressprint=TRUE)
#' pp1 = plot(qfit4, emmval=1, suppressprint=TRUE)
#' pp2 = plot(qfit4, emmval=2, suppressprint=TRUE)
#' pp1 + theme_linedraw() # can use with other ggplot functions
#'
#' # overlay (fussy to work with)
#' #ppom <- ggplot_build(pp0 + pp1[2] + pp2[[2]] + scale_color_discrete(guide="none"))
#' ppom <- ggplot_build(pp0ee + pp1ee[2] + pp2ee[[2]] + scale_color_discrete(guide="none"))
#' ppom$data[[1]]$colour <- ppom$data[[2]]$colour <- "gray40" # emmval = 0 -> dark gray
#' ppom$data[[3]]$colour <- ppom$data[[4]]$colour <- "gray80" # emmval = 1 -> light gray
#' ppom$data[[5]]$colour <- ppom$data[[6]]$colour <- "black" # emmval = 2 -> black
#' xincrement = 0.025
#' ppom$data[[1]]$x <- ppom$data[[2]]$x <- ppom$data[[1]]$x - xincrement
#' ppom$data[[2]]$xmin <- ppom$data[[2]]$xmin - xincrement
#' ppom$data[[2]]$xmax <- ppom$data[[2]]$xmax - xincrement
#' ppom$data[[5]]$x <- ppom$data[[6]]$x <- ppom$data[[5]]$x + xincrement
#' ppom$data[[6]]$xmin <- ppom$data[[6]]$xmin + xincrement
#' ppom$data[[6]]$xmax <- ppom$data[[6]]$xmax + xincrement
#' plot(ggplot_gtable(ppom))
#'
#' \dontrun{
#' library(gtable) # may need to separately install gtable
#' # example with no bootstrapping, adding object from bootstrapped fit
#' pp2 <- plot(qfit, emmval = 1, suppressprint=TRUE)
#' grid.draw(pp2)
#' # insert row on top that is 1/2 height of existing plot
#' pp2b = gtable::gtable_add_rows(pp2, heights=unit(0.5, 'null') , pos = 0)
#' # add plot to that row
#' pp3 = gtable::gtable_add_grob(pp2b, ggplot2::ggplotGrob(pp), t=1, l=1, r=2)
#' grid.draw(pp3)
#' }
plot.qgcompemmfit <- function(x, emmval = NULL, suppressprint=FALSE, geom_only=FALSE, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref=1, alpha=0.05, ...) {
  if (is.null(emmval)) {
    stop("emmval must be specified (level of the modifier for which you would like results)")
  }
  if (!x$bootstrap && !inherits(x, "eeqgcompfit")) {
    zwts <- qgcompint::getstratweights(x, emmval=emmval)
    x$pos.weights <- zwts$pos.weights
    x$neg.weights <- zwts$neg.weights
    x$pos.psi <- zwts$pos.psi
    x$neg.psi <- zwts$neg.psi
    class(x) <- setdiff(class(x), "qgcompemmfit")
    # this should plot from the qgcomp package
    p <- plot(x, suppressprint=suppressprint, ...)
  }

  if (x$bootstrap || inherits(x, "eeqgcompfit")) {
    if (is.null(x$q)) {
      stop("Plotting not enabled when fitting models with q=NULL")
    }

    #stop("not yet implemented")
    # TODO: implement some qgcomp defaults
    # variance based on delta method and knowledge that non-linear
    #functions will always be polynomials in qgcomp
    # default plot for bootstrap results (no weights obtained)
    if(geom_only & !suppressprint){
      message("geomonly=TRUE, suppressing output by setting suppressprint=TRUE")
      suppressprint=TRUE
    } else{

    }
    p <- list()

    if (is.null(x$msmfit$family)) {
      # ZI model
      stop("Not implemented for this model")
      #p <- .plot.boot.zi(p, x, modelband=FALSE, flexfit=FALSE, modelfitline=FALSE, pointwisebars=TRUE, pointwiseref, alpha=0.05)
    } else{
      # Cox model or standard GLM
      if (x$msmfit$family$family=='cox') p <-
          .plot.boot.cox(p, x, modelband=modelband, flexfit=flexfit, modelfitline=modelfitline, pointwisebars=pointwisebars, pointwiseref=pointwiseref, alpha=alpha, ..., emmval=emmval)
      if (x$msmfit$family$family=='gaussian') p <-
          .plot.boot.gaussian(p, x, modelband=modelband, flexfit=flexfit, modelfitline=modelfitline, pointwisebars=pointwisebars, pointwiseref=pointwiseref, alpha=alpha, ..., emmval=emmval)
      if (x$msmfit$family$family=='binomial') p <-
          .plot.boot.binomial(p, x, modelband=modelband, flexfit=flexfit, modelfitline=modelfitline, pointwisebars=pointwisebars, pointwiseref=pointwiseref, alpha=alpha, ..., emmval=emmval)
      if (x$msmfit$family$family=='poisson') p <-
          .plot.boot.poisson(p, x, modelband=modelband, flexfit=flexfit, modelfitline=modelfitline, pointwisebars=pointwisebars, pointwiseref=pointwiseref, alpha=alpha, ..., emmval=emmval)
    }

    p <- c(p, list(
    				scale_fill = scale_fill_grey(name="", start=.9),
			        scale_color = scale_colour_grey(name="", start=0.0, end=0.6),
			        theme = theme_classic()
			      )
		  )
    if (!suppressprint) print(ggplot2::ggplot() + p)
  }
  if (suppressprint & geom_only ) return(p)
  if (suppressprint & !geom_only ) return(ggplot2::ggplot() + p)
}
#
#qgcomp::
#
#
#plot(qfit2, emmval=0)
#plot(qfit2, emmval=1)
#plot(qfit2, emmval=2)
#


