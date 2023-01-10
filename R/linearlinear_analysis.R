#' Analysis: Linear-Linear
#'
#' This function performs linear linear regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param middle	A scalar in [0,1]. This represents the range that the change-point can occur in. 0 means the change-point must occur at the middle of the range of x-values. 1 means that the change-point can occur anywhere along the range of the x-values.
#' @param CI Whether or not a bootstrap confidence interval should be calculated. Defaults to FALSE because the interval takes a non-trivial amount of time to calculate
#' @param bootstrap.samples	 The number of bootstrap samples to take when calculating the CI.
#' @param sig.level	What significance level to use for the confidence intervals.
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.position legend position (\emph{default} is "top")
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param linetype line type
#' @param pointshape format point (default is 21)
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @import purrr
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); breakpoint and the graph using ggplot2 with the equation automatically.
#' @details
#' The linear-linear model is defined by:
#' First curve:
#' \deqn{y = \beta_0 + \beta_1 \times x (x < breakpoint)}
#'
#' Second curve:
#' \deqn{y = \beta_0 + \beta_1 \times breakpoint + w \times x (x > breakpoint)}
#'
#' @export
#' @author Model imported from the SiZer package
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression theory and applications. Journal of the American Statistical Association 101:542-553.
#' @references Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for identifying ecological thresholds. Ecology 84:2034-2041.
#' @seealso \link{quadratic.plateau}, \link{linear.plateau}
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' linear.linear(time,WL)

linear.linear=function (trat,
                        resp,
                        middle = 1,
                        CI = FALSE,
                        bootstrap.samples = 1000,
                        sig.level = 0.05,
                        error = "SE",
                        ylab = "Dependent",
                        xlab = "Independent",
                        theme = theme_classic(),
                        point = "all",
                        width.bar = NA,
                        legend.position = "top",
                        textsize = 12,
                        pointsize = 4.5,
                        linesize = 0.8,
                        linetype=1,
                        pointshape = 21,
                        fillshape = "gray",
                        colorline = "black",
                        round = NA,
                        xname.formula="x",
                        yname.formula="y",
                    comment=NA,
                    fontfamily="sans"){
  requireNamespace("ggplot2")
  requireNamespace("purrr")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  piecewise.linear.simple <- function(x, y, middle=1){
    piecewise.linear.likelihood <- function(alpha, x, y){
      N <- length(x);
      w <- (x-alpha);
      w[w<0] <- 0;
      fit <- stats::lm(y ~ x + w);
      Beta <- stats::coefficients(fit);
      Mu <- Beta[1] + Beta[2]*x + Beta[3]*w;
      SSE <- sum(fit$residuals^2);
      sigma2 <- SSE/N;                    # MLE estimate of sigma^2
      likelihood <- sum( log( stats::dnorm(y, mean=Mu, sd=sqrt(sigma2)) ) );
      return(likelihood);
    }

    r <- range(x);
    offset <- r * (1-middle)/2;
    low <- min(x)  + offset;
    high <- max(x) - offset;
    temp <- stats::optimize(piecewise.linear.likelihood, c(low, high), x=x, y=y, maximum=TRUE);
    return(temp$maximum);
  }
  x=trat
  y=resp
  alpha <- piecewise.linear.simple(x, y, middle)
  w <- x - alpha
  w[w < 0] <- 0
  model <- stats::lm(y ~ x + w)
  out <- NULL
  out$change.point <- alpha
  out$model <- model
  out$x <- seq(min(x), max(x), length = 1000)
  w <- out$x - alpha
  w[w < 0] <- 0
  out$y <- stats::predict(out$model, data.frame(x = out$x,
                                                w = w))
  out$CI <- CI
  class(out) <- "PiecewiseLinear"
  if (CI == TRUE){
    data <- data.frame(x = x, y = y)
    my.cp <- function(data, index) {
      x <- data[index, 1]
      y <- data[index, 2]
      cp <- piecewise.linear.simple(x, y)
      w <- x - cp
      w[w < 0] <- 0
      model <- stats::lm(y ~ x + w)
      out <- c(cp, model$coefficients[2], model$coefficients[3],
               model$coefficients[2] + model$coefficients[3])
      return(out)
    }
    boot.result <- boot::boot(data, my.cp, R = bootstrap.samples)
    out$intervals <- apply(boot.result$t, 2, stats::quantile,
                           probs = c(sig.level/2, 1 - sig.level/2))
    colnames(out$intervals) <- c("Change.Point", "Initial.Slope",
                                 "Slope.Change", "Second.Slope")
    out$CI <- t(out$CI)
  }
  ymean=tapply(y,x,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(x,x,mean)
  mod=out
  breaks=mod$change.point

  if(is.na(round)==TRUE){
  b0=mod$model$coefficients[1]
  b1=mod$model$coefficients[2]
  b11=mod$model$coefficients[3]}

  if(is.na(round)==FALSE){
    b0=round(mod$model$coefficients[1],round)
    b1=round(mod$model$coefficients[2],round)
    b11=round(mod$model$coefficients[3],round)}
  r2=round(summary(mod$model)$r.squared,2)
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e %s %0.3e*%s~(%s~'<'~%0.3e)~%s %0.3e*%s~(%s>%0.3e)~~~R^2==%0.2e",
                   yname.formula,
                   b0,
                   ifelse(b1 >= 0, "+", "-"),
                   abs(b1),
                   xname.formula,
                   xname.formula,
                   breaks,
                   ifelse(b11 >= 0, "+", "-"),
                   abs(b11),
                   xname.formula,
                   xname.formula,
                   breaks,
                   r2)
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  s=equation
  temp1=mod$x
  result=mod$y
  preditos1=data.frame(x=temp1,
                        y=result)
  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  model=mod$model
  result1=predict(mod$model)
  rmse=sqrt(mean((result1-resp)^2))
  if(point=="mean"){
    graph=ggplot(data,aes(x=xmean,y=ymean))
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                                 width=width.bar,
                                                 size=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  graph=graph+theme+
    geom_line(data=preditos1,aes(x=x,y=y,color="black"),size=linesize,lty=linetype)+
    scale_color_manual(name="",values=colorline,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  maximo=breaks
  respmax=predict(mod$model,newdata=data.frame(x=maximo,w=0))
  minimo=NA
  respmin=NA
  aic=AIC(mod$model)
  bic=BIC(mod$model)
  graphs=data.frame("Parameter"=c("Breakpoint",
                                  "Breakpoint Response",
                                  "X Minimum",
                                  "Y Minimum",
                                  "AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(maximo,
                               respmax,
                               minimo,
                               respmin,
                               aic,
                               bic,
                               r2,
                               rmse))
  graficos=list("Coefficients"=summary(mod$model),
                "values"=graphs,
                "plot"=graph)
  graficos
}
