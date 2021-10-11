#' @title Analysis: Plateau-quadratic
#'
#' @name plateau.quadratic
#' @rdname plateau.quadratic
#' @description This function performs the plateau-quadratic regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param width.bar	Bar width
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param x Numeric vector with dependent variable.
#' @param a The plateau value
#' @param breakpoint breakpoint value
#' @param b Linear term
#' @param c Quadratic term
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details
#' The Plateau-quadratic model is defined by:
#'
#' First curve:
#' \deqn{y = \beta_0 + \beta_1 \cdot breakpoint + \beta_2 \cdot breakpoint^2 (x < breakpoint)}
#'
#' Second curve:
#' \deqn{y = \beta_0 + \beta_1 \cdot x + \beta_2 \cdot x^2 (x > breakpoint)}
#'
#' or
#'
#' \deqn{y = a + b(x+breakpoint) + c(x+breakpoint)^2 (x > breakpoint)}
#'
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Miguez, F. (2020). nlraa: nonlinear Regression for Agricultural Applications. R package version 0.65.
#' @references Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression theory and applications. Journal of the American Statistical Association 101:542-553.
#' @references Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for identifying ecological thresholds. Ecology 84:2034-2041.
#' @import minpack.lm
#' @import rcompanion
#' @importFrom broom tidy
#' @importFrom stats nls.control
#' @importFrom stats selfStart
#' @importFrom stats sortedXyData
#' @seealso \link{linear.linear}, \link{linear.plateau}
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' x=time[length(time):1]
#' plateau.quadratic(x,WL)
NULL

#' @rdname plateau.quadratic
#' @export
plateau.quadratic=function(trat,resp,
                           ylab="Dependent",
                           xlab="Independent",
                           theme=theme_classic(),
                           legend.position="top",
                           error="SE",
                           r2="all",
                           point="all",
                           width.bar=NA,
                           scale="none",
                           textsize = 12,
                           pointsize = 4.5,
                           linesize = 0.8,
                           pointshape = 21,
                           round=NA,
                           yname.formula="y",
                           xname.formula="x",
                           comment=NA){
  requireNamespace("minpack.lm")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("dplyr")
  requireNamespace("rcompanion")
  mod=nls(resp~plquadratic(trat,a,breakpoint,b,c))
  model2=summary(mod)
  m.ini <- mean(resp)
  nullfunct <- function(x, m) {
    m}
  data=data.frame(trat,resp)
  null <- nls(resp~ nullfunct(trat, m),
              start = list(m = m.ini),
              trace = FALSE,
              data = data,
              nls.control(maxiter = 1000))
  r2 <- nagelkerke(mod, null)$Pseudo.R.squared.for.model.vs.null[2]

  requireNamespace("drc")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  r2=floor(r2*100)/100

  if(is.na(round)==TRUE){
    a=coef(mod)[1]
    breakpoint=coef(mod)[2]
    b=coef(mod)[3]
    c=coef(mod)[4]}

  if(is.na(round)==FALSE){
    a=round(coef(mod)[1],round)
    breakpoint=round(coef(mod)[2],round)
    b=round(coef(mod)[3],round)
    c=round(coef(mod)[4],round)}

  s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^2~(%s>%e) ~~~~~ italic(R^2) ==  %0.2f",
               yname.formula,
               a,
               ifelse(b >= 0, "+", "-"),
               abs(b),
               xname.formula,
               ifelse(c >= 0, "+", "-"),
               abs(c),
               xname.formula,
               xname.formula,
               breakpoint,
               r2)
  equation=s
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  predesp=predict(mod)
  predobs=resp
  model=mod
  rmse=sqrt(mean((predesp-predobs)^2))

  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  if(point=="mean"){
    graph=ggplot(data,aes(x=xmean,y=ymean))
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                                 width=width.bar,
                                                 size=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  xp=seq(min(trat),max(trat),length=1000)
  yp=predict(mod,newdata=data.frame(trat=xp))
  preditos=data.frame(x=xp,y=yp)
  temp1=xp
  result=yp
  x=xp
  y=yp
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=linesize)+
    scale_color_manual(name="",values="black",label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black"),
          axis.title = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  aic=AIC(model)
  bic=BIC(model)
  ybreakpoint=predict(mod,newdata = data.frame(trat=breakpoint))
  graphs=data.frame("Parameter"=c("Breakpoint",
                                  "Response",
                                  "AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(breakpoint,
                               ybreakpoint,
                               aic,
                               bic,
                               r2,
                               rmse))
  graficos=list("Coefficients quadratic model"=model2,
                "values"=graphs,
                graph)
  print(graficos)

}


pquadInit <- function(mCall, LHS, data, ...){
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Small sample size")}
  xy1 <- xy[1:floor(nrow(xy)/2),]
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  xy2$x2 <- xy2[,"x"] - min(xy2[,"x"])
  fit2 <- stats::lm(xy2[,"y"] ~ xy2[,"x2"] + I(xy2[,"x2"]^2))
  a <- coef(fit2)[1]
  b <- coef(fit2)[2]
  c <- coef(fit2)[3]
  objfun <- function(cfs){
    pred <- pquad(xy[,"x"], a=cfs[1], breakpoint=cfs[2], b=cfs[3], c=cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans}
  op <- try(stats::optim(c(a, mean(xy[,"x"]),b,c), objfun,
                         method = "L-BFGS-B",
                         upper = c(Inf, max(xy[,"x"]), Inf, Inf),
                         lower = c(-Inf, min(xy[,"x"]), -Inf, -Inf)),
            silent = TRUE)

  if(class(op) != "try-error"){
    a <- op$par[1]
    breakpoint <- op$par[2]
    b <- op$par[3]
    c <- op$par[4]}else{
      a <- mean(xy1[,"y"])
      breakpoint <- mean(xy[,"x"])
      b <- b
      c <- c}
  value <- c(a, breakpoint, b, c)
  names(value) <- mCall[c("a","breakpoint","b","c")]
  value}

pquad <- function(x, a, breakpoint, b, c){
  .value <- (x < breakpoint) * a + (x >= breakpoint) * (a + b * (x - breakpoint) + c * (x - breakpoint)^2)
  .exp1 <- 1
  .exp2 <- ifelse(x < breakpoint, 0, -b + -(c * (2 * (x - breakpoint))))
  .exp3 <- ifelse(x < breakpoint, 0, x - breakpoint)
  .exp4 <- ifelse(x < breakpoint, 0, (x - breakpoint)^2)
  .actualArgs <- as.list(match.call()[c("a","breakpoint","b","c")])
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a","breakpoint","b","c")))
    .grad[, "a"] <- .exp1
    .grad[, "breakpoint"] <- .exp2
    .grad[, "b"] <- .exp3
    .grad[, "c"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad}
  .value}

#' @rdname plateau.quadratic
#' @export

plquadratic <- selfStart(model = pquad,
                  initial = pquadInit,
                  parameters = c("a","breakpoint","b","c"))
