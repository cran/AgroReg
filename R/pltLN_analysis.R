#' Analysis: Plateau-Linear
#'
#' This function performs the plateau-linear regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param round round equation
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); breakpoint and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression theory and applications. Journal of the American Statistical Association 101:542-553.
#' @references Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for identifying ecological thresholds. Ecology 84:2034-2041.
#' @seealso \link{quadratic.plateau}, \link{linear.linear}
#' @importFrom dplyr if_else
#' @details
#' The plateau-linear model is defined by:
#' First curve:
#' \deqn{y = \beta_0 + \beta_1 \times breakpoint (x < breakpoint)}
#'
#' Second curve:
#' \deqn{y = \beta_0 + \beta_1 \times x (x > breakpoint)}
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' x=time[length(time):1]
#' plateau.linear(x,WL)

plateau.linear=function(trat,resp,
                        sample.curve=1000,
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
                        fillshape = "gray",
                        colorline = "black",
                        round=NA,
                        xname.formula="x",
                        yname.formula="y",
                        comment=NA,
                        fontfamily="sans"){
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  lp <- function(x, a, b, c) {
    if_else(condition = x < c,
            true = a + b * c,
            false = a + b * x)
  }
  data=data.frame(trat,resp)
  ini_fit <- lm(data = data, formula = resp ~ trat)
  ini_a <- ini_fit$coef[[1]]
  ini_b <- ini_fit$coef[[2]]
  ini_c <- mean(data$trat)

  lp_model <- nlsLM(
    formula = resp ~ lp(trat, a, b, c),
    data = data,
    start = list(a = ini_a, b = ini_c, c = ini_c))
  model2=summary(lp_model)
  breakpoint=model2$coefficients[3,1]
  ybreakpoint=predict(lp_model,newdata = data.frame(trat=breakpoint))
  m.ini <- mean(data$resp)

  nullfunct <- function(x, m) {
    m
  }
  null <- nls(resp~ nullfunct(trat, m),
              data = data,
              start = list(m = m.ini),
              trace = FALSE,
              nls.control(maxiter = 1000))
  r2 <- nagelkerke(lp_model, null)$Pseudo.R.squared.for.model.vs.null[2]
  requireNamespace("drc")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  r2=floor(r2*100)/100

  if(is.na(round)==TRUE){
    coef1=coef(model2)[1]
    coef2=coef(model2)[2]
    coef3=breakpoint}

  if(is.na(round)==FALSE){
    coef1=round(coef(model2)[1],round)
    coef2=round(coef(model2)[2],round)
    coef3=round(breakpoint,round)}

  s <- sprintf("~~~%s == %e %s %e * %s ~(%s~'>'~%e) ~~~~~ italic(R^2) ==  %0.2f",
               yname.formula,
               coef1,
               ifelse(coef2 >= 0, "+", "-"),
               abs(coef2),
               xname.formula,
               xname.formula,
               coef3,
               r2)

  equation=s
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  predesp=predict(lp_model)
  predobs=resp
  model=lp_model
  rmse=sqrt(mean((predesp-predobs)^2))

  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
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
  xp=seq(min(trat),max(trat),length=sample.curve)
  yp=predict(lp_model,newdata=data.frame(trat=xp))
  preditos=data.frame(x=xp,y=yp)
  temp1=xp
  result=yp
  x=xp
  y=yp
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=linesize)+
    scale_color_manual(name="",values=colorline,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  aic=AIC(lp_model)
  bic=BIC(lp_model)
  minimo=NA
  respmin=NA
  graphs=data.frame("Parameter"=c("Breakpoint",
                                  "Breakpoint response",
                                  "X Minimum",
                                  "Y Minimum",
                                  "AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(breakpoint,
                               ybreakpoint,
                               minimo,
                               respmin,
                               aic,
                               bic,
                               r2,
                               rmse))
  graficos=list("Coefficients segmented"=model2,
                "values"=graphs,
                graph)
  print(graficos)

}
