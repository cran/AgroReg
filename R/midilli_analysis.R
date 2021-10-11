#' Analysis: Midilli
#'
#' This function performs Midilli regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param initial Starting estimates
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details
#' The exponential model is defined by:
#' \deqn{y = \alpha \times e^{-\beta \cdot x^n} + \theta \cdot x}
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @export
#'
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' midilli(time,100-WL)

midilli=function(trat,
                 resp,
                 initial=NA,
                 ylab="Dependent",
                 xlab = "Independent",
                 theme = theme_classic(),
                 legend.position = "top",
                 error = "SE",
                 r2 = "all",
                 point = "all",
                 width.bar = NA,
                 scale = "none",
                 textsize = 12,
                 pointsize = 4.5,
                 linesize = 0.8,
                 pointshape = 21,
                 round = NA,
                 yname.formula="y",
                 xname.formula="x",
                 comment=NA){
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)

  if(is.na(initial[1])==TRUE){
  theta.0 <- min(resp) * 0.5
  model.0 <- lm(log(resp - theta.0) ~ trat)
  alpha.0 <- exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  theta.0<-theta.0/mean(trat)
  n=1
  initial <- list(alpha = alpha.0, beta = beta.0, theta = theta.0, n=n)}
  model <- nls(resp ~ alpha * exp(beta * trat^n) + theta*trat, start = initial)
  coef=summary(model)

  if(is.na(round)==TRUE){
    b=coef$coefficients[,1][1]
    d=coef$coefficients[,1][2]
    e=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    b=round(coef$coefficients[,1][1],round)
    d=round(coef$coefficients[,1][2],round)
    e=round(coef$coefficients[,1][3],round)}

  # if(r2=="all"){r2=cor(resp, fitted(model))^2}
  # if(r2=="mean"){r2=cor(ymean, predict(model,newdata=data.frame(trat=unique(trat))))^2}
  if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
  if(r2=="mean"){
    model1 <- nls(ymean ~ alpha * exp(beta * xmean^n) + theta*xmean, start = initial)
    r2=1-deviance(model1)/deviance(lm(ymean~1))}

  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e*e^{%0.3e*%s}+ %s %0.3e ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   b,
                   d,
                   xname.formula,
                   ifelse(e<0,"-","+"),
                   abs(e),
                   r2)
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(model,newdata = data.frame(trat=xp)))
  predesp=predict(model)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  x=preditos$x
  y=preditos$y
  s=equation
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

  graph=graph+theme+geom_line(data=preditos,aes(x=x,
                                                y=y,color="black"),size=linesize)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black"),
          axis.title = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(model,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  aic=AIC(model)
  bic=BIC(model)
  graphs=data.frame("Parameter"=c("X Maximum",
                                  "Y Maximum",
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
  graficos=list("Coefficients"=coef,
                "values"=graphs,
                graph)
  print(graficos)
}