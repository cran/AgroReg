#' Analysis: Mitscherlich
#'
#' This function performs Mitscherlich regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param initial List Initial parameters (A, b, e)
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @details
#' The Mitscherlich model is defined by:
#' \deqn{y = A \times (1-10^{-eb-ex})}
#'
#' where "y" is the yield obtained when "b" units of a nutrient are in the soil and
#' "x" units of it are added as fertilizer, "A" is the maximum yield, and "e" is the
#' proportionality factor, has recently received increasing interest.
#'
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @export
#'
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' mitscherlich(time,WL)

mitscherlich=function(trat,
                     resp,
                     initial=NA,
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
                     round=NA,
                     yname.formula="y",
                     xname.formula="x",
                     comment=NA,
                     fontfamily="sans"){
  requireNamespace("ggplot2")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  if(is.na(initial[1])==TRUE){
    A=max(resp)
    b=0.1/max(resp)
    e=max(resp)/2
    initial=list(A=A,b=b,e=e)}
  model <- nls(resp ~ A*(1-10^(-e*b-b*trat)), start = initial)
  coef=summary(model)

  if(is.na(round)==TRUE){
  A=coef$coefficients[,1][1]
  b=coef$coefficients[,1][2]
  e=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    A=round(coef$coefficients[,1][1],round)
    b=round(coef$coefficients[,1][2],round)
    e=round(coef$coefficients[,1][3],round)}


  if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
  if(r2=="mean"){
    model1 <- nls(ymean ~ A*(1-10^(-e*b-b*xmean)), start = initial)
    r2=1-deviance(model1)/deviance(lm(ymean~1))}

  r2=floor(r2*100)/100
  eb=e*b
  equation=sprintf("~~~%s==%0.3e*(1-10^{%0.3e %s %0.3e*%s}) ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   A,
                   eb,
                   ifelse(b >= 0, "+", "-"),
                   abs(b),
                   xname.formula,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(model,newdata = data.frame(trat=xp)))
  predesp=predict(model)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  x=preditos$x
  y=preditos$y
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
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
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  temp1=seq(min(trat),max(trat),length.out=sample.curve)
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
