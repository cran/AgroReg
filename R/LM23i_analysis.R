#' Analysis: Cubic inverse without beta1
#'
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Degree 3 polynomial inverse model without the beta 1 coefficient.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param legend.position legend position (\emph{default} is "top")
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
#'
#' @details
#' Inverse degree 3 polynomial model without the beta 1 coefficient  is defined by:
#' \deqn{y = \beta_0 + \beta_2\cdot x^{1/2} + \beta_3\cdot x^{1/3}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @keywords regression linear
#' @export
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' LM23i(time, WL)

LM23i=function(trat,
               resp,
               ylab="Dependent",
               error="SE",
               xlab="Independent",
               theme=theme_classic(),
               legend.position="top",
               point="all",
               width.bar=NA,
               scale="none",
               textsize = 12,
               pointsize = 4.5,
               linesize = 0.8,
               pointshape = 21,
               round=NA,
               xname.formula="x",
               yname.formula="y",
               comment=NA){
  requireNamespace("ggplot2")
  requireNamespace("crayon")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  dose=tapply(trat, trat, mean)
  moda=lm(resp~I(trat^(1/2))+I(trat^(1/3)))
  mods=summary(moda)$coefficients
  modm=lm(media~I(dose^(1/2))+I(dose^(1/3)))
  r2=round(summary(modm)$r.squared,2)
  coef1=coef(moda)[1]
  coef2=coef(moda)[2]
  coef3=coef(moda)[3]
  s1 = s =  sprintf("~~~%s == %e %s %e*%s^{1/2} %s %e*%s^{1/3} ~~~~~ italic(R^2) == %0.2f",
               yname.formula,
               coef1,
               ifelse(coef2 >= 0, "+", "-"),
               abs(coef2),
               xname.formula,
               ifelse(coef3 >= 0, "+", "-"),
               abs(coef3),
               xname.formula,
               r2)
  if(is.na(comment)==FALSE){s1=paste(s1,"~\"",comment,"\"")}
  data1=data.frame(trat=unique(trat),
                   media=media,
                   desvio)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(moda,newdata = data.frame(trat=xp)))
  x=preditos$x
  y=preditos$y
  predesp=predict(moda)
  predobs=resp
  if(point=="mean"){
    graph=ggplot(data1,aes(x=trat,y=media))
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=media-desvio,
                                                     ymax=media+desvio),
                                                 width=width.bar,
                                                 size=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  grafico=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,color="black"),size=linesize)+
    scale_color_manual(name="",values=1,label=parse(text = s1))+
    theme(axis.text = element_text(size=textsize,color="black"),
          axis.title = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){grafico=grafico+scale_x_log10()}
  models=mods
  model=moda
  r2=summary(modm)$r.squared
  aic=AIC(moda)
  bic=BIC(moda)
  vif=NA
  predesp=predict(moda)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  temp1=seq(min(trat),max(trat),length.out=5000)
  result=predict(moda,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  cat("\n")
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
  graficos=list("Coefficients"=models,
                "values"=graphs,
                "VIF"=vif,
                grafico)
  print(graficos)
}
