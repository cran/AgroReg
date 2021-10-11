#' Analysis: Thompson
#'
#' @description This function performs Thompson regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
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
#' The logarithmic model is defined by:
#' \deqn{y = \beta_1 ln(\cdot x) + \beta_2 ln(\cdot x)^2}
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @references Sadeghi, E., Haghighi Asl, A., & Movagharnejad, K. (2019). Mathematical modelling of infrared-dried kiwifruit slices under natural and forced convection. Food science & nutrition, 7(11), 3589-3606.
#' @export
#'
#' @examples
#' library(AgroReg)
#' resp=c(10,8,6.8,6,5,4.3,4.1,4.2,4.1)
#' trat=seq(1,9,1)
#' thompson(trat,resp)


thompson = function(trat,
                    resp,
                    ylab = "Dependent",
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
                    comment = NA) {
  requireNamespace("drc")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  model <- lm(resp ~ log(trat)+I(log(trat)^2)-1)
  model1<- lm(ymean ~ log(xmean)+I(log(xmean)^2)-1)
  coef=summary(model)

  if(is.na(round)==TRUE){
    b=coef$coefficients[,1][1]
    d=coef$coefficients[,1][2]}

  if(is.na(round)==FALSE){
    b=round(coef$coefficients[,1][1],round)
    d=round(coef$coefficients[,1][2],round)}
  if(r2=="all"){r2=coef$r.squared}
  if(r2=="mean"){
    coef1=summary(model1)
    r2=coef1$r.squared}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e*ln(%s) %s %0.3e*ln(%s)^2 ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   d,
                   xname.formula,
                   ifelse(b<0,"-","+"),
                   abs(b),
                   xname.formula,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
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
