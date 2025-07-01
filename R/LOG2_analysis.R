#' Analysis: Logarithmic quadratic
#'
#' This function performs logarithmic quadratic regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
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
#' @param linetype line type
#' @param pointshape format point (default is 21)
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @param print.on Print output
#'
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details
#' The logarithmic model is defined by:
#' \deqn{y = \beta_0 + \beta_1 ln(\cdot x) + \beta_2 ln(\cdot x)^2}
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @export
#'
#' @examples
#' library(AgroReg)
#' resp=c(10,8,6.8,6,5,4.3,4.1,4.2,4.1)
#' trat=seq(1,9,1)
#' LOG2(trat,resp)


LOG2=function(trat,
             resp,
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
             linetype=1,
             pointshape = 21,
             fillshape = "gray",
             colorline = "black",
             round=NA,
             xname.formula="x",
             yname.formula="y",
             comment=NA,
             fontfamily="sans",
             print.on=TRUE){
  requireNamespace("drc")
  requireNamespace("ggplot2")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  model <- lm(resp ~ log(trat)+I(log(trat)^2))
  model1<- lm(ymean ~ log(xmean)+I(log(xmean)^2))
  coef=summary(model)

  if(is.na(round)==TRUE){
    b=coef$coefficients[,1][1]
    d=coef$coefficients[,1][2]
    e=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    b=round(coef$coefficients[,1][1],round)
    d=round(coef$coefficients[,1][2],round)
    e=round(coef$coefficients[,1][3],round)}
  if(r2=="all"){r2=coef$r.squared}
  if(r2=="mean"){
    coef1=summary(model1)
    r2=coef1$r.squared}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e %s %0.3e*ln(%s) %s %0.3e*ln(%s)^2 ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   d,
                   ifelse(b<0,"-","+"),
                   abs(b),
                   xname.formula,
                   ifelse(e<0,"-","+"),
                   abs(e),
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
                                                 linewidth=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}

  graph=graph+theme+geom_line(data=preditos,aes(x=x,
                                                y=y,color="black"),linewidth=linesize,lty=linetype)+
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
                "plot"=graph,
                "expression"=s,
                "xaxisp"=temp1,
                "yaxisp"=result,
                "trt"=trat,
                "resp"=resp,
                "desvio"=desvio,
                "model"=model)
  if(print.on==TRUE){print(graficos[1:3])}
  output=graficos
}
