#' Analysis: Steinhart-Hart
#'
#' The Steinhart-Hart model. The Steinhart-Hart equation is a model used to explain the behavior of a semiconductor at different temperatures, however, Zhai et al. (2020) used this model to relate plant density and grain yield.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param initial Starting estimates
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position Legend position (\emph{default} is "top")
#' @param r2 Coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param point Defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	Shape size
#' @param linesize	Line size
#' @param pointshape Format point (default is 21)
#' @param round round equation
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The model function for the Steinhart-Hart model is:
#' \deqn{ y = \frac{1}{A+B \times ln(x)+C \times ln(x)^3}}
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Zhai, L., Li, H., Song, S., Zhai, L., Ming, B., Li, S., ... & Zhang, L. (2021). Intra-specific competition affects the density tolerance and grain yield of maize hybrids. Agronomy Journal, 113(1), 224-23. doi:10.1002/agj2.20438
#' @seealso \link{LL}, \link{CD},\link{GP}
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' SH(trat,resp)

SH=function(trat,
            resp,
            initial=NA,
            sample.curve=1000,
            ylab="Dependent",
            xlab="Independent",
            theme=theme_classic(),
            legend.position="top",
            r2="all",
            error="SE",
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
            yname.formula="y",
            xname.formula="x",
            comment=NA,
            fontfamily="sans"){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  ymean=tapply(resp,trat,mean)
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd

  if(is.na(initial[1])==TRUE){
    ix=1/max(resp)
    initial=list(A=ix,B=ix,C=ix)}
  mod=nls(resp~1/(A+B*log(trat)+C*log(trat)^3),start = initial)
  model=mod
  coef=summary(mod)

  if(is.na(round)==TRUE){
  A=coef$coefficients[,1][1]
  B=coef$coefficients[,1][2]
  C=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    A=round(coef$coefficients[,1][1],round)
    B=round(coef$coefficients[,1][2],round)
    C=round(coef$coefficients[,1][3],round)}


  if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
  if(r2=="mean"){
    model1=nls(ymean~1/(A+B*log(xmean)+C*log(xmean)^3),start = initial)
    r2=1-deviance(model1)/deviance(lm(ymean~1))}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==frac(1, %0.3e %s %0.3e*ln(%s) %s %0.3e*ln(%s)^3) ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   A,
                   ifelse(B >= 0, "+", "-"),
                   abs(B),
                   xname.formula,
                   ifelse(C >= 0, "+", "-"),
                   C,
                   xname.formula,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  x=preditos$x
  y=preditos$y
  s=equation
  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  if(point=="mean"){
    graph=ggplot(data,aes(x=xmean,y=ymean))
    if(error!="FALSE"){graph=graph+
      geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),size=linesize,
                    width=width.bar)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
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
  temp1=seq(min(trat),max(trat),length.out=sample.curve)
  result=predict(mod,newdata = data.frame(trat=temp1),
                 type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  aic=AIC(mod)
  bic=BIC(mod)
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
