#' Analysis: Michaelis-Menten
#'
#' This function performs regression analysis using the Michaelis-Menten model.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param npar Number of parameters (mm2 or mm3)
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param ic Add interval of confidence
#' @param fill.ic Color interval of confidence
#' @param alpha.ic confidence interval transparency level
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
#' @param fontfamily Font family
#'
#' @details
#' The two-parameter Michaelis-Menten model is defined by:
#' \deqn{y = \frac{Vm \times x}{k + x}}
#' The three-parameter Michaelis-Menten model is defined by:
#' \deqn{y = c + \frac{Vm \times x}{k + x}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Gabriel Danilo Shimizu
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @examples
#' data("granada")
#' attach(granada)
#' MM(time,WL)
#' MM(time,WL,npar="mm3")
#' @md

MM=function(trat,
            resp,
            npar="mm2",
            sample.curve=1000,
            error="SE",
            ylab="Dependent",
            xlab="Independent",
            theme=theme_classic(),
            legend.position="top",
            point="all",
            width.bar=NA,
            r2="all",
            ic=FALSE,
            fill.ic="gray70",
            alpha.ic=0.5,
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
  requireNamespace("drc")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd

  if(npar=="mm2"){
    mod=drm(resp~trat,fct=MM.2())
    mod=nls(resp~((trat*Vm)/(k+trat)),start = list(Vm=coef(mod)[1],
                                                   k=coef(mod)[2]))
  model=mod
  coef=summary(mod)

  if(is.na(round)==TRUE){
  d=coef$coefficients[,1][1]
  e=coef$coefficients[,1][2]}

  if(is.na(round)==FALSE){
    d=round(coef$coefficients[,1][1],round)
    e=round(coef$coefficients[,1][2],round)}

  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==frac(%0.3e*%s, %0.3e + %s)~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   d,
                   xname.formula,
                   e,
                   xname.formula,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  x=preditos$x
  y=preditos$y
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  mod=drm(resp~trat,fct=MM.2())}

  if(npar=="mm3"){
    mod=drm(resp~trat,fct=MM.3())
    mod=nls(resp~((trat*Vm)/(k+trat))+c,start = list(c=coef(mod)[1],
                                                    Vm=coef(mod)[2],
                                                    k=coef(mod)[3]))
    model=mod
    coef=summary(mod)
    if(is.na(round)==TRUE){
      c=coef$coefficients[,1][1]
      d=coef$coefficients[,1][2]
      e=coef$coefficients[,1][3]}

    if(is.na(round)==FALSE){
      c=round(coef$coefficients[,1][1],round)
      d=round(coef$coefficients[,1][2],round)
      e=round(coef$coefficients[,1][3],round)}

    if(r2=="all"){r2=cor(resp, fitted(mod))^2}
    if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
    r2=floor(r2*100)/100
    equation=sprintf("~~~%s== %0.3e + frac(%0.3e*%s, %0.3e + %s) ~~~~~ italic(R^2) == %0.2f",
                     yname.formula,
                     c,
                     d,
                     xname.formula,
                     e,
                     xname.formula,
                     r2)
    xp=seq(min(trat),max(trat),length.out = sample.curve)
    preditos=data.frame(x=xp,
                        y=predict(mod,newdata = data.frame(trat=xp)))
    x=preditos$x
    y=preditos$y
    if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
    mod=drm(resp~trat,fct=MM.3())}
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
  if(ic==TRUE){
    pred=data.frame(x=xp,
                    y=predict(mod,interval = "confidence",newdata = data.frame(trat=xp)))
    preditosic=data.frame(x=c(pred$x,pred$x[length(pred$x):1]),
                          y=c(pred$y.Lower,pred$y.Upper[length(pred$x):1]))
    graph=graph+geom_polygon(data=preditosic,aes(x=x,y),fill=fill.ic,alpha=alpha.ic)}

  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,y=y,color="black"), size=linesize)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  temp1=seq(min(trat),max(trat),length.out=sample.curve)
  result=predict(mod,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  aic=AIC(mod)
  bic=BIC(mod)
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
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
