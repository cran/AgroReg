#' Analysis: Von Bertalanffy
#'
#' The Von Bertalanffy model. It's a kind of growth curve for a time series and takes its name from its creator, Ludwig von Bertalanffy. It is a special case of the generalized logistic function. The growth curve (biology) is used to model the average length from age in animals.
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
#' @param linetype line type
#' @param pointshape Format point (default is 21)
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The model function for the von Bertalanffy model is:
#' \deqn{ y = L(1-exp(-k(t-t0)))}
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @examples
#' library(AgroReg)
#' x=seq(1,20)
#' y=c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.91,
#'     0.92, 0.94, 0.96, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00)
#' VB(x,y)

VB=function(trat,
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
            linetype=1,
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
    L=max(resp)
    k=1
    t0=1
    initial=list(L=L,k=k,t0=t0)}
  mod=nls(resp~L*(1-exp(-k*(trat-t0))),start = initial)
  model=mod
  coef=summary(mod)

  if(is.na(round)==TRUE){
    L=coef$coefficients[,1][1]
    k=coef$coefficients[,1][2]
    t0=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    L=round(coef$coefficients[,1][1],round)
    k=round(coef$coefficients[,1][2],round)
    t0=round(coef$coefficients[,1][3],round)}

  if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
  if(r2=="mean"){
    model1=nls(ymean~L*(1-exp(-k*(xmean-t0))),start = initial)
    r2=1-deviance(model1)/deviance(lm(ymean~1))}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e*(1-e^{%s %0.3e*(%s %s %0.3e)}) ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   L,
                   ifelse(k<0,"+","-"),
                   abs(k),
                   xname.formula,
                   ifelse(t0<0,"+","-"),
                   abs(t0),
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
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
                                color="black"),size=linesize,lty=linetype)+
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
