#' Analysis: Valcam
#'
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description This function performs Valcam regression analysis.
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
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param fontfamily Font family
#'
#' @details
#' The Valcam model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x + \beta_2\cdot x^1.5 + \beta_3\cdot x^2}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @references Siqueira, V. C., Resende, O., & Chaves, T. H. (2013). Mathematical modelling of the drying of jatropha fruit: an empirical comparison. Revista Ciencia Agronomica, 44, 278-285.
#' @export
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' valcam(trat,resp)

valcam=function(trat,
                resp,
                error = "SE",
                ylab = "Dependent",
                xlab = "Independent",
                theme = theme_classic(),
                legend.position = "top",
                r2 = "mean",
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
                comment = NA,
                fontfamily="sans") {
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  theta.0 <- min(resp) * 0.5
  model <- lm(resp ~ trat+I(trat^1.5)+I(trat^2))
  coef=summary(model)

  if(is.na(round)==TRUE){
  a=coef$coefficients[,1][1]
  b=coef$coefficients[,1][2]
  c=coef$coefficients[,1][3]
  d=coef$coefficients[,1][4]}

  if(is.na(round)==FALSE){
    a=round(coef$coefficients[,1][1],round)
    b=round(coef$coefficients[,1][2],round)
    c=round(coef$coefficients[,1][3],round)
    d=round(coef$coefficients[,1][4],round)}

  modm=lm(ymean~xmean+I(xmean^1.5)+I(xmean^2))
  if(r2=="all"){r2=round(summary(model)$r.squared,2)}
  if(r2=="mean"){r2=round(summary(modm)$r.squared,2)}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s == %e %s %e * %s %s %e * %s^1.5 %s %0.e * %s^2 ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   coef(modm)[1],
                   ifelse(coef(modm)[2] >= 0, "+", "-"),
                   abs(coef(modm)[2]),
                   xname.formula,
                   ifelse(coef(modm)[3] >= 0, "+", "-"),
                   abs(coef(modm)[3]),
                   xname.formula,
                   ifelse(coef(modm)[4] >= 0, "+", "-"),
                   abs(coef(modm)[4]),
                   xname.formula,
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
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
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
