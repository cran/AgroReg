#' Analysis: Linear, quadratic, quadratic inverse, cubic and quartic
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Linear, quadratic, quadratic inverse, cubic and quartic regression.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param grau degree of the polynomial (0.5, 1, 2, 3 or 4)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ic Add interval of confidence
#' @param fill.ic Color interval of confidence
#' @param alpha.ic confidence interval transparency level
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
#' The linear model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x}
#' The quadratic model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x + \beta_2\cdot x^2}
#' The quadratic inverse model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x + \beta_2\cdot x^{0.5}}
#' The cubic model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x + \beta_2\cdot x^2 + \beta_3\cdot x^3}
#' The quartic model is defined by:
#' \deqn{y = \beta_0 + \beta_1\cdot x + \beta_2\cdot x^2 + \beta_3\cdot x^3+ \beta_4\cdot x^4}
#'
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @keywords regression linear
#' @export
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' LM(trat,resp, grau = 3)

LM=function(trat,
            resp,
            grau=NA,
            ylab="Dependent",
            xlab="Independent",
            error="SE",
            ic=FALSE,
            fill.ic="gray70",
            alpha.ic=0.5,
            point="all",
            theme=theme_classic(),
            legend.position="top",
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
  if(is.na(grau)==TRUE){grau=1}
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  mod=c()
  mod1=c()
  mod2=c()
  mod05=c()

  modm=c()
  mod1m=c()
  mod2m=c()
  mod05m=c()

  text1=c()
  text2=c()
  text3=c()
  text05=c()

  mods=c()
  mod1s=c()
  mod2s=c()
  mod05s=c()

  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  dose=tapply(trat, trat, mean)
  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mod3a=lm(resp~trat+I(trat^2)+I(trat^3)+I(trat^4))
  mod05a=lm(resp~trat+I(trat^0.5))

  mods=summary(moda)$coefficients
  mod1s=summary(mod1a)$coefficients
  mod2s=summary(mod2a)$coefficients
  mod3s=summary(mod3a)$coefficients
  mod05s=summary(mod05a)$coefficients

  modm=lm(media~dose)
  mod1m=lm(media~dose+I(dose^2))
  mod2m=lm(media~dose+I(dose^2)+I(dose^3))
  mod3m=lm(media~dose+I(dose^2)+I(dose^3)+I(dose^4))
  mod05m=lm(media~dose+I(dose^0.5))

  if(grau=="1"){r2=round(summary(modm)$r.squared,2)}
  if(grau=="2"){r2=round(summary(mod1m)$r.squared,2)}
  if(grau=="3"){r2=round(summary(mod2m)$r.squared,2)}
  if(grau=="4"){r2=round(summary(mod3m)$r.squared,2)}
  if(grau=="0.5"){r2=round(summary(mod05m)$r.squared,2)}
  if(grau=="1"){
    if(is.na(round)==TRUE){
    coef1=coef(moda)[1]
    coef2=coef(moda)[2]}
    if(is.na(round)==FALSE){
      coef1=round(coef(moda)[1],round)
      coef2=round(coef(moda)[2],round)}
    s1=s <- sprintf("~~~%s == %e %s %e* %s ~~~~~ italic(R^2) == %0.2f",
                    yname.formula,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.formula,
                    r2)
  if(is.na(comment)==FALSE){s1=paste(s1,"~\"",comment,"\"")}}
  if(grau=="2"){
    if(is.na(round)==TRUE){
      coef1=coef(mod1a)[1]
      coef2=coef(mod1a)[2]
      coef3=coef(mod1a)[3]}
    if(is.na(round)==FALSE){
      coef1=round(coef(mod1a)[1],round)
      coef2=round(coef(mod1a)[2],round)
      coef3=round(coef(mod1a)[3],round)}
    s2=s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^2 ~~~~~ italic(R^2) ==  %0.2f",
                    yname.formula,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.formula,
                    ifelse(coef3 >= 0, "+", "-"),
                    abs(coef3),
                    xname.formula,
                    r2)
  if(is.na(comment)==FALSE){s2=paste(s2,"~\"",comment,"\"")}}
  if(grau=="3"){
    if(is.na(round)==TRUE){
      coef1=coef(mod2a)[1]
      coef2=coef(mod2a)[2]
      coef3=coef(mod2a)[3]
      coef4=coef(mod2a)[4]}
    if(is.na(round)==FALSE){
      coef1=round(coef(mod2a)[1],round)
      coef2=round(coef(mod2a)[2],round)
      coef3=round(coef(mod2a)[3],round)
      coef4=round(coef(mod2a)[4],round)}
    s3=s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^2 %s %0.e * %s^3 ~~~~~ italic(R^2) == %0.2f",
                    yname.formula,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.formula,
                    ifelse(coef3 >= 0, "+", "-"),
                    abs(coef3),
                    xname.formula,
                    ifelse(coef4 >= 0, "+", "-"),
                    abs(coef4),
                    xname.formula,
                    r2)
  if(is.na(comment)==FALSE){s3=paste(s3,"~\"",comment,"\"")}}
  if(grau=="4"){
    if(is.na(round)==TRUE){
      coef1=coef(mod3a)[1]
      coef2=coef(mod3a)[2]
      coef3=coef(mod3a)[3]
      coef4=coef(mod3a)[4]
      coef5=coef(mod3a)[5]}
    if(is.na(round)==FALSE){
      coef1=round(coef(mod3a)[1],round)
      coef2=round(coef(mod3a)[2],round)
      coef3=round(coef(mod3a)[3],round)
      coef4=round(coef(mod3a)[4],round)
      coef5=round(coef(mod3a)[5],round)}
    s4=s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^2 %s %0.e * %s^3 %s %0.e * %s^4 ~~~~~ italic(R^2) == %0.2f",
                    yname.formula,
                    coef1,
                    ifelse(coef2 >= 0, "+", "-"),
                    abs(coef2),
                    xname.formula,
                    ifelse(coef3 >= 0, "+", "-"),
                    abs(coef3),
                    xname.formula,
                    ifelse(coef4 >= 0, "+", "-"),
                    abs(coef4),
                    xname.formula,
                    ifelse(coef5 >= 0, "+", "-"),
                    abs(coef5),
                    xname.formula,
                    r2)
  if(is.na(comment)==FALSE){s4=paste(s4,"~\"",comment,"\"")}}
  if(grau=="0.5"){
    if(is.na(round)==TRUE){
      coef1=coef(mod05a)[1]
      coef2=coef(mod05a)[2]
      coef3=coef(mod05a)[3]}
    if(is.na(round)==FALSE){
      coef1=round(coef(mod05a)[1],round)
      coef2=round(coef(mod05a)[2],round)
      coef3=round(coef(mod05a)[3],round)}
    s05=s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^0.5 ~~~~~ italic(R^2) ==  %0.2f",
                     yname.formula,
                     coef1,
                     ifelse(coef2 >= 0, "+", "-"),
                     abs(coef2),
                     xname.formula,
                     ifelse(coef3 >= 0, "+", "-"),
                     abs(coef3),
                     xname.formula,
                     r2)
  if(is.na(comment)==FALSE){s05=paste(s05,"~\"",comment,"\"")}}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=unique(trat),
                   resp=media,
                   desvio)
  if(point=="mean"){
    grafico=ggplot(data1,aes(x=trat,y=resp))+
      geom_point(aes(fill=as.factor(rep(1,length(resp)))),na.rm=TRUE,
                 size=pointsize,color="black",shape=pointshape)
    if(error!="FALSE"){grafico=grafico+
      geom_errorbar(aes(ymin=resp-desvio, ymax=resp+desvio),
                    width=width.bar,
                    size=linesize)}}
  if(point=="all"){
    grafico=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))+
      geom_point(aes(fill=as.factor(rep(1,length(resp)))),size=pointsize,
                 shape=pointshape,color="black")}

  grafico=grafico+
    theme+ylab(ylab)+xlab(xlab)
  if(grau=="1"){grafico=grafico+geom_smooth(method = "lm",se=ic, fill=fill.ic, alpha=alpha.ic, na.rm=TRUE, formula = y~x,size=linesize,color="black")}
  if(grau=="2"){grafico=grafico+geom_smooth(method = "lm",se=ic, fill=fill.ic, alpha=alpha.ic,na.rm=TRUE, formula = y~x+I(x^2),size=linesize,color="black")}
  if(grau=="3"){grafico=grafico+geom_smooth(method = "lm",se=ic, fill=fill.ic, alpha=alpha.ic,na.rm=TRUE, formula = y~x+I(x^2)+I(x^3),size=linesize,color="black")}
  if(grau=="4"){grafico=grafico+geom_smooth(method = "lm",se=ic, fill=fill.ic, alpha=alpha.ic,na.rm=TRUE, formula = y~x+I(x^2)+I(x^3)+I(x^4),size=linesize,color="black")}
  if(grau=="0.5"){grafico=grafico+geom_smooth(method = "lm",se=ic, fill=fill.ic, alpha=alpha.ic,na.rm=TRUE, formula = y~x+I(x^0.5),size=linesize,color="black")}
  if(grau=="1"){grafico=grafico+
    scale_fill_manual(values="gray",label=c(parse(text=s1)),name="")}
  if(grau=="2"){grafico=grafico+
      scale_fill_manual(values="gray",label=c(parse(text=s2)),name="")}
  if(grau=="3"){grafico=grafico+
      scale_fill_manual(values="gray",label=c(parse(text=s3)),name="")}
  if(grau=="4"){grafico=grafico+
    scale_fill_manual(values="gray",label=c(parse(text=s4)),name="")}
  if(grau=="0.5"){grafico=grafico+
    scale_fill_manual(values="gray",label=c(parse(text=s05)),name="")}

  grafico=grafico+
    theme(text = element_text(size=textsize,color="black"),
          axis.text = element_text(size=textsize,color="black"),
          axis.title = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)
  if(scale=="log"){grafico=grafico+scale_x_log10()}

  moda=lm(resp~trat)
  mod1a=lm(resp~trat+I(trat^2))
  mod2a=lm(resp~trat+I(trat^2)+I(trat^3))
  mod3a=lm(resp~trat+I(trat^2)+I(trat^3)+I(trat^4))
  mod05a=lm(resp~trat+I(trat^0.5))

  if(grau=="1"){
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
  }

  if(grau=="2"){
  models=mod1s
  model=mod1a
  r2=summary(mod1m)$r.squared
  aic=AIC(mod1a)
  bic=BIC(mod1a)
  vif=car::vif(mod1a)
  predesp=predict(mod1a)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))

  temp1=seq(min(trat),max(trat),length.out=5000)
  result=predict(mod1a,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  }


  if(grau=="3"){
  models=mod2s
  model=mod2a
  r2=summary(mod2m)$r.squared
  aic=AIC(mod2a)
  bic=BIC(mod2a)
  vif=car::vif(mod2a)
  predesp=predict(mod2a)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))

  temp1=seq(min(trat),max(trat),length.out=5000)
  result=predict(mod2a,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]}

  if(grau=="4"){
    models=mod3s
    model=mod3a
    r2=summary(mod3m)$r.squared
    aic=AIC(mod3a)
    bic=BIC(mod3a)
    vif=car::vif(mod3a)
    predesp=predict(mod3a)
    predobs=resp
    rmse=sqrt(mean((predesp-predobs)^2))
    temp1=seq(min(trat),max(trat),length.out=5000)
    result=predict(mod3a,newdata = data.frame(trat=temp1),type="response")
    maximo=temp1[which.max(result)]
    respmax=result[which.max(result)]
    minimo=temp1[which.min(result)]
    respmin=result[which.min(result)]}

  if(grau=="0.5"){
    models=mod05s
    model=mod05a
    r2=summary(mod05m)$r.squared
    aic=AIC(mod05a)
    bic=BIC(mod05a)
    vif=car::vif(mod05a)
    predesp=predict(mod05a)
    predobs=resp
    rmse=sqrt(mean((predesp-predobs)^2))
    temp1=seq(min(trat),max(trat),length.out=5000)
    result=predict(mod05a,newdata = data.frame(trat=temp1),type="response")
    maximo=temp1[which.max(result)]
    respmax=result[which.max(result)]
    minimo=temp1[which.min(result)]
    respmin=result[which.min(result)]
    }
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
