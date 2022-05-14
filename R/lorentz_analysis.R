#' Analysis: Lorentz
#'
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param npar number of parameters (lo3 or lo4)
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
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
#' @param fontfamily Font family
#'
#' @importFrom stats var
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The model to the three-parameter Lorentz is:
#' \deqn{y = frac{d}{1+b(x-e)^2}}
#' The model to the three-parameter Lorentz is:
#' \deqn{y = c+frac{d-c}{1+b(x-e)^2}}
#' @export
#' @author Model imported from the aomisc package (Onofri, 2020)
#' @author Gabriel Danilo Shimizu
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @references Onofri A. (2020) The broken bridge between biologists and statisticians: a blog and R package, Statforbiology, IT, web: https://www.statforbiology.com
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' x=time[length(time):1]
#' lorentz(x,WL)

lorentz=function(trat,
                 resp,
                 npar = "lo3",
                 sample.curve=1000,
                 ylab = "Dependent",
                 xlab = "Independent",
                 theme = theme_classic(),
                 error = "SE",
                 legend.position = "top",
                 r2 = "all",
                 point = "all",
                 width.bar = NA,
                 scale = "none",
                 textsize = 12,
                 pointsize = 4.5,
                 linesize = 0.8,
                 pointshape = 21,
                 round = NA,
                 xname.formula = "x",
                 yname.formula = "y",
                 comment = NA,
                 fontfamily = "sans"){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  ymean=tapply(resp,trat,mean)
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)

  DRC.lorentz.3 <- function(){
    lorentz.3.fun <- function(X, b, d, e){
      d / ( 1 + b * (X - e)^2)}
    fct <- function(x, parm){lorentz.3.fun(x, parm[,1], parm[,2], parm[,3])}
    ssfct <- function(data){
      x <- data[, 1]
      y <- data[, 2]
      d <- max(y)
      e <- x[which.max(y)]
      pseudoY <- ( d - y )/ y
      pseudoX <- (x - e)^2
      coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
      b <- coefs[1]
      start <- c(b, d, e)
      return( start )}
    names <- c("b", "d", "e")
    text <- "Three parameters"
    returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
    class(returnList) <- "drcMean"
    invisible(returnList)}


  DRC.lorentz.4 <- function(){
    lorentz.4.fun <- function(X, b, c, d, e){
      c + (d - c) / ( 1 + b * (X - e)^2)}
    fct <- function(x, parm) {
      lorentz.4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])}
    ssfct <- function(data){
      x <- data[, 1]
      y <- data[, 2]
      d <- max(y)
      c <- min(y) * 0.95
      e <- x[which.max(y)]
      pseudoY <- (d - y)/(y - c)
      pseudoX <- (x - e)^2
      coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
      b <- coefs[1]
      start <- c(b, c, d, e)
      return(start)}
    names <- c("b", "c", "d", "e")
    text <- "Four parameters"
    returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
    class(returnList) <- "drcMean"
    invisible(returnList)}
  if(npar=="lo3"){
    mod=drm(resp ~ trat, fct = DRC.lorentz.3())
    mod=nls(resp ~ d/(1+b*(trat-e)^2),
            start = list(b=coef(mod)[1],
                         d=coef(mod)[2],
                         e=coef(mod)[3]))
    model=mod
    coef=summary(mod)
    if(is.na(round)==TRUE){
      b=coef$coefficients[,1][1]
      d=coef$coefficients[,1][2]
      e=coef$coefficients[,1][3]}

    if(is.na(round)==FALSE){
      b=round(coef$coefficients[,1][1],round)
      d=round(coef$coefficients[,1][2],round)
      e=round(coef$coefficients[,1][3],round)}
    if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
    if(r2=="mean"){
      model1 <- nls(ymean ~ d/(1+b*(xmean-e)^2),
                    start = list(b=coef(mod)[1],
                                 d=coef(mod)[2],
                                 e=coef(mod)[3]))
      r2=1-deviance(model1)/deviance(lm(ymean~1))}

    r2=floor(r2*100)/100

    equation=sprintf("~~~%s==frac(%0.3e,1 %s %0.3e(%s %s %0.3e)^2) ~~~~~ italic(R^2) == %0.2f",
                     yname.formula,
                     d,
                     ifelse(b >= 0, "+", "-"),
                     abs(b),
                     xname.formula,
                     ifelse(e <= 0, "+", "-"),
                     abs(e),
                     r2)
    if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}}
  if(npar=="lo4"){
    mod=drm(resp ~ trat, fct = DRC.lorentz.4())
    mod=nls(resp ~ c+(d-c)/(1+b*(trat-e)^2),
            start = list(b=coef(mod)[1],
                         c=coef(mod)[2],
                         d=coef(mod)[3],
                         e=coef(mod)[4]))
    model=mod
    coef=summary(mod)

    if(is.na(round)==TRUE){
      b=coef$coefficients[,1][1]
      c=coef$coefficients[,1][2]
      d=coef$coefficients[,1][3]
      dc=d-c
      e=coef$coefficients[,1][4]}

    if(is.na(round)==FALSE){
      b=round(coef$coefficients[,1][1],round)
      c=round(coef$coefficients[,1][2],round)
      d=round(coef$coefficients[,1][3],round)
      dc=d-c
      e=round(coef$coefficients[,1][4],round)}
    if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
    if(r2=="mean"){
      model1 <- nls(ymean ~ c+(d-c)/(1+b*(xmean-e)^2),
                    start = list(b=coef(mod)[1],
                                 c=coef(mod)[2],
                                 d=coef(mod)[3],
                                 e=coef(mod)[4]))
      r2=1-deviance(model1)/deviance(lm(ymean~1))}

    r2=floor(r2*100)/100

    equation=sprintf("~~~%s==%0.3e + frac(%0.3e,1 %s %0.3e(%s %s %0.3e)^2) ~~~~~ italic(R^2) == %0.2f",
                     yname.formula,
                     c,
                     dc,
                     ifelse(b >= 0, "+", "-"),
                     abs(b),
                     xname.formula,
                     ifelse(e >= 0, "+", "-"),
                     abs(e),
                     r2)
    if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}}
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
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
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                                 width=width.bar,
                                                 size=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
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
  result=predict(mod,newdata = data.frame(trat=temp1),type="response")
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
