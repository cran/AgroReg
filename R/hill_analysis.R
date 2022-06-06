#' Analysis: Hill
#'
#' This function performs regression analysis using the Hill model.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @details
#' The Hill model is defined by:
#' \deqn{y = \frac{a \times x^c}{b+x^c}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Model imported from the aomisc package (Onofri, 2020)
#' @author Gabriel Danilo Shimizu
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @references Onofri A. (2020) The broken bridge between biologists and statisticians: a blog and R package, Statforbiology, IT, web: https://www.statforbiology.com
#' @examples
#' data("granada")
#' attach(granada)
#' hill(time,WL)
#' @md

hill=function(trat,
              resp,
              sample.curve=1000,
              error = "SE",
              ylab = "Dependent",
              xlab = "Independent",
              theme = theme_classic(),
              legend.position = "top",
              point = "all",
              width.bar = NA,
              r2 = "all",
              textsize = 12,
              pointsize = 4.5,
              linesize = 0.8,
              pointshape = 21,
              fillshape = "gray",
              colorline = "black",
              round = NA,
              xname.formula="x",
              yname.formula="y",
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
  DRC.hill<-function(fixed = c(NA, NA, NA), names = c("a", "b", "c")){
    hillCurveMean <- function(predictor, a, b, c) {
      (a * predictor^c)/(b + predictor^c)}
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    fct <- function(x, parm){
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
      (a * x^c)/(b + x^c)}
    ssfct <- function(dataf){
      x <- dataf[, 1]
      y <- dataf[, 2]
      a <- max(y) * 1.05
      pseudoY <-  log(( a - y )/ y)
      pseudoX <- log(x)
      coefs <- coef( lm(pseudoY ~ pseudoX ) )
      b <- exp(coefs[1])
      c <- - coefs[2]
      return(c(a, b, c)[notFixed])}
    pnames <- names[notFixed]
    text <- "Hill function (Morgan-Mercer-Flodin)"
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    class(returnList) <- "drcMean"
    invisible(returnList)}
  mod=drm(resp~trat,fct=DRC.hill())
  mod=nls(resp~(a * trat^c)/(b + trat^c),
          start = list(a=coef(mod)[1],
                       b=coef(mod)[2],
                       c=coef(mod)[3]))
  model=mod
  coef=summary(mod)
  if(is.na(round)==TRUE){
    a=coef$coefficients[,1][1]
    b=coef$coefficients[,1][2]
    c=coef$coefficients[,1][3]}

  if(is.na(round)==FALSE){
    a=round(coef$coefficients[,1][1],round)
    b=round(coef$coefficients[,1][2],round)
    c=round(coef$coefficients[,1][3],round)}
  if(r2=="all"){r2=1-deviance(model)/deviance(lm(resp~1))}
  if(r2=="mean"){
    model1 <- nls(ymean~(a * xmean^c)/(b + xmean^c),
                  start = list(a=coef(mod)[1],
                               b=coef(mod)[2],
                               c=coef(mod)[3]))
    r2=1-deviance(model1)/deviance(lm(ymean~1))}

  r2=floor(r2*100)/100
  equation=sprintf("~~~%s== frac(%0.3e * %s ^ %0.3e, %0.3e + %s ^%0.3e) ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   a,
                   xname.formula,
                   c,
                   b,
                   xname.formula,
                   c,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
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
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}

  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,y=y,color="black"))+
    scale_color_manual(name="",values=colorline,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize, family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
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
