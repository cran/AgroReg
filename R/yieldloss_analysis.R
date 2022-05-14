#' Analysis: Yield-loss
#'
#' This function performs regression analysis using the Yield loss model.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
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
#' @details
#' The Yield Loss model is defined by:
#' \deqn{y = \frac{i \times x}{1+\frac{i}{A} \times x}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Model imported from the aomisc package (Onofri, 2020)
#' @author Gabriel Danilo Shimizu
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @references Onofri A. (2020) The broken bridge between biologists and statisticians: a blog and R package, Statforbiology, IT, web: https://www.statforbiology.com
#' @examples
#' data("granada")
#' attach(granada)
#' yieldloss(time,WL)
#' @md

yieldloss=function(trat,
            resp,
            sample.curve=1000,
            error="SE",
            ylab="Dependent",
            xlab="Independent",
            theme=theme_classic(),
            legend.position="top",
            point="all",
            width.bar=NA,
            r2="all",
            textsize = 12,
            pointsize = 4.5,
            linesize = 0.8,
            pointshape = 21,
            round=NA,
            yname.formula="y",
            xname.formula="x",
            comment=NA,
            scale="none",
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
  DRC.YL <- function(fixed = c(NA, NA), names = c("i", "A")) {
    YL.fun <- function(predictor, i, A) {
      i * predictor/(1 + i/A * predictor)}
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    fct <- function(x, parm) {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      i <- parmMat[, 1]; A <- parmMat[, 2]
      YL.fun(x, i, A)}
    ssfct <- function(dataf) {
      x <- dataf[, 1]
      y <- dataf[, 2]
      pseudoY <- 1 /  y[x > 0]
      pseudoX <- 1 / x [x > 0]
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      A <- 1 / coefs[1]; i <- 1 / coefs[2]
      return(c(i, A)[notFixed])}
    pnames <- names[notFixed]
    text <- "Yield-Loss function (Cousens, 1985)"
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))
    class(returnList) <- "drcMean"
    invisible(returnList)}
    mod=drm(resp~trat,fct=DRC.YL())
    mod=nls(resp~i * trat/(1 + i/A * trat),
            start = list(i=coef(mod)[1],
                         A=coef(mod)[2]))
    model=mod
    coef=summary(mod)
    if(is.na(round)==TRUE){
      b=coef$coefficients[,1][1]
      c=coef$coefficients[,1][2]
      bc=b/c}

    if(is.na(round)==FALSE){
      b=round(coef$coefficients[,1][1],round)
      c=round(coef$coefficients[,1][2],round)
      bc=b/c}

    if(r2=="all"){r2=cor(resp, fitted(mod))^2}
    if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
    r2=floor(r2*100)/100
    equation=sprintf("~~~%s== frac(%0.3e * %s, 1 %s %0.3e * %s) ~~~~~ italic(R^2) == %0.2f",
                     yname.formula,
                     b,
                     xname.formula,
                     ifelse(bc <= 0, "-","+"),
                     bc,
                     xname.formula,
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
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,y=y,color="black"))+
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
  aic=AIC(mod)
  bic=BIC(mod)
  print(coef)
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  graphs=data.frame("Parameter"=c("AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(aic,
                               bic,
                               r2,
                               rmse))
  graficos=list("teste"=graphs,graph)
  print(graficos)
}
