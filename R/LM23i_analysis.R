#' Analysis: Cubic inverse without beta1
#'
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Degree 3 polynomial inverse model without the beta 1 coefficient.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param legend.position legend position (\emph{default} is "top")
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
#' @details
#' Inverse degree 3 polynomial model without the beta 1 coefficient  is defined by:
#' \deqn{y = \beta_0 + \beta_2\cdot x^{1/2} + \beta_3\cdot x^{1/3}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @keywords regression linear
#' @export
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' LM23i(time, WL)

LM23i=function(trat,
               resp,
               sample.curve=1000,
               ylab="Dependent",
               error="SE",
               xlab="Independent",
               theme=theme_classic(),
               legend.position="top",
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
  requireNamespace("ggplot2")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  dose=tapply(trat, trat, mean)
  moda=lm(resp~I(trat^(1/2))+I(trat^(1/3)))
  mods=summary(moda)$coefficients
  modm=lm(media~I(dose^(1/2))+I(dose^(1/3)))
  if(r2=="mean"){r2=round(summary(modm)$r.squared,2)}
  if(r2=="all"){r2=round(summary(moda)$r.squared,2)}
  coef1=coef(moda)[1]
  coef2=coef(moda)[2]
  coef3=coef(moda)[3]
  s1 = s =  sprintf("~~~%s == %e %s %e*%s^{1/2} %s %e*%s^{1/3} ~~~~~ italic(R^2) == %0.2f",
               yname.formula,
               coef1,
               ifelse(coef2 >= 0, "+", "-"),
               abs(coef2),
               xname.formula,
               ifelse(coef3 >= 0, "+", "-"),
               abs(coef3),
               xname.formula,
               r2)
  if(is.na(comment)==FALSE){s1=paste(s1,"~\"",comment,"\"")}
  data1=data.frame(trat=unique(trat),
                   media=media,
                   resp=media,
                   desvio)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(moda,newdata = data.frame(trat=xp)))
  x=preditos$x
  y=preditos$y
  predesp=predict(moda)
  predobs=resp
  if(point=="mean"){
    graph=ggplot(data1,aes(x=trat,y=media))
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=media-desvio,
                                                     ymax=media+desvio),
                                                 width=width.bar,
                                                 linewidth=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  grafico=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,color="black"),linewidth=linesize,lty=linetype)+
    scale_color_manual(name="",values=colorline,label=parse(text = s1))+
    theme(axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){grafico=grafico+scale_x_log10()}
  models=mods
  model=moda
  r2=summary(modm)$r.squared
  aic=AIC(moda)
  bic=BIC(moda)
  vif.test=function (mod){
    if (any(is.na(coef(mod))))
      stop("there are aliased coefficients in the model")
    v <- vcov(mod)
    assign <- attr(model.matrix(mod), "assign")
    if (names(coefficients(mod)[1]) == "(Intercept)") {
      v <- v[-1, -1]
      assign <- assign[-1]}
    else warning("No intercept: vifs may not be sensible.")
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2)
      stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    for (term in 1:n.terms) {
      subs <- which(assign == term)
      result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs,
                                                                         -subs]))/detR
      result[term, 2] <- length(subs)}
    if (all(result[, 2] == 1))
      result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result}
  vif=vif.test(moda)
  predesp=predict(moda)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  temp1=seq(min(trat),max(trat),length.out=sample.curve)
  result=predict(moda,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
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
                "plot"=grafico,
                "VIF"=vif,
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
