#' Analysis: Asymptotic Regression
#'
#' This function performs Asymptotic regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param width.bar	Bar width
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param comment Add text after equation
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @seealso \link{biexponential}, \link{exponential_neg}
#' @details
#' The Asymptotic model is defined by:
#' \deqn{y=Asym+(R0-Asym) \times e^{-e^{lrc} \times x}}
#' Asym: a numeric parameter representing the horizontal asymptote on the right side.
#' R0	a numeric parameter representing the response when input is zero.
#' lrc a numeric parameter representing the natural logarithm of the rate constant.
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley \& Sons (p. 330).
#' @export
#' @importFrom stats SSasymp
#'
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' asymptotic(time,100-WL)

asymptotic =function(trat,
                     resp,
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
                     pointshape = 21,
                     comment=NA){
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  model <- nls(resp ~ SSasymp(trat, Asym, R0, lrc))
  coef=summary(model)
  Asym=coef$coefficients[,1][1]
  R0=coef$coefficients[,1][2]
  lrc=coef$coefficients[,1][3]
  if(r2=="all"){r2=cor(resp, fitted(model))^2}
  if(r2=="mean"){r2=cor(ymean, predict(model,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  R0Asym=R0-Asym
  equation=sprintf("~~~y==%0.3e %s %0.3e*e^{-e^{%0.3e}*x} ~~~~~ italic(R^2) == %0.2f",
                   Asym,
                   ifelse(R0Asym >= 0, "+", "-"),
                   abs(R0Asym),
                   lrc,r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
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
                                                 width=0.5)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  # if(point=="ic"){
  #   ic1=confint(mod)
  #   lower=ic1[1,1]+(ic1[2,1]-ic1[1,1])*exp(-exp(ic1[3,1])*xp)
  #   upper=ic1[1,2]+(ic1[2,2]-ic1[1,2])*exp(-exp(ic1[3,2])*xp)
  #   graph=ggplot(data,aes(x=xmean,y=ymean))
  #   graph=graph+geom_polygon(data=data.frame(x=c(xp,xp[length(xp):1]),
  #                                y=c(lower,upper[length(xp):1])),
  #                aes(x=x,y=y),fill="gray")
  #   graph=graph+
  #     geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")}
  # if(point=="all_ic"){
  #   ic1=confint(mod)
  #   lower=ic1[1,1]+(ic1[2,1]-ic1[1,1])*exp(-exp(ic1[3,1])*xp)
  #   upper=ic1[1,2]+(ic1[2,2]-ic1[1,2])*exp(-exp(ic1[3,2])*xp)
  #   graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
  #   graph=graph+geom_polygon(data=data.frame(x=c(xp,xp[length(xp):1]),
  #                                            y=c(lower,upper[length(xp):1])),
  #                            aes(x=x,y=y),fill="gray")
  #   graph=graph+
  #     geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}

  graph=graph+theme+geom_line(data=preditos,aes(x=x,
                                                y=y,color="black"),size=linesize)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=textsize),
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
