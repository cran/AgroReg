#' Analysis: Weibull
#'
#' @description  The w3' and 'w4' logistical models provide Weibull. This model was extracted from the 'drc' package.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param npar Number of model parameters (\emph{default} is  w3)
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position Legend position (\emph{default} is "top")
#' @param r2 Coefficient of determination of the mean or all values (\emph{default} is all)
#' @param ic Add interval of confidence
#' @param fill.ic Color interval of confidence
#' @param alpha.ic confidence interval transparency level
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
#' @param print.on Print output
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The three-parameter Weibull model is given by the expression
#' \deqn{y = d\exp(-\exp(b(\log(x)-e)))}
#' Fixing the lower limit at 0 yields the four-parameter model
#' \deqn{y = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e)))))}
#' @export
#' @author Model imported from the drc package (Ritz et al., 2016)
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @references Ritz, C.; Strebig, J.C. and Ritz, M.C. Package ‘drc’. Creative Commons: Mountain View, CA, USA, 2016.
#' @seealso \link{LL}, \link{CD},\link{GP}
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' weibull(trat,resp)

weibull=function(trat,
            resp,
            npar="w3",
            sample.curve=1000,
            ylab="Dependent",
            xlab="Independent",
            theme=theme_classic(),
            legend.position="top",
            r2="all",
            ic=FALSE,
            fill.ic="gray70",
            alpha.ic=0.5,
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
            fontfamily="sans",
            print.on=TRUE){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  ymean=tapply(resp,trat,mean)
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd
  if(npar=="w3"){mod=drm(resp~trat,fct=w3())
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

  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s==%0.3e * e^{-e^{%0.3e*(log(%s)-log(%0.3e))}} ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   d,
                   b,
                   xname.formula,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  }
  if(npar=="w4"){mod=drm(resp~trat,fct=w4())
  model=mod
  coef=summary(mod)

  if(is.na(round)==TRUE){
    b=coef$coefficients[,1][1]
    c=coef$coefficients[,1][2]
    d=coef$coefficients[,1][3]
    e=coef$coefficients[,1][4]
    cd=d-c}

  if(is.na(round)==FALSE){
    b=round(coef$coefficients[,1][1],round)
    c=round(coef$coefficients[,1][2],round)
    d=round(coef$coefficients[,1][3],round)
    e=round(coef$coefficients[,1][4],round)
    cd=d-c}
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~%s == %0.3e %s %0.3e * e^{-e^{%0.3e*(log(%s)-log(%0.3e))}} ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   c,
                   ifelse(cd<0,"-","+"),
                   abs(cd),
                   b,
                   xname.formula,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))}
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
      geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),linewidth=linesize,
                    width=width.bar)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(ic==TRUE){
    pred=data.frame(x=xp,
                    y=predict(mod,interval = "confidence",newdata = data.frame(trat=xp)))
    preditosic=data.frame(x=c(pred$x,pred$x[length(pred$x):1]),
                          y=c(pred$y.Lower,pred$y.Upper[length(pred$x):1]))
    graph=graph+geom_polygon(data=preditosic,aes(x=x,y),fill=fill.ic,alpha=alpha.ic)}

  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),linewidth=linesize,lty=linetype)+
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
                "plot"=graph,
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
