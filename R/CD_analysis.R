#' Analysis: Logistic regression Cedergreen-Ritz-Streibig model
#'
#' The 'CRS.4' and 'CRS.5' logistical models provide Brain-Cousens modified logistical models to describe u-shaped hormesis. This model was extracted from the 'drc' package.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param npar Number of model parameters
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
#' @param comment Add text after equation
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The four-parameter model is given by the expression:
#'
#' \deqn{f(x) = 0 + \frac{d-0+f \exp(-1/x)}{1+\exp(b(\log(x)-\log(e)))}}
#'
#' while the five-parameter is:
#'
#' \deqn{f(x) = c + \frac{d-c+f \exp(-1/x)}{1+\exp(b(\log(x)-\log(e)))}}
#' @seealso \link{LL}, \link{BC}, \link{GP}
#' @export
#' @author Model imported from the drc package (Ritz et al., 2016)
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley \& Sons (p. 330).
#' @references Ritz, C.; Strebig, J.C.; Ritz, M.C. Package 'drc'. Creative Commons: Mountain View, CA, USA, 2016.
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' CD(trat,resp)

CD=function(trat,
            resp,
            npar="CRS.4",
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
  requireNamespace("ggplot2")
  requireNamespace("drc")
  requireNamespace("crayon")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd

  if(npar=="CRS.4"){mod=drm(resp~trat,fct=CRS.4a())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  d=coef$coefficients[,1][2]
  e=coef$coefficients[,1][3]
  f=coef$coefficients[,1][4]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y==frac(%0.3e %s %0.3e*exp(-1/x), 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
                   d,
                   ifelse(f >= 0, "+", "-"),
                   abs(f),
                   b,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  }
  if(npar=="CRS.5"){mod=drm(resp~trat,fct=CRS.5a())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  c=coef$coefficients[,1][2]
  d=coef$coefficients[,1][3]
  e=coef$coefficients[,1][4]
  f=coef$coefficients[,1][5]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y == %0.3e + frac(%0.3e %s %0.3e %s %0.3e*exp(-1/x), 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
                   c,
                   d,
                   ifelse(c >= 0, "+", "-"),
                   abs(c),
                   ifelse(f >= 0, "+", "-"),
                   abs(f),
                   b,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
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
  graph=ggplot(data,aes(x=xmean,y=ymean))
  if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                               width=width.bar)}
  graph=graph+
    geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")+
    theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=linesize)+
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
