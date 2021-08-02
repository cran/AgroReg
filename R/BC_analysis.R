#' Analysis: Logistic regression Brain-Cousens hormesis models
#'
#' The 'BC.4' and 'BC.5' logistical models provide Brain-Cousens' modified logistical models to describe u-shaped hormesis. This model was extracted from the 'drc' package.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param npar Number of model parameters (\emph{default} is  BC.4)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details The model function for the Brain-Cousens model (Brain and Cousens, 1989) is
#' \deqn{ f(x, b,c,d,e,f) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
#' and it is a five-parameter model, obtained by extending the four-parameter log-logistic model (LL.4 to take into account inverse u-shaped hormesis effects.
#' Fixing the lower limit at 0 yields the four-parameter model
#' \deqn{f(x) = 0 + \frac{d-0+fx}{1+\exp(b(\log(x)-\log(e)))}}
#' used by van Ewijk and Hoekstra (1993).
#' @export
#' @import ggplot2
#' @import drc
#' @importFrom crayon green
#' @importFrom stats coef
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom stats cor
#' @importFrom stats predict
#' @importFrom stats fitted
#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom car vif
#' @importFrom stats glm
#' @importFrom stats loess
#' @importFrom stats nls
#' @author Model imported from the drc package (Ritz et al., 2016)
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley \& Sons (p. 330).
#' @references Ritz, C.; STREBIG, J.C. and RITZ, M.C. Package ‘drc’. Creative Commons: Mountain View, CA, USA, 2016.
#' @seealso \link{LL}, \link{CD},\link{GP}
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' BC(trat,resp)

BC=function(trat,
            resp,
            npar="BC.4",
            ylab="Dependent",
            xlab="Independent",
            theme=theme_classic(),
            legend.position="top",
            r2="all",
            error="SE",
            point="all",
            scale="none"){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  requireNamespace("crayon")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd
  if(npar=="BC.4"){mod=drm(resp~trat,fct=BC.4())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  d=coef$coefficients[,1][2]
  e=coef$coefficients[,1][3]
  f=coef$coefficients[,1][4]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y==frac(%0.3e %s %0.3e*x, 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
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
  if(npar=="BC.5"){mod=drm(resp~trat,fct=BC.5())
  coef=summary(mod)
  b=coef$coefficients[,1][1]
  c=coef$coefficients[,1][2]
  d=coef$coefficients[,1][3]
  e=coef$coefficients[,1][4]
  f=coef$coefficients[,1][5]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y == %0.3e + frac(%0.3e %s %0.3e %s %0.3e*x, 1+e^(%0.3e*(log(x)-log(%0.3e)))) ~~~~~ italic(R^2) == %0.2f",
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
                                                 width=0.5)}
    graph=graph+
      geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=4.5,shape=21,fill="gray")}
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=0.8)+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  temp1=seq(min(trat),max(trat),length.out=10000)
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
