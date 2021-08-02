#' Analysis: Michaelis-Menten Model
#'
#' This function performs regression analysis using the Michaelis-Menten model.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @details
#' The Michaelis-Menten model is defined by:
#' \deqn{f(x, (VM,k)) = c + \frac{Vm*x}{k+x}}
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Gabriel Danilo Shimizu
#' @references Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
#' @examples
#' data("granada")
#' attach(granada)
#' MM(time,WL)
#' @md

MM=function(trat,
            resp,
            error="SE",
            ylab="ylab",
            xlab="xlab",
            theme=theme_classic(),
            legend.position="top",
            point="all",
            r2="all"){
  requireNamespace("ggplot2")
  requireNamespace("drc")
  requireNamespace("crayon")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  xmean=tapply(trat,trat,mean)
  desvio=ysd
  mod=nls(resp~SSmicmen(trat, Vm, K))
  coef=summary(mod)
  d=coef$coefficients[,1][1]
  e=coef$coefficients[,1][2]
  if(r2=="all"){r2=cor(resp, fitted(mod))^2}
  if(r2=="mean"){r2=cor(ymean, predict(mod,newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100
  equation=sprintf("~~~y==frac(%0.3e*x, %0.3e+x)~~~~~ italic(R^2) == %0.2f",
                   d,
                   e,
                   r2)
  xp=seq(min(trat),max(trat),length.out = 1000)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
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
    geom_line(data=preditos,aes(x=x,y=y,color="black"))+
    scale_color_manual(name="",values=1,label=parse(text = equation))+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  temp1=seq(min(trat),max(trat),length.out=10000)
  result=predict(mod,newdata = data.frame(trat=temp1),type="response")
  aic=AIC(mod)
  bic=BIC(mod)
  print(coef)
  predesp=predict(mod)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))
  cat("\n===================================\n")
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
