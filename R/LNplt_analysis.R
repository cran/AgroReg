#' Analysis: Linear-Plateau Regression
#'
#' This function performs the quadratic-plateau regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); breakpoint and the graph using ggplot2 with the equation automatically.
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression theory and applications. Journal of the American Statistical Association 101:542-553.
#' @references Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for identifying ecological thresholds. Ecology 84:2034-2041.
#' @seealso \link{quadratic.plateau}, \link{linear.linear}
#' @importFrom dplyr if_else
#' @details
#' The linear-plateau model is defined by:
#' First curve:
#' \deqn{f(x) = b0+b1*x (x<breakpoint)}
#'
#' Second curve:
#' \deqn{f(x) = b0+b1*breakpoint (x>breakpoint)}
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' linear.plateau(time,WL)

linear.plateau=function(trat,resp,
                      ylab="Dependent",
                      xlab="Independent",
                      theme=theme_classic(),
                      legend.position="top",
                      error="SE",
                      r2="all",
                      point="all",
                      scale="none"){
  lp <- function(x, a, b, c) {
    if_else(condition = x < c,
            true = a + b * x,
            false = a + b * c)
  }
  data=data.frame(trat,resp)
  ini_fit <- lm(data = data, formula = resp ~ trat)
  ini_a <- ini_fit$coef[[1]]
  ini_b <- ini_fit$coef[[2]]
  ini_c <- mean(data$trat)

    lp_model <- nlsLM(
      formula = resp ~ lp(trat, a, b, c),
      data = data,
      start = list(a = ini_a, b = ini_c, c = ini_c))
  model2=summary(lp_model)
  breakpoint=model2$coefficients[3,1]
  ybreakpoint=predict(lp_model,newdata = data.frame(trat=breakpoint))
  m.ini <- mean(data$resp)

  nullfunct <- function(x, m) {
    m
  }
  null <- nls(resp~ nullfunct(trat, m),
              data = data,
              start = list(m = m.ini),
              trace = FALSE,
              nls.control(maxiter = 1000))
  r2 <- nagelkerke(lp_model, null)$Pseudo.R.squared.for.model.vs.null[2]
  requireNamespace("drc")
  requireNamespace("crayon")
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  r2=floor(r2*100)/100
  s <- sprintf("~~~y == %e %s %e * x ~(x<%e) ~~~~~ italic(R^2) ==  %0.2f",
               coef(model2)[1],
               ifelse(coef(model2)[2] >= 0, "+", "-"),
               abs(coef(model2)[2]),
               breakpoint,
               r2)

  equation=s
  predesp=predict(lp_model)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))

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
  xp=seq(min(trat),max(trat),length=1000)
  yp=predict(lp_model,newdata=data.frame(trat=xp))
  preditos=data.frame(x=xp,y=yp)
  temp1=xp
  result=yp
  x=xp
  y=yp
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=0.8)+
    scale_color_manual(name="",values="black",label=parse(text = equation))+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=12),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  aic=AIC(lp_model)
  bic=BIC(lp_model)
  graphs=data.frame("Parameter"=c("Breakpoint",
                                  "Response",
                                  "AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(breakpoint,
                               ybreakpoint,
                               aic,
                               bic,
                               r2,
                               rmse))
  graficos=list("Coefficients segmented"=model2,
                "values"=graphs,
                graph)
  print(graficos)

}
