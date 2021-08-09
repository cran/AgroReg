#' Analysis: Quadratic-plateau Regression
#'
#' This function performs the quadratic-plateau regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position legend position (\emph{default} is "top")
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param width.bar	Bar width
#' @param r2 coefficient of determination of the mean or all values (\emph{default} is all)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param scale Sets x scale (\emph{default} is none, can be "log")
#' @param textsize Font size
#' @param pointsize	shape size
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param comment Add text after equation
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details
#' The linear-linear model is defined by:
#'
#' First curve:
#' \deqn{f(x) = \beta_0 + \beta_1 \cdot x + \beta_2 \cdot x^2 (x < breakpoint)}
#'
#' Second curve:
#' \deqn{f(x) = \beta_0 + \beta_1 \cdot breakpoint + \beta_2 \cdot breakpoint^2 (x > breakpoint)}
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Chiu, G. S., R. Lockhart, and R. Routledge. 2006. Bent-cable regression theory and applications. Journal of the American Statistical Association 101:542-553.
#' @references Toms, J. D., and M. L. Lesperance. 2003. Piecewise regression: a tool for identifying ecological thresholds. Ecology 84:2034-2041.
#' @import minpack.lm
#' @import rcompanion
#' @importFrom broom tidy
#' @importFrom stats nls.control
#' @seealso \link{linear.linear}, \link{linear.plateau}
#' @examples
#' library(AgroReg)
#' data("granada")
#' attach(granada)
#' quadratic.plateau(time,WL)

quadratic.plateau=function(trat,resp,
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
  requireNamespace("minpack.lm")
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("dplyr")
  requireNamespace("rcompanion")
  qp <- function(x, b0, b1, b2) {
    jp <- -0.5 * b1 / b2
    if_else(
      condition = x < jp,
      true  = b0 + (b1 * x) + (b2 * x * x),
      false = b0 + (b1 * jp) + (b2 * jp * jp)
    )
  }
  qp_jp <- function(x, b0, b1, jp) {
    b2 <- -0.5 * b1 / jp
    if_else(
      condition = x < jp,
      true  = b0 + (b1 * x) + (b2 * x * x),
      false = b0 + (b1 * jp) + (b2 * jp * jp)
    )
  }
  data=data.frame(trat,resp)
  start <- lm(resp ~ poly(trat, 2, raw = TRUE),
              data = data)
  start_b0 <- start$coef[[1]]
  start_b1 <- start$coef[[2]]
  start_b2 <- start$coef[[3]]
  start_jp <- mean(data$resp)
  try(corr_model <- minpack.lm::nlsLM(
    formula = resp ~ qp(trat, b0, b1, b2),
    data = data,
    start = list(b0 = start_b0,
                 b1 = start_b1,
                 b2 = start_b2)
  ))
  model1=summary(corr_model)
  try(corr_model_jp <- minpack.lm::nlsLM(
    formula = resp ~ qp_jp(trat, b0, b1, jp),
    data = data,
    start = list(b0 = start_b0,
                 b1 = start_b1,
                 jp = start_jp)))
  model2=summary(corr_model_jp)
  breakpoint=model2$coefficients[3,1]
  ybreakpoint=predict(corr_model,newdata = data.frame(trat=breakpoint))
  nullfunct <- function(x, m) {
    m
  }

  m_ini <- mean(data$resp)

  null <- nls(resp ~ nullfunct(trat, m),
              data = data,
              start = list(m = m_ini),
              trace = FALSE)
  model_error <- round(summary(corr_model)$sigma, 2)
  r2 <- rcompanion::nagelkerke(corr_model, null)$Pseudo.R.squared.for.model.vs.null[2]
  b0 <- tidy(corr_model)$estimate[1]
  b1 <- tidy(corr_model)$estimate[2]
  b2 <- tidy(corr_model)$estimate[3]
  cc <- tidy(corr_model_jp)$estimate[3]
  plateau <- round(b0 + (b1 * cc) + (b2 * cc * cc), 1)

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
  s <- sprintf("~~~y == %e %s %e * x %s %e * x^2~(x<%e) ~~~~~ italic(R^2) ==  %0.2f",
                                coef(model1)[1],
                                ifelse(coef(model1)[2] >= 0, "+", "-"),
                                abs(coef(model1)[2]),
                                ifelse(coef(model1)[3] >= 0, "+", "-"),
                                abs(coef(model1)[3]),
                                breakpoint,
                                r2)
  equation=s
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  predesp=predict(corr_model)
  predobs=resp
  rmse=sqrt(mean((predesp-predobs)^2))

  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  if(point=="mean"){
    graph=ggplot(data,aes(x=xmean,y=ymean))
    if(error!="FALSE"){graph=graph+geom_errorbar(aes(ymin=ymean-ysd,ymax=ymean+ysd),
                                                 width=width.bar)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  xp=seq(min(trat),max(trat),length=1000)
  yp=predict(corr_model,newdata=data.frame(trat=xp))
  preditos=data.frame(x=xp,y=yp)
  temp1=xp
  result=yp
  x=xp
  y=yp
  graph=graph+theme+
    geom_line(data=preditos,aes(x=x,
                                y=y,
                                color="black"),size=linesize)+
    scale_color_manual(name="",values="black",label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text = element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  aic=AIC(corr_model_jp)
  bic=BIC(corr_model_jp)
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
  graficos=list("Coefficients quadratic model"=model1,
                "Coefficients segmented"=model2,
                "values"=graphs,
                graph)
  print(graficos)

}
