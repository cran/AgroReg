#' Analysis: Quadratic-plateau
#'
#' This function performs the quadratic-plateau regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
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
#' @param colorline Color lines
#' @param fillshape Fill shape
#' @param round round equation
#' @param xname.formula Name of x in the equation
#' @param yname.formula Name of y in the equation
#' @param comment Add text after equation
#' @param fontfamily Font family
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @details
#' The quadratic-plateau model is defined by:
#'
#' First curve:
#' \deqn{y = \beta_0 + \beta_1 \cdot x + \beta_2 \cdot x^2 (x < breakpoint)}
#'
#' Second curve:
#' \deqn{y = \beta_0 + \beta_1 \cdot breakpoint + \beta_2 \cdot breakpoint^2 (x > breakpoint)}
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
                           sample.curve=1000,
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
                      fillshape = "gray",
                      colorline = "black",
                      round=NA,
                      yname.formula="y",
                      xname.formula="x",
                      comment=NA,
                      fontfamily="sans"){
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
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  r2=floor(r2*100)/100

  if(is.na(round)==TRUE){
  coef1=coef(model1)[1]
  coef2=coef(model1)[2]
  coef3=coef(model1)[3]
  coef4=breakpoint}

  if(is.na(round)==FALSE){
    coef1=round(coef(model1)[1],round)
    coef2=round(coef(model1)[2],round)
    coef3=round(coef(model1)[3],round)
    coef4=round(breakpoint,round)}

  s <- sprintf("~~~%s == %e %s %e * %s %s %e * %s^2~(%s^'<'^%e) ~~~~~ italic(R^2) ==  %0.2f",
               yname.formula,
               coef1,
               ifelse(coef2 >= 0, "+", "-"),
               abs(coef2),
               xname.formula,
               ifelse(coef3 >= 0, "+", "-"),
               abs(coef3),
               xname.formula,
               xname.formula,
               coef4,
               r2)
  equation=s
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  predesp=predict(corr_model)
  predobs=resp
  model=corr_model
  rmse=sqrt(mean((predesp-predobs)^2))

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
  xp=seq(min(trat),max(trat),length=sample.curve)
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
  minimo=NA
  respmin=NA
  aic=AIC(corr_model_jp)
  bic=BIC(corr_model_jp)
  graphs=data.frame("Parameter"=c("Breakpoint",
                                  "Breakpoint Response",
                                  "X Minimum",
                                  "Y Minimum",
                                  "AIC",
                                  "BIC",
                                  "r-squared",
                                  "RMSE"),
                    "values"=c(breakpoint,
                               ybreakpoint,
                               minimo,
                               respmin,
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
