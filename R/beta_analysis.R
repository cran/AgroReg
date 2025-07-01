#' Analysis: Beta
#'
#' This function performs beta regression analysis.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param legend.position Legend position (\emph{default} is "top")
#' @param r2 Coefficient of determination of the mean or all values (\emph{default} is all)
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
#' @details
#' The beta model is defined by:
#' \deqn{Y = d \times \{(\frac{X-X_b}{X_o-X_b})(\frac{X_c-X}{X_c-X_o})^{\frac{X_c-X_o}{X_o-X_b}}\}^b}
#' @author Model imported from the aomisc package (Andrea Onofri)
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @references Onofri, A., 2020. The broken bridge between biologists and statisticians: a blog and R package. Statforbiology. http://www.statforbiology.com/tags/aomisc/
#' @export
#' @examples
#' library(AgroReg)
#' X <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' Y <- c(0, 0, 0, 7.7, 12.3, 19.7, 22.4, 20.3, 6.6, 0, 0)
#' beta_reg(X,Y)


beta_reg=function(trat,
                  resp,
                  sample.curve=1000,
                  ylab = "Dependent",
                  xlab = "Independent",
                  theme = theme_classic(),
                  legend.position = "top",
                  error = "SE",
                  r2 = "all",
                  point = "all",
                  width.bar = NA,
                  scale = "none",
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
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  beta.fun <- function(X, b, d, Xb, Xo, Xc){
    .expr1 <-  (X - Xb)/(Xo - Xb)
    .expr2 <- (Xc - X)/(Xc - Xo)
    .expr3 <- (Xc - Xo)/(Xo - Xb)
    ifelse(X > Xb & X < Xc, d * (.expr1*.expr2^.expr3)^b, 0)}
  DRC.beta <- function(){
    fct <- function(x, parm) {
      beta.fun(x, parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
    }
    ssfct <- function(data){
      x <- data[, 1]
      y <- data[, 2]

      d <- max(y)
      Xo <- x[which.max(y)]
      firstidx <- min( which(y !=0) )
      Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
      secidx <- max( which(y !=0) )
      Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
      c(1, d, Xb, Xo, Xc)
    }
    names <- c("b", "d", "Xb", "Xo", "Xc")
    text <- "Beta function"
    returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
    class(returnList) <- "drcMean"
    invisible(returnList)}
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  model <- drm(resp~trat,fct=DRC.beta())
  coef=summary(model)

  if(is.na(round)==TRUE){
  b=coef$coefficients[,1][1]
  d=coef$coefficients[,1][2]
  Xb=coef$coefficients[,1][3]
  Xo=coef$coefficients[,1][4]
  Xc=coef$coefficients[,1][5]}

  if(is.na(round)==FALSE){
    b=round(coef$coefficients[,1][1],round)
    d=round(coef$coefficients[,1][2],round)
    Xb=round(coef$coefficients[,1][3],round)
    Xo=round(coef$coefficients[,1][4],round)
    Xc=round(coef$coefficients[,1][5],round)}

  if(r2=="all"){r2=cor(resp, fitted(model))^2}
  if(r2=="mean"){r2=cor(ymean, predict(model,
                                       newdata=data.frame(trat=unique(trat))))^2}
  r2=floor(r2*100)/100

  xoxb=Xo-Xb
  xcxo=Xc-Xo
  xcxoxoxb=(Xc-Xo)/(Xo-Xb)
  equation=sprintf("~~~%s==%0.3e*bgroup(\"{\",group(\"(\",frac(%s %s %0.3e, %0.3e),\")\")*group(\"(\",frac(%0.3e-%s, %0.3e),\")\")^%0.3e,\"}\")^%0.3e ~~~~~ italic(R^2) == %0.2f",
                   yname.formula,
                   d,
                   xname.formula,
                   ifelse(Xb >= 0, "+", "-"),
                   abs(Xb),
                   xoxb,
                   Xc,
                   xname.formula,
                   xcxo,
                   xcxoxoxb,
                   b,
                   r2)
  if(is.na(comment)==FALSE){equation=paste(equation,"~\"",comment,"\"")}
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(model,newdata = data.frame(trat=xp)))
  predesp=predict(model)
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
                                                 width=width.bar,
                                                 linewidth=linesize)}
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}
  if(point=="all"){
    graph=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    graph=graph+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill=fillshape)}

  graph=graph+theme+geom_line(data=preditos,aes(x=x,
                                                y=y,color="black"),linewidth=linesize,lty=linetype)+
    scale_color_manual(name="",values=colorline,label=parse(text = equation))+
    theme(axis.text = element_text(size=textsize,color="black", family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text = element_text(size=textsize, family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)+
    ylab(ylab)+xlab(xlab)
  if(scale=="log"){graph=graph+scale_x_log10()}
  temp1=seq(min(trat),max(trat),length.out=sample.curve)
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
