#' Analysis: loess regression (degree 0, 1 or 2)
#'
#' Fit a polynomial surface determined by one or more numerical predictors, using local fitting.
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param degree Degree polynomial (0,1 or 2)
#' @param sample.curve Provide the number of observations to simulate curvature (default is 1000)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_bw())
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
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
#' @param fontfamily Font family
#' @return The function returns a list containing the loess regression and graph using ggplot2.
#' @seealso \link{loess}
#' @export
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' loessreg(trat,resp)

loessreg=function(trat,
                  resp,
                  degree=2,
                  sample.curve=1000,
                  ylab="Dependent",
                  xlab="Independent",
                  theme=theme_classic(),
                  legend.position="top",
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
                  fontfamily="sans"){
  requireNamespace("ggplot2")
  ymean=tapply(resp,trat,mean)
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  if(error=="SE"){ysd=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){ysd=tapply(resp,trat,sd)}
  if(error=="FALSE"){ysd=0}
  desvio=ysd
  xmean=tapply(trat,trat,mean)
  mod=loess(resp~trat,degree = degree)
  xp=seq(min(trat),max(trat),length.out = sample.curve)
  preditos=data.frame(x=xp,
                      y=predict(mod,newdata = data.frame(trat=xp)))
  x=preditos$x
  y=preditos$y
  data=data.frame(xmean,ymean)
  data1=data.frame(trat=xmean,resp=ymean)
  s="~~~ Loess~regression"
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
  graph=graph+theme+
    geom_line(data=preditos,aes(x=preditos$x,
                                y=preditos$y,color="black"),size=linesize,lty=linetype)+
    scale_color_manual(name="",values=colorline,label="Loess regression")+
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
  result=predict(mod,newdata = data.frame(trat=temp1),type="response")
  maximo=temp1[which.max(result)]
  respmax=result[which.max(result)]
  minimo=temp1[which.min(result)]
  respmin=result[which.min(result)]
  graphs=data.frame("Parameter"=c("X Maximum",
                                  "Y Maximum",
                                  "X Minimum",
                                  "Y Minimum"),
                    "values"=c(maximo,
                               respmax,
                               minimo,
                               respmin))
  graficos=list("test"=graphs,graph)
  print(graficos)
}
