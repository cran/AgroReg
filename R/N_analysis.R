#' Analysis: Graph for not significant trend
#' @author Gabriel Danilo Shimizu
#' @author Leandro Simoes Azeredo Goncalves
#' @description Graph for non-significant trend. Can be used within the multicurve command
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param ylab Dependent variable name (Accepts the \emph{expression}() function)
#' @param xlab Independent variable name (Accepts the \emph{expression}() function)
#' @param error Error bar (It can be SE - \emph{default}, SD or FALSE)
#' @param width.bar	Bar width
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.position legend position (\emph{default} is "top")
#' @param legend.text legend text
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param textsize Font size
#' @param pointsize	shape size
#' @param pointshape format point (default is 21)
#' @return The function returns an exploratory graph of segments
#' @keywords non-significant
#' @export
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' Nreg(trat,resp)

Nreg=function(trat,
              resp,
              ylab="Dependent",
              xlab="Independent",
              error="SE",
              theme=theme_classic(),
              legend.position="top",
              legend.text="not~significant",
              width.bar=NA,
              point="all",
              textsize = 12,
              pointsize = 4.5,
              pointshape = 21){
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("ggplot2")
  dados=data.frame(trat,resp)
  medias=c()
  dose=tapply(trat, trat, mean)
  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=as.numeric(as.character(names(media))),
                   resp=media,
                   desvio)
  temp1=dose
  result=media
  s=legend.text
  if(point=="mean"){
    grafico=ggplot(data1,aes(x=media,y=desvio))
    if(error!="FALSE"){grafico=grafico+
      geom_errorbar(aes(ymin=media-desvio,ymax=media+desvio),
                                                 width=width.bar)}
    grafico=grafico+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(point=="all"){
    grafico=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    grafico=grafico+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}

  grafico=grafico+theme+ylab(ylab)+xlab(xlab)+
    scale_color_manual(values="black",label=c(parse(text=s)),name="")+
    theme(text = element_text(size=textsize,color="black"),
          axis.text = element_text(size=textsize,color="black"),
          axis.title = element_text(size=textsize,color="black"),
          legend.position = legend.position,
          legend.text=element_text(size=textsize),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)
  graficos=list(teste="not significant",grafico)
  print(graficos)
}
