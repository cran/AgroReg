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
#' @param legend.add.mean Add average in legend
#' @param legend.add.mean.name Add media name
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param textsize Font size
#' @param pointsize	shape size
#' @param add.line Add line
#' @param add.line.mean Add line mean
#' @param linesize	line size
#' @param pointshape format point (default is 21)
#' @param fontfamily Font family
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
              legend.add.mean=TRUE,
              legend.add.mean.name="hat(y)",
              width.bar=NA,
              point="all",
              textsize = 12,
              add.line=FALSE,
              add.line.mean=FALSE,
              linesize=0.8,
              pointsize = 4.5,
              pointshape = 21,
              fontfamily="sans"){
  if(is.na(width.bar)==TRUE){width.bar=0.01*mean(trat)}
  requireNamespace("ggplot2")
  dados=data.frame(trat,resp)
  dose=tapply(trat, trat, mean)
  media=tapply(resp, trat, mean)
  if(error=="SE"){desvio=tapply(resp,trat,sd)/sqrt(tapply(resp,trat,length))}
  if(error=="SD"){desvio=tapply(resp,trat,sd)}
  if(error=="FALSE"){desvio=0}
  data1=data.frame(trat,resp)
  data1=data.frame(trat=dose,
                   resp=media,
                   desvio)
  temp1=dose
  result=media
  s=legend.text
  if(legend.add.mean==TRUE){s=paste(legend.text,"~","(",legend.add.mean.name,"==",mean(media),")")}
  if(point=="mean"){
    grafico=ggplot(data1,aes(x=trat,y=resp))
    if(add.line==TRUE){grafico=grafico+geom_line(size=linesize)}
    if(add.line.mean==TRUE){grafico=grafico+geom_hline(yintercept = mean(resp,na.rm=TRUE),lty=2)}
    if(error!="FALSE"){grafico=grafico+
      geom_errorbar(aes(ymin=resp-desvio,ymax=resp+desvio),
                                                 width=width.bar,
                    size=linesize)}
    grafico=grafico+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")
    }
  if(point=="all"){
    grafico=ggplot(data.frame(trat,resp),aes(x=trat,y=resp))
    if(add.line==TRUE){grafico=grafico+stat_summary(geom="line",fun = "mean",
                                                    size=linesize)}
    if(add.line.mean==TRUE){grafico=grafico+geom_hline(yintercept = mean(resp,na.rm=TRUE),lty=2)}
    grafico=grafico+
      geom_point(aes(color="black"),size=pointsize,shape=pointshape,fill="gray")}
  if(add.line.mean==TRUE){result=mean(media)}
  grafico=grafico+theme+ylab(ylab)+xlab(xlab)+
    scale_color_manual(values="black",label=c(parse(text=s)),name="")+
    theme(text = element_text(size=textsize,color="black"),
          axis.text = element_text(size=textsize,color="black",family = fontfamily),
          axis.title = element_text(size=textsize,color="black",family = fontfamily),
          legend.position = legend.position,
          legend.text=element_text(size=textsize,family = fontfamily),
          legend.direction = "vertical",
          legend.text.align = 0,
          legend.justification = 0)
  graficos=list(teste="not significant",grafico)
  print(graficos)
}
