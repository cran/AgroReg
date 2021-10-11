#' Merge multiple curves into a single graph
#' @param plots list with objects of type analysis.
#' @param theme ggplot2 theme (\emph{default} is theme_classi())
#' @param legend.title caption title
#' @param trat name of the curves
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
#' @param gray gray scale (\emph{default} is FALSE)
#' @param widthbar bar width (\emph{default} is 0.3)
#' @param pointsize shape size
#' @param linesize line size
#' @param textsize Font size
#' @return The function returns a graph joining the outputs of the functions LM_model, LL_model, BC_model, CD_model, loess_model, normal_model, piecewise_model and N_model
#' @author Gabriel Danilo Shimizu
#' @export
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' a=LM(trat,resp)
#' b=LL(trat,resp,npar = "LL.3")
#' c=BC(trat,resp, npar = "BC.4")
#' d=CD(trat,resp, npar = "CRS.4")
#' plot_arrange(list(a,b,c,d))

plot_arrange=function(plots,
                      point="mean",
                      theme = theme_classic(),
                      legend.title = NULL,
                      legend.position = "top",
                      trat = NA,
                      gray = FALSE,
                      ylab = "Dependent",
                      xlab = "Independent",
                      widthbar = 0,
                      pointsize = 4.5,
                      linesize = 0.8,
                      textsize = 12) {
  requireNamespace("ggplot2")
  equation=1:length(plots)
  grafico=ggplot()
  if(gray==FALSE & point=="mean"){
  for(i in 1:length(plots)){
    equation[[i]]=plots[[i]][[]]$plot$s
    x=plots[[i]][[]]$plot$temp1
    y=plots[[i]][[]]$plot$result
    data=data.frame(x,y,color=as.factor(i))
    pontosx=plots[[i]][[]]$plot$data1$trat
    pontosy=plots[[i]][[]]$plot$data1$resp
    desvio=plots[[i]][[]]$plot$desvio
    pontos=data.frame(x=pontosx,
                      y=pontosy,
                      desvio=desvio,
                      color=as.factor(i))
    color=pontos$color
    grafico=grafico+
      geom_errorbar(data=pontos,
                    aes(x=x,
                        y=y,
                        ymin=y-desvio,
                        ymax=y+desvio,
                        color=color,
                        group=color),width=widthbar, size=linesize)+
      geom_point(data=pontos,aes(x=x,y=y,
                                 color=color,
                                 group=color),size=pointsize)+
      geom_line(data=data,aes(x=x,
                              y=y,
                              color=color,
                              group=color),size=linesize)
  }
  texto=parse(text=paste(trat,"~",unlist(equation)))
  grafico=grafico+
    scale_color_discrete(label=texto)+
    theme+labs(color=legend.title)+
    theme(axis.text = element_text(size=12,color="black"),
          legend.position = legend.position,
          legend.justification='left',
          legend.direction = "vertical",
          legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  if(gray==TRUE & point=="mean"){
    for(i in 1:length(plots)){
      equation[[i]]=plots[[i]][[]]$plot$s
      x=plots[[i]][[]]$plot$temp1
      y=plots[[i]][[]]$plot$result
      data=data.frame(x,y,color=as.factor(i))
      pontosx=plots[[i]][[]]$plot$data1$trat
      pontosy=plots[[i]][[]]$plot$data1$resp
      desvio=plots[[i]][[]]$plot$desvio
      pontos=data.frame(x=pontosx,y=pontosy,desvio=desvio,color=as.factor(i))
      grafico=grafico+
        geom_errorbar(data=pontos,
                      aes(x=x,
                          y=y,
                          ymin=y-desvio,
                          ymax=y+desvio),width=widthbar, size=linesize)+
        geom_point(data=pontos,aes(x=x,
                                   y=y,
                                   pch=color,
                                   group=color),
                   size=pointsize,fill="gray")+
        geom_line(data=data,aes(x=x,
                                y=y,
                                lty=color,
                                group=color),size=linesize)
    }
    texto=parse(text=paste(trat,"~",unlist(equation)))
    grafico=grafico+
      scale_linetype_discrete(label=texto)+
      scale_shape_discrete(label=texto)+
      theme+labs(lty=legend.title,shape=legend.title)+
      theme(axis.text = element_text(size=textsize,color="black"),
            axis.title = element_text(size=textsize,color="black"),
            legend.position = legend.position,
            legend.text=element_text(size=textsize),
            legend.justification='left',
            legend.direction = "vertical",
            legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  if(gray==FALSE & point=="all"){
    for(i in 1:length(plots)){
      equation[[i]]=plots[[i]][[]]$plot$s
      x=plots[[i]][[]]$plot$temp1
      y=plots[[i]][[]]$plot$result
      data=data.frame(x,y,color=as.factor(i))
      pontosx=plots[[i]][[]]$plot$trat
      pontosy=plots[[i]][[]]$plot$resp
      pontos=data.frame(x=pontosx,
                        y=pontosy,
                        color=as.factor(i))
      color=pontos$color
      grafico=grafico+
        geom_point(data=pontos,aes(x=x,y=y,
                                   color=color,
                                   group=color),size=pointsize)+
        geom_line(data=data,aes(x=x,
                                y=y,
                                color=color,
                                group=color),size=linesize)}
    texto=parse(text=paste(trat,"~",unlist(equation)))
    grafico=grafico+
      scale_color_discrete(label=texto)+
      theme+labs(color=legend.title)+
      theme(axis.text = element_text(size=12,color="black"),
            legend.position = legend.position,
            legend.justification='left',
            legend.direction = "vertical",
            legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  if(gray==TRUE){
    for(i in 1:length(plots)){
      equation[[i]]=plots[[i]][[]]$plot$s
      x=plots[[i]][[]]$plot$temp1
      y=plots[[i]][[]]$plot$result
      data=data.frame(x,y,color=as.factor(i))
      pontosx=plots[[i]][[]]$plot$trat
      pontosy=plots[[i]][[]]$plot$resp
      pontos=data.frame(x=pontosx,
                        y=pontosy,
                        color=as.factor(i))
      grafico=grafico+
        geom_point(data=pontos,aes(x=x,
                                   y=y,
                                   pch=color,
                                   group=color),
                   size=pointsize,fill="gray")+
        geom_line(data=data,aes(x=x,
                                y=y,
                                lty=color,
                                group=color),size=linesize)
    }
    texto=parse(text=paste(trat,"~",unlist(equation)))
    grafico=grafico+
      scale_linetype_discrete(label=texto)+
      scale_shape_discrete(label=texto)+
      theme+labs(lty=legend.title,shape=legend.title)+
      theme(axis.text = element_text(size=textsize,color="black"),
            axis.title = element_text(size=textsize,color="black"),
            legend.position = legend.position,
            legend.text=element_text(size=textsize),
            legend.justification='left',
            legend.direction = "vertical",
            legend.text.align = 0)+ylab(ylab)+xlab(xlab)}
  print(grafico)
}
