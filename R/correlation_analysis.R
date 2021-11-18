#' Graph: Plot correlation
#'
#' @description Correlation analysis function (Pearson or Spearman)
#' @param x Numeric vector with independent variable
#' @param y Numeric vector with dependent variable
#' @param method Method correlation (\emph{default} is Pearson)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param textsize Axis text size
#' @param pointsize Point size
#' @param pointshape shape format
#' @param linesize	line size
#' @param ic Add interval of confidence
#' @param fill.ic Color interval of confidence
#' @param alpha.ic confidence interval transparency level
#' @param title title
#' @param fontfamily Font family
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @return The function returns a graph for correlation
#' @importFrom stats cor
#' @importFrom stats cor.test
#' @export
#' @examples
#' data("aristolochia")
#' with(aristolochia, correlation(trat,resp))

correlation = function(x,
                       y,
                       method = "pearson",
                       ylab = "Dependent",
                       xlab = "Independent",
                       theme = theme_classic(),
                       textsize = 12,
                       pointsize = 5,
                       pointshape = 21,
                       linesize = 0.8,
                       fill.ic = "gray70",
                       alpha.ic = 0.5,
                       ic = TRUE,
                       title = NA,
                       fontfamily = "sans") {
  if(is.na(title)==TRUE){
    if(method=="pearson"){title="Pearson correlation"}
    if(method=="spearman"){title="Spearman correlation"}
  }
  if(method=="pearson"){corre=cor(y,x)
  pvalor=cor.test(y,x,method="pearson",exact=FALSE)$p.value}
  if(method=="spearman"){corre=cor(y,x,method = "spearman")
  pvalor=cor.test(y,x,method="spearman",exact=FALSE)$p.value}
  requireNamespace("ggplot2")
  data=data.frame(y,x)
  ggplot(data,aes(y=y,x=x))+
    geom_point(shape=pointshape,fill="gray",color="black",size=pointsize)+
    geom_smooth(
      color = "black",
      se = ic,
      size = linesize,
      fill = fill.ic,
      alpha = alpha.ic,
      na.rm = TRUE,
      method = "lm")+
    theme+labs(title=title,x=xlab,y=ylab)+
    theme(axis.text = element_text(size=textsize,
                                   color="black",family = fontfamily),
          axis.title = element_text(family = fontfamily))+
    annotate(geom = "text",x = -Inf,y=Inf,
             label=paste("R = ", round(corre,2), ", p-value =",
                         format(round(pvalor,3),scientific = TRUE),sep = ""),
             hjust = -0.1, vjust = 1.1)
}
