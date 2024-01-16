#' Utils: Adjust y scale
#'
#' Adjust y scale for chart or charts
#'
#' @param plots Object of analysis or plot_arrange
#' @param scale y-axis scale (use vector)
#' @param limits limits in y-axis (use vector)
#'
#' @return Returns the scaled graph
#' @export
#'
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' a=LM(trat,resp)
#' b=LL(trat,resp,npar = "LL.3")
#' a=plot_arrange(list(a,b),gray = TRUE)
#' adjust_scale_y(a,scale = seq(0,100,10),limits = c(0,100))

adjust_scale_y=function(plots,
                        scale="default",
                        limits="default"){
  if(length(plots)==3 | length(plots)==3){plots=plots[[3]]}else{plots=plots}
  requireNamespace("ggplot2")
  if(limits[1]=="default"){limits=c(min(plots$plot$data$y),
                                    max(plots$plot$data$y))}
  if(scale[1]!="default"){plots=plots+
    scale_y_continuous(breaks=scale,limits = limits)}
  plots
}
