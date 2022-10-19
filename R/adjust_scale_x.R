#' Utils: Adjust x scale
#'
#' Adjust x scale for chart or charts
#'
#' @param plots Object of analysis or plot_arrange
#' @param scale x-axis scale (use vector)
#' @param limits limits in x-axis (use vector)
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
#' adjust_scale_x(a,scale = seq(10,40,5),limits = c(10,40))
adjust_scale_x=function(plots,
                        scale="default",
                        limits="default"){
  if(length(plots)==9){plots=plots}else{plots=plots[[3]]}
  requireNamespace("ggplot2")
  if(limits[1]=="default"){limits=c(min(plots$plot$data$x),
                                    max(plots$plot$data$x))}
  if(scale[1]!="default"){plots=plots+
    scale_x_continuous(breaks=scale,limits = limits)}
  plots
}
