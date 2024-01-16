#' Utils: Adjust y and x scale
#'
#' Adjust y and x scale for chart or charts
#'
#' @param plots Object of analysis or plot_arrange
#' @param scale.x x-axis scale (use vector)
#' @param limits.x limits in x-axis (use vector)
#' @param scale.y y-axis scale (use vector)
#' @param limits.y limits in y-axis (use vector)
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
#' adjust_scale(a,scale.y = seq(0,100,10),limits.y = c(0,100))

adjust_scale=function(plots,
                      scale.x="default",
                      limits.x="default",
                      scale.y="default",
                      limits.y="default"){
  if(length(plots)==3 | length(plots)==3){plots=plots[[3]]}else{plots=plots}
  requireNamespace("ggplot2")

  if(limits.y[1]=="default"){limits.y=c(min(plots$plot$data$y),
                                    max(plots$plot$data$y))}
  if(scale.y[1]!="default"){plots=plots+
    scale_y_continuous(breaks=scale.y,limits = limits.y)}

  if(limits.x[1]=="default"){limits.x=c(min(plots$plot$data$x),
                                      max(plots$plot$data$x))}
  if(scale.x[1]!="default"){plots=plots+
    scale_x_continuous(breaks=scale.x,limits = limits.x)}
  plots
}
