#' Change the colors of a graph from the plot_arrange function
#' @param graphs object from a plot_arrange function
#' @param color color curve and point
#' @return The function changes the colors of a graph coming from the plot_arrange function
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
#' graph=plot_arrange(list(a,b,c,d))
#' coloredit_arrange(graph,color=c("black","red","blue","green"))

coloredit_arrange = function(graphs, color = NA) {
  if (is.na(color[1]) == TRUE) {
    graphs = graphs + scale_color_discrete(labels = parse(text = graphs$plot$equation))
  }
  if (is.na(color[1]) == FALSE) {
    graphs = graphs + scale_color_manual(values = color,
                                         labels = parse(text = graphs$plot$equation))
  }
  graphs
}
