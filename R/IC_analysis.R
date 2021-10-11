#' Analysis: Interval of confidence
#'
#' @description Interval of confidence in model regression
#' @param model Object analysis
#' @author Gabriel Danilo Shimizu
#' @return Return in the interval of confidence
#' @importFrom stats confint
#' @export
#' @examples
#' data("granada")
#' attach(granada)
#' a=LM(time, WL)
#' interval.confidence(a)

interval.confidence=function(model){
  modelo=model[[]]$plot$model
  confint(modelo)}
