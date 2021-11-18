#' Analysis: Extract models
#'
#' @description This function allows extracting the model (type="model") or residuals (type="resids"). The model class depends on the function and can be (lm, drm or nls). This function also allows you to perform graphical analysis of residuals (type="residplot"), graphical analysis of standardized residuals (type="stdresidplot"), graph of theoretical quantiles (type="qqplot").
#' @param model Object returned from an analysis function
#' @param type output type
#' @export
#' @return Returns an object of class drm, lm or nls (type="model"), or vector of residuals (type="resids"), or graph of the residuals (type="residplot", type="stdresidplot", type=" qqplot").
#' @importFrom stats resid
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @examples
#'
#' data("aristolochia")
#' attach(aristolochia)
#' a=linear.linear(trat,resp,point = "mean")
#' extract.model(a,type = "qqplot")

extract.model=function(model,type="model"){
  if(type=="model"){results=model[[]]$plot$model}
  if(type=="resid"){results=resid(model[[]]$plot$model)}
  if(type=="residgraph"){
    resids=resid(model[[]]$plot$model)
    data=data.frame(ID=1:length(resids),resids)
    requireNamespace("ggplot2")
    results=ggplot(data,aes(x=ID,y=resids))+
      geom_hline(yintercept = 0,lty=2, size=1)+
      geom_point(shape=21,color="black",fill="lightblue",size=4.5)+
      geom_text(aes(label=ID),color="red")+
      theme_classic()+
      ylab("Residuals")+
      # scale_x_continuous(breaks=seq(1,length(resids)))+
      theme(axis.text = element_text(size=12))}
  if(type=="stdresidgraph"){
    resids=resid(model[[]]$plot$model)
    resids=resids/sd(resids)
    ID=1:length(resids)
    data=data.frame(ID=ID,resids)
    requireNamespace("ggplot2")
    results=ggplot(data,aes(x=ID,y=resids))+
      geom_hline(yintercept = c(3,0,-3),lty=2,
                 color=c("red","black","red"),
                 size=1)+
      geom_point(shape=21,color="black",fill="lightblue",size=4.5)+
      geom_text(aes(label=ID),color="red")+
      theme_classic()+
      ylab("Standart residuals")+
      #scale_x_continuous(breaks=seq(1,length(resids)))+
      theme(axis.text = element_text(size=12))}
  if(type=="qqplot"){
    yres=sort(resid(model[[]]$plot$model))
    distribution = qnorm
    probs=c(0.25,0.75)
    stopifnot(length(probs) == 2, is.function(distribution))
    y <- quantile(yres, c(0.25,0.75), names = FALSE, type = 7, na.rm = TRUE)
    x <- distribution(probs)
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    data=data.frame(ID=1:length(yres),yres)
    requireNamespace("ggplot2")
    results=ggplot(data,aes(sample=yres))+
      stat_qq(alpha=0.5,size=4.5,shape=21,fill="blue",color="black")+
      geom_abline(slope = slope,
                  intercept = int,
                  color="black",
                  size=1,
                  lty=2)+
      theme_classic()+
      ylab("Sample Quantiles")+
      xlab("Theoretical Quantiles")+
      theme(axis.text = element_text(size=12))}
  results}


