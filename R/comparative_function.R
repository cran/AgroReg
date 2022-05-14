#' Analysis: Comparative models
#'
#' This function allows the construction of a table and/or graph with the statistical parameters to choose the model from the analysis functions.
#' @param models List with objects of type analysis
#' @param names_model Names of the models
#' @param plot Plot in the parameters
#' @param round.label Round label plot
#' @return Returns a table and/or graph with the statistical parameters for choosing the model.
#' @author Gabriel Danilo Shimizu
#' @export
#' @examples
#' library(AgroReg)
#' data(granada)
#' attach(granada)
#' a=LM(time,WL)
#' b=LL(time,WL)
#' c=BC(time,WL)
#' d=weibull(time,WL)
#' comparative_model(models=list(a,b,c,d),names_model=c("LM","LL","BC","Weibull"))
#' @importFrom egg ggarrange

comparative_model=function(models,
                           names_model=NA,
                           plot=FALSE,
                           round.label=2){
  tabela=matrix(rep(NA,length(models)*4),ncol=4)
  for(i in 1:length(models)){
    tabela[i,]=c(models[[i]]$values[5:8,2])}
  tabela=data.frame(tabela)
  if(is.na(names_model[1])==TRUE){rownames(tabela)=
    paste("Model",1:length(models))}else{
    rownames(tabela)=names_model}
  colnames(tabela)=c("AIC","BIC","R2","RMSE")
  requireNamespace("ggplot2")
  if(plot==TRUE){
    modelo=rownames(tabela)
    R2=tabela$R2
    RMSE=tabela$RMSE
    a=ggplot(tabela,aes(y=modelo,x=AIC))+
      geom_col(aes(fill=modelo),color="black",
               show.legend = FALSE,size=1)+
      geom_label(aes(x=(AIC/2),label=round(AIC,round.label)),show.legend = FALSE)+
      theme_bw()+labs(y="Models")+
      theme(axis.text = element_text(size=12,color="black"),
            axis.title = element_text(size=12))
    b=ggplot(tabela,aes(y=modelo,x=BIC))+
        geom_col(aes(fill=modelo),color="black",
                 show.legend = FALSE,size=1)+
      geom_label(aes(x=(BIC/2),label=round(BIC,round.label)),
                 show.legend = FALSE)+
      theme_bw()+labs(y="Models")+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=12,color="black"),
            axis.title.x = element_text(size=12))
    c=ggplot(tabela,aes(y=modelo,x=R2))+
        geom_col(aes(fill=modelo),color="black",
                 show.legend = FALSE,size=1)+
        theme_bw()+labs(y="Models",x=expression(R^2))+
      geom_label(aes(x=(R2/2),label=round(R2,2)),show.legend = FALSE)+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=12,color="black"),
            axis.title.x = element_text(size=12))
    d=ggplot(tabela,aes(y=modelo,x=RMSE))+
        geom_col(aes(fill=modelo),color="black",
                 show.legend = FALSE,size=1)+
        theme_bw()+labs(y="Models")+
      geom_label(aes(x=(RMSE/2),label=round(RMSE,round.label)),show.legend = FALSE)+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=12,color="black"),
            axis.title.x = element_text(size=12))
    egg::ggarrange(a,b,c,d,nrow=1)}
  tabela}
