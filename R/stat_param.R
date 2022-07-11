#' Analysis: Other statistical parameters
#'
#' This function calculates other statistical parameters such as Mean (Bias) Error, Relative Mean (Bias) Error, Mean Absolute Error, Relative Mean Absolute Error, Root Mean Square Error, Relative Root Mean Square Error, Modeling Efficiency, Standard deviation of differences, Coefficient of Residual Mass.
#' @param models List with objects of type analysis
#' @param names_model Names of the models
#' @param round Round numbers
#' @return Returns a table with the statistical parameters for choosing the model.
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
#' stat_param(models=list(a,b,c,d))

stat_param<-function(models,
                     names_model=NA,
                     round=3){
  mbe=c()
  rmbe=c()
  mae=c()
  rmae=c()
  se=c()
  mse=c()
  rmse=c()
  rrmse=c()
  me=c()
  sd=c()
  crm=c()
  AC=c()
  ACU=c()
  ACs=c()
  for(i in 1:length(models)){
  mod=models[[i]][[3]]$plot$model
  resp=models[[i]][[3]]$plot$resp
  mbe[i]=mean((predict(mod)-resp))
  rmbe[i]=mean((predict(mod)-resp))/mean(resp)*100
  mae[i]=1/length(resp)*sum(abs(predict(mod)-resp))
  rmae[i]=(1/length(resp)*sum(abs(predict(mod)-resp)))/mean(resp)*100
  se[i]=sum((predict(mod)-resp)^2)
  mse[i]=se[i]/length(resp)
  rmse[i]=sqrt(mse[i])
  rrmse[i]=rmse[i]/mean(resp)*100
  me[i]=1 - (sum((predict(mod)-resp)^2)/sum((resp - mean(resp))^2))
  sd[i]=sd(predict(mod)-resp)
  crm[i] <- (mean(resp - predict(mod), na.rm = T))/mean(resp)
  respe=predict(mod)
  dif <- respe-resp
  mdif <- mean(dif, na.rm = T)
  sdif <- dif^2
  msdif <- mean(sdif, na.rm = T)
  mmea <- mean(resp, na.rm = T)
  mcal <- mean(respe, na.rm = T)
  SSD <- sum(sdif)
  SPOD <- sum((abs(mcal-mmea)+abs(respe-mcal))*(abs(mcal-mmea)+abs(resp-mmea)))
  b <- sqrt((sum((resp - mean(resp))^2))/(sum((respe - mean(respe))^2)))
  if(!is.na(cor(respe,resp)) & cor(respe, resp) < 0){b <- -b}
  a <- mmea - b * mcal
  dY <- a + b * respe
  dX <- -a/b + (1/b) * resp
  SPDu <- sum((abs(respe - dX)) * (abs(resp - dY)))
  SPDs <- SSD - SPDu
  AC[i] <- 1 - (SSD/SPOD)
  ACU[i] <- 1 - (SPDu/SPOD)
  ACs[i] <- 1 - (SPDs/SPOD)}
  tabela=rbind(round(mbe,round),
               round(rmbe,round),
               round(mae,round),
               round(rmae,round),
               round(se,round),
               round(mse,round),
               round(rmse,round),
               round(rrmse,round),
               round(me,round),
               round(sd,round),
               round(crm,round),
               round(AC,round),
               round(ACU,round),
               round(ACs,round))
  if(is.na(names_model[1])==TRUE){colnames(tabela)=paste("Model",1:length(models))}
  rownames(tabela)=c("Mean Bias Error",
                     "Relative Mean Bias Error",
                     "Mean Absolute Error",
                     "Relative Mean Absolute Error",
                     "Squared error",
                     "Mean squared error",
                     "Root Mean Square Error",
                     "Relative Root Mean Square Error",
                     "Modelling Efficiency",
                     "Standard deviation of differences",
                     "Coefficient of Residual Mass",
                     "Agreement Coefficient",
                     "Unsystematic Agreement Coefficient",
                     "Systematic Agreement Coefficient")
  tabela}
