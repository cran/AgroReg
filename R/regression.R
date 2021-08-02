#' Analysis: Regression linear or nolinear
#'
#' @description This function is a simplification of all the analysis functions present in the package.
#' @export
#' @param trat Numeric vector with dependent variable.
#' @param resp Numeric vector with independent variable.
#' @param model model regression (\emph{default} is LM1)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param xlab treatments name (Accepts the \emph{expression}() function)
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param legend.position legend position (\emph{default} is c(0.3,0.8))
#' @param point defines whether you want to plot all points ("all") or only the mean ("mean")
#' @details To change the regression model, change the "model" argument to:
#'
#' * **N:** Graph for not significant trend.
#' * **LM0.5:** Quadratic inverse regression.
#' * **LM1:** Linear regression.
#' * **LM2:** Quadratic regression.
#' * **LM3:** Cubic regression.
#' * **LM4:** Quartic regression.
#' * **LL3:** Three-parameter logistics.
#' * **LL4:** Four-parameter logistics.
#' * **BC4:** Logistic regression Brain-Cousens with four parameter.
#' * **BC5:** Logistic regression Brain-Cousens with five parameter.
#' * **CD4:** Logistic regression Cedergreen-Ritz-Streibig with four parameter.
#' * **CD5:** Logistic regression Cedergreen-Ritz-Streibig with fice parameter.
#' * **loess:** Loess non-parametric regression.
#' * **gaussian:** Gaussian model.
#' * **MM:** Michaelis-Menten.
#' * **linear.linear:** Linear-Linear regression.
#' * **linear.plateau:** Linear-plateau regression.
#' * **quadratic.plateau:** Quadratic-plateau regression.
#' * **exponential:** Exponential regression.
#' * **exponential_neg:** Exponential negative regression.
#' * **biexponential:** Biexponential regression.
#' * **log:** Logarithmic regression.
#' * **GP:** Gompertz regression.
#' @return The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
#' @md
#' @examples
#' library(AgroReg)
#' data("aristolochia")
#' attach(aristolochia)
#' regression(trat, resp)

regression=function(trat,
                    resp,
                    model="LM1",
                    ylab="Dependent",
                    xlab="Independent",
                    theme=theme_classic(),
                    legend.position="top",
                    point="all"){
  if(model=="N"){a=Nreg(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="LM0.5"){a=LM(trat, resp, grau=0.5, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM1"){a=LM(trat, resp, grau=1, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM2"){a=LM(trat, resp, grau=2, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM3"){a=LM(trat, resp, grau=3, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM4"){a=LM(trat, resp, grau=4, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LL3"){a=LL(trat, resp, npar="LL.3", ylab, xlab, theme, legend.position,point = point)}
  if(model=="LL4"){a=LL(trat, resp, npar="LL.4", ylab, xlab, theme, legend.position,point = point)}
  if(model=="BC4"){a=BC(trat, resp, npar="BC.4", ylab, xlab, theme, legend.position,point = point)}
  if(model=="BC5"){a=BC(trat, resp, npar="BC.5", ylab, xlab, theme, legend.position,point = point)}
  if(model=="CD4"){a=CD(trat, resp, npar="CD.4", ylab, xlab, theme, legend.position,point = point)}
  if(model=="CD5"){a=CD(trat, resp, npar="CD.5", ylab, xlab, theme, legend.position,point = point)}
  if(model=="loess"){a=loessreg(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="gaussian"){a=gaussianreg(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="linear.linear"){a=linear.linear(trat, resp, ylab, xlab, theme,legend.position,point = point)}
  if(model=="MM"){a=MM(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="linear.plateau"){a=linear.plateau(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="quadratic.plateau"){a=quadratic.plateau(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="exponential"){a=exponential(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="exponential_neg"){a=exponential_neg(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="biexponential"){a=biexponential(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="log"){a=LOG(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  if(model=="GP"){a=GP(trat, resp, ylab, xlab, theme, legend.position,point = point)}
  a
}
