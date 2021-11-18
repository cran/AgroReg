#' Analysis: Regression linear or nonlinear
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
#' 1. **N:** Graph for not significant trend.
#' 2. **loess0:** Loess non-parametric degree 0
#' 3. **loess1:** Loess non-parametric degree 1
#' 4. **loess2:** Loess non-parametric degree 2
#' 5. **LM0.5:** Quadratic inverse
#' 6. **LM1:** Linear regression.
#' 7. **LM2:** Quadratic
#' 8. **LM3:** Cubic
#' 9. **LM4:** Quartic
#' 10. **LM0.5_i:** Quadratic inverse without intercept.
#' 11. **LM1_i:** Linear without intercept.
#' 12. **LM2_i:** Quadratic regression without intercept.
#' 13. **LM3_i:** Cubic without intercept.
#' 14. **LM4_i:** Quartic without intercept.
#' 15. **LM13:** Cubic without beta2
#' 16. **LM13i:** Cubic inverse without beta2
#' 17. **LM23:** Cubic without beta1
#' 18. **LM23i:** Cubic inverse without beta2
#' 19. **LM2i3:** Cubic without beta1, with inverse beta3
#' 20. **valcam:** Valcam
#' 21. **L3:** Three-parameter logistics.
#' 22. **L4:** Four-parameter logistics.
#' 23. **L5:** Five-parameter logistics.
#' 24. **LL3:** Three-parameter log-logistics.
#' 25. **LL4:** Four-parameter log-logistics.
#' 26. **LL5:** Five-parameter log-logistics.
#' 27. **BC4:** Brain-Cousens with four parameter.
#' 28. **BC5:** Brain-Cousens with five parameter.
#' 29. **CD4:** Cedergreen-Ritz-Streibig with four parameter.
#' 30. **CD5:** Cedergreen-Ritz-Streibig with five parameter.
#' 31. **weibull3:** Weibull with three parameter.
#' 32. **weibull4:** Weibull with four parameter.
#' 33. **GP2:** Gompertz with two parameter.
#' 34. **GP3:** Gompertz with three parameter.
#' 35. **GP4:** Gompertz with four parameter.
#' 36. **VB:** Von Bertalanffy
#' 37. **lo3:** Lorentz with three parameter
#' 38. **lo4:** Lorentz with four parameter
#' 39. **beta:** Beta
#' 40. **gaussian3:** Analogous to the Gaussian model/Bragg with three parameters.
#' 41. **gaussian4:** Analogous to the Gaussian model/Bragg with four parameters.
#' 42. **linear.linear:** Linear-linear
#' 43. **linear.plateau:** Linear-plateau
#' 44. **quadratic.plateau:** Quadratic-plateau
#' 45. **plateau.linear:** Plateau-linear
#' 46. **plateau.quadratic:** Plateau-Quadratic
#' 47. **log:** Logarithmic
#' 48. **log2:** Logarithmic quadratic
#' 49. **thompson:** Thompson
#' 50. **asymptotic:** Exponential
#' 51. **asymptotic_neg:** Exponential negative
#' 52. **asymptotic_i:** Exponential without intercept.
#' 53. **asymptotic_ineg:** Exponential negative without intercept.
#' 54. **biexponential:** Biexponential
#' 55. **mitscherlich:** Mitscherlich
#' 56. **yieldloss:** Yield-loss
#' 57. **hill:** Hill
#' 58. **MM2:** Michaelis-Menten with two parameter.
#' 59. **MM3:** Michaelis-Menten with three parameter.
#' 60. **SH:** Steinhart-Hart
#' 61. **page:** Page
#' 62. **newton:** Newton
#' 63. **potential:** Potential
#' 64. **midilli:** Midilli
#' 65. **midillim:** Modified Midilli
#' 66. **AM:** Avhad and Marchetti
#' 67. **peleg:** Peleg
#' 68. **VG:** Vega-Galvez

#
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

  if(model=="LM0.5"){a=LM(trat, resp, degree=0.5, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM1"){a=LM(trat, resp, degree=1, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM2"){a=LM(trat, resp, degree=2, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM3"){a=LM(trat, resp, degree=3, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM4"){a=LM(trat, resp, degree=4, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}

  if(model=="LM0.5_i"){a=LM_i(trat, resp, degree=0.5, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM1_i"){a=LM_i(trat, resp, degree=1, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM2_i"){a=LM_i(trat, resp, degree=2, ylab = ylab, xlab = xlab, theme = theme, legend.position = legend.position,point = point)}
  if(model=="LM3_i"){a=LM_i(trat, resp, degree=3, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM4_i"){a=LM_i(trat, resp, degree=4, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}

  if(model=="LM13"){a=LM13(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM13i"){a=LM13i(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM23"){a=LM23(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM23i"){a=LM23i(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}
  if(model=="LM2i3"){a=LM2i3(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position =  legend.position,point = point)}

  if(model=="L3"){a=logistic(trat, resp, npar="L.3", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="L4"){a=logistic(trat, resp, npar="L.4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="L5"){a=logistic(trat, resp, npar="L.5", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="LL3"){a=LL(trat, resp, npar="LL.3", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="LL4"){a=LL(trat, resp, npar="LL.4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="LL5"){a=LL(trat, resp, npar="LL.5", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="BC4"){a=BC(trat, resp, npar="BC.4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="BC5"){a=BC(trat, resp, npar="BC.5", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="CD4"){a=CD(trat, resp, npar="CD.4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="CD5"){a=CD(trat, resp, npar="CD.5", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="weibull3"){a=weibull(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="weibull4"){a=weibull(trat, resp, npar="w4",ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="GP2"){a=GP(trat, resp, npar = "g2", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="GP3"){a=GP(trat, resp, npar = "g3", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="GP4"){a=GP(trat, resp, npar = "g4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="lo3"){a=lorentz(trat, resp, npar = "lo3", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="lo4"){a=lorentz(trat, resp, npar = "lo4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="beta"){a=beta_reg(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="gaussian3"){a=gaussianreg(trat, resp, npar="g3", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="gaussia4"){a=gaussianreg(trat, resp,npar = "g4", ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="linear.linear"){a=linear.linear(trat, resp, ylab = ylab, xlab = xlab, theme = theme,legend.position,point = point)}
  if(model=="linear.plateau"){a=linear.plateau(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="quadratic.plateau"){a=quadratic.plateau(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="plateau.linear"){a=plateau.linear(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="plateau.quadratic"){a=plateau.quadratic(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="log"){a=LOG(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="log2"){a=LOG2(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="thompson"){a=thompson(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="asymptotic"){a=asymptotic(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="asymptotic_neg"){a=asymptotic_neg(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="asymptotic_i"){a=asymptotic_i(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="asymptotic_ineg"){a=asymptotic_ineg(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="biexponential"){a=biexponential(trat, resp,ylab =  ylab, xlab = xlab, theme = theme, legend.position,point = point)}

  if(model=="mitscherlich"){a=mitscherlich(trat, resp, ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="MM2"){a=MM(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="MM3"){a=MM(trat, resp, npar="mm3",ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="loess0"){a=loessreg(trat, resp, degree = 0, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="loess1"){a=loessreg(trat, resp, degree = 1, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="loess2"){a=loessreg(trat, resp, degree = 2, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="SH"){a=SH(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="page"){a=PAGE(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="newton"){a=newton(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="valcam"){a=valcam(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="potential"){a=potential(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="midilli"){a=midilli(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="midillim"){a=midillim(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="hill"){a=hill(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="AM"){a=AM(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="yieldloss"){a=yieldloss(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="VB"){a=VB(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="peleg"){a=peleg(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  if(model=="VG"){a=VG(trat, resp, ylab = ylab, xlab = xlab, theme = theme, legend.position,point = point)}
  a
}
