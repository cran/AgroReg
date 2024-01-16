# AgroReg 1.2.9

* Fix `regression` bug for Cedergreen-Ritz-Streibig model.

* The third element of the output has been named "plot" and function outputs to models no longer return a print.

* The shiny apps are now implemented on the server of the State University of Londrina.

# AgroReg 1.2.7

* Added `linetype` argument to functions

* `adjusted_scale_x`, `adjusted_scale_y` and `adjusted_scale` functions have been added

* Fix `quadratic.plateau` bug  

* The `plot_arrange` function bug when adding a number greater than 10 curves has been fixed.

# AgroReg 1.2.5

* Fixed the bug of the `LM` function, where when the quantitative factor is negative and when using the midpoint, the curve is generated in reverse order.

* The `plot_arrange` function has been improved

* Added colorline and fillshape arguments

# AgroReg 1.2.4

* Fix `stat_param` bug in me, mse and rmse

* The function `comparative_model` now returns a data.frame object

* Fix `regression` bug  

# AgroReg 1.2.3

* Fix breakpoint response bug in `linear.linear` function

* Fix regression "N" model bug  

# AgroReg 1.2.1

* Bug fix of LM23i, LM2i3, LM13, LM23, LM13i functions when object is grouped by plot_arrange function

* Change argument name from LM and LM_i and from "grau" to degree

* Nreg now accepts average line when grouped in plot_arrange

* Fixed the plot_arrange function curve grouping bug when gray=TRUE and point="mean" 

* `correlation` function was added

# AgroReg 1.2.0

* The Michaelis-Menten function now has a choice of two (mm2) or three parameters (mm3) 

* Fixed decimal place issue when options(OutDec=",") for LM function. 

* Correction of the equation for the `mitscherlich` function 

* Modification of the beta_reg function example 

* The `plot_arrange` function now accepts plotting all points by setting argument point="all". 

* `plateau.linear` function was added

* `plateau.quadratic` function was added

* `yieldloss` function was added

* `lorentz` function was added

* `LM3` function was added

* `LM13` function was added

* `hill` function was added

* `AM` function was added

* `GP` function was added

* `LM_i` function was added

* `asymptotic_i` function was added

* `asymptotic_ineg` function was added

* `LOG2` function was added

* `midilli` function was added

* `midillim` function was added

* `newton` function was added

* `PAGE` function was added

* `peleg` function was added

* `potential` function was added

* `SH` function was added

* `thompson` function was added

* `valcam` function was added

* `VG` function was added

* `weibull` function was added

* `comparative_model` function was added

* `AM` function was added

* The `exponential` and `exponential_neg` function has been changed to `asymptotic` and `asymptotic_neg` 


# AgroReg 1.1.0

* The `width.bar`, `textsize`, `pointsize`, `linesize`, `pointshape`,`comment` arguments have been added to the parsing functions. 

* `asymptotic` function was added 

* `mitscherlich` function was added 

* `biexponential` function bug fix

* Correction of `exponential` function equation
