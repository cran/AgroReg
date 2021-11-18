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
