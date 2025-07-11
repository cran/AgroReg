% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CD_analysis.R
\name{CD}
\alias{CD}
\title{Analysis: Cedergreen-Ritz-Streibig}
\usage{
CD(
  trat,
  resp,
  npar = "CRS.4",
  sample.curve = 1000,
  ylab = "Dependent",
  xlab = "Independent",
  theme = theme_classic(),
  legend.position = "top",
  error = "SE",
  r2 = "all",
  ic = FALSE,
  fill.ic = "gray70",
  alpha.ic = 0.5,
  point = "all",
  width.bar = NA,
  scale = "none",
  textsize = 12,
  pointsize = 4.5,
  linesize = 0.8,
  linetype = 1,
  pointshape = 21,
  fillshape = "gray",
  colorline = "black",
  round = NA,
  xname.formula = "x",
  yname.formula = "y",
  comment = NA,
  fontfamily = "sans",
  print.on = TRUE
)
}
\arguments{
\item{trat}{Numeric vector with dependent variable.}

\item{resp}{Numeric vector with independent variable.}

\item{npar}{Number of model parameters}

\item{sample.curve}{Provide the number of observations to simulate curvature (default is 1000)}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{treatments name (Accepts the \emph{expression}() function)}

\item{theme}{ggplot2 theme (\emph{default} is theme_classic())}

\item{legend.position}{legend position (\emph{default} is "top")}

\item{error}{Error bar (It can be SE - \emph{default}, SD or FALSE)}

\item{r2}{coefficient of determination of the mean or all values (\emph{default} is all)}

\item{ic}{Add interval of confidence}

\item{fill.ic}{Color interval of confidence}

\item{alpha.ic}{confidence interval transparency level}

\item{point}{defines whether you want to plot all points ("all") or only the mean ("mean")}

\item{width.bar}{Bar width}

\item{scale}{Sets x scale (\emph{default} is none, can be "log")}

\item{textsize}{Font size}

\item{pointsize}{shape size}

\item{linesize}{line size}

\item{linetype}{line type}

\item{pointshape}{format point (default is 21)}

\item{fillshape}{Fill shape}

\item{colorline}{Color lines}

\item{round}{round equation}

\item{xname.formula}{Name of x in the equation}

\item{yname.formula}{Name of y in the equation}

\item{comment}{Add text after equation}

\item{fontfamily}{Font family}

\item{print.on}{Print output}
}
\value{
The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
}
\description{
The 'CRS.4' and 'CRS.5' logistical models provide Brain-Cousens modified logistical models to describe u-shaped hormesis. This model was extracted from the 'drc' package.
}
\details{
The four-parameter model is given by the expression:

\deqn{y = 0 + \frac{d-0+f \exp(-1/x)}{1+\exp(b(\log(x)-\log(e)))}}

while the five-parameter is:

\deqn{y = c + \frac{d-c+f \exp(-1/x)}{1+\exp(b(\log(x)-\log(e)))}}
}
\examples{
library(AgroReg)
data("aristolochia")
attach(aristolochia)
CD(trat,resp)
}
\references{
Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).

Ritz, C.; Strebig, J.C.; Ritz, M.C. Package 'drc'. Creative Commons: Mountain View, CA, USA, 2016.
}
\seealso{
\link{LL}, \link{BC}, \link{GP}
}
\author{
Model imported from the drc package (Ritz et al., 2016)

Gabriel Danilo Shimizu

Leandro Simoes Azeredo Goncalves
}
