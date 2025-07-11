% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asymptotic_sup_analysis.R
\name{asymptotic_neg}
\alias{asymptotic_neg}
\title{Analysis: Asymptotic or Exponential Negative}
\usage{
asymptotic_neg(
  trat,
  resp,
  sample.curve = 1000,
  ylab = "Dependent",
  xlab = "Independent",
  theme = theme_classic(),
  legend.position = "top",
  error = "SE",
  r2 = "all",
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

\item{sample.curve}{Provide the number of observations to simulate curvature (default is 1000)}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{treatments name (Accepts the \emph{expression}() function)}

\item{theme}{ggplot2 theme (\emph{default} is theme_bw())}

\item{legend.position}{legend position (\emph{default} is "top")}

\item{error}{Error bar (It can be SE - \emph{default}, SD or FALSE)}

\item{r2}{coefficient of determination of the mean or all values (\emph{default} is all)}

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

\item{print.on}{Print Output}
}
\value{
The function returns a list containing the coefficients and their respective values of p; statistical parameters such as AIC, BIC, pseudo-R2, RMSE (root mean square error); largest and smallest estimated value and the graph using ggplot2 with the equation automatically.
}
\description{
This function performs asymptotic regression analysis.
}
\details{
The asymptotic model is defined by:
\deqn{y = -\alpha \times e^{-\beta \cdot x}+\theta}
}
\examples{
library(AgroReg)
data("granada")
attach(granada)
asymptotic_neg(time,WL)
}
\references{
Seber, G. A. F. and Wild, C. J (1989) Nonlinear Regression, New York: Wiley & Sons (p. 330).
}
\author{
Gabriel Danilo Shimizu

Leandro Simoes Azeredo Goncalves
}
