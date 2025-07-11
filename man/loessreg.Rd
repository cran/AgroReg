% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loess_analysis.R
\name{loessreg}
\alias{loessreg}
\title{Analysis: loess regression (degree 0, 1 or 2)}
\usage{
loessreg(
  trat,
  resp,
  degree = 2,
  sample.curve = 1000,
  ylab = "Dependent",
  xlab = "Independent",
  theme = theme_classic(),
  legend.position = "top",
  error = "SE",
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
  fontfamily = "sans",
  print.on = TRUE
)
}
\arguments{
\item{trat}{Numeric vector with dependent variable.}

\item{resp}{Numeric vector with independent variable.}

\item{degree}{Degree polynomial (0,1 or 2)}

\item{sample.curve}{Provide the number of observations to simulate curvature (default is 1000)}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{treatments name (Accepts the \emph{expression}() function)}

\item{theme}{ggplot2 theme (\emph{default} is theme_bw())}

\item{legend.position}{legend position (\emph{default} is c(0.3,0.8))}

\item{error}{Error bar (It can be SE - \emph{default}, SD or FALSE)}

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

\item{fontfamily}{Font family}

\item{print.on}{Print output}
}
\value{
The function returns a list containing the loess regression and graph using ggplot2.
}
\description{
Fit a polynomial surface determined by one or more numerical predictors, using local fitting.
}
\examples{
library(AgroReg)
data("aristolochia")
attach(aristolochia)
loessreg(trat,resp)
}
\seealso{
\link{loess}
}
\author{
Gabriel Danilo Shimizu

Leandro Simoes Azeredo Goncalves
}
