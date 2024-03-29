% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color_edition.R
\name{coloredit_arrange}
\alias{coloredit_arrange}
\title{Change the colors of a graph from the plot_arrange function}
\usage{
coloredit_arrange(graphs, color = NA)
}
\arguments{
\item{graphs}{object from a plot_arrange function}

\item{color}{color curve and point}
}
\value{
The function changes the colors of a graph coming from the plot_arrange function
}
\description{
Change the colors of a graph from the plot_arrange function
}
\examples{
library(AgroReg)
data("aristolochia")
attach(aristolochia)
graph1=LM(trat,resp)
graph2=LL(trat,resp,npar = "LL.3")
graph=plot_arrange(list(graph1,graph2))
coloredit_arrange(graph,color=c("red","blue"))
}
\author{
Gabriel Danilo Shimizu
}
