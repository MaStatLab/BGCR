% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BGCR.R
\name{plot.BGCR}
\alias{plot.BGCR}
\title{Visualizing the BGCR results on the phylogenetic tree}
\usage{
\method{plot}{BGCR}(x, BF = FALSE, ind = FALSE, cex = 0.3,
  main = "PMAP", legend = TRUE, ...)
}
\arguments{
\item{x}{a BGCR object returned by the \code{"BGCR"} function.}

\item{BF}{logicals. If true, the bayes factors are plotted at each internal node.}

\item{ind}{logicals. If true, plot the PMAPs under BCR, the model without the graphical structure.}

\item{cex}{a positive value controls the size of the node label.}

\item{main}{a string contains the title of the plot.}

\item{legend}{logicals. If true, plot the legend.}

\item{...}{further arguments to be passed to plot or to plot.BGCR.}
}
\description{
This function plots the BGCR results on the phylogenetic tree.
}
