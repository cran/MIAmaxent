% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotResp.R, R/plotResp2.R
\name{plotResp}
\alias{plotResp}
\alias{plotResp2}
\title{Plot model response.}
\usage{
plotResp(model, transformations, EV, logscale = FALSE, ...)

plotResp2(model, transformations, EV, logscale = FALSE, ...)
}
\arguments{
\item{model}{The model for which the response is to be plotted. This may be
the object returned by \code{\link{chooseModel}}, or the 'selectedmodel'
returned by \code{\link{selectEV}}.}

\item{transformations}{Transformation functions used to create the derived
variables in the model. I.e. the 'transformations' returned by
\code{\link{deriveVars}}. Equivalently, the full file pathway of the
'transformations.Rdata' file saved as a result of \code{\link{deriveVars}}.}

\item{EV}{Character. Name of the explanatory variable for which the response
curve is to be plotted. Interaction terms not allowed.}

\item{logscale}{Logical. Plot the common logarithm of PRO rather than PRO
itself.}

\item{...}{Arguments to be passed to \code{plot} or \code{barplot} to control
the appearance of the plot. For example: \itemize{ \item \code{lwd} for
line width \item \code{cex.main} for size of plot title \item \code{space}
for space between bars }}
}
\description{
Plots the response of a given model over any of the explanatory variables
(EVs) included in that model. For categorical variables, a bar plot is
returned rather than a line plot. Single-effect response curves show the
response of a model containing the explanatory variable of interest only,
while marginal effect response curves show the response of the model when all
other explanatory variables are held constant at their mean values (cf.
\code{plotResp}, \code{plotResp2}).
}
\section{Functions}{
\itemize{
\item \code{plotResp()}: Plot single-effect model response.

\item \code{plotResp2()}: Plot marginal-effect model response.

}}
\examples{
\dontrun{
# From vignette:
plotResp(grasslandmodel, grasslandDVs$transformations, "pr.bygall")
plotResp(grasslandmodel, grasslandDVs$transformations, "geolmja1")

plotResp2(grasslandmodel, grasslandDVs$transformations, "pr.bygall")
}

}
