\name{EffToxSelect}
\docType{class}
\alias{EffToxSelect}
\alias{plot, EffToxSelect-method}
\alias{superiority, EffToxSelect-method}

\title{EffToxSelect Class and Methods}

\description{
An R6 class for an EffTox dose selection analysis.
}

\section{Constructor}{
\preformatted{EffToxSelect$new(design, level, utility, samples, mpsrf)}
\subsection{Arguments}{
\describe{
\item{\code{design}}{\emph{\link{EffToxDesign}} upon which the simulations were based.}

\item{\code{level}}{selected dosing level.}

\item{\code{utility}}{vector of dose utility values.}

\item{\code{samples}}{MCMC sampling results from a call to \code{\link[rstan:stanmodel-method-sampling]{sampling}}.}

\item{\code{mpsrf}}{Gelman and Brooks multivariate potential scale reduction factor convergence assessment diagnostic.}
}
}
\subsection{Value}{
An \emph{EffToxSelect} object.
}
}

\section{Methods}{
\describe{

\item{\code{$plot(par = c("utility", "eff", "tox"),
                  type = c("density", "curve"), prob = 0.95)}}{
Posterior plots.
\describe{
\item{\code{par}}{character string indicating the model parameters to plot.}
\item{\code{type}}{character string indicating the type of plot.}
\item{\code{prob}}{credible interval probability for curve plots.}
}
}

\item{\code{$dose}}{
Selected dose.
}

\item{\code{$level}}{
Selected dosing level.
}

\item{\code{$superiority()}}{
Dose superiority matrix such that the value in row i and col j is the posterior probability of dose i exceeding dose j.
}

\item{\code{$utility}}{
Vector of dose utility values.
}

}
}

\seealso{
\code{\link{EffToxDesign}}
}

\keyword{classes}
