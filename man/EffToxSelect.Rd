\name{EffToxSelect}
\alias{EffToxSelect}
\alias{density}
\alias{superiority}

\title{EffToxSelect Class and Methods}

\description{
An R6 class for an EffTox dose selection analysis.
}

\usage{
EffToxSelect$new(design, level, utility, samples, mpsrf)
}

\arguments{
\item{design}{\emph{\link{EffToxDesign}} upon which the simulations were based.}

\item{level}{selected dosing level.}

\item{utility}{vector of dose utility values.}

\item{samples}{MCMC sampling results from a call to \code{\link[rstan:stanmodel-method-sampling]{sampling}}.}

\item{mpsrf}{Gelman and Brooks multivariate potential scale reduction factor convergence assessment diagnostic.}
}

\value{
An \emph{EffToxSelect} object.
}

\section{Fields and Methods}{
\describe{

\item{\code{$density()}}{
Dose utility density plots.
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
