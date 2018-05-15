\name{EffToxSim}
\alias{EffToxSim}
\alias{summary}

\title{EffToxSim Class and Methods}

\description{
An R6 class for simulated EffTox clinical trials.
}

\usage{
EffToxSim$new(design, eff, tox, outcomes)
}

\arguments{
\item{design}{\emph{\link{EffToxDesign}} upon which the simulations were based.}

\item{eff}{vector of true efficacy rates in the simulations.}

\item{tox}{vector of true toxicity rates in the simulations.}

\item{outcomes}{list of efficacies, toxicities, doses given, and doses selecting in the simulated trials.}
}

\value{
An \emph{EffToxSim} object.
}

\section{Methods}{
\describe{

\item{\code{$summary()}}{
Simulation-based summary of the operating characteristics for the trial design.
}

}
}

\seealso{
\code{\link{EffToxDesign}}
}

\keyword{classes}