\name{TwoPhaseSim}
\alias{TwoPhaseSim}
\alias{summary}

\title{TwoPhaseSim Class and Methods}

\description{
An R6 class for simulated phase I 3+3 design to establish maximum tolerable dose (MTD)
followed by a phase II expansion of the MTD cohort to assess efficacy.
}

\usage{
TwoPhaseSim$new(design, eff, tox, outcomes)
}

\arguments{
\item{design}{\emph{\link{TwoPhaseDesign}} upon which the simulations were based.}

\item{eff}{vector of true efficacy rates in the simulations.}

\item{tox}{vector of true toxicity rates in the simulations.}

\item{outcomes}{list of simulated dose-specific selections probabilities in the first and second trial phases and the average number of patients.}
}

\value{
An \emph{TwoPhaseSim} object.
}

\section{Methods}{
\describe{

\item{\code{$summary()}}{
Simulation-based summary of the operating characteristics for the trial design.
}

}
}

\seealso{
\code{\link{TwoPhaseDesign}}
}

\keyword{classes}