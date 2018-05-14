\name{TwoPhaseDesign}
\alias{TwoPhaseDesign}
\alias{simulate}

\title{TwoPhaseDesign Class and Methods}

\description{
An R6 class for phase I 3+3 designs to establish maximum tolerable dose (MTD)
followed by phase II expansions of the MTD cohort to assess efficacy.
}

\usage{
TwoPhaseDesign$new(doses, pi0, pi1, alpha, beta, starting_dose = doses[1],
                   n = NULL, r = NULL)
}

\arguments{
\item{doses}{vector of positive and strictly increasing numeric dose values under consideration in the trial.}

\item{pi0}{upper bound on the true efficacy rate in the null hypothesis.}

\item{pi1}{efficacy rate at which the phase II trial is powered.}

\item{alpha}{type I error (false positive) rate of the phase II trial.}

\item{beta}{type II error (1 - power) rate of the phase II trial.}

\item{starting_dose}{starting dose value to be given to the first cohort of patients.}

\item{n}{total number of patients to include in the phase II trial, as an alterantive to pi0, pi1, alpha, and beta.}

\item{r}{number of responses above which the treatment will be considered efficacious, to be specified along with n.}
}

\value{
An \emph{EffToxDesign} object.
}

\section{Methods}{
\describe{

\item{\code{$simulate(true_eff, true_tox)}}{
Simulate a two phase I-II trial design.
\describe{
\item{\code{true_eff}}{vector of true efficacy rates for the K doses.}
\item{\code{true_tox}}{vector of true toxicity rates for the K doses.}
}
\subsection{Value}{
An \emph{\link{TwoPhaseSim}} object of simulated trials.
}
}

}
}

\examples{
design <- TwoPhaseDesign$new(1:5, 0.3, 0.50, 0.10, 0.15)
sim <- design$simulate(true_eff = c(0.3, 0.375, 0.45, 0.525, 0.6),
                       true_tox = c(0.1, 0.1, 0.15, 0.2, 0.3))
sim$summary()
}

\keyword{classes}