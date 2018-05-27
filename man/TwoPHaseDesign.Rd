\name{TwoPhaseDesign}
\docType{class}
\alias{TwoPhaseDesign}
\alias{simulate, TwoPhaseDesign-method}

\title{TwoPhaseDesign Class and Methods}

\description{
An R6 class for phase I 3+3 designs to establish maximum tolerable dose (MTD)
followed by phase II expansions of the MTD cohort to assess efficacy.
}

\section{Constructor}{
\preformatted{TwoPhaseDesign$new(doses, pi0, pi1, alpha, beta, starting_dose = doses[1],
                   n = NULL, r = NULL)}
\subsection{Arguments}{
\describe{
\item{\code{doses}}{vector of positive and strictly increasing numeric dose values under consideration in the trial.}

\item{\code{pi0}}{upper bound on the true efficacy rate in the null hypothesis.}

\item{\code{pi1}}{efficacy rate at which the phase II trial is powered.}

\item{\code{alpha}}{type I error (false positive) rate of the phase II trial.}

\item{\code{beta}}{type II error (1 - power) rate of the phase II trial.}

\item{\code{starting_dose}}{starting dose value to be given to the first cohort of patients.}

\item{\code{n}}{total number of patients to include in the phase II trial, as an alternative to pi0, pi1, alpha, and beta.}

\item{\code{r}}{number of responses above which the treatment will be considered efficacious, to be specified along with n.}
}
}
\subsection{Value}{
A \emph{TwoPhaseDesign} object.
}
}

\section{Methods}{
\describe{

\item{\code{$simulate(eff, tox)}}{
Simulate a two phase I-II trial design.
\describe{
\item{\code{eff}}{vector of true efficacy rates for the K doses.}
\item{\code{tox}}{vector of true toxicity rates for the K doses.}
}
\subsection{Value}{
A \emph{\link{TwoPhaseSim}} object of simulated trials.
}
}

}
}

\seealso{
\code{\link{TwoPhaseSim}}
}

\examples{
design <- TwoPhaseDesign$new(1:5, 0.3, 0.50, 0.10, 0.15)
sim <- design$simulate(eff = c(0.3, 0.375, 0.45, 0.525, 0.6),
                       tox = c(0.1, 0.1, 0.15, 0.2, 0.3))
sim$summary()
sim$barplot()
sim$barplot(stats = c("E(n_tox)", "E(n)"))
}

\keyword{classes}