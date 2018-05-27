\name{TwoPhaseSim}
\docType{class}
\alias{TwoPhaseSim}
\alias{summary, TwoPhaseSim-method}

\title{TwoPhaseSim Class and Methods}

\description{
An R6 class for simulated phase I 3+3 design to establish maximum tolerable dose (MTD)
followed by a phase II expansion of the MTD cohort to assess efficacy.
}

\section{Constructor}{
\preformatted{TwoPhaseSim$new(design, eff, tox, outcomes)}
\subsection{Arguments}{
\describe{
\item{\code{design}}{\emph{\link{TwoPhaseDesign}} upon which the simulations were based.}

\item{\code{eff}}{vector of true efficacy rates in the simulations.}

\item{\code{tox}}{vector of true toxicity rates in the simulations.}

\item{\code{outcomes}}{list of simulated dose-specific selections probabilities in the first and second trial phases and the average number of patients.}
}
}
\subsection{Value}{
A \emph{TwoPhaseSim} object.
}
}

\section{Methods}{
\describe{

\item{\code{$barplot(stats = "selected")}}{
Bar plots of simulation summary statistics.
\describe{
\item{\code{stats}}{character string or vector of summary statistics to plot as returned by the \code{summary} method function.}
}
}

\item{\code{$summary()}}{
Simulation summary statistics of the operating characteristics for the trial design.
}

}
}

\seealso{
\code{\link{TwoPhaseDesign}}
}

\keyword{classes}