\name{EffToxSim}
\docType{class}
\alias{EffToxSim}
\alias{summary, EffToxSim-method}

\title{EffToxSim Class and Methods}

\description{
An R6 class for simulated EffTox clinical trials.
}

\section{Constructor}{
\preformatted{EffToxSim$new(design, eff, tox, outcomes)
}
\subsection{Arguments}{
\describe{
\item{\code{design}}{\emph{\link{EffToxDesign}} upon which the simulations were based.}

\item{\code{eff}}{vector of true efficacy rates in the simulations.}

\item{\code{tox}}{vector of true toxicity rates in the simulations.}

\item{\code{outcomes}}{list of efficacies, toxicities, doses given, and doses selecting in the simulated trials.}
}
}
\subsection{Value}{
An \emph{EffToxSim} object.
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
\code{\link{EffToxDesign}}
}

\keyword{classes}