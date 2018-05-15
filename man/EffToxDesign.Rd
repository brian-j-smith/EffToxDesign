\name{EffToxDesign}
\alias{EffToxDesign}
\alias{add}
\alias{drop}
\alias{keep}
\alias{reset}
\alias{decision}
\alias{simulate}

\title{EffToxDesign Class and Methods}

\description{
An R6 class for EffTox clinical trial designs.
}

\usage{
EffToxDesign$new(doses, piE, pEL, piT, pTL, pi1E, pi2T, pi3E, pi3T,
                 muE, betaE1, betaE2, muT, betaT1, psi, cohort_sizes,
                 starting_dose = doses[1], run_in = 0)
}

\arguments{
\item{doses}{vector of positive and strictly increasing numeric dose values under consideration in the trial.}

\item{piE}{fixed lower limit for which dose \eqn{d}{d} is considered to have
acceptable efficacy if \eqn{\Pr[\pi_E(d) > \pi_E | \text{data}] \ge p_{E,L}}{Pr[piE(d) > piE | data] >= pEL}.}

\item{pEL}{small threshold on the probability that the efficacy rate is above the fixed \code{piE} limit.}

\item{piT}{fixed upper limit for which dose \eqn{d}{d} is considered to have acceptable toxicity if \eqn{\Pr[\pi_T(d) < \pi_T | \text{data}] \ge p_{T,L}}{Pr[piT(d) < piT | data] >= pTL}.}

\item{pTL}{small threshold on the probability that the toxicity rate is below the fixed \code{piT} limit.}

\item{pi1E}{first efficacy target rate \eqn{(\pi_{1,E}, 0)}{(pi1E, 0)} defining the efficacy-toxicity trade-off contour.}

\item{pi2T}{second toxicity target rate \eqn{(0, \pi_{2,T})}{(0, pi2T)}.}

\item{pi3E}{third efficacy target rate in the pair \eqn{(\pi_{3,E}, \pi_{3,T})}{(pi3E, pi3T)}.}

\item{pi3T}{third toxicity target rate in the pair \eqn{(\pi_{3,E}, \pi_{3,T})}{(pi3E, pi3T)}.}

\item{p}{shape parameter of the contour.}

\item{muE}{vector of prior mean and standard deviation for the intercept in the efficacy logit model.}

\item{betaE1}{vector of prior mean and standard deviation for the linear effect.}

\item{betaE2}{vector of prior mean and standard deviation for the quadratic effect.}

\item{muT}{vector of prior mean and standard deviation for the intercept in the toxicity logit model.}

\item{betaT1}{vector of prior mean and standard deviation for the linear effect.}

\item{psi_mean}{vector of prior mean and standard deviation for the efficacy-toxicity association parameter.}

\item{starting_dose}{starting dose value to be given to the first cohort of patients.}

\item{run_in}{number of patients to enroll before checking for acceptable toxicity and efficacy.}

\item{cohort_sizes}{vector giving the number of patients treated at each successive selected dose.}
}

\value{
An \emph{EffToxDesign} object.
}

\section{Methods}{
\describe{

\item{\code{$add(yE, yT, doses, levels = NULL)}}{
Add one or more patient outcomes.
\describe{
\item{\code{yE}}{vector of 0-1 efficacy indicators for the patients.}
\item{\code{yT}}{vector of 0-1 toxicity indicators for the patients.}
\item{\code{doses}}{vector of doses given to the patients.}
\item{\code{levels}}{vector of dosing levels (1:K), as an alterenative to doses.}
}
}

\item{\code{$drop(n)}}{
Drop a terminal number of patient outcomes.
\describe{
\item{\code{n}}{number of outcomes to drop.}
}
}

\item{\code{$keep(n)}}{
Keep an initial number of patient outcomes.
\describe{
\item{\code{n}}{number of outcomes to keep.}
}
}

\item{\code{$reset()}}{
Reset the design by removing any existing patient outcomes.
}

\item{\code{$contour(eff = NULL, tox = NULL, num_levels = 10, num_points = 101)}}{
Plot utility contours.
\describe{
\item{\code{eff}}{optional vector of efficacy probabilities to include as points in the plot.}
\item{\code{tox}}{optional vector of toxicity probabilities to include as points in the plot.}
\item{\code{num_levels}}{number of contour levels.}
\item{\code{num_points}}{number of points at which to evaluate the contour levels.}
}
}

\item{\code{$ESS()}}{
Effective sample sizes of the efficacy and toxicity logit model priors.
\subsection{Value}{
Resulting ESS estimates in a vector of components named \dQuote{eff} and \dQuote{tox}.
}
}

\item{\code{$decision(mcmcdiag = FALSE, ...)}}{
Identify an optimal dose based on an EffTox statistical analysis of the current patient outcomes.
\describe{
\item{\code{mcmcdiag}}{logical indicating whether to return convergence diagnostics for the MCMC output upon which analysis results are based.}
\item{\code{...}}{parameters passed to \code{\link[rstan:stanmodel-method-sampling]{sampling}}.}
}
\subsection{Value}{
A list with the following components:
\describe{
\item{\code{dose}}{optimal dose.}
\item{\code{level}}{optimal dosing level (1:K).}
\item{\code{samples}}{MCMC samples for the analysis.}
\item{\code{mpsrf}}{Gelman and Brooks multivariate potential scale reduction factor MCMC convergence diagnostic if \code{mcmcdiag = TRUE} or \code{NA} otherwise.}
}
}
}

\item{\code{$dtp(n = 1, ...)}}{
Calculate dose-transition pathways.
\describe{
\item{\code{n}}{number of future cohorts over which to calculate pathways.}
\item{\code{...}}{parameters passed to \code{\link[rstan:stanmodel-method-sampling]{sampling}}.}
}
}

\item{\code{$simulate(n, eff, tox, seed = sample.int(.Machine$integer.max, 1), mcmcdiag = FALSE, ...)}}{
Simulate an EffTox trial design.
\describe{
\item{\code{n}}{number of trials to simulate.}
\item{\code{eff}}{vector of true efficacy rates for the K doses.}
\item{\code{tox}}{vector of true toxicity rates for the K doses.}
\item{\code{seed}}{integer value for random number generation.}
\item{\code{mcmcdiag}}{logical indicating whether to return convergence diagnostics for the MCMC output upon which analysis results are based.}
\item{\code{...}}{parameters passed to \code{\link[rstan:stanmodel-method-sampling]{sampling}}.}
}
\subsection{Value}{
An \emph{\link{EffToxSim}} object of simulated trials.
}
}

}
}

\references{
Brock K, Billingham L, Copland M, Siddique S, Sirovica M, Yap C (2017). Implementing the EffTox dose-finding design in the Matchpoint trial. BMC Medical Research Methodology, 17(112):1-15.

Thall PF, Cook JD (2004). Dose-finding based on efficacy-toxcity trade-offs. Biometrics, 60:684-693.

Thall PF, Herrick RC, Nguyen HQ, Venier JJ, Norris JC (2014). Effective sample size for computing prior hyperparameters in Bayesian phase I-II dose-finding. Clinical Trials, 11(6):657-666.
}

\seealso{
\code{\link{EffToxSim}}
}

\examples{
design <- EffToxDesign$new(
    doses = c(1, 2, 4, 6.6, 10),

    piE = 0.5,
    pEL = 0.1,
    piT = 0.3,
    pTL = 0.1,

    pi1E = 0.5,
    pi2T = 0.65,
    pi3E = 0.7,
    pi3T = 0.25,
    
    thetaE_mean = c(0.7367, 3.4181, 0),
    thetaE_sd = c(2.5423, 2.4406, 0.2),
    thetaT_mean = c(-7.9593, 1.5482),
    thetaT_sd = c(3.5487, 3.5018),
    psi_mean = 0,
    psi_sd = 1,
    
    cohort_sizes = rep(3, 13)
)
design$contour()
design$ESS()
design$add(yE = c(0, 1, 0, 1, 0, 1),
           yT = c(0, 0, 0, 1, 0, 0),
           doses = c(1, 1, 1, 2, 2, 2))
design$decision()

\dontrun{
design$reset()
design$dtp()
sim <- design$simulate(n = 10,
                       eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
                       tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       chains = 2)
sim$summary()}
}

\keyword{classes}