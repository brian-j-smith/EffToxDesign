# Part of the EffToxDesign package for estimating model parameters
# Copyright (C) 2018 Brian J Smith
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' @title Calculate the p-index for EffTox utility contours
#'
#' @description
#' Calculate the p-index for EffTox utility contours so that the neutral utility
#' contour intersects the following points in the
#' Prob(Efficacy) - Prob(Toxicity) plane:
#' (\code{pi1E}, 0), (1, \code{pi2T}) and (\code{pi3E}, \code{pi3T})
#'
#' @param pi1E Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param pi2T Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#' @param pi3E Efficacy probability of an equi-utility third point
#' @param pi3T Toxicity probability of an equi-utility third point
#'
#' @return The p-index
#'
#' @examples
#' efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
#'
#' @noRd
#'
efftox_solve_p <- function(pi1E, pi2T, pi3E, pi3T) {
  # Calculate p for the efficacy/toxicity contours that will intersect points
  # (pi1E, 0), (eff.star, tox.star), and (1, pi2T)
  
  .objective = function(p, pi1E, pi2T, pi3E, pi3T) {
    a <- ((1 - pi3E) / (1 - pi1E))
    b <- pi3T / pi2T
    return(a^p + b^p - 1)
  }
  
  rt <- stats::uniroot(.objective, interval = c(0, 100),
                       pi1E = pi1E, pi2T = pi2T, pi3E = pi3E,
                       pi3T = pi3T)
  return(rt$root)
}


#' @title Estimate EffTox logit model coefficients
#'
#' @description
#' Estimate regression coefficients for the efficacy or toxicity logit models in
#' an EffTox trial design, given a set of probabilities and doses.
#'
#' @param probs vector of probabilites at which to estimate coefficients.
#' @param doses vector of corresponding doses.
#' @param degree degree of the polynomial dose effect in the logit model.
#'
#' @return A vector of coefficients.
#'
#' @seealso
#' \code{\link{EffToxDesign}}
#'
efftox_theta <- function(probs, doses, degree = 1) {
  x <- log(doses) - mean(log(doses))
  y <- log(probs / (1 - probs))
  structure(lm(y ~ poly(x, degree, raw = TRUE))$coef, names = NULL)
}


#' @title Get parameters to run the EffTox demo
#'
#' @description Get parameters to run the EffTox demo. These match those used
#' to demonstrate EffTox in Thall et al. 2014.
#'
#' @return a \code{list} of parameters, described in \code{efftox_params}
#'
#' @examples
#' design <- efftox_parameters_demo()
#' names(design)
#' design$doses == c(1, 2, 4, 6.6, 10)
#'
#' @seealso
#' \code{\link{efftox_params}}
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
#'
#' @noRd
#'
efftox_parameters_demo <- function() {
  # Demonstration from 'Effective sample size for computing prior
  # hyperparameters in Bayesian phase I-II dose-finding', Thall et al., 2014
  EffToxDesign$new(
    doses = c(1, 2, 4, 6.6, 10),

    piT = 0.3,
    pTL = 0.1,
    piE = 0.5,
    pEL = 0.1,

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
}


#' @title Get the utility of efficacy & toxicity probability pairs
#'
#' @description Get the utility of efficacy & toxicity probability pairs
#'
#' @param p p-index of EffTox utility contours. Use \code{efftox_solve_p}
#' @param pi1E Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param pi2T Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#' @param prob_eff Probability of efficacy; number between 0 and 1
#' @param prob_tox Probability of toxicity; number between 0 and 1
#'
#' @return Utility value(s)
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = 0.7, prob_tox = 0.25)
#' round(u, 4) == 0
#'
#' u <- efftox_utility(p, 0.5, 0.65, prob_eff = c(0.6, 0.7, 0.8), prob_tox = c(0.1, 0.2, 0.3))
#' round(u, 2) == c(0.04, 0.08, 0.12)
#'
#' @seealso \code{\link{efftox_solve_p}}
#'
#' @noRd
#'
efftox_utility <- function(p, pi1E, pi2T, prob_eff, prob_tox) {
  a <- ((1 - prob_eff) / (1 - pi1E))
  b <- prob_tox / pi2T
  r <- (a^p + b^p)^(1/p)
  return(1 - r)
}


#' @noRd
#' @title Parse a string of EffTox outcomes to binary vector notation.
#'
#' @description Parse a string of EffTox outcomes to the binary vector notation
#' required by Stan for model invocation. The outcome string describes the doses
#' given and outcomes observed. The format of the string is described in Brock
#' et al. (2017). The letters E, T, N and B are used to represents patients that
#' experienced (E)fficacy only, (T)oxicity only, (B)oth efficacy and toxicity,
#' and (B)oth. These letters are concatenated after numerical dose-levels to
#' convey the outcomes of cohorts of patients. For instance, \code{2ETB}
#' represents a cohort of three patients that were treated at dose-level 2, and
#' experienced efficacy, toxicity and both events, respectively. The results of
#' cohorts are separated by spaces. Thus, \code{2ETB 1NN} extends our previous
#' example, where the next cohort of two were treated at dose-level 1 and both
#' patients experienced neither efficacy nor toxicity. See examples.
#'
#' We present the notation in the EffTox setting but it is applicable in
#' general seamless phase I/II dose-finding scenarios.
#'
#' @param outcome_string character string, conveying doses given and outcomes
#' observed.
#'
#' @return a list with elements \code{yE}, \code{yT}, \code{levels} and
#' \code{n}. These elements are congruent with the same in
#' \code{efftox_params}.
#'
#' @examples
#' x = efftox_parse_outcomes('1NNE 2EEN 3TBB')
#' x$n == 9
#' x$yE == c(0, 0, 1, 1, 1, 0, 0, 1, 1)
#' sum(x$yT) == 3
#'
#' @references Brock et al. (submitted 2017), Implementing the EffTox
#' Dose-Finding Design in the Matchpoint Trial.
#'
#' @noRd
#'
efftox_parse_outcomes <- function(outcome_string) {
  regex = '\\w*\\d+[ETBN]+\\w*'
  outcomes = outcome_string
  m = gregexpr(pattern = regex, text = outcomes)
  cohorts = regmatches(outcomes, m)
  if(length(cohorts) <= 0) stop(paste('Matching', outcomes,
                                      'to EffTox outcomes failed'))
  cohorts = cohorts[[1]]
  levels = c(); yE = c(); yT = c();
  n = 0
  for(cohort in cohorts) {
    # print(cohort)
    mc = regexpr('(\\d+)', cohort)
    dl_start = mc[1]
    dl_end = mc[1] + attr(mc,"match.length") - 1
    dl = substr(cohort, dl_start, dl_end)
    dl = as.integer(dl)
    outcomes = substr(cohort, dl_end + 1, nchar(cohort))
    outcomes = strsplit(outcomes, "")[[1]]
    these_doses = rep(dl, length(outcomes))
    levels = c(levels, these_doses)
    these_eff = as.integer((outcomes == 'E') | (outcomes == 'B'))
    yE =c(yE, these_eff)
    these_tox = as.integer((outcomes == 'T') | (outcomes == 'B'))
    yT = c(yT, these_tox)
    n = n + length(outcomes)
  }
  return(list(
    levels = levels, yE = yE, yT = yT, n = n
  ))
}
