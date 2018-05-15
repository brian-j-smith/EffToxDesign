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
#' @export
#'
#' @examples
#' efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' @references Thall, Herrick, Nguyen, Venier & Norris. 2014, Effective sample
#' size for computing prior hyperparameters in Bayesian phase I-II dose-finding
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
#' @export
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
#' @export
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
#' @export
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
efftox_utility <- function(p, pi1E, pi2T, prob_eff, prob_tox) {
  a <- ((1 - prob_eff) / (1 - pi1E))
  b <- prob_tox / pi2T
  r = (a^p + b^p) ^ (1/p)
  return(1 - r)
}


#' @title Process RStan samples from an EffTox model
#'
#' @description Process RStan samples from an EffTox model to make inferences
#' about dose-acceptability, dose-utility and which dose should be recommended
#' next.
#'
#' @param design An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model.
#' @return An instance of \code{\link{efftox_analysis}}.
#' @export
#'
#' @examples
#' design <- efftox_parameters_demo()
#' design$n <- 3
#' design$yE <- c(0, 1, 1)
#' design$yT <- c(0, 0, 1)
#' design$levels <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = as.list(design))
#' decision <- efftox_process(design, fit)
#' decision$recommended_dose == 3
#' @seealso
#' \code{\link{efftox_params}}
#'
#' \code{\link{efftox_parameters_demo}}
#'
efftox_process <- function(design, fit) {
  # Posterior mean estimates
  prob_eff <- colMeans(rstan::extract(fit, 'prob_eff')[[1]])
  prob_acc_eff <- colMeans(rstan::extract(fit, 'prob_acc_eff')[[1]])
  prob_tox <- colMeans(rstan::extract(fit, 'prob_tox')[[1]])
  prob_acc_tox <- colMeans(rstan::extract(fit, 'prob_acc_tox')[[1]])
  post_utility <- colMeans(rstan::extract(fit, 'utility')[[1]])
  # Derived quantities
  utility = efftox_utility(design$p, design$pi1E, design$pi2T,
                           prob_eff, prob_tox)
  # Dose admissibility
  dose_indices <- seq(design$doses)
  lowest <- min(design$levels)
  highest <- max(design$levels)
  in_range <- sapply(dose_indices,
                     function(x) (x >= lowest - 1) & (x <= highest + 1))
  acceptable <- (prob_acc_eff > design$pEL) & (prob_acc_tox > design$pTL) & in_range
  if(sum(acceptable) > 0) {
    recommended_dose <- which.max(ifelse(acceptable, utility, NA))  # 2
  } else {
    recommended_dose <- NA
  }
  
  l <- list(dose_indices = dose_indices, recommended_dose = recommended_dose,
            prob_eff = prob_eff, prob_tox = prob_tox,
            prob_acc_eff = prob_acc_eff, prob_acc_tox = prob_acc_tox,
            utility = utility, post_utility = post_utility,
            acceptable = acceptable)
  class(l) <- "efftox_analysis"
  return(l)
}


#' @title EffTox analysis to data.frame
#'
#' @description Convenient function to turn an \code{\link{efftox_analysis}}
#' into a \code{data.frame}.
#'
#' @param x An \code{\link{efftox_analysis}}
#'
#' @return a \code{data.frame}
#' @export
#'
#' @examples
#' design <- efftox_parameters_demo()
#' design$n <- 3
#' design$yE <- c(0, 1, 1)
#' design$yT <- c(0, 0, 1)
#' design$levels <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = as.list(design))
#' decision <- efftox_process(design, fit)
#' df <- efftox_analysis_to_df(decision)
#' round(df$Utility, 2) == c(-0.64, 0.04, 0.24, -0.05, -0.19)
#'
efftox_analysis_to_df <- function(x) {
  df <- data.frame(DoseLevel = factor(x$dose_indices),
                   ProbEff = x$prob_eff, ProbTox = x$prob_tox,
                   ProbAccEff = x$prob_acc_eff, ProbAccTox = x$prob_acc_tox,
                   Utility = x$utility, Acceptable = x$acceptable)
  return(df)
}


#' @title Run EffTox simulations
#'
#' @description Run EffTox simulations for assumed true efficacy and toxicity
#'
#' @param design An instance of \code{\link{efftox_params}}, a list of EffTox
#' parameters. An example is yielded by \code{\link{efftox_parameters_demo}}.
#' @param num_sims integer, number of simulated iterations
#' @param first_dose integer, the dose-level to give to patient 1, e.g. 1 for
#' the lowest dose.
#' @param true_eff the true probabilities of efficacy at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param true_tox the true probabilities of toxicity at the doses under
#' investigation; a vector of numbers between 0 and 1.
#' @param cohort_sizes a vector of integer cohort sizes. A dose decision is made
#' when each cohort is completed and the next cohort is treated at the
#' recommended dose. To conduct a trial using at most 20 patients, where dose is
#' re-evaluated after every second patient, use \code{rep(2, 10)}. To conduct a
#' trial of 8 patients where dose is re-evaluated after each single patient, use
#' \code{rep(1, 8)}. Cohort size need not be uniform. E.g.
#' \code{c(rep(1, 5), rep(3, 10))} represents a trial where the dose is
#' re-evaluated after each patient for the first 5 patients, and then after
#' every third patient for a further 30 patients.
#' @param ... Extra parameters provided via the ellipsis are passed to
#' \code{stan::sampling}
#'
#' @return A list with named elements \code{recommended_dose},
#' \code{efficacies}, \code{toxicities}, and \code{doses_given}.
#' @export
#'
#' @examples
#' design <- efftox_parameters_demo()
#' set.seed(123)
#' # Let's say we want to use only 2 chains. Extra args are passed to stan
#' \dontrun{
#' sims <- efftox_simulate(design, num_sims = 10, first_dose = 1,
#'                         true_eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
#'                         true_tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
#'                         cohort_sizes = rep(3, 13),
#'                         chains = 2)
#' table(sims$recommended_dose) / length(sims$recommended_dose)
#' table(unlist(sims$doses_given)) / length(unlist(sims$doses_given))
#' table(unlist(sims$doses_given)) / length(sims$recommended_dose)
#' }
#' # In real life, we would run thousands of iterations, not 10.
#' # This is an example.
#'
efftox_simulate <- function(design, num_sims, first_dose, true_eff, true_tox,
                            cohort_sizes, ...) {
  
  recommended_dose <- integer(length = num_sims)
  efficacies <- list()
  toxicities <- list()
  doses_given <- list()
  
  for(i in 1:num_sims) {
    print(paste('Starting iteration', i))
    design <- as.list(design)
    dose <- first_dose
    for(cohort_size in cohort_sizes) {
      prob_eff <- true_eff[dose]
      prob_tox <- true_tox[dose]
      # Simulate new efficacy events
      new_eff <- stats::rbinom(n = cohort_size, size = 1, prob = prob_eff)
      new_tox <- stats::rbinom(n = cohort_size, size = 1, prob = prob_tox)
      # And append to trial data
      design$yE <- c(design$yE, new_eff)
      design$yT <- c(design$yT, new_tox)
      # Also reflect doses delivered
      design$levels <- c(design$levels, rep(dose, cohort_size))
      design$n <- design$n + cohort_size
      samp <- rstan::sampling(stanmodels$EffTox, data = design, ...)
      l <- efftox_process(design, samp)
      # Select a dose?
      if(sum(l$acceptable) > 0) {
        # Select dose
        dose = which.max(ifelse(l$acceptable, l$utility, NA))
      } else {
        dose <- NA
        break()
      }
    }
    recommended_dose[i] = dose
    efficacies[[i]] = design$yE
    toxicities[[i]] = design$yT
    doses_given[[i]] = design$levels
  }
  
  return(list(recommended_dose = recommended_dose,
              efficacies = efficacies,
              toxicities = toxicities,
              doses_given = doses_given))
}


#' @title Get the probability of toxicity for probability-of-efficacy and utility pairs
#'
#' @description Get the probability of toxicity for probability-of-efficacy and utility pairs
#'
#' @param eff Probability of efficacy; number between 0 and 1
#' @param util Utility score; number
#' @param p p-index of EffTox utility contours. Use \code{efftox_solve_p}
#' @param pi1E Efficacy probability required when toxicity is impossible;
#' a number between 0 and 1
#' @param pi2T Toxicity probability permitted when efficacy is guaranteed;
#' a number between 0 and 1
#'
#' @return Probability(s) of toxicity
#' @export
#'
#' @examples
#' p <- efftox_solve_p(0.5, 0.65, 0.7, 0.25)
#'
#' prob_tox <- efftox_get_tox(0.7, 0, p, pi1E = 0.5, pi2T = 0.65)
#' round(prob_tox, 2) == 0.25
#'
#' prob_tox <- efftox_get_tox(0.7, seq(-0.5, 0.25, by = 0.25), p, pi1E = 0.5, pi2T = 0.65)
#' round(prob_tox, 2) == c(0.57, 0.41, 0.25, 0.09)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), 0.25, p, pi1E = 0.5, pi2T = 0.65)
#' round(prob_tox, 2) == c(NaN, 0.09, 0.22)
#'
#' prob_tox <- efftox_get_tox(c(0.5, 0.7, 0.8), c(-1, 0, 1), p, pi1E = 0.5, pi2T = 0.65)
#' round(prob_tox, 2) == c(0.63, 0.25, NaN)
#'
#' @note Various ways of vectorising the function are demonstrated in the examples
#'
#' @seealso \code{\link{efftox_solve_p}}
#'
efftox_get_tox <- function(eff, util, p, pi1E, pi2T) {
  
  a = ((1 - eff) / (1 - pi1E))
  return(pi2T * ((1 - util)^p - a^p)^(1 / p))
}


#' @title Plot densities of EffTox dose utilities
#'
#' @description Plot densities of EffTox dose utilities. Optionally plot only a
#' subset of the doses by specifying the \code{doses} parameter. This function
#' requires ggplot2 be installed.
#'
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model.
#' @param doses optional, vector of integer dose-levels to plot. E.g. to plot
#' only dose-levels 1, 2 & 3 (and suppress the plotting of any other doses), use
#' \code{doses = 1:3}
#'
#' @return an instance of \code{ggplot}. Omit assignment to just view the plot.
#' @export
#'
#' @note This function requires that ggplot2 be installed.
#'
#' @examples
#' design <- efftox_parameters_demo()
#' design$n <- 3
#' design$yE <- c(0, 1, 1)
#' design$yT <- c(0, 0, 1)
#' design$levels <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = as.list(design))
#' efftox_utility_density_plot(fit) + ggplot2::ggtitle('My doses')  # Too busy?
#' # Specify subset of doses to make plot less cluttered
#' efftox_utility_density_plot(fit, doses = 1:3) + ggplot2::ggtitle('My doses')
#'
efftox_utility_density_plot <- function(fit, doses = NULL) {
  if(!('ggplot2' %in% utils::installed.packages()))
    stop('THis function requires ggplot2 be installed.')
  
  u <- rstan::extract(fit, par = 'utility')[[1]]
  df <- data.frame(Utility = as.numeric(u),
                   D = rep(1:5, each = nrow(u))
  )
  df$Dose = factor(df$D)
  if(!is.null(doses))
    df = df[df$D %in% doses, ]
  Dose <- Utility <- NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Utility, group = Dose,
                                        colour = Dose)) +
    ggplot2::geom_density()
  return(p)
}


#' @title Get dose-superiority matrix in EffTox
#'
#' @description Get a dose-superiority matrix from an EffTox dose analysis.
#' EffTox seeks to choose the dose with the highest utility, thus superiority
#' is inferred by posterior utility. The item in row i, col j is the posterior
#' probability that the utility of dose j exceeds that of dose i.
#'
#' @param fit An instance of \code{rstan::stanmodel}, derived by sampling an
#' EffTox model.
#'
#' @return n by n matrix, where n is number of doses under investigation.
#' The item in row i, col j is the posterior probability that the utility of
#' dose j exceeds that of dose i.
#' @export
#'
#' @examples
#' design <- efftox_parameters_demo()
#' design$n <- 3
#' design$yE <- c(0, 1, 1)
#' design$yT <- c(0, 0, 1)
#' design$levels <- c(1, 2, 3)
#' fit <- rstan::sampling(stanmodels$EffTox, data = as.list(design))
#' sup_mat <- efftox_superiority(fit)
#'
efftox_superiority <- function(fit) {
  u <- rstan::extract(fit, par = 'utility')[[1]]
  superiority_mat <- sapply(1:ncol(u), function(i) sapply(1:ncol(u), function(j)
    mean(u[ , i] > u[ , j])))
  diag(superiority_mat) <- NA
  dimnames(superiority_mat) = list(1:ncol(u), 1:ncol(u))
  return(superiority_mat)
}


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
#' @export
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
