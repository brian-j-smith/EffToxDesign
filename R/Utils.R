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

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) Rcpp::loadModule(m, what = TRUE)
}


.nextcohorts <- function(design, ...) {
  n <- design$n
  if(n == 0) {
    starting_cohort <- 1
    starting_level <- design$starting_level
    sizes <- design$cohort_sizes
  } else {
    cum_cohort_sizes <- cumsum(design$cohort_sizes)
    starting_cohort <- min(which(n <= cum_cohort_sizes))
    if(cum_cohort_sizes[starting_cohort] == n) {
      starting_cohort <- starting_cohort + 1
      starting_level <- design$decision(...)$level
    } else {
      starting_level <- tail(design$levels, 1)
    }   
    sizes <- c(cum_cohort_sizes[starting_cohort] - n,
               tail(design$cohort_sizes, -starting_cohort))
  }
  list(starting_cohort = starting_cohort, starting_level = starting_level,
       sizes = sizes)
}


.simsummary <- function(eff, tox, selected, given, n_eff, n_tox, n) {
  x <- rbind("Pr(eff)" = eff,
             "Pr(tox)" = tox,
             "selected" = selected,
             "given" = given,
             "E(eff)" = n_eff,
             "E(tox)" = n_tox,
             "E(n)" = n)
  cbind(x, "all" = c(NA, NA, rowSums(x[-(1:2), ])))
}
