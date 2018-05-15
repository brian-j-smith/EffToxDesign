EffToxDesign <- R6Class("EffToxDesign",
  public = list(
    doses = NULL,
    K = NULL,

    piE = NULL,
    pEL = NULL,
    piT = NULL,
    pTL = NULL,
    
    pi1E = NULL,
    pi2T = NULL,
    pi3E = NULL,
    pi3T = NULL,
    p = NULL,
    
    muE_mean = NULL,
    muE_sd = NULL,
    betaE1_mean = NULL,
    betaE1_sd = NULL,
    betaE2_mean = NULL,
    betaE2_sd = NULL,
    muT_mean = NULL,
    muT_sd = NULL,
    betaT1_mean = NULL,
    betaT1_sd = NULL,
    psi_mean = NULL,
    psi_sd = NULL,
    
    starting_level = NULL,
    burn_in = NULL,
    cohort_sizes = NULL,
    
    yE = numeric(0),
    yT = numeric(0),
    levels = numeric(0),
    n = 0,
    
    
    initialize = function(doses, piE, pEL, piT, pTL, pi1E, pi2T, pi3E, pi3T,
                          thetaE_mean, thetaE_sd, thetaT_mean, thetaT_sd,
                          psi_mean, psi_sd, cohort_sizes,
                          starting_dose = doses[1], burn_in = 0) {
      self$doses <- doses
      self$K <- length(doses)
      self$piE <- piE
      self$pEL <- pEL
      self$piT <- piT
      self$pTL <- pTL
      self$pi1E <- pi1E
      self$pi2T <- pi2T
      self$pi3E <- pi3E
      self$pi3T <- pi3T
      self$p <- efftox_solve_p(pi1E, pi2T, pi3E, pi3T)
      self$muE_mean <- thetaE_mean[1]
      self$muE_sd <- thetaE_sd[1]
      self$betaE1_mean <- thetaE_mean[2]
      self$betaE1_sd <- thetaE_sd[2]
      self$betaE2_mean <- thetaE_mean[3]
      self$betaE2_sd <- thetaE_sd[3]
      self$muT_mean <- thetaT_mean[1]
      self$muT_sd <- thetaT_sd[1]
      self$betaT1_mean <- thetaT_mean[2]
      self$betaT1_sd <- thetaT_sd[2]
      self$psi_mean <- psi_mean
      self$psi_sd <- psi_sd
      
      self$starting_level <- which(starting_dose == doses)
      if(length(self$starting_level) != 1) stop("Invalid starting dose value")
      
      self$burn_in <- burn_in
      self$cohort_sizes <- cohort_sizes
      invisible(self)
    },
    
    
    add = function(yE, yT, doses, levels = NULL) {
      if(is.null(levels)) {
        levels <- match(doses, self$doses)
        if(any(is.na(levels))) stop("Invalid dose value")
      }
      stopifnot(length(yE) == length(levels) && length(yT) == length(levels))
      self$yE <- c(self$yE, yE)
      self$yT <- c(self$yT, yT)
      self$levels <- c(self$levels, levels)
      self$n <- length(self$levels)
      invisible(self)
    },
    
    
    drop = function(n) {
      self$keep(self$n - n)
      invisible(self)
    },
    
    
    keep = function(n) {
      n <- max(n, 0)
      self$yE <- head(self$yE, n)
      self$yT <- head(self$yT, n)
      self$levels <- head(self$levels, n)
      self$n <- length(self$levels)
      invisible(self)
    },
    
    
    reset = function() {
      self$keep(0)
      invisible(self)
    },
    
    
    ESS = function() {
      
      .ESS <- function(mu, sigma) {
        X <- momentsLogitnorm(mu, sigma)
        a <- X["mean"] * (X["mean"] * (1 - X["mean"]) / X["var"] - 1)
        b <- a / X["mean"] * (1 - X["mean"])
        a + b
      }
      
      x <- log(self$doses) - mean(log(self$doses))
      etaE_mean <- self$muE_mean + self$betaE1_mean * x + self$betaE2_mean * x^2
      etaE_sd <- sqrt(self$muE_sd^2 + (self$betaE1_sd * x)^2 +
                        (self$betaE2_sd * x^2)^2)
      betaT1_trunc_mean <- etruncnorm(0, Inf, self$betaT1_mean, self$betaT1_sd)
      betaT1_trunc_var <- vtruncnorm(0, Inf, self$betaT1_mean, self$betaT1_sd)
      etaT_mean <- self$muT_mean + betaT1_trunc_mean * x
      etaT_sd <- sqrt(self$muT_sd^2 + betaT1_trunc_var * x^2)
      c("eff" = mean(mapply(.ESS, etaE_mean, etaE_sd)),
        "tox" = mean(mapply(.ESS, etaT_mean, etaT_sd)))
    },
    
    
    decision = function(mcmcdiag = FALSE, ...) {
      samples <- sampling(stanmodels$EffTox, data = as.list(self), ...)
      
      # MCMC convergence diagnostic
      mpsrf <- NA
      if(mcmcdiag && (d2 <- dim(samples)[2]) > 1) {
        theta <- extract(
          samples,
          par = c("muE", "betaE1", "betaE2", "muT", "betaT1", "psi"),
          permuted = FALSE
        )
        x <- as.mcmc.list(lapply(1:d2, function(i) as.mcmc(theta[,i,])))
        try(mpsrf <- gelman.diag(x)$mpsrf, silent = TRUE)
      }
      
      # Acceptable dosing levels
      acc_efftox <- if(self$n > self$burn_in) {
        prob_acc_eff <- colMeans(extract(samples, "prob_acc_eff")[[1]])
        prob_acc_tox <- colMeans(extract(samples, "prob_acc_tox")[[1]])
        (prob_acc_eff > self$pEL) & (prob_acc_tox > self$pTL)
      } else {
        TRUE
      }
      dosing_levels <- 1:self$K
      acc_range <- (dosing_levels >= min(self$levels) - 1) &
                   (dosing_levels <= max(self$levels) + 1)
      acceptable <- acc_efftox & acc_range
      
      # Recommended level
      utility <- rep(NA, self$K)
      if(any(acceptable)) {
        prob_eff <- colMeans(extract(samples, "prob_eff")[[1]][, acceptable, drop = FALSE])
        prob_tox <- colMeans(extract(samples, "prob_tox")[[1]][, acceptable, drop = FALSE])
        utility[acceptable] <- efftox_utility(self$p, self$pi1E, self$pi2T,
                                              prob_eff, prob_tox)
      }
      level <- which.max(utility)

      list(dose = self$doses[level], level = level, utilty = utility,
           samples = samples, mpsrf = mpsrf)
    },
    
    
    dtp = function(n = 1, ...) {
      n0 <- self$n
      next_cohorts <- .nextcohorts(self)
      level <- next_cohorts$starting_level
      next_cohort_sizes <- head(next_cohorts$sizes, n)
      
      # Construct all possible outcome combinations
      cohort_paths <- next_cohort_sizes %>%
        lapply(function(size) {
          combinations(4, size, c("E", "T", "N", "B"), repeats.allowed = TRUE) %>%
            apply(1, paste0, collapse="")
        }) %>%
        expand.grid(stringsAsFactors = FALSE)
      
      # Record of dose levels given
      levels_given = matrix(NA, nrow(cohort_paths), ncol(cohort_paths))
      
      cache <- list()
      for(i in 1:nrow(cohort_paths)) {
        dtp <- ""
        cohort_level <- level
        for(j in 1:ncol(cohort_paths)) {
          dtp <- paste0(dtp, " ", cohort_level, cohort_paths[i, j])
          if(dtp %in% names(cache)) {
            cohort_level <- cache[[dtp]]
          } else {
            these_outcomes <- efftox_parse_outcomes(dtp)
            self$add(yE = these_outcomes$yE,
                     yT = these_outcomes$yT,
                     levels = these_outcomes$levels)
            cohort_level <- self$decision(...)$level
            cache[[dtp]] <- cohort_level
            self$keep(n0)
          }
          if(length(cohort_level)) levels_given[i, j] <- cohort_level
        }
      }
      df <- data.frame(L0 = rep(level, nrow(cohort_paths)))
      for(k in 1:ncol(cohort_paths)) {
        df[, paste0("C", k - 1)] = cohort_paths[, k]
        df[, paste0("L", k)] = levels_given[, k]
      }
      return(df)
    },
    
    
    simulate = function(num_sims, true_eff, true_tox,
                        seed = sample.int(.Machine$integer.max, 1),
                        mcmcdiag = FALSE, ...) {
      set.seed(seed)
      seeds <- sample.int(.Machine$integer.max, num_sims)
      outcomes <- foreach(i = 1:num_sims, .export = "self") %dopar% {
        print(paste("EffToxDesign Trial Simulation", i))
        set.seed(seeds[i])
        n0 <- self$n
        next_cohorts <- .nextcohorts(self)
        level <- next_cohorts$starting_level
        next_cohort_sizes <- next_cohorts$sizes
        mpsrf <- NULL
        for(cohort_size in next_cohort_sizes) {
          yE <- rbinom(cohort_size, 1, true_eff[level])
          yT <- rbinom(cohort_size, 1, true_tox[level])
          self$add(yE = yE, yT = yT, levels = rep(level, cohort_size))
          decision <- self$decision(mcmcdiag = mcmcdiag, ...)
          level <- decision$level
          mpsrf <- c(mpsrf, decision$mpsrf)
          if(length(level) == 0) break
        }
        outcome <- list(dose = decision$dose, yE = self$yE, yT = self$yT,
                        doses_given = self$doses[self$levels], mpsrf = mpsrf)
        self$keep(n0)
        return(outcome)
      }
      EffToxSim$new(design = self$clone(), true_eff = true_eff,
                    true_tox = true_tox, outcomes = outcomes)
    }
  )
)


EffToxSim <- R6Class("EffToxSim",
  public = list(
    design = NULL,
    truth = NULL,
    outcomes = NULL,
   
   
    initialize = function(design, true_eff, true_tox, outcomes) {
      self$design <- design
      self$truth <- data.frame(eff = true_eff, tox = true_tox)
      self$outcomes <- outcomes
      invisible(self)
    },
   
   
    summary = function() {
     
      .extract <- function(name) {
        unlist(lapply(self$outcomes, function(outcome) outcome[[name]]))
      }
     
      num_sims <- length(self$outcomes)
      doses_selected <- factor(.extract("dose"), self$design$doses)
      doses_given <- factor(.extract("doses_given"), self$design$doses)
      .simsummary(self$truth$eff,
                  self$truth$tox,
                  table(doses_selected) / num_sims,
                  prop.table(table(doses_given)),
                  tapply(.extract("yE"), doses_given, sum) / num_sims,
                  tapply(.extract("yT"), doses_given, sum) / num_sims,
                  table(doses_given) / num_sims)
    }
  )
)
