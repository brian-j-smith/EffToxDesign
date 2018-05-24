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
      stopifnot(all(doses > 0))
      stopifnot(all(diff(doses) > 0))
      self$doses <- doses
      self$K <- length(doses)
      
      self$piE <- piE
      self$pEL <- pEL
      self$piT <- piT
      self$pTL <- pTL
      
      stopifnot(pi3E > pi1E)
      stopifnot(pi3T < pi2T)
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
    
    
    contour = function(eff = NULL, tox = NULL, bins = 8, n = 25) {
      stopifnot(length(eff) == length(tox))
      
      contours <- expand.grid(eff_vals = seq(0, 1, length = n),
                              tox_vals = seq(0, 1, length = n))
      contours$util_vals <- efftox_utility(self$p, self$pi1E, self$pi2T,
                                           contours$eff_vals, contours$tox_vals)

      target_points <- data.frame(
        eff = c(self$pi1E, 1, self$pi3E),
        tox = c(0, self$pi2T, self$pi3T)
      )
      
      plt <- ggplot() +
        geom_contour(data = contours,
                     aes(x = eff_vals, y = tox_vals, z = util_vals),
                     bins = bins) +
        geom_point(data = target_points, aes(x = eff, y = tox),
                   col = 'red', size = 3) +
        xlim(0, 1) +
        ylim(0, 1) +
        xlab("Pr(Efficacy)") +
        ylab("Pr(Toxicity)")
      
      if(length(eff)) {
        user_points <- data.frame(eff = eff, tox = tox, dl = seq(eff))
        plt <- plt + geom_text(data = user_points,
                               aes(x = eff, y = tox, label = dl),
                               col = 'red', size = 4)
      }
      
      return(plt)
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
    
    
    select = function(mcmcdiag = FALSE, ...) {
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
      dosing_levels <- 1:self$K
      acc_efftox <- if(self$n >= self$burn_in) {
        acc_eff <- colMeans(extract(samples, "prob_acc_eff")[[1]]) > self$pEL
        acc_tox <- colMeans(extract(samples, "prob_acc_tox")[[1]]) > self$pTL
        acc_lowest_untried <- (dosing_levels == max(self$levels) + 1) & acc_tox
        (acc_eff & acc_tox) | acc_lowest_untried
      } else {
        TRUE
      }
      acc_range <- (dosing_levels >= min(self$levels) - 1) &
                   (dosing_levels <= max(self$levels) + 1)
      acceptable <- acc_efftox & acc_range
      
      # Recommended level
      utility <- rep(NA, self$K)
      if(any(acceptable)) {
        subsamples <- extract(samples, "utility")[[1]]
        utility[acceptable] <- colMeans(subsamples[, acceptable, drop = FALSE])
      }
      level <- which.max(utility)

      EffToxSelect$new(design = self, level = level, utility = utility,
                       samples = samples, mpsrf = mpsrf)
    },
    
    
    dtp = function(n = 1, ...) {
      num_keep <- self$n
      next_cohorts <- .nextcohorts(self, ...)

      # Construct all possible outcome combinations
      cohort_paths <- next_cohorts$sizes %>%
        head(n) %>%
        lapply(function(size) {
          combinations(4, size, c("E", "T", "N", "B"),
                       repeats.allowed = TRUE) %>%
            apply(1, paste0, collapse="")
        }) %>%
        expand.grid(stringsAsFactors = FALSE)
      
      # Dose-transition pathways
      dtp <- matrix(NA_character_, nrow(cohort_paths), ncol(cohort_paths)) %>%
        data.frame(stringsAsFactors = FALSE)

      starting_level <- next_cohorts$starting_level
      if(length(starting_level) == 0) starting_level <- NA
      cohort_levels <- rep(starting_level, nrow(cohort_paths))      

      cache <- list()
      for(j in 1:ncol(cohort_paths)) {
        for(i in which(!is.na(cohort_levels))) {
          dtp[i, j] <- paste0(cohort_levels[i], cohort_paths[i, j])
          path <- paste(dtp[i, 1:j], collapse = " ")
          if(path %in% names(cache)) {
            cohort_levels[i] <- cache[[path]]
          } else {
            outcomes <- efftox_parse_outcomes(path)
            self$add(yE = outcomes$yE, yT = outcomes$yT, levels = outcomes$levels)
            selected <- self$select(...)
            cohort_levels[i] <-
              ifelse(length(selected$level), selected$level, NA)
            cache[[path]] <- cohort_levels[i]
            self$keep(num_keep)
          }
        }
      }
      dtp <- cbind(dtp, cohort_levels)
      colnames(dtp) <- paste0("cohort",
                              next_cohorts$starting_cohort + seq(dtp) - 1)
      return(dtp)
    },
    
    
    simulate = function(n, eff, tox, seed = sample.int(.Machine$integer.max, 1),
                        mcmcdiag = FALSE, ...) {
      set.seed(seed)
      seeds <- sample.int(.Machine$integer.max, n)
      outcomes <- foreach(i = 1:n, .export = "self") %dopar% {
        print(paste("EffToxDesign Trial Simulation", i))
        set.seed(seeds[i])
        num_keep <- self$n
        next_cohorts <- .nextcohorts(self)
        cohort_level <- next_cohorts$starting_level
        next_cohort_sizes <- next_cohorts$sizes
        mpsrf <- NULL
        for(cohort_size in next_cohort_sizes) {
          yE <- rbinom(cohort_size, 1, eff[cohort_level])
          yT <- rbinom(cohort_size, 1, tox[cohort_level])
          self$add(yE = yE, yT = yT, levels = rep(cohort_level, cohort_size))
          selected <- self$select(mcmcdiag = mcmcdiag, ...)
          cohort_level <- selected$level
          mpsrf <- c(mpsrf, selected$mpsrf)
          if(length(cohort_level) == 0) break
        }
        outcome <- list(dose = selected$dose, yE = self$yE, yT = self$yT,
                        doses_given = self$doses[self$levels], mpsrf = mpsrf)
        self$keep(num_keep)
        return(outcome)
      }
      EffToxSim$new(design = self$clone(), eff = eff, tox = tox,
                    outcomes = outcomes)
    }
  )
)
