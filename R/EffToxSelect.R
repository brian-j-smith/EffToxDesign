EffToxSelect <- R6Class("EffToxSelect",
  public = list(
    design = NULL,
    dose = NULL,
    level = NULL,
    utility = NULL,
    samples = NULL,
    mpsrf = NULL,
  
  
    initialize = function(design, level, utility, samples, mpsrf) {
      self$design <- design$clone()
      self$dose <- design$doses[level]
      self$level <- level
      self$utility <- utility
      self$samples <- samples
      self$mpsrf <- mpsrf
      invisible(self)
    },
  
  
    density = function(par = "utility") {
      switch(match.arg(par, c("utility", "eff", "tox")),
             "utility" = {
               par <- "utility"
               lab <- "Utility"
               lim <- rep(NA_real_, 2)
             },
             "eff" = {
               par <- "prob_eff"
               lab <- "Pr(Efficacy)"
               lim <- c(0, 1)
             },
             "tox" = {
               par <- "prob_tox"
               lab <- "Pr(Toxicity)"
               lim = c(0, 1)
             }
      )
      x <- extract(self$samples, par = par)[[1]]
      df <- data.frame(x = as.numeric(x))
      df$Dose <- factor(rep(self$design$doses, each = nrow(x)))
      ggplot(df, aes(x = x, group = Dose, color = Dose)) +
        geom_density() +
        ylab("Density") +
        xlab(lab) +
        xlim(lim[1], lim[2])
    },
  
  
    superiority = function() {
      utility <- extract(self$samples, par = "utility")[[1]]
      n <- ncol(utility)
      expand.grid(1:n, 1:n) %>%
        apply(1, function(ind) mean(utility[,ind[1]] > utility[,ind[2]])) %>%
        matrix(n, n, dimnames = list(self$design$doses, self$design$doses))
    }
  )
)
