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
    
    
    plot = function(par = c("utility", "eff", "tox"),
                    type = c("density", "curve"), prob = 0.95) {
      par = match.arg(par)
      type = match.arg(type)
      
      switch(par,
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
      
      if(type == "density") {
        df <- data.frame(x = as.numeric(x))
        df$Dose <- factor(rep(self$design$doses, each = nrow(x)))
        ggplot(df, aes(x = x, group = Dose, color = Dose)) +
          geom_density() +
          ylab("Density") +
          xlab(lab) +
          xlim(lim[1], lim[2])
      } else if(type == "curve") {
        alpha <- (1 - prob) / 2
        df <- data.frame(
          Dose = self$design$doses,
          Mean = apply(x, 2, mean),
          Lower = apply(x, 2, quantile, probs = alpha),
          Upper = apply(x, 2, quantile, probs = 1 - alpha)
        )
        interp <- function(x, y) {
          fit <- spline(x, y, method = "natural", n = 101)
          data.frame(
            x = fit$x,
            y = pmin(pmax(fit$y, lim[1], na.rm = TRUE), lim[2], na.rm = TRUE)
          )
        }
        ggplot(mapping = aes(x, y)) +
          geom_line(data = interp(df$Dose, df$Mean)) +
          geom_line(data = interp(df$Dose, df$Lower),
                    linetype = "dashed", color = "red") +
          geom_line(data = interp(df$Dose, df$Upper),
                    linetype = "dashed", color = "red") +
          ylab(lab) +
          xlab("Dose") +
          ylim(lim[1], lim[2])
      }
    },
  
  
    density = function(par = "utility") {
      warning("density() is deprecated; use plot() instead")
      self$plot(par)
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
