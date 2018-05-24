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
  
  
    density = function() {
      utility <- extract(self$samples, par = "utility")[[1]]
      df <- data.frame(
        Utility = as.numeric(utility),
        Dose = factor(rep(self$design$doses, each = nrow(utility)))
      )
      ggplot(df, aes(x = Utility, group = Dose, color = Dose)) +
        geom_density()
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
