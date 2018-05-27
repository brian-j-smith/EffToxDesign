TwoPhaseSim <- R6Class("TwoPhaseSim",
  public = list(
    design = NULL,
    truth = NULL,
    outcomes = NULL,
   
   
    initialize = function(design, eff, tox, outcomes) {
      self$design <- design$clone()
      self$truth <- data.frame(eff = eff, tox = tox)
      self$outcomes <- outcomes
      invisible(self)
    },
   
   
    barplot = function(stats = "selected") .simbarplot(self, stats),
    
    
    summary = function() {
      .simsummary(self$truth$eff,
                  self$truth$tox,
                  self$outcomes$stage1_selected * self$outcomes$stage2_selected,
                  prop.table(self$outcomes$n_average),
                  self$outcomes$n_average * self$truth$eff,
                  self$outcomes$n_average * self$truth$tox,
                  self$outcomes$n_average)
    }
  )
)
