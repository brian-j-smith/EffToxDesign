EffToxSim <- R6Class("EffToxSim",
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
     
      .extract <- function(name) {
        unlist(lapply(self$outcomes, function(outcome) outcome[[name]]))
      }
     
      n <- length(self$outcomes)
      doses_selected <- factor(.extract("dose"), self$design$doses)
      doses_given <- factor(.extract("doses_given"), self$design$doses)
      .simsummary(self$truth$eff,
                  self$truth$tox,
                  table(doses_selected) / n,
                  prop.table(table(doses_given)),
                  tapply(.extract("yE"), doses_given, sum) / n,
                  tapply(.extract("yT"), doses_given, sum) / n,
                  table(doses_given) / n)
    }
  )
)
