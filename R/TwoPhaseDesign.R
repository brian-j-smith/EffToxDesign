TwoPhaseDesign <- R6Class("TwoPhaseDesign",
  public = list(
    doses = NULL,
    starting_level = NULL,
    n = NULL,
    r = NULL,

    
    initialize = function(doses, pi0, pi1, alpha, beta,
                          starting_dose = doses[1], n = NULL, r = NULL) {
      self$doses <- doses
      self$starting_level <- which(starting_dose == doses)
      if(length(self$starting_level) != 1) stop("Invalid starting dose value")
      
      if(is.null(n) || is.null(r)) {
        design2 <- ph2single(pu = pi0, pa = pi1, ep1 = alpha, ep2 = beta,
                             nsoln = 1)
        n <- design2[1, "n"]
        r <- design2[1, "r"]
      }
      self$n <- n
      self$r <- r
      invisible(self)
    },
    
    
    simulate = function(true_eff, true_tox) {
      design1 <- threep3(true_tox, start = self$starting_level)
      mtd <- tapply(design1$prob, design1$mtd, sum)[-1]
      n_average <- design1$n.average + mtd * (self$n - 6)
      TwoPhaseSim$new(
        design = self$clone(),
        true_eff = true_eff,
        true_tox = true_tox,
        outcomes = list(stage1_selected = mtd,
                        stage2_selected = 1 - pbinom(self$r, self$n, true_eff),
                        n_average = n_average)
      )
    }
  )
)


TwoPhaseSim <- R6Class("TwoPhaseSim",
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