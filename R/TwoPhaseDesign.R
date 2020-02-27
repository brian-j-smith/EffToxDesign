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
    
    
    simulate = function(eff, tox) {
      design1 <- threep3(tox, threep3.start = self$starting_level)
      mtd <- tapply(design1$prob, design1$mtd, sum)[-1]
      n_average <- design1$n.average + mtd * (self$n - 6)
      TwoPhaseSim$new(
        design = self,
        eff = eff,
        tox = tox,
        outcomes = list(stage1_selected = mtd,
                        stage2_selected = 1 - pbinom(self$r, self$n, eff),
                        n_average = n_average)
      )
    }
  )
)
