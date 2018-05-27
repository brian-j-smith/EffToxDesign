# EffToxDesign: EffTox Phase I/II Clinical Trial Design

R package implementation of an adaptive Bayesian design for dose-finding in combined phase I/II clinical trials based on trade-offs between the probabilities of treatment efficacy and toxicity, as developed by Thall and Cook (2004, Biometrics 60, 684-693).


## Requirements

Use of the package requires that a C++ compiler be installed, such as the one provided by the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) software collection available for Windows.


## Installation

The package can be installed directly from GitHub.  First, install the devtools package in R if not done so already.

```R
install.packages("devtools")
```

Then, execute the following command.

```R
devtools::install_github("brian-j-smith/EffToxDesign")
```

## Usage

```R
library(EffToxDesign)
```

## Example

### Trial Design Parameters

An ``EffToxDesign`` R6 class and method functions are provided for EffTrox trial design.  Designs are specified with the ``new`` method function.

```R
doses <- c(1, 2, 4, 6.6, 10)
design <- EffToxDesign$new(
  # Treatment doses
  doses = doses,

  # Acceptable efficacy and toxicity parameters
  piE = 0.5,
  pEL = 0.1,
  piT = 0.3,
  pTL = 0.1,

  # Points defining the utility function
  pi1E = 0.5,
  pi2T = 0.65,
  pi3E = 0.7,
  pi3T = 0.25,
  
  # Logit model prior hyperparameters
  thetaE_mean = efftox_theta(c(0.2, 0.4, 0.6, 0.8, 0.90), doses, 2),
  thetaE_sd = c(2.5423, 2.4406, 0.2),
  thetaT_mean = efftox_theta(c(0.02, 0.04, 0.06, 0.08, 0.10), doses),
  thetaT_sd = c(3.5487, 3.5018),
  psi_mean = 0,
  psi_sd = 1,
  
  # Starting dose
  starting_dose = doses[2],
  
  # Cohort sizes
  cohort_sizes = rep(3, 13)
)
```

Once created, characteristics of the design can be obtained, including effective samples sizes of the prior specifications and utility contours.

```R
design$ESS()

      eff       tox 
0.8305102 0.6750093 

design$contour()
```

### Simulation Studies

Trials can be simulated for a design in order to estimate its operating characteristics.  With use of the *doParallel* package, simulations can be run in parallel to shorten runtimes.

```R
library(doParallel)
registerDoParallel(cores = 6)

sim <- design$simulate(n = 10,
                       eff = c(0.20, 0.40, 0.60, 0.80, 0.90),
                       tox = c(0.05, 0.10, 0.15, 0.20, 0.40),
                       seed = 123)
sim$summary()

                 1         2          4       6.6        10  all
Pr(eff)  0.2000000 0.4000000  0.6000000 0.8000000 0.9000000   NA
Pr(tox)  0.0500000 0.1000000  0.1500000 0.2000000 0.4000000   NA
selected 0.0000000 0.1000000  0.4000000 0.2000000 0.3000000  1.0
given    0.1230769 0.1692308  0.4076923 0.1538462 0.1461538  1.0
E(n_eff) 1.3000000 2.6000000  9.5000000 4.9000000 5.1000000 23.4
E(n_tox) 0.1000000 0.8000000  2.9000000 1.3000000 1.5000000  6.6
E(n)     4.8000000 6.6000000 15.9000000 6.0000000 5.7000000 39.0

sim$barplot()
sim$barplot(stat = c("E(n_tox)", "E(n)"))
```

### Dose-Transition Pathways

Dose-transition pathways show all possible outcomes (E = efficacy only, T = toxicity only, B = both, N = neither) and the doses that would result from n = 1 or more subsequent cohorts.

```R
design$dtp(n = 1, seed = 123)

   cohort1 cohort2
1     2BBB       1
2     2BBE       1
3     2BBN       1
4     2BBT       1
5     2BEE       1
6     2BEN       1
7     2BET       1
8     2BNN       1
9     2BNT       1
10    2BTT       1
11    2EEE       2
12    2EEN       3
13    2EET       1
14    2ENN       3
15    2ENT       1
16    2ETT       1
17    2NNN       3
18    2NNT       3
19    2NTT      NA
20    2TTT      NA
```

### Trial Conduct

Efficacy and toxicity events can be added to the design as they are observed, and dose selection performed on the events.

```R
design$add(yE = c(0, 1, 0), yT = c(0, 0, 0), doses = c(2, 2, 2))
selected <- design$select(seed = 123)
selected$dose
selected$density()
selected$superiority()
```

Events can be removed from the design altogether with the ``reset`` method or selectively with the ``drop`` and ``keep`` methods (see documenation).

```R
design$reset()
```
