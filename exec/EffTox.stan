// Dose-finding based on efficacy-toxicity trade-offs (Thall, Cook 2004)
// Effective sample size for computing prior hyperparameters in Bayesian
//   phase I-II dose-finding (Thall, Herrick, Nguyen, Venier, Norris 2014)

functions {
  real log_joint_pdf(real[] doses, real[] doses_sq,
                     int n, int[] yE, int[] yT, int[] levels,
                     real muE, real betaE1, real betaE2, real muT, real betaT1,
                     real psi) {
    real p;
    p = 0;
    for(j in 1:n) {
      real prob_eff;
      real prob_tox;
      real p_j;
      prob_eff = inv_logit(muE + betaE1 * doses[levels[j]] +
                           betaE2 * doses_sq[levels[j]]);
      prob_tox = inv_logit(muT + betaT1 * doses[levels[j]]);
      p_j = prob_eff^yE[j] * (1 - prob_eff)^(1 - yE[j]) * prob_tox^yT[j] *
              (1 - prob_tox)^(1 - yT[j]) + (-1)^(yE[j] + yT[j]) * prob_eff *
              prob_tox * (1 - prob_eff) * (1 - prob_tox) *
              (exp(psi) - 1) / (exp(psi) + 1);
      p = p + log(p_j);
    }
    return p;
  }
}

data {
  // Hyperparameters
  real muE_mean;
  real<lower=0> muE_sd;
  real betaE1_mean;
  real<lower=0> betaE1_sd;
  real betaE2_mean;
  real<lower=0> betaE2_sd;
  real muT_mean;
  real<lower=0> muT_sd;
  real betaT1_mean;
  real<lower=0> betaT1_sd;
  real psi_mean;
  real<lower=0> psi_sd;
  // Fixed trial parameters
  int<lower=1> K;
  real doses[K];
  // L^p norm for the efficacy-toxicity indifference contours
  real p;
  // Minimum required Pr(Efficacy) when Pr(Toxicity) = 0
  real pi1E;
  // Maximum permissable Pr(Toxicity) when Pr(Efficacy) = 1
  real pi2T;
  // A dose is acceptable if prob(eff) exceeds this hurdle
  real piE;
  // A dose is acceptable if prob(eff) is less than this hurdle
  real piT;
  // Observed trial outcomes
  int<lower=0> n;
  int<lower=0, upper=1> yE[n];
  int<lower=0, upper=1> yT[n];
  int<lower=1, upper=K> levels[n];
}

transformed data {
  real doses_sq[K];
  doses_sq = square(doses);
}

parameters {
  // Coefficients in efficacy logit model
  real muE;
  real betaE1;
  real betaE2;
  // Coefficients in toxicity logit model
  real muT;
  real<lower=0> betaT1;
  // Association paramater
  real psi;
}

transformed parameters {
  // Posterior probability of efficacy at doses i = 1, ..., K
  real<lower=0, upper=1> prob_eff[K];
  // Posterior probability of toxicity at doses i = 1, ..., K
  real<lower=0, upper=1> prob_tox[K];
  // Probability efficacy is acceptable at doses i = 1, ..., K
  real<lower=0, upper=1> prob_acc_eff[K];
  // Probability toxicity is acceptable at doses i = 1, ..., K
  real<lower=0, upper=1> prob_acc_tox[K];
  real utility[K];
  for(i in 1:K) {
    real r_to_the_p;
    prob_eff[i] = inv_logit(muE + betaE1 * doses[i] + betaE2 * doses_sq[i]);
    prob_tox[i] = inv_logit(muT + betaT1 * doses[i]);
    prob_acc_eff[i] = int_step(prob_eff[i] - piE);
    prob_acc_tox[i] = int_step(piT - prob_tox[i]);
    r_to_the_p = ((1 - prob_eff[i]) / (1 - pi1E))^p + (prob_tox[i] / pi2T)^p;
    utility[i] = 1 - r_to_the_p^(1/p);
  }
}

model {
  target += normal_lpdf(muE | muE_mean, muE_sd);
  target += normal_lpdf(betaE1 | betaE1_mean, betaE1_sd);
  target += normal_lpdf(betaE2 | betaE2_mean, betaE2_sd);
  target += normal_lpdf(muT | muT_mean, muT_sd);
  target += normal_lpdf(betaT1 | betaT1_mean, betaT1_sd);
  target += normal_lpdf(psi | psi_mean, psi_sd);
  target += log_joint_pdf(doses, doses_sq, n, yE, yT, levels,
                          muE, betaE1, betaE2, muT, betaT1, psi);
}
