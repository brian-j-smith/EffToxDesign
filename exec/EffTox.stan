// Dose-Finding Based on Efficacy-Toxicity Trade-Offs, by Thall, Cook
// Effective sample size for computing prior hyperparameters in Bayesian phase I-II dose-finding,
//  by Thall, Herrick, Nguyen, Venier, Norris

functions {
  real log_joint_pdf(real[] coded_doses, real[] coded_doses_squ,
                     int n, int[] yE, int[] yT, int[] levels,
                     real muE, real betaE1, real betaE2, real muT, real betaT1,
                     real psi) {
    real p;
    p = 0;
    for(j in 1:n) {
      real prob_eff;
      real prob_tox;
      real p_j;
      prob_eff = inv_logit(muE + betaE1 * coded_doses[levels[j]] +
                           betaE2 * coded_doses_squ[levels[j]]);
      prob_tox = inv_logit(muT + betaT1 * coded_doses[levels[j]]);
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
  real<lower=0> doses[K]; // Doses under investigation, e.g. 10, 20, 30 for 10mg, 20mg, 30mg
  real p;  // The p of the L^p norm used to model the efficacy-toxicity indifference contours.
           // See Efficacy-Toxicity trade-offs based on L-p norms: Technical Report UTMDABTR-003-06
  real pi1E; // Minimum required Pr(Efficacy) when Pr(Toxicity) = 0
  real pi2T; // Maximum permissable Pr(Toxicity) when Pr(Efficacy) = 1
  real piE; // A dose is acceptable if prob(eff) exceeds this hurdle...
  real piT; //  ... and prob(tox) is less than this hurdle
  // Observed trial outcomes
  int<lower=0> n;
  int<lower=0, upper=1> yE[n]; // Binary efficacy event for patients j=1,..,n
  int<lower=0, upper=1> yT[n]; // Binary toxicity event for patients j=1,..,n
  int<lower=1, upper=K> levels[n];  // Dose-levels given for patients j=1,..,n.
                                    // Dose-levels are 1-based indices of real_doses
                                    // E.g. 1 means 1st dose in real_doses was given
}

transformed data {
  // Thall & Cook transform the actual doses by logging and centralising:
  real coded_doses[K];
  real coded_doses_squ[K]; // The square of coded_doses
  real mean_log_doses; // Variable created for convenience
  mean_log_doses = 0.0;
  for(i in 1:K) mean_log_doses += log(doses[i]);
  mean_log_doses = mean_log_doses / K;
  for(i in 1:K) {
    coded_doses[i] = log(doses[i]) - mean_log_doses;
    coded_doses_squ[i] = coded_doses[i]^2;
  }
}

parameters {
  // Coefficients in efficacy logit model:
  real muE;
  real betaE1;
  real betaE2;
  // Coefficients in toxicity logit model:
  real muT;
  real<lower=0> betaT1;
  // Association:
  real psi;
}

transformed parameters {
  real<lower=0, upper=1> prob_eff[K]; // Posterior probability of efficacy at doses i=1,...,K
  real<lower=0, upper=1> prob_tox[K]; // Posterior probability of toxicity at doses i=1,...,K
  real<lower=0, upper=1> prob_acc_eff[K]; // Probability efficacy is acceptable at doses i=1,...,K
  real<lower=0, upper=1> prob_acc_tox[K]; // Probability toxicity is acceptable at doses i=1,...,K
  real utility[K]; // Posterior utility of doses i=1,...,K
  // Calculate the utility of each dose using the method described in
  // "Efficacy-Toxicity trade-offs based on L-p norms: Technical Report UTMDABTR-003-06", John Cook
  for(i in 1:K) {
    real r_to_the_p; // Convenience variable, as in (2) of Cook.
    prob_eff[i] = inv_logit(muE + betaE1 * coded_doses[i] +
                            betaE2 * coded_doses_squ[i]);
    prob_tox[i] = inv_logit(muT + betaT1 * coded_doses[i]);
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
  target += log_joint_pdf(coded_doses, coded_doses_squ, n, yE, yT, levels,
                          muE, betaE1, betaE2, muT, betaT1, psi);
}
