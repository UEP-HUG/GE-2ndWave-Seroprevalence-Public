//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic using data from lab validation study.
// used for analyses of Specchio and Corona Immunitas in Geneva, Switzerland

// We include household as random effects

// Hierarchical model for specificity and sensitivity following Gelman & Carpenter 2020
// For code see https://github.com/bob-carpenter/diagnostic-testing/blob/master/src/specificity-santa-clara/stan/santa-clara-hierarchical.stan
// and https://github.com/bob-carpenter/diagnostic-testing/blob/master/src/specificity-santa-clara/R/rstan-santa-clara.R
// For initial discussion see https://statmodeling.stat.columbia.edu/2020/05/01/simple-bayesian-analysis-inference-of-coronavirus-infection-rate-from-the-stanford-study-in-santa-clara-county/
// and for final paper see http://www.stat.columbia.edu/~gelman/research/published/specificity.pdf

// We have input data from a validation set (which determines the diagnostic performance) and the survey (from which we'd like to estimate seropos).
data {
    int<lower=1> N_survey; //numbber of participants in the survey
    int<lower=1> H; //number of households in the survey
    int<lower=0> survey_pos[N_survey]; //observation positive or negative
    int<lower=1, upper=H> hh[N_survey]; //household of observation
    int<lower=1> p_vars; //number of variables to adjust for
    matrix[N_survey, p_vars] X; //covariate model matrix (age, sex, week in these analyses)

    int<lower=0> J_studies; // number of studies with validation data
    int<lower=0> N_pos_control[J_studies]; //number of positive controls in the validation data
    int<lower=0> control_tp[J_studies]; // number of true positive tests in the validation data
    int<lower=0> N_neg_control[J_studies]; // number of negative controls in the validation data
    int<lower=0> control_fp[J_studies];// number of false positives by the diagnostic test in the validation study
}

parameters {
    vector[p_vars] beta; // fixed regression coefficients
    real<lower=0> sigma_h; // variability of household random effect
    vector[H] eta_h; // standard normals for the household random effect

    real mu_logit_spec; // mean for hierarchical specificity
    real mu_logit_sens; // mean for hierarchical sensitivity
    real<lower = 0> sigma_logit_spec; // variance for hierarchical specificity
    real<lower = 0> sigma_logit_sens; // variance for hierarchical sensitivity
    vector<offset = mu_logit_spec, multiplier = sigma_logit_spec>[J_studies] logit_spec; // normal dist on log odds scale for hierarchical specificity
    vector<offset = mu_logit_sens, multiplier = sigma_logit_sens>[J_studies] logit_sens; // normal dist on log odds scale for hierarchical sensitivity
}

transformed parameters {
  vector<lower=0, upper=1>[N_survey] p; // probability of seropositivity for an observation

  vector[J_studies] spec = inv_logit(logit_spec); // unwind the logistic/log odds transform (I think so its normal value is used in the model??)
  vector[J_studies] sens = inv_logit(logit_sens);

  p = inv_logit(X * beta + sigma_h*eta_h[hh]);
}

//  We observe 'survey_pos' cases as a bernoulli distribution based on the
//  survey with observations coming as the sum of the true
//  positive rate p*sens and false negative rate (1-p)*(1-spec).
model {
    target+= bernoulli_lpmf(survey_pos | p*sens[1]+(1-p)*(1-spec[1]));
    target+= binomial_lpmf(control_tp | N_pos_control, sens);
    target+= binomial_lpmf(control_fp | N_neg_control, 1-spec);
    target+= normal_lpdf(eta_h | 0, 1);
    target+= normal_lpdf(beta | 0, 1); // priors for coefficients

    target+= normal_lpdf(logit_spec | mu_logit_spec, sigma_logit_spec); // hierarchical prior for spec
    target+= normal_lpdf(logit_sens | mu_logit_sens, sigma_logit_sens); // hierarchical prior for sens
    target+= normal_lpdf(sigma_logit_spec | 0, 1);// weak truncated half-normal prior for spec variance
    target+= normal_lpdf(sigma_logit_sens | 0, 1);// weak truncated half-normal prior for sens variance
    target+= normal_lpdf(mu_logit_spec | 4, 2);  // weak prior on mean of distribution of spec
    target+= normal_lpdf(mu_logit_sens | 4, 2);  // weak prior on mean of distribution of sens
}

generated quantities {
  vector[N_survey] log_lik;

    for(i in 1:N_survey){
      log_lik[i] = bernoulli_logit_lpmf(survey_pos[i] | p[i]*sens+(1-p[i])*(1-spec));
    }

}
