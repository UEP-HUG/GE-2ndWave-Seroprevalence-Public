# Preamble ----------------------------------------------------------------
library(tidyverse)
library(rstan)
library(uuid) # for generating a new Universally Unique Identifier
library(parallel)
library(foreach)

source("scripts_specchio/R/functions_Stan_hierarchicalSpecSens.R")

# use a session ID for similar filenames from same code run
session_ID = substr(uuid::UUIDgenerate(), 1, 8)

## Stan control settings that can be changed for testing
# (for production use 4 chains and at least 1500 iter, 250 warmup)
n_chains = 4
n_iter = 1500
n_warmup = 250
random_seed = 1543

with_households = TRUE
input_dat = read_csv("data/example_dat4stan.csv")

if (with_households) {
  Stan_script = "scripts_specchio/Stan/specchio-seroprev-hh-sens-spec-hierarchical.stan"
} else {
  Stan_script = "scripts_specchio/Stan/specchio-seroprev-basic-sens-spec-hierarchical.stan"
  input_dat = input_dat %>% filter(!family)
}

## Stan settings, don't need to change
options(mc.cores = 4)
p_delta = 0.99
n_treedepth = 20
rstan_options(auto_write = TRUE)

# Define model --------------------------------------------------------------
model = "Sex + age_cat"

## Define desired age cuts
age_cuts = c(0, 6, 12, 18, 25, 35, 50, 65, 75, 105)

age_ref = "[25,35)"

sex_ref = "Female"

ID_note = "_HH_silviaAges"

session_ID = paste0(session_ID, ID_note)

# specchio data -------------------------------------------
input_dat = input_dat %>%
  mutate(
    age_cat = factor(cut(age, age_cuts, right = FALSE, include.lowest = TRUE)),
    Sex = fct_recode(sex, "Female" = "f", "Male" = "m"),
    UEP_result = fct_recode(
      UEP_result,
      "0" = "neg", "1" = "pos"
    ) %>%
      as.character() %>%
      as.integer()
  ) %>%
  group_by(hh_id) %>%
  mutate(
    hh_obs = n(),
    hh_inf = sum(UEP_result) - UEP_result,
    other_inf = hh_inf > 0
  ) %>%
  ungroup() %>%
  droplevels()

# some fct_relevelling
input_dat = input_dat %>%
  mutate(
    Sex = fct_relevel(Sex, ref = "Female"),
    age_cat = fct_relevel(age_cat, ref = age_ref)
  ) %>%
  group_by(hh_id, hh_obs) %>%
  mutate(tot_inf = sum(UEP_result)) %>%
  arrange(desc(hh_obs), tot_inf, hh_id)

# Control data ------------------------------------------------------------

## number of positive controls from validation data
pos_control = c(172, 1423)
## number of negative controls
neg_control = c(185, 5991)
## number of true positives for cases
control_tp = c(159, 1406)
## number of false positives for controls (1-specificity)
control_fp = c(1, 1)

# Geneva population ages data ---------------------------------------------------------------
population_data = getGE_age_cats(age_cuts)

population_data = population_data %>%
  mutate(
    Sex = fct_relevel(Sex, ref = "Female"),
    age_cat = fct_relevel(age_cat, ref = age_ref)
  )

population_data = population_data %>%
  filter(age_cat %in% levels(input_dat$age_cat),
         Sex %in% levels(input_dat$Sex)) %>%
  droplevels()

# Where to save output --------------------------------------------------------------
if (!dir.exists("output/results")) {
  dir.create("output/results")
}

output_result_list_filename = paste0(
  "output/results/",
  session_ID,
  "_results-list_",
  format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
  ".rds"
)

# Run models --------------------------------------------------------------
if (with_households) {
  ## model the overall seroprevalence with household effect
  fun = run_analysis_stan_with_hh
} else {
  ## model the overall seroprevalence (only index participants - no household effect)
  fun = run_analysis_stan_no_hh
}

calc_seropos = fun(
  model_script = Stan_script,
  dat = input_dat,
  coef_eqn = model,
  pos_control = pos_control,
  neg_control = neg_control,
  control_tp = control_tp,
  control_fp = control_fp,
  n_cores = getOption("mc.cores"),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  Stan_control = list(
    adapt_delta = p_delta,
    max_treedepth = n_treedepth
  ),
  seed = random_seed,
  pop_age_cats = population_data,
  session_ID = session_ID,
  sex_ref = sex_ref,
  age_ref = age_ref
)

# Calculate seroprevalence probabilities --------------------------------------

my_stanfit = calc_seropos$stan_posterior

if (with_households) {
  seroprev_probs = calc_integrated_seroprev_probs(my_stanfit,
    coef_eqn = model,
    pop_age_cats = population_data,
    n_cores = getOption("mc.cores"),
    session_ID = session_ID
  )
} else {
  seroprev_probs = calc_seroprev_probs(my_stanfit,
    coef_eqn = model,
    pop_age_cats = population_data,
    n_cores = getOption("mc.cores"),
    session_ID = session_ID
  )
}

pop_cats_probs = seroprev_probs$pop_cats_probs
subset_estimates = compute_weighted_estimates(pop_cats_probs)

all_results_list = c(
  model_matrix = list(calc_seropos$model_matrix),
  seroprev_probs,
  subset_estimates = list(subset_estimates)
)
saveRDS(all_results_list, output_result_list_filename)
cat("\n-------\n", output_result_list_filename, " saved at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n\n", sep = "")

# Save final output table ----------------------------------------------------
final_table = make_final_table(pop_cats_probs, subset_estimates, input_dat, age_ref, sex_ref)
res_table_file = paste0("output/results/", session_ID, "_results_table", ".csv")
final_table %>%
  write_csv(res_table_file)

cat("\n---- Done ----\n", res_table_file, " written at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n-------\n", sep = "")

# library(shinystan)
# launch_shinystan(my_stanfit)
