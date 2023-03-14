## Compartment indexes: age groups (i), vaccine class (j)
# Note k index is used for shape parameter of erlang distributed rates
n_groups <- user()
n_vacc_classes <- user()


## Definition of time-step and output as "time"
steps_per_day <- user(integer = TRUE)
dt <- 1 / steps_per_day
initial(time) <- 0
update(time) <- (step + 1) * dt


## Core state variables (model compartments)
update(S[, ]) <- new_S[i, j]
update(E[, , ]) <- new_E[i, j, k]
update(I_A[, ]) <- new_I_A[i, j]
update(I_P[, ]) <- new_I_P[i, j]
update(I_C_1[, ]) <- new_I_C_1[i, j]
update(I_C_2[, ]) <- new_I_C_2[i, j]
update(H_R_conf[, , ]) <- new_H_R_conf[i, j, k]
update(H_R_unconf[, , ]) <- new_H_R_unconf[i, j, k]
update(H_D_conf[, , ]) <- new_H_D_conf[i, j, k]
update(H_D_unconf[, , ]) <- new_H_D_unconf[i, j, k]
update(R[, ]) <- new_R[i, j]
update(D_hosp[, ]) <- D_hosp[i, j] + delta_D_hosp[i, j]
update(D_non_hosp[, ]) <- D_non_hosp[i, j] + n_I_C_2_to_G_D[i, j]


## Parallel flows state variables
update(T_sero_pre[, , ]) <- new_T_sero_pre[i, j, k]
update(T_sero_pos[, , ]) <- new_T_sero_pos[i, j, k]
update(T_sero_neg[, ]) <- new_T_sero_neg[i, j]

update(T_PCR_pre[, , ]) <- new_T_PCR_pre[i, j, k]
update(T_PCR_pos[, , ]) <- new_T_PCR_pos[i, j, k]
update(T_PCR_neg[, , ]) <- new_T_PCR_neg[i, j, k]


## Work out new individuals entering compartments
new_S[, ] <- S[i, j] + n_R_progress[i, j] + n_infected_to_S[i,j] -
  n_S_progress[i, j] - n_S_next_vacc_class[i, j] +
  (if (j == 1) n_S_next_vacc_class[i, n_vacc_classes] else
    n_S_next_vacc_class[i, j - 1])
new_E[, , ] <- E[i, j, k] +
  (if (k == 1) n_S_progress[i, j] else n_E_progress[i, j, k - 1]) -
  n_E_progress[i, j, k] -
  n_E_next_vacc_class[i, j, k] +
  (if (j == 1) n_E_next_vacc_class[i, n_vacc_classes, k] else
    n_E_next_vacc_class[i, j - 1, k])
new_I_A[, ] <- I_A[i, j] + n_EI_A[i, j] -
  n_I_A_progress[i, j] - n_I_A_next_vacc_class[i, j] +
  (if (j == 1) n_I_A_next_vacc_class[i, n_vacc_classes] else
    n_I_A_next_vacc_class[i, j -1])
new_I_P[, ] <- I_P[i, j] + n_EI_P[i, j] -
  n_I_P_progress[i, j] - n_I_P_next_vacc_class[i, j] +
  if (j == 1) n_I_P_next_vacc_class[i, n_vacc_classes] else
    n_I_P_next_vacc_class[i, j - 1]
new_I_C_1[, ] <- I_C_1[i, j] + n_I_P_progress[i, j] -
  n_I_C_1_progress[i, j]
new_I_C_2[, ] <- I_C_2[i, j] + n_I_C_1_progress[i, j] -
  n_I_C_2_progress[i, j]

aux_H_R_conf[, , ] <- H_R_conf[i, j, k] +
  (if (k > 1) n_H_R_conf_progress[i, j, k - 1] else 0) -
  n_H_R_conf_progress[i, j, k]
aux_H_R_unconf[, , ] <- H_R_unconf[i, j, k] +
  (if (k > 1) n_H_R_unconf_progress[i, j, k - 1] else 0) -
  n_H_R_unconf_progress[i, j, k]
new_H_R_conf[, , ] <-
  aux_H_R_conf[i, j, k] + n_H_R_unconf_to_conf[i, j, k] +
  (if (k == 1) n_I_C_2_to_H_R_conf[i, j] else 0)
new_H_R_unconf[, , ] <-
  aux_H_R_unconf[i, j, k] - n_H_R_unconf_to_conf[i, j, k] +
  (if (k == 1) n_I_C_2_to_H_R[i, j] - n_I_C_2_to_H_R_conf[i, j] else 0)

aux_H_D_conf[, , ] <- H_D_conf[i, j, k] +
  (if (k > 1) n_H_D_conf_progress[i, j, k - 1] else 0) -
  n_H_D_conf_progress[i, j, k]
aux_H_D_unconf[, , ] <- H_D_unconf[i, j, k] +
  (if (k > 1) n_H_D_unconf_progress[i, j, k - 1] else 0) -
  n_H_D_unconf_progress[i, j, k]
new_H_D_conf[, , ] <-
  aux_H_D_conf[i, j, k] + n_H_D_unconf_to_conf[i, j, k] +
  (if (k == 1) n_I_C_2_to_H_D_conf[i, j] else 0)
new_H_D_unconf[, , ] <-
  aux_H_D_unconf[i, j, k] - n_H_D_unconf_to_conf[i, j, k] +
  (if (k == 1) n_I_C_2_to_H_D[i, j] - n_I_C_2_to_H_D_conf[i, j] else 0)

new_R[, ] <- R[i, j] -
  n_R_progress[i, j] - n_R_next_vacc_class[i, j] +
  n_infected_to_R[i, j] +
  (if (j == 1) n_R_next_vacc_class[i, n_vacc_classes] else
    n_R_next_vacc_class[i, j - 1])
delta_D_hosp[, ] <- sum(n_H_D_unconf_progress[i, j, ]) +
  sum(n_H_D_conf_progress[i, j, ])


## Binomial draws for numbers who will progress
n_S_progress[, ] <- rbinom(S[i, ], p_SE[i, ])
n_S_next_vacc_class[, ] <- rbinom(S[i, j] - n_S_progress[i, j],
                                  p_S_next_vacc_class[i, j])
n_E_progress[, , ] <- rbinom(E[i, j, k], p_E_progress)
n_E_next_vacc_class[, , ] <- rbinom(E[i, j, k] - n_E_progress[i, j, k],
                                    p_E_next_vacc_class[i, j])
n_EI_A[, ] <- rbinom(n_E_progress[i, j, k_E], 1 - p_C[i, j])
n_EI_P[, ] <- n_E_progress[i, j, k_E] - n_EI_A[i, j]

n_I_A_progress[, ] <- rbinom(I_A[i, j], p_I_A_progress)
n_I_A_next_vacc_class[, ] <- rbinom(I_A[i, j] - n_I_A_progress[i, j],
                                    p_I_A_next_vacc_class[i, j])
n_I_P_progress[, ] <- rbinom(I_P[i, j], p_I_P_progress[j])
n_I_P_next_vacc_class[, ] <- rbinom(I_P[i, j] - n_I_P_progress[i, j],
                                    p_I_P_next_vacc_class[i, j])

n_I_C_1_progress[, ] <- rbinom(I_C_1[i, j], p_I_C_1_progress[j])
n_I_C_2_progress[, ] <- rbinom(I_C_2[i, j], p_I_C_2_progress[j])

n_I_C_2_to_RS[, ] <- rbinom(n_I_C_2_progress[i, j], 1 - p_H[i, j])
n_I_C_2_to_G_D[, ] <- rbinom(n_I_C_2_progress[i, j] - n_I_C_2_to_RS[i, j],
                           p_G_D[i, j])
n_I_C_2_to_hosp[, ] <- n_I_C_2_progress[i, j] - n_I_C_2_to_RS[i, j] -
  n_I_C_2_to_G_D[i, j]

n_I_C_2_to_H_R[, ] <- rbinom(n_I_C_2_to_hosp[i, j], 1 - p_H_D[i, j])
n_I_C_2_to_H_D[, ] <- n_I_C_2_to_hosp[i, j] - n_I_C_2_to_H_R[i, j]

n_I_C_2_to_H_R_conf[, ] <- rbinom(n_I_C_2_to_H_R[i, j], p_star[i])
n_I_C_2_to_H_D_conf[, ] <- rbinom(n_I_C_2_to_H_D[i, j], p_star[i])
n_I_C_2_to_H_R_unconf[, ] <- n_I_C_2_to_H_R[i, j] - n_I_C_2_to_H_R_conf[i, j]
n_I_C_2_to_H_D_unconf[, ] <- n_I_C_2_to_H_D[i, j] - n_I_C_2_to_H_D_conf[i, j]

n_H_R_unconf_to_conf[, , ] <- rbinom(aux_H_R_unconf[i, j, k], p_test)
n_H_D_unconf_to_conf[, , ] <- rbinom(aux_H_D_unconf[i, j, k], p_test)
n_H_R_conf_progress[, , ] <- rbinom(H_R_conf[i, j, k], p_H_R_progress[j])
n_H_R_unconf_progress[, , ] <- rbinom(H_R_unconf[i, j, k], p_H_R_progress[j])
n_H_D_conf_progress[, , ] <- rbinom(H_D_conf[i, j, k], p_H_D_progress[j])
n_H_D_unconf_progress[, , ] <- rbinom(H_D_unconf[i, j, k], p_H_D_progress[j])

n_infection_end[, ] <- n_I_A_progress[i, j] + n_I_C_2_to_RS[i, j] +
  n_H_R_conf_progress[i, j, k_H_R] + n_H_R_unconf_progress[i, j, k_H_R]
n_infected_to_R[, ] <- rbinom(n_infection_end[i, j], p_R[i, j])
n_infected_to_S[, ] <- n_infection_end[i, j] - n_infected_to_R[i, j]

n_R_progress[, ] <- rbinom(R[i, j], p_R_progress)
n_R_next_vacc_class[, ] <- rbinom(R[i, j] - n_R_progress[i, j],
                                  p_R_next_vacc_class[i, j])


## Parallel serology flows
# Note flow is split between those who will and won't seroconvert
n_com_to_T_sero_pre[, ] <- n_E_progress[i, j, k_E]

new_T_sero_pre[, , ] <- T_sero_pre[i, j, k] -
  n_T_sero_pre_progress[i, j, k] +
  (if (k == 1) n_com_to_T_sero_pre[i, j] else
    n_T_sero_pre_progress[i, j, k - 1])

n_T_sero_pre_to_T_sero_pos[, ] <-
  rbinom(n_T_sero_pre_progress[i, j, k_sero_pre], p_sero_pos[i])

new_T_sero_pos[, , ] <- T_sero_pos[i, j, k] -
  n_T_sero_pos_progress[i, j, k] +
  (if (k == 1)  n_T_sero_pre_to_T_sero_pos[i, j] else
    n_T_sero_pos_progress[i, j, k - 1])

new_T_sero_neg[, ] <- T_sero_neg[i, j] +
  n_T_sero_pre_progress[i, j, k_sero_pre] -
  n_T_sero_pre_to_T_sero_pos[i, j] +
  n_T_sero_pos_progress[i, j, k_sero_pos]

n_T_sero_pre_progress[, , ] <- rbinom(T_sero_pre[i, j, k],
                                      p_T_sero_pre_progress)
n_T_sero_pos_progress[, , ] <- rbinom(T_sero_pos[i, j, k],
                                      p_T_sero_pos_progress)


## Clinical progression probabilities
p_SE[, ] <- 1 - exp(-lambda_susc[i, j] * dt) # age/vacc dependent
p_E_progress <- 1 - exp(-gamma_E * dt)
p_I_A_progress <- 1 - exp(-gamma_A * dt)
p_I_P_progress <- 1 - exp(-gamma_P * dt)
p_I_C_1_progress <- 1 - exp(-gamma_C_1 * dt)
p_I_C_2_progress <- 1 - exp(-gamma_C_2 * dt)
p_H_R_progress <- 1 - exp(-gamma_H_R * dt)
p_H_D_progress <- 1 - exp(-gamma_H_D * dt)
p_R_progress <- 1 - exp(-waning_rate * dt)
p_test <- 1 - exp(-gamma_U * dt)

p_sero_pos[] <- user()
p_T_sero_pre_progress <- 1 - exp(-gamma_sero_pre * dt)
p_T_sero_pos_progress <- 1 - exp(-gamma_sero_pos * dt)


## Vaccination model
p_S_next_vacc_class[, ] <- vaccine_probability[i, j]
p_E_next_vacc_class[, ] <- vaccine_probability[i, j]
p_I_A_next_vacc_class[, ] <- vaccine_probability[i, j]
p_I_P_next_vacc_class[, ] <- vaccine_probability[i, j]
p_R_next_vacc_class[, ] <- vaccine_probability[i, j]


## Time-varying probabilities
p_C[, ] <- if (as.integer(step) >= n_p_C_steps)
  min(p_C_step[n_p_C_steps, i] * rel_p_sympt[i, j], as.numeric(1)) else
    min(p_C_step[step + 1, i] * rel_p_sympt[i, j], as.numeric(1))

p_H[, ] <- if (as.integer(step) >= n_p_H_steps)
  min(p_H_step[n_p_H_steps, i] * rel_p_hosp_if_sympt[i, j], as.numeric(1)) else
    min(p_H_step[step + 1, i] * rel_p_hosp_if_sympt[i, j], as.numeric(1))

p_H_D[, ] <- if (as.integer(step) >= n_p_H_steps)
  min(p_H_D_step[n_p_H_D_steps, i] * rel_p_H_D[i, j], as.numeric(1)) else
    min(p_H_D_step[step + 1, i] * rel_p_H_D[i, j], as.numeric(1))

p_star[] <- if (as.integer(step) >= n_p_star_steps)
  p_star_step[n_p_star_steps, i] else p_star_step[step + 1, i]



####### TODO: need to define serology flow progression probabilities
# p_T_PCR_pre_progress <- 1 - exp(-gamma_PCR_pre * dt)
# p_T_PCR_pos_progress <- 1 - exp(-gamma_PCR_pos * dt)
##########


## Seeding model for first epidemic wave
seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
  seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)

seed_age_band <- user()
seed_step_start <- user()
seed_value[] <- user()
dim(seed_value) <- user()

n_S_progress[seed_age_band, 1] <- n_S_progress[seed_age_band, 1] +
  min(S[seed_age_band, 1], seed)


## TODO: Force of infection
m[, ] <- user()
temp[] <- I[i]
s_ij[, ] <- m[i, j] * temp[j] # note temp[j] not temp[i]
lambda[] <- beta * sum(s_ij[i, ])

beta_step[] <- user()
dim(beta_step) <- user()
beta <- if (step >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]

# Output betas - useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta




## Space observation equations
update(new_I[]) <- if (step %% steps_per_day == 0) n_EI[i] else
  new_I[i] + n_EI[i]

delta_case_under_15 <- sum(new_I[1:3]) * rho
delta_case_15_19 <- new_I[4] * rho
delta_case_20_29 <- sum(new_I[5:6]) * rho
delta_case_30_39 <- sum(new_I[7:8]) * rho
delta_case_40_49 <- sum(new_I[9:10]) * rho
delta_case_50_plus <- sum(new_I[11:n_groups]) * rho

update(cases_under_15) <- if (step %% steps_per_day == 0)
  delta_case_under_15 else delta_case_under_15 + cases_under_15
update(cases_15_19) <- if (step %% steps_per_day == 0)
  delta_case_15_19 else delta_case_15_19 + cases_15_19
update(cases_20_29) <- if (step %% steps_per_day == 0)
  delta_case_20_29 else delta_case_20_29 + cases_20_29
update(cases_30_39) <- if (step %% steps_per_day == 0)
  delta_case_30_39 else delta_case_30_39 + cases_30_39
update(cases_40_49) <- if (step %% steps_per_day == 0)
  delta_case_40_49 else delta_case_40_49 + cases_40_49
update(cases_50_plus) <- if (step %% steps_per_day == 0)
  delta_case_50_plus else delta_case_50_plus + cases_50_plus


### Initial states are all zeroed ###
# S and E will be initialiased based on the seeding model
initial(S[]) <- 0
initial(E1[]) <- 0
initial(E2[]) <- 0
initial(I[]) <- 0
initial(R[]) <- 0
initial(D[]) <- 0

initial(new_I[]) <- 0
initial(cases_under_15) <- 0
initial(cases_15_19) <- 0
initial(cases_20_29) <- 0
initial(cases_30_39) <- 0
initial(cases_40_49) <- 0
initial(cases_50_plus) <- 0


### User defined parameters ###
gamma_E <- user()
gamma_I <- user()
gamma_R <- user()
rho <- user()
gamma_D[] <- user() # age-specific IFR


## Dimensions of arrays and vectors ###
dim(m) <- c(n_groups, n_groups)

dim(S) <- n_groups
dim(p_SE) <- n_groups
dim(n_SE) <- n_groups
dim(delta_S) <- n_groups

dim(E1) <- c(n_groups)
dim(delta_E1) <- c(n_groups)
dim(n_EE) <- c(n_groups)

dim(E2) <- c(n_groups)
dim(delta_E2) <- c(n_groups)
dim(n_EI) <- c(n_groups)

dim(I) <- c(n_groups)
dim(delta_I) <- c(n_groups)
dim(n_IR) <- c(n_groups)

dim(R) <- c(n_groups)
dim(delta_R) <- c(n_groups)
dim(n_RS) <- c(n_groups)

dim(D) <- c(n_groups)
dim(p_ID) <- c(n_groups)
dim(delta_D) <- c(n_groups)
dim(n_ID) <- c(n_groups)

dim(lambda) <- n_groups
dim(s_ij) <- c(n_groups, n_groups)
dim(temp) <- n_groups

dim(gamma_D) <- n_groups

dim(new_I) <- n_groups
