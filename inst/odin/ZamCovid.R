## Compartment indexes: age groups (i), vaccine (j)
n_groups <- user()
# n_vacc_classes <- user()


## Definition of time-step and output as "time"
steps_per_day <- user(integer = TRUE)
dt <- 1 / steps_per_day
initial(time) <- 0
update(time) <- (step + 1) * dt


## State variables
update(S[]) <- S[i] + delta_S[i]
update(E1[]) <- E1[i] + delta_E1[i]
update(E2[]) <- E2[i] + delta_E2[i]
update(I[]) <- I[i] + delta_I[i]
update(R[]) <- R[i] + delta_R[i]
update(D[]) <- D[i] + delta_D[i]


## Progression probabilities
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to I - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_E * dt) # progression through latent period
p_IR <- 1 - exp(-gamma_I * dt) # progression through infectious period
p_ID[] <- 1 - exp(-gamma_D[i] * dt) # I to D - age dependent
p_RS <- 1 - exp(-gamma_R * dt) # waning immunity


## Binomial draws for numbers who will be changing compartments
n_SE[] <- rbinom(S[i], p_SE[i])
n_EE[] <- rbinom(E1[i], p_EE)
n_EI[] <- rbinom(E2[i], p_EI)
n_IR[] <- rbinom(I[i], p_IR)
n_ID[] <- rbinom(D[i], p_ID[i])
n_RS[] <- rbinom(R[i], p_RS)


## Vectors handling movement of individuals across compartments
delta_S[] <- n_RS[i] - n_SE[i]
delta_E1[] <- n_SE[i] - n_EE[i]
delta_E2[] <- n_EE[i] - n_EI[i]
delta_I[] <- n_EI[i] - n_IR[i] - n_ID[i]
delta_R[] <- n_IR[i] - n_RS[i]
delta_D[] <- n_ID[i]


## Seeding model for first epidemic wave
seed_step_end <- seed_step_start + length(seed_value)
seed_rate <- if (step >= seed_step_start && step < seed_step_end)
  seed_value[as.integer(step - seed_step_start + 1)] else 0
seed <- rpois(seed_rate)

seed_age_band <- user()
seed_step_start <- user()
seed_value[] <- user()
dim(seed_value) <- user()

n_SE[seed_age_band] <- n_SE[seed_age_band] + min(S[seed_age_band], seed)


## Force of infection
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


## Observation equations and states
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
