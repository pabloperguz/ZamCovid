## Definition of time-step and output as "time"
steps_per_day <- user(integer = TRUE)
dt <- 1 / steps_per_day
initial(time) <- 0
update(time) <- (step + 1) * dt


### State model equations and transition probabilities

# Compartments
update(S[]) <- S[i] + delta_S[i]
update(E1[]) <- E1[i] + delta_E1[i]
update(E2[]) <- E2[i] + delta_E2[i]
update(I[]) <- I[i] + delta_I[i]
update(R[]) <- R[i] + delta_R[i]

# Transition probabilities
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to I - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_E * dt) # progression through latent period
p_IR <- 1 - exp(-gamma_I * dt) # progression through infectious period
p_RS <- 1 - exp(-eta_R * dt) # waning immunity

# Binomial draws for numbers who will be changing compartments
n_SE[] <- rbinom(S[i], p_SE[i])
n_EE[] <- rbinom(E1[i], p_EE)
n_EI[] <- rbinom(E2[i], p_EI)
n_IR[] <- rbinom(I[i], p_IR)
n_RS[] <- rbinom(R[i], p_RS)

# Movement of individuals across compartments
delta_S[] <- n_RS[i] - n_SE[i]
delta_E1[] <- n_SE[i] - n_EE[i]
delta_E2[] <- n_EE[i] - n_EI[i]
delta_I[] <- n_EI[i] - n_IR[i]
delta_R[] <- n_IR[i] - n_RS[i]


### Observation equations and states

update(new_I[]) <- if (step %% steps_per_day == 0) n_EI[i] else
  new_I[i] + n_EI[i]

delta_case_under_15 <- sum(new_I[1:3]) * rho
delta_case_15_19 <- new_I[4] * rho
delta_case_20_29 <- sum(new_I[5:6]) * rho
delta_case_30_39 <- sum(new_I[7:8]) * rho
delta_case_40_49 <- sum(new_I[9:10]) * rho
delta_case_50_plus <- sum(new_I[11:N_age]) * rho

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



### Force of infection ###

# Lambda as function of contact matrix and beta
m[, ] <- user()
temp[] <- I[i]
s_ij[, ] <- m[i, j] * temp[j] # note temp[j] not temp[i]
lambda[] <- beta * sum(s_ij[i, ])

# Time-varying beta
beta_step[] <- user()
dim(beta_step) <- user()
beta <- if (step >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]

# Output betas - useful for debugging
initial(beta_out) <- beta_step[1]
update(beta_out) <- beta


### Initial states ###
initial(S[]) <- S_ini[i]
initial(E1[]) <- E_ini[i]
initial(E2[]) <- E2_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- R_ini[i]

initial(new_I[]) <- 0
initial(cases_under_15) <- 0
initial(cases_15_19) <- 0
initial(cases_20_29) <- 0
initial(cases_30_39) <- 0
initial(cases_40_49) <- 0
initial(cases_50_plus) <- 0

# Initial vectors
S_ini[] <- user()
E_ini[] <- user()
E2_ini[] <- user()
I_ini[] <- user()
R_ini[] <- user()

### User defined parameters ###
gamma_E <- user()
gamma_I <- user()
eta_R <- user()
rho <- user()

## Dimensions of arrays and vectors ###

# Number of age (i), vaccine (j) and strain (k) classes
N_age <- user()

dim(m) <- c(N_age, N_age)

dim(S) <- N_age
dim(S_ini) <- N_age
dim(p_SE) <- N_age
dim(n_SE) <- N_age
dim(delta_S) <- N_age

dim(E1) <- c(N_age)
dim(E_ini) <- c(N_age)
dim(E2_ini) <- c(N_age)
dim(delta_E1) <- c(N_age)
dim(n_EE) <- c(N_age)

dim(E2) <- c(N_age)
dim(delta_E2) <- c(N_age)
dim(n_EI) <- c(N_age)

dim(I_ini) <- c(N_age)
dim(I) <- c(N_age)
dim(delta_I) <- c(N_age)
dim(n_IR) <- c(N_age)

dim(R) <- c(N_age)
dim(R_ini) <- c(N_age)
dim(delta_R) <- c(N_age)
dim(n_RS) <- c(N_age)

dim(lambda) <- N_age
dim(s_ij) <- c(N_age,N_age)
dim(temp) <- N_age

dim(new_I) <- N_age
