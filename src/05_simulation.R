# Install required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("data.table", "dplyr", "optimx", "Rcpp", "purrr", "svDialogs", "futile.logger")


# Load the necessary functions for the estimation of the simulated data
Rcpp::sourceCpp(file.path("src", "loglikelihoodnew.cpp"))
source(file.path("src", "load_functions.R"))


# Define the number of states - here we are examining the 4 state model
nstates <- 4


# Load the estimated parameters from the real dataset
load(file.path("parameter estimations", "vij.RData"))
load(file.path("parameter estimations", "sij.RData"))
load(file.path("parameter estimations", "aij.RData"))
load(file.path("parameter estimations", "bij.RData"))
params <- c(vij, sij, aij, bij)


# Scaling parameter to control the initial parameters for the estimation - always set to 1
parscale <- 1


# How many samples of simulated data do you want to produce? - It is advised to have a minimum number of 500 simulations, for 
# statistically significant results. However, if you want to quickly check the code, you can run it for 1 or 2 dataset simulations.
n_bootstrap <- as.numeric(
  dlg_input("Please enter the number of simulated data: ")$res)

    # Check for potential invalid bootstrapping number
    if (n_bootstrap == "" || is.na(n_bootstrap) || n_bootstrap < 0 || n_bootstrap == 0 ) {
      cat("No valid number entered.\n")
    } else {
      cat("You entered:", n_bootstrap, "\n")
    }


# --------------------------------- BOOTSTRAP ---------------------------------
# -----------------------------------------------------------------------------

# Create a list to save the estimation results
estimation_list <- list()


# Calculate the initial distribution of the 4 states in the real dataset for sampling the 
# first patient -  For confidentiality reasons, the initial state distribution is not 
# calculated here. Instead, I will display only the final probabilities
initial_dist <- c(0.04, 0.02, 0.78, 0.16)


# Create the loop for the bootstrapping
n_estim <- 1
while(n_estim <= n_bootstrap){
  
 # Start the dataset simulation
 flog.info(paste0('Bootstrapping sample ', n_estim))
 results_table <- data.table(PATIENT = integer(),
                             state = integer(),
                             obstime = numeric(),
                             deltaobstime = numeric()
 )
 
 
 
 # Calculate the number of patients from the real dataset and produce the same 
 # number of simulated patients. For confidentiality reasons, the number of patients
 # from the real dataset is not calculated here. Instead, I will display only the 
 # final number of total patients.
 n_simulations <- 5932
 death <- 0
 
 # Start the loop for the data simulation
 for (w in 1:n_simulations) {
   if (w == n_simulations) flog.info(paste0('The simulation of the dataset is finished.'))
   
   # Set the initial time and current state of the patient
   t <- 0
   current_state <- sample(1:nstates, 1, prob = initial_dist)
   states <- current_state
   
   sojourn_times <- c()
   observation_times <- t
   states[1] <- current_state
   C <- 12
   iter <- 1
   
   while (iter == 1 || sum(sojourn_times) < C) {
     current_state <- states[iter]
     x <- sojourn_times[iter - 1]
     t <- observation_times[iter]
     iter <- iter + 1
     
     # Calculate the probability matrix Pij based on the current state of the patient
     P <- matrix(0, nstates, nstates)
     for (i in 1:nstates){
       for (j in 1:nstates){
         if (i != j){
           P[i, j] = cpp_p(i, j, aij, bij, t, nstates) 
         }
       }
     }  
     
     
     # Sample the next_state from the probability matrix
     next_state <- sample(nstates, 1, prob = P[current_state, ])
     states[iter] <- next_state
     
     # Sample the sojourn time (x) from weibull distribution
     sojourn_times[iter - 1] <- rweibull(1, vij[states[iter - 1], states[iter]], 
                                         sij[states[iter - 1], states[iter]])
     
     # Calculate the observation_times (t)
     observation_times[iter] <- sum(sojourn_times[1:(iter-1)])
   }
  
   # Append the results of the current simulation to the data table
   results_table <- rbind(results_table, data.table(PATIENT = w, 
                                                    state = states, 
                                                    obstime = observation_times,
                                                    deltaobstime = c(sojourn_times, NA))
   )
 }
 
 
 
 
 # Prepare the results_table for the simulated data parameter estimation
 results_table <- results_table %>% 
   .[, state_prev := shift(state, type = 'lag'), by = PATIENT] %>% 
   .[, state_next := shift(state, type = 'lead'), by = PATIENT] %>% 
   .[, death := 0, by = PATIENT] %>% 
   .[, rowpos := .I]
 
 results_table_c <- results_table[PATIENT %in% results_table[, .N, PATIENT][N > 1, PATIENT]]
 
 
 # ------------------------------ Start the estimation ------------------------------
 
 # Define the initial parameters to start the estimation - They should be the estimated
 # parameters from the real dataset.
 vij_s <- vij
 sij_s <- sij
 aij_s <- aij
 bij_s <- bij
 
 # Put the parameters into one vector for the optimization
 params <- c(vij_s, sij_s, aij_s, bij_s)
 
 # Start the estimation of the parameters. Optimize the loglikelihood function
 result <- optimx(
   par = params, 
   fn = loglikelihood_fast,
   gr = my_gradient,
   method = "BFGS",
   obstimes = results_table_c,
   control = list(
     fnscale = -1,
     trace = 1,
     REPORT = 50,
     maxit = 2000,
     kkt = FALSE
   )
 )
 
 # Put the estimated parameters in a vector and then put them in the matrix that they belong
 params_b <- result[1:(4*nstates**2)] %>% map(1) %>% unlist() %>% c()
 vij_b <- est_vij(params = params_b, nstates = nstates, parscale = parscale)
 sij_b <- est_sij(params = params_b, nstates = nstates, parscale = parscale)
 aij_b <- est_aij(params = params_b, nstates = nstates, parscale = parscale)
 bij_b <- est_bij(params = params_b, nstates = nstates, parscale = parscale)
 
 # Put the estimated parameters from the simulated data in one vector
 params_new <- c(vij_b,sij_b,aij_b,bij_b)
 
 # Store the results of each iteration in the bootstrap estimation list
 estimation_list[[n_estim]] <- params_new
 
 # Start the next iteration
 n_estim <- n_estim + 1
}

# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
 
# Turn the bootstrap estimation list into a data table
boot_sample_logit <- as.data.table(transpose(estimation_list))
estim_list <- copy(estimation_list)

# Calculate some descriptive statistics for the estimated parameters from the simulated data
statistics <- list()

for(q in 1:(4*nstates^2)){
  statistics[[q]] <- t.test(as.numeric(boot_sample_logit[[q]]))
}


# ---------------- Display the bootstrap statistics for each parameter -----------------

# Initialize vectors to hold extracted values
t_values <- p_values <- means <- ci_lower <- ci_upper <- numeric(4*nstates^2)

# Extract values from each t.test object
for (i in 1:(4*nstates^2)) {
  result <- statistics[[i]]
  t_values[i] <- result$statistic
  p_values[i] <- result$p.value
  means[i] <- result$estimate
  ci_lower[i] <- result$conf.int[1]
  ci_upper[i] <- result$conf.int[2]
}


# Convert the metrics into arrays for display (order of parameters: vij, sij, aij, bij)
t_array <- array(t_values, dim = c(nstates, nstates, nstates))
p_array <- array(p_values, dim = c(nstates, nstates, nstates))
mean_array <- array(means, dim = c(nstates, nstates, nstates))
ci_array <- array(paste0("[",round(ci_lower,3),", ", round(ci_upper,3), "]"), 
                  dim = c(nstates, nstates, nstates))

