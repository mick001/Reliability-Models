rm(list=ls())

#-------------------------------------------------------------------------------
# Boolean table of the system

# Parallel system
bool_table <- matrix(c(0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1), nrow=4, ncol=3, byrow = T)
# Syst with 5 comp
#load(file="C:\\users\\michy\\desktop\\reliability-experiment\\bool-table.Rdata")

#-------------------------------------------------------------------------------

# Time horizon considered (in years) for reliability calculation
T_HORIZON <- 450
# Time mission considered in years
T_M <- 200
# Time for calculating conditional reliability in years
T_COND <- 20
# Number of Monte Carlo iterations
N_MC_SIM <- 5000
# Number of system components
N_COMP <- ncol(bool_table) - 1
# For reproducibility
set.seed(2561)
# System state matrix
system_state_matrix <- matlab::zeros(N_MC_SIM, T_HORIZON + 1)
# Fix equal to 1 element k, 1 of the state matrix for each row since R(0)=1
system_state_matrix[1:N_MC_SIM, 1] <- 1
# All possible existing combinations of system state
ALL_POSSIBLE_STATES <- bool_table[, 1:N_COMP]
# Working and not working states: working statets are those with 1
SYS_STATES <- bool_table[, N_COMP + 1]

#-------------------------------------------------------------------------------
# Failure rate distribution modelling and dependence structure

# Failure rate of each component [years^(-1)]
LAMBDA = 0.01
# Dependence structure. Can try different copulas...
COPULA_TYPE <- copula::tCopula(param=-0.5, df=1,dim=N_COMP)#copula::normalCopula(param=-0.9, dim = N_COMP)#copula::indepCopula(dim = N_COMP)
# Can try different margins for instance weibull
MARGINS <- rep("exp", N_COMP)
PARAM_MARGINS <- rep(list(rate=LAMBDA), 2)

# Multivariate distribution from copula
my_dist <- copula::mvdc(COPULA_TYPE,
                        margins = MARGINS,
                        paramMargins = PARAM_MARGINS)
# Sample N_MC_SIM time to failure from the distribution
sim_ttf <- copula::rMvdc(N_MC_SIM, my_dist)

#PerformanceAnalytics::chart.Correlation(sim_ttf)

#-------------------------------------------------------------------------------
# Monte Carlo iterations

# Monte Carlo iterations outer loop
for(k in 1:N_MC_SIM)
{
    # Components state matrix: Ncomp x (T_horizon + 1)
    components_state_matrix = matlab::zeros(N_COMP, T_HORIZON + 1)
    #% Fix equal to 1 element (k, 1) of the state matrix for each row k
    #% since reliability is 1 at beginning of the period
    components_state_matrix[1:N_COMP, 1] = 1

    #% Sample tf for each component
    for(i in 1:N_COMP)
    {
        # Sample tf (tf is in years)
        #tf = rexp(1, rate=LAMBDA)
        tf = sim_ttf[k, i]
        #% From year 1 to ceil(tf) component has worked, set 1.
        #% Note that index is shifted! Year 1 is at index 2.
        components_state_matrix[i, 2:min(ceiling(tf), T_HORIZON + 1)] = 1
        #% From year ceil(tf) to year T_m component has not worked, set 0
        components_state_matrix[i, min(ceiling(tf) + 1, T_HORIZON + 1):(T_HORIZON + 1)] = 0
    }
    
    # Determine system history.
    for(i in 2:(T_HORIZON + 1))
    {
    
        # State system at time i is equal to column i of
        # components_state_matrix
        current_system_state = components_state_matrix[1:N_COMP, i]
        
        # Find if current_system_state is a working state
        for(j in 1:2^N_COMP)
        {
            # This finds the corresponding state from all possible states
            index_sum = sum(current_system_state == ALL_POSSIBLE_STATES[j, 1:N_COMP])
            # Once found the corresponding state we exit the loop
            if(index_sum == N_COMP)
            {
                break
            }
        }
    
        # Corresponding state is the j-th element of SYS_STATES
        system_status = SYS_STATES[j]
    
        # If system_status is 0 system is not working
        if(system_status == 0)
        {
            system_state_matrix[k, i:(T_HORIZON+1)] = 0 
            break
        }
    
        # Else the system is still working
        system_state_matrix[k, i] = 1
    }
}

#-------------------------------------------------------------------------------
# Calculate quantities of interest

# System reliability R(t)
R <- colSums(system_state_matrix)/N_MC_SIM

# MTTF: integrating with unit increment
MTTF <- sum(R)

# Reliability at T_M
R_TM <- R[T_M + 1]

# Conditional reliability at T = 20 years -> R(t + T)/R(T)
cond_R <- R[(T_COND*2 + 1):(T_HORIZON + 1)]/R[T_COND + 1]
cond_t_axis <- T_COND:(T_HORIZON - T_COND)

# Residual reliability
residual_reliability <- mean(c(R[floor(MTTF) + 1], R[ceiling(MTTF) + 1]))

# Print results

# Print MTTF
paste('MTTF', MTTF, 'years')

# Reliability at mission time T_M
paste('Reliability at mission time', R_TM)

# Residual reliability
paste('Residual reliability', residual_reliability)

# Estimated probability that system survives 20 years
paste('Estimated reliability at 20 years', R[21])

#-------------------------------------------------------------------------------
# Plot results

# Plot of reliability
plot(0:T_HORIZON, R, type="l", xlab = 'Time [years]', ylab = 'R(t)', main='System Reliability', xlim=c(0 , T_HORIZON), ylim=c(0, 1))

# Plot conditional reliability P(tf > t + T | tf > T)
lines(cond_t_axis, cond_R, lwd=1.5, col='red')

# Plot a 1 indicator for reference %refline(0, 1)
abline(v=T_COND, col='blue', lwd= 0.5)

# Reliability in conditions of indepence betweeen the failure rates: analytical expression
lines(0:T_HORIZON, 2*exp(-LAMBDA*0:T_HORIZON) - exp(-2*LAMBDA*0:T_HORIZON), col="green")

legend("topright", legend=c("R(t) MC", "R(t) | T_M -- MC", "Ref T_M", "R(t) analytic"),
       col=c("black", "red", "blue", "green"), lty=rep(1, 4), cex=0.8)
