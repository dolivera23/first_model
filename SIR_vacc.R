np<- user(3)               #  Number of populations to be modelled 

# -------- Differential equations for the SIR model -------- #
dim(S) <- np
dim(I) <- np
dim(R) <- np
deriv(S[1:np]) <- - lambda[i]*S[i] - delta[i]* S[i] 
deriv(I[1:np]) <- lambda[i]*S[i] - rec[i]*I[i]
deriv(R[1:np]) <- rec[i]*I[i] + delta[i] * S[i]

# -------- Initial conditions -------- #

initial(S[1:np]) <- N[i] - I0[i] - R0[i]
initial(I[1:np]) <- I0[i]
initial(R[1:np]) <- R0[i]

# -------- Outbreak size --------#
dim(It) <- np
deriv(It[1:np]) <- lambda[i]*S[i]
initial (It[1:np]) <- I0[i]

# -------- Parameters -------- #

# User parameters
N[] <- user()               #  Population on each patch 
rec[] <-user()              #  Recovery rate on each patch 
ro[] <- user()              #  Reproductive number on each patch  
cov[] <- user()             #  Vaccination coverage on each patch 


# Force of infection
dim(beta) <- np
beta[] <- ro[i] * rec[i]     #  Individual rate of infection on each patch   
alpha[,] <- user()           #  Strength of transmission between populations (patches)


dim(rows) <- c(np,np)
rows[,] <- alpha[i,j]*I[j]
lambda[] <- (beta[i] / N[i]) * sum(rows[i,])


# Vaccination rate 
vacc_time <- user()
dim(delta) <- np 
delta[] <- if(R0[i] < 1 && t >= vacc_time && t < vacc_time +2 ) log(1/(1-cov[i]))  else  0




# Initial conditions
I0[] <- user()
R0[] <- user()

# -------- Dimensions -------- #
dim(rec) <- np
dim(cov) <- np
dim(N) <- np
dim (alpha) <- c(np,np)
dim(ro) <-np
dim(lambda) <- np
dim (I0) <- np
dim(R0) <- np


# -------- Outputs -------- #

#output(lambda[]) <- lambda

config(base) <- "sir_mult"