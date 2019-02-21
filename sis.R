np<- user(3)               #  Number of populations to be modelled 

# -------- Differential equations for the SIS model -------- #
dim(S) <- np
dim(I) <- np

dim(beta) <- np
 
deriv(S[1:np]) <- - lambda[i]*S[i] +rec[i]*I[i]
deriv(I[1:np]) <- lambda[i]*S[i] - rec[i]*I[i]


# -------- Initial conditions -------- #

initial(S[1:np]) <- N[i] - I0[i]
initial(I[1:np]) <- I0[i]


# -------- Parameters -------- #


N[] <- user()                 #  Population on each patch 
rec[] <-user()                #  Recovery rate on each patch 
ro[] <- user()                #  Reproductive number on each patch  
beta[] <- ro[i]* rec [i]    #  Individual rate of infection on each patch   
alpha[,] <-user()           #  Strength of transmission between populations (patches)


# Force of infection
dim(rows) <- c(np,np)
rows[,] <- alpha[i,j]*I[j]
lambda[] <- (beta[i] /N[i]) * sum(rows[i,])
# Initial conditions
I0[] <- user()


# -------- Dimensions -------- #
dim(N) <- np
dim(rec) <- np
dim (alpha) <- c(np,np)
dim(ro) <-np
dim(lambda) <- np
dim (I0) <- np

# -------- Outputs -------- #

#output(lambda[]) <- lambda

config(base) <- "sis_mult"





