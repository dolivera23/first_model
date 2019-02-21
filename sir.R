
# Differential equations for the SIR model

deriv(S) <- - lambda*S
deriv(I) <- lambda*S- delta*I
deriv(R) <- delta*I

# Initial conditions

initial(S) <- N - I0
initial(I) <- I0
initial(R) <- 0

# Parameters 

N <- 1e6
delta <-user()
ro <- user()
lambda <- ro* delta
I0 <- user(1)

config(base) <- "sir"