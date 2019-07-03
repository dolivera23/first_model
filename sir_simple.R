
# Differential equations for the SIR model with vaccination
# Vaccination is estimated as  vaccinated people starting in the recovered compartmente 

deriv(S) <- - lambda*S*I/N
deriv(I) <- lambda*S*I/N- rec*I
deriv(R) <- rec*I

# Initial conditions

initial(S) <- N - I0 - R0
initial(I) <- I0
initial(R) <- R0


# -------- Outbreak size --------#

deriv(It) <- lambda*S*I/N
initial (It) <- I0

# Parameters 

N <- user()
rec <-user()
ro <- user()
lambda <- ro* rec
I0 <- user(1)
R0 <- user()

config(base) <- "sir"