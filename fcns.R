###  This function creates the model and returns a data frame with the model outputs.
##   @param path  Path to the model to be run
##   @param tt    Vector of the time points where the model is going to be run

runModel <- function (path, tt, np, N, rec, ro, I0, alpha, ...) {
  ##-------- Packages required -------- ##
  
  library("odin")
  library("here")
  library("ggplot2")
  library("reshape2")
  
  # Path to  model
  path_mult <- here::here(path)
  
  # Model generator
  model_generator <- odin::odin(path_mult)
  mod <-
    model_generator(
      np = np,
      N = N,
      rec = rec,
      ro = ro,
      alpha = alpha,
      I0 = I0,
      ...
    )
  
  # Running the model
  out_model <- mod$run(tt)
  out_model <- data.frame(out_model)
  
  return (out_model)
  
}

### This function estimates the basic reproductive number for a metapopulation with three patches
#   and with a SIR and SIS dynamics. It uses the Next Generation Matrix
#   @param alpha I a matrix 3x3 with the interaction between patches
#   @param rec Recovery rates for each path
#   @param ro Basic reproductive of each patch if it were isolated

R0 <- function (alpha, rec, ro) {
  # Rate of infection
  beta <- ro * rec
  
  # Transmission matrix
  firstRow <-
    c(beta[1] * alpha[1, 1] / rec[1], beta[1] * alpha[1, 2] / rec[2], beta[1] * alpha[1, 3] /
        rec[3])
  secondRow <-
    c(beta[2] * alpha[2, 1] / rec[1], beta[2] * alpha[2, 2] / rec[2], beta[2] * alpha[2, 3] /
        rec[3])
  thirdRow <-
    c(beta[3] * alpha[3, 1] / rec[1], beta[3] * alpha[3, 2] / rec[2], beta[3] * alpha[3, 3] /
        rec[3])
  
  T <-
    as.matrix(data.frame(firstRow, secondRow, thirdRow), byrow = TRUE)
  
  
  x <- eigen(T)
  
  ro <- max(x$values)
  
  return(ro)
}

### This function estimates the basic reproductive number for a metapopulation with three patches
#   and with a SIR and SIS dynamics. It uses the Next Generation Matrix
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param cov    Vaccine coverage at each path

R0_vacc  <- function (alpha, rec, ro, cov) {
  # Rate of infection
  beta <- ro * rec
  
  # Transmission matrix
  firstRow <-
    c((1 - cov[1]) * beta[1] * alpha[1, 1] / rec[1],
      (1 - cov[1]) * beta[1] * alpha[1, 2] / rec[2],
      (1 - cov[1]) * beta[1] * alpha[1, 3] / rec[3]
    )
  secondRow <-
    c((1 - cov[2]) * beta[2] * alpha[2, 1] / rec[1],
      (1 - cov[2]) * beta[2] * alpha[2, 2] / rec[2],
      (1 - cov[2]) * beta[2] * alpha[2, 3] / rec[3]
    )
  thirdRow <-
    c((1 - cov[3]) * beta[3] * alpha[3, 1] / rec[1],
      (1 - cov[3]) * beta[3] * alpha[3, 2] / rec[2],
      (1 - cov[3]) * beta[3] * alpha[3, 3] / rec[3]
    )
  
  T <-
    as.matrix(data.frame(firstRow, secondRow, thirdRow), byrow = TRUE)
  
  
  x <- eigen(T)
  
  ro <- max(x$values)
  
  return(ro)
}


### Function that returns the infectous equations for 3 interconected patches with SIS disease dynamics.
#   @param  I      Vector with the infectous population (Variable to be solved)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated


SS_SIS <- function(I, alpha, rec, ro, N) {
  beta <- ro * rec
  
  dIL <-
    beta[1] / N[1] * (I[1] + alpha[1, 2] * I[2] + alpha[1, 3] * I[3]) * (N[1] -
                                                                           I[1]) - rec[1] * I[1]
  
  dIM <-
    beta[2] / N[2] * (I[1] * alpha[2, 1] + I[2] + alpha[2, 3] * I[3]) * (N[2] -
                                                                           I[2]) - rec[2] * I[2]
  
  dIH <-
    beta[3] / N[3] * (I[1] * alpha[3, 1] + I[2] * alpha[3, 2] + I[3]) * (N[3] -
                                                                           I[3]) - rec[3] * I[3]
  
  c(dIL = dIL, dIM = dIM, dIH = dIH)
  
}

### Function that returns the infectous equations for 3 interconected patches with SIS disease dynamics with vaccination
#   @param  I      Vector with the infectous population (Variable to be solved)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param cov    Vaccine coverage at each path

SS_SIS_Vacc <- function(I, alpha, rec, ro, cov, N) {
  beta <- ro * rec
  
  dIL <-
    beta[1] / N[1] * (I[1] + alpha[1, 2] * I[2] + alpha[1, 3] * I[3]) * (N[1] *
                                                                           (1 - cov[1]) - I[1]) - rec[1] * I[1]
  
  dIM <-
    beta[2] / N[2] * (I[1] * alpha[2, 1] + I[2] + alpha[2, 3] * I[3]) * (N[2] *
                                                                           (1 - cov[2]) - I[2]) - rec[2] * I[2]
  
  dIH <-
    beta[3] / N[3] * (I[1] * alpha[3, 1] + I[2] * alpha[3, 2] + I[3]) * (N[3] *
                                                                           (1 - cov[3]) - I[3]) - rec[3] * I[3]
  
  c(dIL = dIL, dIM = dIM, dIH = dIH)
  
}



##################################################################################################
#                                   Optimisation functions                                      #
#################################################################################################

### Objective function to minimise the basic reproductive number (Ro) for 3 patch system
#   and with a SIR and SIS dynamics and vaccinations. It uses the Next Generation Matrix
#   @param cov     Variables used in the minimisation (vaccine coverage in each path)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param N      Population size of each path
#   @param Vmax   Number of avaialbe vaccines

f_min  <- function (cov, alpha, rec, ro, N, Vmax) {
  # Rate of infection
  beta <- ro * rec
  
  
  # Transmission matrix
  firstRow <-
    c((1 - cov[1]) * beta[1] * alpha[1, 1] / rec[1],
      (1 - cov[1]) * beta[1] * alpha[1, 2] / rec[2],
      (1 - cov[1]) * beta[1] * alpha[1, 3] / rec[3]
    )
  secondRow <-
    c((1 - cov[2]) * beta[2] * alpha[2, 1] / rec[1],
      (1 - cov[2]) * beta[2] * alpha[2, 2] / rec[2],
      (1 - cov[2]) * beta[2] * alpha[2, 3] / rec[3]
    )
  thirdRow <-
    c((1 - cov[3]) * beta[3] * alpha[3, 1] / rec[1],
      (1 - cov[3]) * beta[3] * alpha[3, 2] / rec[2],
      (1 - cov[3]) * beta[3] * alpha[3, 3] / rec[3]
    )
  
  T <-
    as.matrix(data.frame(firstRow, secondRow, thirdRow), byrow = TRUE)
  
  
  x <- eigen(T)
  
  ro <- max(x$values)
  
  return (ro)
  
  #return((1-ro)^2)
}

### Objective function to maximise the number of cases averted with vaccination for
#   3 patch system with SIS dynamics
#   @param cov     Variables used in the minimisation (vaccine coverage in each path)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param N      Population size of each path
#   @param Vmax   Number of avaialbe vaccines


sis_min <- function(cov, alpha, rec, ro, N, Vmax) {
  
  library ("rootSolve")
  
  ss <- multiroot(
    f = SS_SIS,
    start = c(1e6, 1e6, 1e6),
    alpha = alpha,
    ro = ro,
    rec = rec,
    N = N
  )
  
  
  ss_vacc <-  multiroot(
    f = SS_SIS_Vacc,
    start = c(1e6, 1e6, 1e6),
    alpha = alpha,
    ro = ro,
    rec = rec,
    N = N,
    cov = cov
  )
  
  cases_averted <- -sum(ss$root- ss_vacc$root)
  
  return(cases_averted)
  
}

### Objective function to maximise the number of cases averted with vaccination for
#   3 patch system with SIR dynamics
#   @param cov     Variables used in the minimisation (vaccine coverage in each path)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param N      Population size of each path
#   @param Vmax   Number of avaialbe vaccines


sir_min <- function(cov, alpha, rec, ro, N, Vmax) {
  
  npop <- 3
  I0 <- c(1,1,1)
  
  sir_mult <- runModel("sir_mult.R", tt, alpha=alpha,rec=rec,ro=ro,N=N, np=npop, I0=I0 )
  outbreakSize <- as.numeric(tail(sir_mult, 1)[1,11:13])
  
  sir_vacc <- runModel("SIR_vacc.R", tt,np=npop, N=N, rec=rec, ro=ro,alpha=alpha,I0=I0, cov=cov, vacc_time=1, R0=R0 )
  outbreakSize_vacc <- as.numeric(tail(sir_vacc, 1)[1,11:13])
  
 
 cases_averted <- -100*(sum(outbreakSize -outbreakSize_vacc))/sum(outbreakSize) 
  
  return (cases_averted)
  
} 

##  Inequality contraints for the minimisation process.
#   The number of vaccines used can not excede the number of vaccine available
#   @param cov     Variables used in the minimisation (vaccine coverage in each path)
#   @param alpha   Matrix 3x3 with the interaction between patches
#   @param rec     Recovery rates for each path
#   @param ro     Basic reproductive of each patch if it were isolated
#   @param N      Population size of each path
#   @param Vmax   Number of avaialbe vaccines

g_min <- function (cov, alpha, rec, ro, N, Vmax) {
  return (cov[1] * N[1] + cov[2] * N[2] + cov[3] * N[3] - Vmax)
}
