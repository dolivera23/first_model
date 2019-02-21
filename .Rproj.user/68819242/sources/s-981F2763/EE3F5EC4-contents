# Load file with all the functions
source ("fcns.R")


# General parameters for all the models
npop <- 3
N <- c(1e6, 1e6, 1e6)
rec <- c(0.5, 0.5, 0.5)
ro <- c(1.2, 5, 15)
alpha <-
  matrix(
    c(1, 0.7, 0.3, 0.7, 1, 0.7, 0.3, 0.7, 1),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
I0 <- c(1, 1, 1)

#Time
tt <- seq(0, 100, length.out = 201)

##################################################################################################
#                                   SIR MODELS                                                   #
#################################################################################################


## ----- Basic SIR model  ---->

sir_mult <-   runModel(
    "sir_mult.R",
    tt,
    np = npop,
    N = N,
    rec = rec,
    ro = ro,
    alpha = alpha,
    I0 = I0
  )

outbreakSize <- as.numeric(tail(sir_mult, 1)[1, 11:13])

# Plotting the model
sir_mult_long = melt(sir_mult, id.vars = "t")
ggplot(sir_mult_long, aes(x = t, y = value, col = variable)) + geom_line() + coord_cartesian(xlim =
                                                                                               c(0, 20))


## --------- Vaccination SIR model ---->

cov <- 1 - 1 / ro
R0 <- cov * N

sir_vacc <- runModel(
    "SIR_vacc.R",
    tt,
    np = npop,
    N = N,
    rec = rec,
    ro = ro,
    alpha = alpha,
    I0 = I0,
    cov = cov,
    vacc_time = 1,
    R0 = R0
  )

# Plotting the model
sir_vacc_long = melt(sir_vacc, id.vars = "t")
ggplot(sir_vacc_long, aes(x = t, y = value, col = variable)) + geom_line()#+ coord_cartesian(xlim=c(0, 20))


#  ------------ Different coverages:

cov_range <- seq(0, 1, length.out = 100)
cases_averted <-
  data.frame(matrix(1, nrow = length(cov_range), ncol = 3))
colnames (cases_averted) <- c("Low", "Medium", "High")

for (i in c(1:length(cov_range))) {
  cov <- c(cov_range[i], cov_range[i], cov_range[i])
  R0 <- cov * N
  
  #tryCatch({
    print(i)
    
    
    sir_vacc <-
      runModel(
        "SIR_vacc.R",
        tt,
        np = npop,
        N = N,
        rec = rec,
        ro = ro,
        alpha = alpha,
        I0 = I0,
        cov = cov,
        vacc_time = 1,
        R0 = R0
      )
    
    
    
 # }, error = function(e) {
  #  cat("ERROR :", conditionMessage(e), "\n")
 # })
  
  outbreakSize_vacc <- as.numeric(tail(sir_vacc, 1)[1, 11:13])
  
  cases_averted[i, ] <-
    100 * (outbreakSize - outbreakSize_vacc) / outbreakSize
  
}

cases_averted [, 4] <- cov_range
sir_averted = melt(cases_averted, id.vars = "V4")
# Plotting the model
ggplot(sir_averted, aes(x = V4, y = value, col = variable)) + geom_line() + xlab(("Coverage")) + ylab ("% Cases averted") + labs(colour = "Patch")

##################################################################################################
#                                   SIS MODELS                                                   #
#################################################################################################


## -------- Basic SIS  ------->

#Steady state
library ("rootSolve")
ss <-
  multiroot(
    f = SS_SIS,
    start = c(1e6, 1e6, 1e6),
    alpha = alpha,
    ro = ro,
    rec = rec,
    N = N
  )


# Dinamic model
sis_mult <-
  runModel(
    "sis.R",
    tt,
    np = npop,
    N = N,
    rec = rec,
    ro = ro,
    alpha = alpha,
    I0 = I0
  )

# Plotting the dynamic  model
sis_long = melt(sis_mult, id.vars = "t")
ggplot(sis_long, aes(x = t, y = value, col = variable)) + geom_line()

## --------- Vaccination SIS model ------->

cov <- 1 - 1 / ro

# Solving steady state
ss_vacc <-
  multiroot(
    f = SS_SIS_Vacc,
    start = c(1e6, 1e6, 1e6),
    alpha = alpha,
    ro = ro,
    rec = rec,
    N = N,
    cov = cov
  )

# Dynamic model
I0_SIS <- as.numeric(tail(sis_mult, 1)[1, 5:7])
sis_vacc <-
  runModel(
    "SIS_vacc.R",
    tt,
    np = npop,
    N = N,
    rec = rec,
    ro = ro,
    alpha = alpha,
    I0 = I0_SIS,
    cov = cov,
    vacc_time = 20
  )

# Plotting the model
sis_vacc_long = melt(sis_vacc, id.vars = "t")
ggplot(sis_vacc_long, aes(x = t, y = value, col = variable)) + geom_line() +
  coord_cartesian(xlim = c(19, 27))


# -------------  Different coverages ----->

cov_range <- seq(0, 1, length.out = 100)
cases_averted_sis <-
  data.frame(matrix(0, nrow = length(cov_range), ncol = 3))
colnames (cases_averted_sis) <- c("Low", "Medium", "High")

for (i in c(1:length(cov_range))) {
  cov <- c(cov_range[i], cov_range[i], cov_range[i])
  
  ss_vacc <-
    multiroot(
      f = SS_SIS_Vacc,
      start = c(1e6, 1e6, 1e6),
      alpha = alpha,
      ro = ro,
      rec = rec,
      N = N,
      cov = cov
    )
  
  I_vacc <- ss_vacc$root
  
  cases_averted_sis [i, ] <- 100 * (ss$root - ss_vacc$root) / ss$root
}

cases_averted_sis [, 4] <- cov_range
sis_averted = melt(cases_averted_sis, id.vars = "V4")
# Plotting the model
ggplot(sis_averted, aes(x = V4, y = value, col = variable)) + geom_line() + xlab(("Coverage")) + ylab ("% Cases averted") + labs(colour = "Patch")#+ coord_cartesian(xlim=c(0, 0.9))


##################################################################################################
#                                          R0                                                   #
#################################################################################################

cov_range <- seq(0, 1, length.out = 100)
ro_vec <- data.frame(matrix(0, nrow = length(cov_range), ncol = 1))


for (i in c(1:length(cov_range))) {
  cov <- c(cov_range[i], cov_range[i], cov_range[i])
  
  ro_vec[i, 1] <- R0_vacc(alpha, rec, ro, cov)
  
}

ro_vec [, 2] <- cov_range
colnames(ro_vec) <- c("ro", "cov")

# Plotting the model
ggplot(ro_vec, aes(x = cov, y = ro)) + geom_line() + xlab(("Coverage")) + ylab ("Ro") #+ coord_cartesian(xlim=c(0, 0.9))





##################################################################################################
#                                    Minimisation                                               #
#################################################################################################
library("nloptr")

# RO minimization for both SIS and SIR models

# NLOPT_LN_COBYLA - NLOPT_GN_ISRES - NLOPT_GN_DIRECT - NLOPT_GN_DIRECT_L

cov_0 <- c(0.5,0.5,0.5)
opts <- list("algorithm" = "NLOPT_GN_DIRECT_L", "xtol_rel" = 1.0e-8)

R0_min <- nloptr(
    x0 = cov_0,
    eval_f = f_min,
    eval_g_ineq = g_min,
    lb = c(0, 0, 0),
    ub = c(1, 1, 1),
    opts = opts,
    alpha = alpha,
    rec = rec,
    ro = ro,
    N = N,
    Vmax = 1.5e6
  )

min_sis <- 
  nloptr(
    x0 = cov_0,
    eval_f = sis_min ,
    eval_g_ineq = g_min,
    lb = c(0, 0, 0),
    ub = c(1, 1, 1),
    opts = opts,
    alpha = alpha,
    rec = rec,
    ro = ro,
    N = N,
    Vmax = 1.5e6
  )

min_sir <-
  nloptr(
    x0 = cov_0,
    eval_f = sir_min ,
    eval_g_ineq = g_min,
    lb = c(0, 0, 0),
    ub = c(1, 1, 1),
    opts = opts,
    alpha = alpha,
    rec = rec,
    ro = ro,
    N = N,
    Vmax = 1.5e6
  )

