#### Manual screening for optimal coverage 
source ("fcns.R")
library ("rootSolve")


CL <- seq(0,1,0.05)
CM <- seq(0,1,0.05)
CH <- seq(0,1,0.05)

C <- expand.grid(CL,CM,CH)
C[,4:6] <- 0
colnames(C) <- c("CL","CM","CH","Ro", "SIS","SIR")

# General parameters for all the models
npop <- 3
N <- c(1e6, 1e6, 1e6)
rec <- c(0.5, 0.5, 0.5)
ro <- c(1.2, 5, 15)
alpha <- matrix(c(0.95, 0.02, 0.03, 0.02, 0.96, 0.02, 0.03, 0.02, 0.95), nrow = 3, ncol = 3, byrow = TRUE)
I0 <- c(1, 1, 1)

#Time
tt <- seq(0, 100, length.out = 201)

ss <- multiroot(f = SS_SIS, start = c(1e6, 1e6, 1e6), alpha = alpha, ro = ro, rec = rec, N = N)

sir_mult <- runModel("sir_mult.R", tt, alpha=alpha,rec=rec,ro=ro,N=N, np=npop, I0=I0 )
outbreakSize <- as.numeric(tail(sir_mult, 1)[1,11:13])


for (i in 1:length(C[,1])){
  
cov <- as.numeric(C[i,1:3])
R0 <- cov * N

ss_vacc <-  multiroot(f = SS_SIS_Vacc,start = c(1e6, 1e6, 1e6), alpha = alpha,ro = ro,rec = rec,N = N,cov = cov)
cases_averted_SIS <- (sum((ss$root- ss_vacc$root))/sum(ss$root))


sir_vacc <- runModel("SIR_vacc.R", tt,np=npop, N=N, rec=rec, ro=ro,alpha=alpha,I0=I0, cov=cov, vacc_time=1, R0=R0 )
outbreakSize_vacc <- as.numeric(tail(sir_vacc, 1)[1,11:13])

cases_averted_SIR <- (sum(outbreakSize - outbreakSize_vacc))/sum(outbreakSize) 

C[i,4] <- R0_vacc(alpha, rec, ro, cov)
C[i,5] <- cases_averted_SIS
C[i,6] <- cases_averted_SIR 
  
}

Vmax <- 2e6
sub <- subset(C,(C[,1]*N[1]+ C[,2]*N[2]+ C[,3]*N[3])<Vmax)

newdata <- sub[order(sub[,4]),] 


