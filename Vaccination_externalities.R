# Load function file and required packages
source ("./fcns.R")
library ("rootSolve")
library("ggpubr")


## ----- Parameters ----- ##
Vmax<- 3e6
vaccine <- seq(0,Vmax,100)

npop <- 3
N <- c(1e6, 1e6, 1e6)
rec <- c(0.5, 0.5, 0.5)
ro <- c(1.5,5,15)
alpha <-  matrix(c(0.95, 0.02, 0.03, 0.02, 0.96, 0.02, 0.03, 0.02, 0.95), nrow = 3,  ncol = 3,byrow = TRUE)
I0 <- c(1, 1, 1)


#Time
tt <- seq(0, 100, length.out = 201)


## ----  Output data ---- ##
vacc_out <- data.frame(matrix(0,nrow=length(vaccine),ncol=9))
colnames(vacc_out) <- c ("Vaccine", "SIS_L", "SIS_M","SIS_H", "SIS", "SIR_L","SIR_M","SIR_H", "SIR")
vacc_out[,1] <- vaccine

pv_out <- data.frame(matrix(0,nrow=length(vaccine),ncol=9))
colnames(pv_out) <-c ("Vaccine","SIS_L", "SIS_M","SIS_H", "SIS", "SIR_L","SIR_M","SIR_H", "SIR")
pv_out[,1] <-vaccine


## ---- Model results without interventions --- ##

# SIS model
ss <- multiroot(f = SS_SIS, start = c(1e6, 1e6, 1e6), alpha = alpha, ro = ro, rec = rec,N = N)
endemic <- ss$root 

# SIR model 
sir_mult <- runModel("sir_mult.R",tt,np = npop,N = N,rec = rec,ro = ro,alpha = alpha,I0 = I0)
outbreakSize <- tail(sir_mult, 1)[1, 11:13]
colnames(outbreakSize) <- c("Low transmission", "Medium transmission", "High transmission")
outbreakSize <- as.numeric (outbreakSize)

## ---- Vaccination screening ---- ##

for (i in 1:length(vaccine)){
  
  cov <- vaccine[i]/(3*N)        # Assuming equal distribution of vaccines so far 
  R0 <- cov * N                  # Initial population of recovered people for the SIR model
  
  
  # SIS model: 
  ss_vacc <-  multiroot(f = SS_SIS_Vacc, start = c(1e6, 1e6, 1e6),alpha = alpha, ro = ro,rec = rec,N = N,cov = cov)
  vacc_out[i,c(2,3,4)] <- (endemic - ss_vacc$root)
  vacc_out[i,5] <- sum((endemic - ss_vacc$root))
  
  pv_sis <- (ss_vacc$root - I0)/(N-I0)
  pv_out[i,c(2,3,4)] <- pv_sis
  pv_out[i,5] <- (sum(ss_vacc$root - I0)) / (sum((N-I0)))
  
  # SIR model: 
  sir_vacc <- runModel("SIR_vacc.R", tt,np = npop,N = N, rec = rec,ro = ro, alpha = alpha, I0 = I0,cov = cov, vacc_time = 1,R0 = R0)
  outbreakSize_vacc <- as.numeric(tail(sir_vacc, 1)[1, 11:13])
  
  vacc_out[i,c(6,7,8)] <- (outbreakSize - outbreakSize_vacc)
  vacc_out[i,9] <- sum(outbreakSize - outbreakSize_vacc)
  
  pv_sir <- as.numeric((outbreakSize_vacc - I0) / (sir_vacc[1,c(2,3,4)]))
  pv_out[i,c(6,7,8)] <- pv_sir
  pv_out[i,9] <- (sum(outbreakSize_vacc - I0)) / sum(sir_vacc[1,c(2,3,4)])
  
 
  
}

# Marginal Effect of Vaccination (MEV): Number of cases prevented per each vaccine added  

MEV <- vacc_out

for (i in 2:length(vaccine)){
  
  MEV[i, c(2:9)] <- (vacc_out[i, c(2:9)] - vacc_out[i-1, c(2:9)]) / (vacc_out[i,1] - vacc_out[i-1,1])
  
}

## ---- Plotting ---- ##

plot_externalities <- function (db, columns,rows, grouping, xlb, ylb){
  
  subset_db <- db[rows ,columns]
  subset_long <- melt(subset_db, id.vars = grouping)
  
  p <- ggplot(subset_long, aes(x = subset_long[,1], y = value , col = variable)) + geom_line()+ xlab(xlb) + ylab(ylb) 
  
  return (p)
  
}

p_MEV <- plot_externalities(MEV, c(1,5,9), 2:length(vaccine), "Vaccine", "Vaccinations", "MEV")
p_MEV_SIS <- plot_externalities (MEV,c(1,2,3,4), 2:length(vaccine),"Vaccine", "Vaccinations", "MEV" )
p_MEV_SIR <- plot_externalities (MEV,c(1,6,7,8), 2:length(vaccine), "Vaccine", "Vaccinations", "MEV" )

ggarrange(ggarrange(p_MEV_SIS, p_MEV_SIR, ncol = 2, labels = c("A", "B")), p_MEV ,nrow = 2)

p_pv <- plot_externalities(pv_out, c(1,5,9), 1:length(vaccine),"Vaccine", "Vaccinations", "P(v)")
p_pv_SIS <- plot_externalities(pv_out, c(1,2,3,4), 1:length(vaccine),"Vaccine", "Vaccinations", "P(v)")
p_pv_SIR <- plot_externalities(pv_out, c(1,6,7,8),2:length(vaccine), "Vaccine", "Vaccinations", "P(v)")

ggarrange(ggarrange(p_pv_SIS, p_pv_SIR, ncol = 2, labels = c("A", "B")), p_pv,nrow = 2)




