##### --------------- One variable equation -----#######

path_logistic <- system.file("examples/logistic.R", package = "odin")

## This is the code with the equations that it is on the file logistics.R
# deriv(N) <- r * N * (1 - N / K)
# initial(N) <- N0
# 
# N0 <- 1
# K <- 100
# r <- 0.5

#Compiles de model
generator <- odin::odin(path_logistic, verbose = FALSE)
# Function that generates an instance of the model
mod <- generator()
mod
# Calculates de diffrential equations with parameters being time and initial conditions
mod$deriv(0, mod$init)
# Shows the content of the model 
mod$contents()

tt <- seq(0, 30, length.out = 101)

#How to run the model
y <- mod$run(tt)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")

# The model can include arbitrary intial conditions 
y2 <- mod$run(tt, 50)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")
lines(y2, col = "red") 




# Generate the same model but user needs to specify parameters
generator <- odin::odin({
  deriv(N) <- r * N * (1 - N / K)
  initial(N) <- N0
  
  N0 <- user(1)
  K <- user(100)
  r <- user()
}, verbose = FALSE)

mod <- generator(r = 1)

y3 <- mod$run(tt)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")
lines(y3, col = "red")


mod$set_user(r = 0.25, K = 75, N0 = 10)
y4 <- mod$run(tt)
plot(y, xlab = "Time", ylab = "N", las = 1, main = "")
lines(y3, col = "red")
lines(y4, col = "blue")


##### --------------- More than one variable equations -----#######


path_lorenz <- system.file("examples/lorenz.R", package = "odin")


# This is what it is on the lorenz.R file:
# deriv(y1) <- sigma * (y2 - y1)
# deriv(y2) <- R * y1 - y2 - y1 * y3
# deriv(y3) <- -b * y3 + y1 * y2
# 
# initial(y1) <- 10.0
# initial(y2) <- 1.0
# initial(y3) <- 1.0
# 
# ## These are the classical parameters:
# sigma <- 10
# R <- 28
# b <-  8 / 3


generator <- odin::odin(path_lorenz, verbose = FALSE)
mod <- generator()

tt <- seq(0, 100, length.out = 20000)
system.time(y <- mod$run(tt))
pairs(y[, -1L], panel = lines, lwd = .2, col = "#00000055")

##### --------------- Delay models -----#######

#delay function uses ylag(t) = y(t - tau). tau can be set up as default or input by the user

gen <- odin::odin({
  ylag <- delay(y, tau)
  initial(y) <- 0.5
  deriv(y) <- 0.2 * ylag * 1 / (1 + ylag^10) - 0.1 * y
  tau <- user(10)
  output(ylag) <- ylag
})

dde <- gen()
t <- seq(0, 300, length.out = 301)
y1 <- dde$run(t)
plot(y1, ylab = "y", mfrow = c(1, 2), which = 1)
plot(y1[, -1L], xlab = "y", ylab = "ylag", mfrow = NULL, type = "l")

#User defines tau
dde$set_user(tau = 20)
y2 <- dde$run(t)
plot(y2, ylab = "y", mfrow = c(1, 2), which = 1)
plot(y2[, -1L], xlab = "y", ylab = "ylag", mfrow = NULL, type = "l")

##### ---------------Models with array equationes  -----#######
gen <- odin::odin({
  deriv(y[]) <- r[i] * y[i] * (1 - sum(ay[i, ]))
  initial(y[]) <- y0[i]
  
  y0[] <- user()
  r[] <- user()
  a[,] <- user()
  ay[,] <- a[i, j] * y[j]
  
  dim(r) <- user()
  n_spp <- length(r)
  
  dim(y) <- n_spp
  dim(y0) <- n_spp
  dim(a) <- c(n_spp, n_spp)
  dim(ay) <- c(n_spp, n_spp)
  
  config(base) <- "lv4"
})

# Parameters and initial conditions
pars <- list(r = c(1.00, 0.72, 1.53, 1.27),
             a = rbind(c(1.00, 1.09, 1.52, 0.00),
                       c(0.00, 1.00, 0.44, 1.36),
                       c(2.33, 0.00, 1.00, 0.47),
                       c(1.21, 0.51, 0.35, 1.00)),
             y0 = c(0.3013, 0.4586, 0.1307, 0.3557))

mod <- gen(user = pars)

t <- seq(0, 2000, length.out = 10001)
y <- mod$run(t)
pairs(y[, -1], panel = lines, col = "#00000055", lwd = 0.2)


## --------------- Debugging ------------------------###

#Generate the model with safe option as TRUE, This will check at every array access (both read and write) that you are not trying to read outside the bounds of the array.

gen <- odin::odin({
  deriv(y[]) <- r[i] * y[i]
  initial(y[]) <- 1
  r[] <- user()
  
  len <- user()
  dim(r) <- user()
  dim(y) <- len
}, safe = TRUE, verbose = FALSE)

####------------- Using interpolation in the model -------------

flux_model <- odin::odin({
  deriv(C) <- flux - kk * C
  initial(C) <- C0
  flux <- interpolate(flux_t, flux_y, "linear")  # flux will be an interpolation of two vectors provided by the user
  C0 <- user()
  kk <- user()
  output(deposition) <- kk * C
  ## Fair bit of boilerplate here that may be removed in future
  ## versions:
  flux_t[] <- user()
  flux_y[] <- user()
  dim(flux_t) <- user()
  dim(flux_y) <- user()
})

flux_t <- c(1, 11, 21, 41, 73, 83, 93, 103, 113, 123, 133, 143, 153,
            163, 173, 183, 194, 204, 214, 224, 234, 244, 254, 264,
            274, 284, 294, 304, 315, 325, 335, 345, 355, 365)
flux_y <- c(0.654, 0.167, 0.06, 0.07, 0.277, 0.186, 0.14, 0.255, 0.231,
            0.309, 1.127, 1.923, 1.091, 1.001, 1.691, 1.404, 1.226, 0.767,
            0.893, 0.737, 0.772, 0.726, 0.624, 0.439, 0.168, 0.28, 0.202,
            0.193, 0.286, 0.599, 1.889, 0.996, 0.681, 1.135)
plot(flux_t, flux_y, type = "l", ylab = "Flux", xlab = "Time")
