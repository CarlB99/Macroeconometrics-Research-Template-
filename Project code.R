library(dplyr)
library(tidyr)
library(ggplot2)
library(ecb)
library(gridExtra)
library(bsvars)
library(lubridate)
library(tseries)
library(zoo)
library(tempdisagg)
library(progress)
library(mvtnorm)
library(MASS)
library(HDInterval)
library(mgcv)

rm(list = ls())

# Downloading Data
start_date_m <- "2003-01" 
end_date_m <- "2023-12"

# Define the datasets
series_info <- list(
  ip = "STS.M.I8.Y.PROD.NS0010.4.000",
  ciss = "CISS.M.U2.Z0Z.4F.EC.SOV_EW.IDX",
  hicp = "ICP.M.U2.Y.000000.3.INX",
  un = "LFSI.M.I9.S.UNEHRT.TOTAL0.15_74.T",
  m2 = "BSI.M.U2.Y.V.M20.X.1.U2.2300.Z01.E"
)

# Making function to collect data
prepare_data <- function(series_id, name, start_date, end_date) {
  get_data(series_id, filter = list(startPeriod = start_date, endPeriod = end_date)) |>
    transmute(obstime, !!name := obsvalue)
}

# Collecting data
datasets <- lapply(names(series_info), function(name) {
  prepare_data(series_info[[name]], name, start_date_m, end_date_m)
})

# Collecting daily data
mro <- get_data("FM.B.U2.EUR.4F.KR.MRR_FR.LEV") |>
  mutate(obstime = as.Date(obstime))

# Creating a full sequence of dates for the range and merging with original data
mro <- data.frame(obstime = seq(from = as.Date("2003-01-01"), to = as.Date("2023-12-31"), by = "day")) |>
  left_join(mro, by = "obstime") |>
  mutate(obsvalue = na.locf(obsvalue, na.rm = FALSE)) |>
  mutate(obsvalue = na.locf(obsvalue, fromLast = TRUE, na.rm = FALSE)) |>
  dplyr::select(obstime, obsvalue)

# Transforming to monthly data
mro <- mro |>
  group_by(obstime = as.yearmon(obstime)) |>
  summarize(mro = last(obsvalue)) |>
  mutate(obstime = as.Date(obstime))

# Collecting quarterly data
consumption <- get_data("MNA.Q.Y.I9.W0.S1M.S1.D.P31._Z._Z._T.EUR.V.N") |>
  dplyr::select(obstime, obsvalue) |>
  slice(-1:-32)  # Removes the first years
consumption_ts <- ts(consumption$obsvalue, start = c(2003, 1), frequency = 4)
monthly_consumption <- td(consumption_ts ~1, to = "monthly", method = "chow-lin-maxlog", conversion = "average")
start_date <- as.Date("2003-01-01")
monthly_dates <- seq(from = start_date, by = "month", length = length(monthly_consumption$values))
consumption <- data.frame(obstime = monthly_dates, consumption = log(monthly_consumption$values) )

# Merge all datasets by 'obstime'
Data_function <- Reduce(function(x, y) left_join(x, y, by = "obstime"), datasets)

Data <- Data_function |>
  mutate(across(c(2, 4, 6), ~log(.)))

# Convert 'obstime' to Date format
Data$obstime <- as.Date(paste0(Data$obstime, "-01"))

# Appending data
Data <- Data |>
  left_join(mro, by = "obstime") |>
  left_join(consumption, by = "obstime")

# Define the plot configurations
plot_settings <- list(
  list(name = "ip", title = "Industrial production", ylab = "LOG(IP)"),
  list(name = "hicp", title = "Inflation", ylab = "LOG(HICP)"),
  list(name = "ciss", title = "Financial Stress", ylab = "CISS"),
  list(name = "un", title = "Unemployment", ylab = "UN"),
  list(name = "m2", title = "M2 - Money Stock", ylab = "LOG(M2)"),
  list(name = "mro", title = "ECB interest rate", ylab ="MRO"),
  list(name = "consumption", title = "Private Consumption Euro Area", ylab="LOG(C)")
)

# Create the plots
plots <- lapply(plot_settings, function(setting) {
  ggplot(Data, aes_string(x = "obstime", y = setting$name)) +
    geom_line() +
    labs(x = "", y = setting$ylab, title = setting$title) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill=NA),  
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      plot.background = element_blank(),  
      plot.title = element_text(hjust = 0.5)  
    )
})

# Arrange all plots in a grid
grid.arrange(grobs = plots, ncol = 2)

# Change to time series
Data <- as.ts(Data)


# ACF and PACF
par(mfrow = c(3, 3))

variable_names <- c("ip", "hicp", "ciss", "un", "m2", "mro", "consumption")

# Loop through the columns and create ACF plots
for (i in seq_along(variable_names)) {
  # Compute and plot ACF
  acf(Data[, i], main = paste("ACF for", variable_names[i]))
}

par(mfrow = c(3, 3))
for (i in seq_along(variable_names)) {
  # Compute and plot PACF
  pacf(Data[, i], main = paste("PACF for", variable_names[i]))
}


for (i in seq_along(variable_names)) {
  variable_name <- paste("adf", variable_names[i], sep = "_")
  assign(variable_name, adf.test(Data[, i], k=12, alternative = "stationary"))
}

adf_results <- data.frame(
  Variable = c("IP", "HICP", "CISS", "UN", "M2", "MRO", "Consumption"),
  ADF_Statistic = c(adf_hicp$statistic, adf_ip$statistic, adf_ciss$statistic, adf_un$statistic, adf_m2$statistic, adf_mro$statistic, adf_consumption$statistic),
  P_Value = c(adf_hicp$p.value, adf_ip$p.value, adf_ciss$p.value, adf_un$p.value, adf_m2$p.value, adf_mro$p.value, adf_consumption$p.value),
  lags = rep(12, 7)
)

knitr::kable(adf_results, digits = 3, align = 'c')

### Baseline model BVAR estimation:

# Setting specifications
N = ncol(Data[ , -1])
p = 12
K = 1+N*p
S = 10000
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -1], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}

# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1 
)

# BVAR function

BVAR = function(Y,X,priors,S){
  
  # normal-inverse Wishart posterior parameters
  V.bar.inv   = t(X)%*%X + diag(1/diag(priors$V.prior))
  V.bar       = solve(V.bar.inv)
  A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(priors$V.prior))%*%priors$A.prior)
  nu.bar      = nrow(Y) + priors$nu.prior
  S.bar       = priors$S.prior + t(Y)%*%Y + t(priors$A.prior)%*%diag(1/diag(priors$V.prior))%*%priors$A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
  S.bar.inv   = solve(S.bar)
  
  #posterior draws
  Sigma.posterior   = rWishart(sum(S), df=nu.bar, Sigma=S.bar.inv)
  Sigma.posterior   = apply(Sigma.posterior,3,solve)
  Sigma.posterior   = array(Sigma.posterior,c(N,N,sum(S)))
  A.posterior       = array(rnorm(prod(c(dim(A.bar),sum(S)))),c(dim(A.bar),sum(S)))
  L                 = t(chol(V.bar))
  
  for (s in 1:sum(S)){
    A.posterior[,,s]= A.bar + L%*%A.posterior[,,s]%*%chol(Sigma.posterior[,,s])
  }
  
  posterior = list(
    Sigma.posterior   = Sigma.posterior,
    A.posterior       = A.posterior
  )
  return(posterior)
}

# Applying BVAR function
posterior.draws = BVAR(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws$A.posterior, 1:2, mean),3)

par(mfrow=c(2,2))
hist(posterior.draws$Sigma.posterior[2,1,1:sum(S)], col="#379683", xlab="Value of Sigma", ylab="Frequency", main ="Histogram of draws of Sigma")
plot.ts(posterior.draws$Sigma.posterior[2,1,1:sum(S)],col="#379683", xlab="Iteration s", ylab=NULL, main="Plot of posterior draws of sigma")
hist(posterior.draws$Sigma.posterior[1,1,1:sum(S)], col="#379683", xlab="Value of A", ylab="Frequency", main="Histogram of draws of A")
plot.ts(posterior.draws$Sigma.posterior[1,1,1:sum(S)],col="#379683", xlab="Iteration s", ylab=NULL, main="Plot of posterior draws of A")


### Showing baseline model
# Specifications
# p = 1
# N = 2
# K = 1+N*p
# S = 1000
# 
# # generate random walk processes
# rw.1    = cumsum(rnorm(1000,0,1))
# rw.2    = cumsum(rnorm(1000,0,1))
# y       = matrix(cbind(rw.1,rw.2),nrow=1000,ncol=N)
# 
# # Set X and Y matrices
# Y       = ts(y[2:nrow(y),])
# X       = matrix(1,nrow(Y),1)
# X       = cbind(X,y[1:nrow(y)-p,])
# 
# # MLE
# A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
# Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/nrow(Y)
# 
# # Minnesota prior
# kappa.1           = 0.02^2
# kappa.2           = 100
# A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
# A.prior[2:(N+1),] = diag(N)
# 
# priors = list(
#   A.prior     = A.prior,
#   V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
#   S.prior     = diag(diag(Sigma.hat)),
#   nu.prior    = N+1
# )
# 
# # Demonstration
# posterior.draws = BVAR(Y=Y, X=X, priors=priors, S=S)
# round(apply(posterior.draws$Sigma.posterior, 1:2, mean),3)
# round(apply(posterior.draws$A.posterior, 1:2, mean),3)




### Extended model BVAR estimation

# Setting specifications
N = ncol(Data[ , -1])
p = 12
K = 1+N*p
S = 10000
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -1], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}

# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1 
)

m_a = 1
s2  = 0.1

a.posterior.store = numeric(S)
Sigma.posterior   = rWishart(1, df=priors$nu.prior, Sigma=priors$S.prior)[,,1]
V.bar             = priors$V.prior
A.posterior       = priors$A.prior

# Posterior of underline a:

compute_sbar2 = function(s2, Sigma, V) {
  sbar2 = 1/s2
  for (i in 1:nrow(Sigma)) {
    sbar2 = sbar2 + 1/(Sigma[i,i]*V[i+1,i+1])
  }
  return(1/sbar2)
}

compute_a_posterior = function(sbar2, Sigma,V,A){
  a_posterior = m_a/sbar2
  n=nrow(Sigma)
  A_1 = A[2:nrow(A),]
  for(i in 1:n) {
    a_posterior = a_posterior + A_1[i,i]/(Sigma[i,i]*V[i+1,i+1])
  }
  return(sbar2*a_posterior)
}

# BVAR function
BVAR.extension = function(Y,X,priors,S){
  
  A.posterior        = array(0.0001, dim = c(K,N,sum(S)))
  Sigma.posterior    = array(0.0001,dim=c(N,N,sum(S)))
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = sum(S), style = 3)
  
  for (s in 1:sum(S)){
    sbar2                = compute_sbar2(s2, Sigma.posterior[,,s],V.bar)
    a_posterior          = compute_a_posterior(sbar2, Sigma.posterior[,,s],V.bar, A.posterior[,,s])
    a.posterior          = rnorm(1, a_posterior, sqrt(sbar2))
    a.posterior.store[s] = a.posterior
    
    # normal-inverse Wishart posterior parameters
    V.bar.inv   = t(X)%*%X + diag(1/diag(priors$V.prior))
    V.bar       = solve(V.bar.inv)
    A.bar       = V.bar%*%(t(X)%*%Y + diag(1/diag(priors$V.prior))%*%priors$A.prior*a.posterior)
    nu.bar      = nrow(Y) + priors$nu.prior
    S.bar       = priors$S.prior + t(Y)%*%Y + t(priors$A.prior*a.posterior)%*%diag(1/diag(priors$V.prior))%*%priors$A.prior*a.posterior - t(A.bar)%*%V.bar.inv%*%A.bar
    S.bar.inv   = solve(S.bar)
    
    #posterior draws
    Sigma.posterior.dist   = rWishart(1, df=nu.bar, Sigma=S.bar.inv)
    Sigma.draw             = apply(Sigma.posterior.dist,3,solve)
    Sigma.posterior[,,s]   = Sigma.draw
    A.posterior[,,s]       = array(rnorm(prod(c(dim(A.bar),1))),c(dim(A.bar),1))
    L                      = t(chol(V.bar))
    A.posterior[,,s]       = A.bar + L%*%A.posterior[,,s]%*%chol(Sigma.posterior[,,s])
    
    
    
    # Update progress bar
    setTxtProgressBar(pb, s)
  }
  # Close progress bar
  close(pb)
  
  posterior = list(
    Sigma.posterior   = Sigma.posterior,
    A.posterior       = A.posterior
  )
  return(posterior)
}

# Applying BVAR function
posterior.draws.extension = BVAR.extension(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws.extension$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws.extension$A.posterior, 1:2, mean),3)

par(mfrow=c(2,2))
hist(posterior.draws.extension$Sigma.posterior[2,1,1:sum(S)], col="#379683", xlab="Value of Sigma", ylab="Frequency", main ="Histogram of draws of Sigma")
plot.ts(posterior.draws.extension$Sigma.posterior[2,1,1:sum(S)],col="#379683", xlab="Iteration s", ylab=NULL, main="Plot of posterior draws of sigma")
hist(posterior.draws.extension$Sigma.posterior[1,1,1:sum(S)], col="#379683", xlab="Value of A", ylab="Frequency", main="Histogram of draws of A")
plot.ts(posterior.draws.extension$Sigma.posterior[1,1,1:sum(S)],col="#379683", xlab="Iteration s", ylab=NULL, main="Plot of posterior draws of A")


## Showing extended model
#Specifications
p = 1
N = 2
K = 1+N*p
S = 1000

# generate random walk processes
rw.1    = cumsum(rnorm(1000,0,1))
rw.2    = cumsum(rnorm(1000,0,1))
y       = matrix(cbind(rw.1,rw.2),nrow=1000,ncol=N)

# Set X and Y matrices
Y       = ts(y[2:nrow(y),])
X       = matrix(1,nrow(Y),1)
X       = cbind(X,y[1:nrow(y)-p,])

# MLE
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/nrow(Y)

# Minnesota prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1
)

# Demonstration
posterior.draws.extension.sim = BVAR.extension(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws.extension.sim$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws.extension.sim$A.posterior, 1:2, mean),3)




### Stochastic volatility

SVcommon.Gibbs.iteration = function(aux, priors){
  # A single iteration of the Gibbs sampler for the SV component
  #
  # aux is a list containing:
  #   Y - a TxN matrix
  #   X - a TxK matrix
  #   H - a Tx1 matrix
  #   h0 - a scalar
  #   sigma.v2 - a scalar
  #   s - a Tx1 matrix
  #   A - a KxN matrix
  #   Sigma - an NxN matrix
  #   sigma2 - a Tx1 matrix
  #
  # priors is a list containing:
  #   h0.v - a positive scalar
  #   h0.m - a scalar
  #   sigmav.s - a positive scalar
  #   sigmav.nu - a positive scalar
  #   HH - a TxT matrix
  
  T             = dim(aux$Y)[1]
  N             = dim(aux$Y)[2]
  alpha.st      = c(1.92677,1.34744,0.73504,0.02266,0-0.85173,-1.97278,-3.46788,-5.55246,-8.68384,-14.65000)
  sigma.st      = c(0.11265,0.17788,0.26768,0.40611,0.62699,0.98583,1.57469,2.54498,4.16591,7.33342)
  pi.st         = c(0.00609,0.04775,0.13057,0.20674,0.22715,0.18842,0.12047,0.05591,0.01575,0.00115)
  
  Lambda        = solve(chol(aux$Sigma))
  Z             = rowSums( ( aux$Y - aux$X %*% aux$A ) %*% Lambda ) / sqrt(N)
  Y.tilde       = as.vector(log((Z + 0.0000001)^2))
  Ytilde.alpha  = as.matrix(Y.tilde - alpha.st[as.vector(aux$s)])
  
  # sampling initial condition
  ############################################################
  V.h0.bar      = 1/((1 / priors$h0.v) + (1 / aux$sigma.v2))
  m.h0.bar      = V.h0.bar*((priors$h0.m / priors$h0.v) + (aux$H[1] / aux$sigma.v2))
  h0.draw       = rnorm(1, mean = m.h0.bar, sd = sqrt(V.h0.bar))
  aux$h0        = h0.draw
  
  # sampling sigma.v2
  ############################################################
  sigma.v2.s    = priors$sigmav.s + sum(c(aux$H[1] - aux$h0, diff(aux$H))^2)
  sigma.v2.draw = sigma.v2.s / rchisq(1, priors$sigmav.nu + T)
  aux$sigma.v2  = sigma.v2.draw
  
  # sampling auxiliary states
  ############################################################
  Pr.tmp        = simplify2array(lapply(1:10,function(x){
    dnorm(Y.tilde, mean = as.vector(aux$H + alpha.st[x]), sd = sqrt(sigma.st[x]), log = TRUE) + log(pi.st[x])
  }))
  Pr            = t(apply(Pr.tmp, 1, function(x){exp(x - max(x)) / sum(exp(x - max(x)))}))
  s.cum         = t(apply(Pr, 1, cumsum))
  r             = matrix(rep(runif(T), 10), ncol = 10)
  ss            = apply(s.cum < r, 1, sum) + 1
  aux$s         = as.matrix(ss)
  
  
  # sampling log-volatilities using functions for tridiagonal precision matrix
  ############################################################
  Sigma.s.inv   = diag(1 / sigma.st[as.vector(aux$s)])
  D.inv         = Sigma.s.inv + (1 / aux$sigma.v2) * priors$HH
  b             = as.matrix(Ytilde.alpha / sigma.st[as.vector(aux$s)] + (aux$h0/aux$sigma.v2)*diag(T)[,1])
  lead.diag     = diag(D.inv)
  sub.diag      = mgcv::sdiag(D.inv, -1)
  D.chol        = mgcv::trichol(ld = lead.diag, sd = sub.diag)
  D.L           = diag(D.chol$ld)
  mgcv::sdiag(D.L,-1) = D.chol$sd
  x             = as.matrix(rnorm(T))
  a             = forwardsolve(D.L, b)
  draw          = backsolve(t(D.L), a + x)
  aux$H         = as.matrix(draw)
  aux$sigma2    = as.matrix(exp(draw))
  
  return(aux)
}


# Setting specifications
N = ncol(Data[ , -1])
p = 12
K = 1+N*p
S = 10000
h = 24
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -1], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
T       = nrow(Y)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}


# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)
H                 = diag(T)
sdiag(H,-1)       = -1
HH                = 2*diag(T)
sdiag(HH,-1)      = -1
sdiag(HH,1)       = -1

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1,
  
  # New priors based on lectures
  h0.v        = 1,
  h0.m        = 0,
  sigmav.s    = 1,
  sigmav.nu   = 1, 
  HH          = HH 
)

# BVAR function

BVAR.SV = function(Y,X,priors,S){

  aux <- list(
    Y = Y, 
    X = X,  
    H = matrix(1,T,1), 
    h0 = 0, 
    sigma.v2 = 1,
    s = matrix(1,T,1),
    A = matrix(0, K, N), 
    Sigma = diag(diag(matrix(1, N, N))),
    sigma2 = matrix(1, T, 1) 
  )
  
  A.posterior        = array(NA, dim = c(K,N,sum(S)))
  Sigma.posterior    = array(NA,dim=c(N,N,sum(S)))
  sigma2.posterior    = matrix(NA, nrow(Y), sum(S)) 

  # Create progress bar
  pb <- txtProgressBar(min = 0, max = sum(S), style = 3)
    
  for (s in 1:sum(S)){
    # normal-inverse Wishart posterior parameters
    V.bar.inv   = t(X)%*%diag(1/as.vector(aux$sigma2))%*%X + diag(1/diag(priors$V.prior))
    V.bar       = solve(V.bar.inv)
    A.bar       = V.bar%*%(t(X)%*%diag(1/as.vector(aux$sigma2))%*%Y + diag(1/diag(priors$V.prior))%*%priors$A.prior)
    nu.bar      = nrow(Y) + priors$nu.prior
    S.bar       = priors$S.prior + t(Y)%*%diag(1/as.vector(aux$sigma2))%*%Y + t(priors$A.prior)%*%diag(1/diag(priors$V.prior))%*%priors$A.prior - t(A.bar)%*%V.bar.inv%*%A.bar
    S.bar.inv   = solve(S.bar)
    
    #posterior draws
    Sigma.posterior.dist   = rWishart(1, df=nu.bar, Sigma=S.bar.inv)
    Sigma.draw             = apply(Sigma.posterior.dist,3,solve)
    Sigma.posterior[,,s]   = Sigma.draw
    A.posterior[,,s]       = array(rnorm(prod(c(dim(A.bar),1))),c(dim(A.bar),1))
    L                      = t(chol(V.bar))
    A.posterior[,,s]       = A.bar + L%*%A.posterior[,,s]%*%chol(Sigma.posterior[,,s])

    
    # Draw using stochastic volatility Gibbs common sampler
    aux                   = SVcommon.Gibbs.iteration(aux, priors)
    sigma2.posterior[,s]  = aux$sigma2
    
    # Update progress bar
    setTxtProgressBar(pb, s)
  }
  # Close progress bar
  close(pb)
  
  posterior = list(
    Sigma.posterior   = Sigma.posterior,
    A.posterior       = A.posterior,
    Sigma2.posterior  = sigma2.posterior
  )
  return(posterior)
}

# Applying BVAR function
posterior.draws.SV = BVAR.SV(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws.SV$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws.SV$A.posterior, 1:2, mean),3)
round(apply(posterior.draws.SV$Sigma2.posterior, 1, mean),3)


### Extended BVAR with stochastic volatility

BVAR.SV.extension = function(Y,X,priors,S){
  
  aux <- list(
    Y = Y, 
    X = X,  
    H = matrix(1,T,1), 
    h0 = 0, 
    sigma.v2 = 1,
    s = matrix(1,T,1),
    A = matrix(0, K, N), 
    Sigma = diag(diag(matrix(1, N, N))),
    sigma2 = matrix(1, T, 1) 
  )
  
  A.posterior        = array(0.1, dim = c(K,N,sum(S)))
  Sigma.posterior    = array(0.1,dim=c(N,N,sum(S)))
  sigma2.posterior    = matrix(NA, nrow(Y), sum(S)) 
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = sum(S), style = 3)
  
  for (s in 1:sum(S)){
    sbar2                = compute_sbar2(s2, Sigma.posterior[,,s],V.bar)
    a_posterior          = compute_a_posterior(sbar2, Sigma.posterior[,,s],V.bar, A.posterior[,,s])
    a.posterior          = rnorm(1, a_posterior, sqrt(sbar2))
    a.posterior.store[s] = a.posterior
    
    # normal-inverse Wishart posterior parameters
    V.bar.inv   = t(X)%*%diag(1/as.vector(aux$sigma2))%*%X + diag(1/diag(priors$V.prior))
    V.bar       = solve(V.bar.inv)
    A.bar       = V.bar%*%(t(X)%*%diag(1/as.vector(aux$sigma2))%*%Y + diag(1/diag(priors$V.prior))%*%priors$A.prior*a.posterior)
    nu.bar      = nrow(Y) + priors$nu.prior
    S.bar       = priors$S.prior + t(Y)%*%diag(1/as.vector(aux$sigma2))%*%Y + t(priors$A.prior*a.posterior)%*%diag(1/diag(priors$V.prior))%*%priors$A.prior*a.posterior - t(A.bar)%*%V.bar.inv%*%A.bar
    S.bar.inv   = solve(S.bar)
    
    #posterior draws
    Sigma.posterior.dist   = rWishart(1, df=nu.bar, Sigma=S.bar.inv)
    Sigma.draw             = apply(Sigma.posterior.dist,3,solve)
    Sigma.posterior[,,s]   = Sigma.draw
    A.posterior[,,s]       = array(rnorm(prod(c(dim(A.bar),1))),c(dim(A.bar),1))
    L                      = t(chol(V.bar))
    A.posterior[,,s]       = A.bar + L%*%A.posterior[,,s]%*%chol(Sigma.posterior[,,s])
    
    
    # Draw using stochastic volatility Gibbs common sampler
    aux                   = SVcommon.Gibbs.iteration(aux, priors)
    sigma2.posterior[,s]  = aux$sigma2
    
    # Update progress bar
    setTxtProgressBar(pb, s)
  }
  # Close progress bar
  close(pb)
  
  posterior = list(
    Sigma.posterior   = Sigma.posterior,
    A.posterior       = A.posterior,
    Sigma2.posterior  = sigma2.posterior
  )
  return(posterior)
}

# Applying BVAR function
posterior.draws.SV.extension = BVAR.SV.extension(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws.SV.extension$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws.SV.extension$A.posterior, 1:2, mean),3)
round(apply(posterior.draws.SV.extension$Sigma2.posterior, 1, mean),3)








### Baseline model forecast

# simulate draws from the predictive density
h = 36 # 3 year ahead forecast
S = 10000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws$A.posterior[,,s], sigma=posterior.draws$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)

### Extended model forecast

# simulate draws from the predictive density
h = 24 # 3 year ahead forecast
S = 1000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.extension$A.posterior[,,s], sigma=posterior.draws.extension$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)



### SV model forecast
# simulate draws from the predictive density
h = 36 # 3 year ahead forecast
S = 10000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.SV$A.posterior[,,s], sigma=posterior.draws.SV$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }

}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)

### SV.extension model forecast
# simulate draws from the predictive density
h = 24 # 3 year ahead forecast
S = 10000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.SV.extension$A.posterior[,,s], sigma=posterior.draws.SV.extension$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
  
}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)



### No Financial Variable models

## No financial variable baseline model
# Setting specifications
N = ncol(Data[ , -c(1,3)])
p = 12
K = 1+N*p
S = 1000
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -c(1,3)], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}

# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1 
)

# Applying BVAR function
posterior.draws.extension.no_finance = BVAR.extension(Y=Y, X=X, priors=priors, S=S)

# simulate draws from the predictive density
h = 24 # 3 year ahead forecast
S = 1000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.extension.no_finance$A.posterior[,,s], sigma=posterior.draws.extension.no_finance$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
}

gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)

## No financial variable extended model
# Setting specifications
N = ncol(Data[ , -c(1,3)])
p = 12
K = 1+N*p
S = 1000
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -c(1,3)], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}

# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1 
)

# Applying BVAR function
posterior.draws.extension.no_finance = BVAR.extension(Y=Y, X=X, priors=priors, S=S)

# simulate draws from the predictive density
h = 24 # 3 year ahead forecast
S = 1000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.extension.no_finance$A.posterior[,,s], sigma=posterior.draws.extension.no_finance$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
}

gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)





## No financial variable SV model
# Setting specifications
N = ncol(Data[ , -c(1,3)])
p = 12
K = 1+N*p
S = 10000
set.seed(1)

# Initializing X and Y matrices
y       = ts(Data[ , -c(1,3)], start=c(2003,1), frequency=12)
Y       = ts(y[13:nrow(y),], start=c(2004,1), frequency=12)
T       = nrow(Y)
X       = matrix(1,nrow(Y),1)
for (i in 1:p){
  X     = cbind(X,y[13:nrow(y)-i,])
}

# Maximum Likelihood Estimator
A.hat       = solve(t(X)%*%X)%*%t(X)%*%Y
Sigma.hat   = t(Y-X%*%A.hat)%*%(Y-X%*%A.hat)/T

# Setting Minnesota Prior
kappa.1           = 0.02^2
kappa.2           = 100
A.prior           = matrix(0,nrow(A.hat),ncol(A.hat))
A.prior[2:(N+1),] = diag(N)
H                 = diag(T)
sdiag(H,-1)       = -1
HH                = 2*diag(T)
sdiag(HH,-1)      = -1
sdiag(HH,1)       = -1

priors = list(
  A.prior     = A.prior,
  V.prior     = diag(c(kappa.2,kappa.1*((1:p)^(-2))%x%rep(1,N))),
  S.prior     = diag(diag(Sigma.hat)),
  nu.prior    = N+1,
  
  # New priors based on lectures
  h0.v        = 1,
  h0.m        = 0,
  sigmav.s    = 1,
  sigmav.nu   = 1, 
  HH          = HH 
)


# Applying BVAR function
posterior.draws.SV.no_finance = BVAR.SV(Y=Y, X=X, priors=priors, S=S)

h = 24 # 3 year ahead forecast
S = 10000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.SV.no_finance$A.posterior[,,s], sigma=posterior.draws.SV.no_finance$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
  
}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)

# Applying BVAR function
posterior.draws.SV.extension.no_finance = BVAR.SV(Y=Y, X=X, priors=priors, S=S)

h = 24 # 3 year ahead forecast
S = 10000
Y.h         = array(NA,c(h,N,S))

# sampling predictive density
for (s in 1:S){
  x.Ti        = Y[(nrow(Y)-h+1):nrow(Y),]
  x.Ti        = x.Ti[p:1,]
  for (i in 1:h){
    x.T         = c(1,as.vector(t(x.Ti)))
    Y.h[i,,s]   = rmvnorm(1, mean = x.T%*%posterior.draws.SV.no_finance$A.posterior[,,s], sigma=posterior.draws.SV.no_finance$Sigma.posterior[,,s])
    x.Ti        = rbind(Y.h[i,,s],x.Ti[1:(p-1),])
  }
  
}


gdp.point.f    = apply(Y.h[,1,],1,mean) 
gdp.interval.f = apply(Y.h[,1,],1,hdi,credMass=0.68)
gdp.range      = range(y[,1],gdp.interval.f)

blue      = "#05386B"
green     = "#379683"
green.rgb = col2rgb(green)
shade     = rgb(green.rgb[1],green.rgb[2],green.rgb[3], alpha=120, maxColorValue=255)


par(mfrow=c(1,1), mar=rep(3,4),cex.axis=1.5)
plot(1:(length(y[,1])+h),c(y[,1],gdp.point.f), type="l", ylim=gdp.range, axes=FALSE, xlab="", ylab="", lwd=2, col=green)
axis(1,c(1,61,205, nrow(y),nrow(y)+h),c("2003-01","2008-01","2020-01","2023-12",""), col=blue)
axis(2,c(gdp.range[1],mean(gdp.range),gdp.range[2]),c("","GDP",""), col=blue)
abline(v=nrow(y), col="black")
text(x=253, y=4.8, srt=90, "2024-01")
abline(v=nrow(y)+12, col="black")
text(x=265, y=4.8, srt=90, "2025-01")
abline(v=nrow(y)+24, col="black")
text(x=277, y=4.8, srt=90, "2026-01")
polygon(c(length(y[,1]):(length(y[,1])+h),(length(y[,1]):(length(y[,1])+h))[25:1]),
        c(y[230,1],gdp.interval.f[1,],gdp.interval.f[2,24:1],y[230,1]),
        col=shade, border=blue)

