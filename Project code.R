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

rm(list = ls())

# Downloading Data
start_date_m <- "2003-01" 
end_date_m <- "2023-12"

# Define the datasets
series_info <- list(
  hicp = "ICP.M.U2.Y.000000.3.INX",
  ciss = "CISS.M.U2.Z0Z.4F.EC.SOV_EW.IDX",
  ip = "STS.M.I8.Y.PROD.NS0010.4.000",
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
  select(obstime, obsvalue)

# Transforming to monthly data
mro <- mro |>
  group_by(obstime = as.yearmon(obstime)) |>
  summarize(mro = last(obsvalue)) |>
  mutate(obstime = as.Date(obstime))

# Collecting quarterly data
consumption <- get_data("MNA.Q.Y.I9.W0.S1M.S1.D.P31._Z._Z._T.EUR.V.N") |>
  select(obstime, obsvalue) |>
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
  list(name = "hicp", title = "Inflation", ylab = "LOG(HICP)"),
  list(name = "ip", title = "Industrial production", ylab = "LOG(IP)"),
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

variable_names <- c("hicp", "ip", "ciss", "un", "m2", "mro", "consumption")

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
  Variable = c("HICP", "IP", "CISS", "UN", "M2", "MRO", "Consumption"),
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
S = 1000
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

### Showing baseline model
# Specifications
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
posterior.draws = BVAR(Y=Y, X=X, priors=priors, S=S)
round(apply(posterior.draws$Sigma.posterior, 1:2, mean),3)
round(apply(posterior.draws$A.posterior, 1:2, mean),3)




### Extended model BVAR estimation: Not filled out as I really do not know how to proceed.



