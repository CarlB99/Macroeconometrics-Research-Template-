library(dplyr)
library(tidyr)
library(ggplot2)
library(ecb)
library(gridExtra)
library(bsvars)
library(lubridate)
library(tseries)

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
  m1 = "BSI.M.U2.Y.V.M10.X.1.U2.2300.Z01.E"
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

# Merge all datasets by 'obstime'
Data <- Reduce(function(x, y) left_join(x, y, by = "obstime"), datasets)

Data <- Data |>
  mutate(across(c(2, 4, 6), ~log(.)))

# Convert 'obstime' to Date format
Data$obstime <- as.Date(paste0(Data$obstime, "-01"))

# Define the plot configurations
plot_settings <- list(
  list(name = "hicp", title = "Inflation", ylab = "LOG(HICP)"),
  list(name = "ip", title = "Industrial production", ylab = "LOG(IP)"),
  list(name = "ciss", title = "Financial Stress", ylab = "CISS"),
  list(name = "un", title = "Unemployment", ylab = "UN"),
  list(name = "m1", title = "M1 - Money Stock", ylab = "M1")
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
par(mfrow = c(2, 3))

variable_names <- c("hicp", "ip", "ciss", "un", "m1")

# Loop through the columns and create ACF plots
for (i in seq_along(variable_names)) {
  # Compute and plot ACF
  acf(Data[, i], main = paste("ACF for", variable_names[i]))
}

par(mfrow = c(2, 3))
for (i in seq_along(variable_names)) {
  # Compute and plot PACF
  pacf(Data[, i], main = paste("PACF for", variable_names[i]))
}


for (i in seq_along(variable_names)) {
  variable_name <- paste("adf", variable_names[i], sep = "_")
  assign(variable_name, adf.test(Data[, i], k=12, alternative = "stationary"))
}

adf_results <- data.frame(
  Variable = c("HICP", "IP", "CISS", "UN", "M1"),
  ADF_Statistic = c(adf_hicp$statistic, adf_ip$statistic, adf_ciss$statistic, adf_un$statistic, adf_m1$statistic),
  P_Value = c(adf_hicp$p.value, adf_ip$p.value, adf_ciss$p.value, adf_un$p.value, adf_m1$p.value),
  lags = rep(12, 5)
)

knitr::kable(adf_results, digits = 3, align = 'c')






kpss_hicp <- kpss.test(Data[,1], null = "Trend")
kpss_ip <- kpss.test(Data[,2], null = "Trend")
kpss_ciss <- kpss.test(Data[,3], null = "Trend")
kpss_un <- kpss.test(Data[,4], null = "Trend")
kpss_m1 <- kpss.test(Data[,5], null = "Trend")

kpss_results <- data.frame(
  Variable = c("HICP", "IP", "CISS", "UN", "M1"),
  KPSS_Statistic = c(kpss_hicp$statistic, kpss_ip$statistic, kpss_ciss$statistic, kpss_un$statistic, kpss_m1$statistic),
  P_Value = c(kpss_hicp$p.value, kpss_ip$p.value, kpss_ciss$p.value, kpss_un$p.value, kpss_m1$p.value)
)

knitr::kable(kpss_results, digits = 3, align = 'c')
