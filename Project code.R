library(dplyr)
library(tidyr)
library(ggplot2)
library(ecb)
library(gridExtra)
library(tseries)
library(bsvars)

rm(list = ls())

# Downloading Data
start_date_m <- "2003-01" 
end_date_m <- "2023-12"


HICP <- get_data("ICP.M.U2.Y.000000.3.INX", filter = list(startPeriod = start_date_m, endPeriod = end_date_m)) |>
  transmute(obstime, hicp = obsvalue)

CISS <- get_data("CISS.M.U2.Z0Z.4F.EC.SOV_EW.IDX", filter = list(startPeriod = start_date_m, endPeriod = end_date_m)) |>
  transmute(obstime, ciss = obsvalue)

IP <- get_data("STS.M.I8.Y.PROD.NS0010.4.000", filter = list(startPeriod = start_date_m, endPeriod = end_date_m)) |>
  transmute(obstime, ip = obsvalue)

UN <- get_data("LFSI.M.I9.S.UNEHRT.TOTAL0.15_74.T", filter = list(startPeriod = start_date_m, endPeriod = end_date_m)) |>
  transmute(obstime, un = obsvalue)

Data <- HICP |>
  left_join(IP, by = "obstime") |>
  left_join(CISS, by = "obstime", suffix = c(".HICP", ".CISS")) |>
  left_join(UN, by = "obstime", suffix = c("", ".UN"))

Data <- Data |>
  mutate(across(-c(1, 4), ~log(.)))

Data <- Data |> 
  mutate(obstime = as.Date(paste0(obstime, "-01")))


# Plotting the variables
p1 <- ggplot(Data, aes(x = obstime, y = hicp)) +
  geom_line() +
  labs(x = "", y = "HICP", title = "Inflation") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )

p1 <- ggplot(Data, aes(x = obstime, y = hicp)) +
  geom_line() +
  labs(x = "", y = "hicp", title = "Inflation") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )

p2 <- ggplot(Data, aes(x = obstime, y = ip)) +
  geom_line() +
  labs(x = "", y = "ip", title = "Industrial production") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )

p3 <- ggplot(Data, aes(x = obstime, y = ciss)) +
  geom_line() +
  labs(x = "", y = "ciss", title = "Financial Stress") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )

p4 <- ggplot(Data, aes(x = obstime, y = un)) +
  geom_line() +
  labs(x = "", y = "un", title = "Unemployment") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.background = element_blank(),  
    plot.title = element_text(hjust = 0.5)  
  )


grid.arrange(p1, p2, p3, p4, ncol = 2)

# ACF + PACF
Data <- as.ts(Data)

par(mfrow = c(2, 2))

acf(Data[,1], main = "ACF for HICP")
acf(Data[,2], main = "ACF for IP")
acf(Data[,3], main = "ACF for CISS")
acf(Data[,4], main = "ACF for UN")

par(mfrow = c(2, 2))

pacf(Data[,1], main = "ACF for HICP")
pacf(Data[,2], main = "ACF for IP")
pacf(Data[,3], main = "ACF for CISS")
pacf(Data[,4], main = "ACF for UN")


# ADF Test

adf_hicp <- adf.test(Data[,1], k=12, alternative = "stationary")
adf_ip <- adf.test(Data[,2], k=12, alternative = "stationary")
adf_ciss <- adf.test(Data[,3], k=12, alternative = "stationary")
adf_un <- adf.test(Data[,4], k=12, alternative = "stationary")

adf_results <- data.frame(
  Variable = c("HICP", "IP", "CISS", "UN"),
  ADF_Statistic = c(adf_hicp$statistic, adf_ip$statistic, adf_ciss$statistic, adf_un$statistic),
  P_Value = c(adf_hicp$p.value, adf_ip$p.value, adf_ciss$p.value, adf_un$p.value),
  Test_Critical_Values = I(list(adf_hicp$cval, adf_ip$cval, adf_ciss$cval, adf_un$cval))
)

Data_diff <- as.data.frame(lapply(Data, diff))

# Run ADF test on the first differences
adf_hicp_diff <- adf.test(Data_diff[,1], k=12, alternative = "stationary")
adf_ip_diff <- adf.test(Data_diff[,2], k=12, alternative = "stationary")
adf_ciss_diff <- adf.test(Data_diff[,3], k=12, alternative = "stationary")
adf_un_diff <- adf.test(Data_diff[,4], k=12, alternative = "stationary")

# Collect the ADF test results for the first differences in a data frame
adf_results_diff <- data.frame(
  Variable = c("HICP", "IP", "CISS", "UN"),
  ADF_Statistic = c(adf_hicp_diff$statistic, adf_ip_diff$statistic, adf_ciss_diff$statistic, adf_un_diff$statistic),
  P_Value = c(adf_hicp_diff$p.value, adf_ip_diff$p.value, adf_ciss_diff$p.value, adf_un_diff$p.value),
  Test_Critical_Values = I(list(adf_hicp_diff$cval, adf_ip_diff$cval, adf_ciss_diff$cval, adf_un_diff$cval))
)

print(adf_results_diff)
