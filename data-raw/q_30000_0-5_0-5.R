## code to prepare `q_30000_0-5_0-5` dataset goes here
dir = getwd()
setwd(dir)
load("q_30000_0.5_0.5_.Rdata")
q_30000_0_5_0_5 = q_muscle
usethis::use_data(q_30000_0_5_0_5, overwrite = TRUE)
