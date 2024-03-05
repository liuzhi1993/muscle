## code to prepare `q_10000` dataset goes here
dir = getwd()
setwd(dir)
q_10000 = readRDS("q_MUSCLE_10000.rds")
usethis::use_data(q_10000, overwrite = TRUE)
