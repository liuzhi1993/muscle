## code to prepare `q_5000` dataset goes here
dir = getwd()
setwd(dir)
q_5000 = readRDS("q_MUSCLE_5000.rds")
usethis::use_data(q_5000, overwrite = TRUE)
