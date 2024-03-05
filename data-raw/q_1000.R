## code to prepare `q_1000` dataset goes here
dir = getwd()
setwd(dir)
q_1000 = readRDS("q_MUSCLE_1000.rds")
usethis::use_data(q_1000, overwrite = TRUE)
