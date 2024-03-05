## code to prepare `q_20000` dataset goes here
dir = getwd()
setwd(dir)
q_20000 = readRDS("q_MUSCLE_20000.rds")
usethis::use_data(q_20000, overwrite = TRUE)
