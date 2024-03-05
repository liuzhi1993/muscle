## code to prepare `q_30000` dataset goes here
dir = getwd()
setwd(dir)
q_30000 = readRDS("q_MUSCLE_30000.rds")
usethis::use_data(q_30000, overwrite = TRUE)
