library(devtools)
install_github("liuzhi1993/muscle")

library(muscle)
packageVersion("muscle")
library(tictoc)

# generate data
n = 2048
alpha = 0.1
df = 3
#change points of the block signal
blocks <- rep(c(0, 14.64, -3.66, 7.32, -7.32, 10.98, -4.39, 3.29, 19.03, 7.68,
                15.37, 0), times =c(204, 62, 41, 164, 40, 308, 82, 430, 225,
                41, 61,390))
blocks.cpt <- c(205, 267, 308, 472, 512, 820, 902,
                1332, 1557, 1598, 1659)

signal = blocks
signal.cpt = sort(c(blocks.cpt,390,666,1445))
sd = c(8,0.5,4,1)
set.seed(96)
noise = c(rt(390,df)*sd[1],rt(278,df)*sd[2],rt(779,df)*sd[3],rt(601,df)*sd[4])/sqrt(df/(df-2))
Y = signal + noise

## MUSCLE beta = 0.5, PRT
beta = 0.5
tic()
q_muscle = simulQuantile_MUSCLE(n, alpha = alpha, beta = beta)
toc()

tic()
reg_muscle_prt_05 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PRT")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.5, PRT")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_prt_05),type = "s",lwd = 2, col = "red")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )


## MUSCLE beta = 0.5, PST
tic()
reg_muscle_pst_05 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PST")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.5, PST")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_pst_05),type = "s",lwd = 2, col = "blue")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.5, original
tic()
reg_muscle_original_05 = MUSCLE(Y, q_muscle, beta = beta, implementation = "original")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.5, original")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_original_05),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.25, PRT
beta = 0.25
tic()
q_muscle = simulQuantile_MUSCLE(n, alpha = alpha, beta = beta)
toc()

tic()
reg_muscle_prt_025 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PRT")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.25, PRT")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_prt_025),type = "s",lwd = 2, col = "red")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.25, PST
tic()
reg_muscle_pst_025 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PST")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.25, PST")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_pst_025),type = "s",lwd = 2, col = "blue")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.25, original
tic()
reg_muscle_original_025 = MUSCLE(Y, q_muscle, beta = beta, implementation = "original")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.25, original")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_original_025),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.75, PRT
beta = 0.75
tic()
q_muscle = simulQuantile_MUSCLE(n, alpha = alpha, beta = beta)
toc()

tic()
reg_muscle_prt_075 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PRT")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.75, PRT")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_prt_075),type = "s",lwd = 2, col = "red")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.75, PST
tic()
reg_muscle_pst_075 = MUSCLE(Y, q_muscle, beta = beta, implementation = "PST")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.75, PST")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_pst_075),type = "s",lwd = 2, col = "blue")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE beta = 0.75, original
tic()
reg_muscle_original_075 = MUSCLE(Y, q_muscle, beta = beta, implementation = "original")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE beta=0.75, original")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle_original_075),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## M-MUSCLE with PST
beta_vec = c(0.25,0.5,0.75)
tic()
q_mmuscle = simulQuantile_MMUSCLE(n, alpha = alpha, beta_vec = beta_vec)
toc()

tic()
reg_mmuscle_pst = MMUSCLE(Y, q_mmuscle, beta_vec = beta_vec, implementation = "PST")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "M-MUSCLE, PST")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(eval_StepFun_M(reg_mmuscle_pst,2),type = "s",lwd = 2, col = "red")
lines(eval_StepFun_M(reg_mmuscle_pst,1),type = "s",lwd = 2, col = "blue")
lines(eval_StepFun_M(reg_mmuscle_pst,3),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## M-MUSCLE with PRT
tic()
reg_mmuscle_prt = MMUSCLE(Y, q_mmuscle, beta_vec = beta_vec, implementation = "PRT")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "M-MUSCLE, PRT")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(eval_StepFun_M(reg_mmuscle_prt,2),type = "s",lwd = 2, col = "red")
lines(eval_StepFun_M(reg_mmuscle_prt,1),type = "s",lwd = 2, col = "blue")
lines(eval_StepFun_M(reg_mmuscle_prt,3),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## M-MUSCLE with original
tic()
reg_mmuscle_original = MMUSCLE(Y, q_mmuscle, beta_vec = beta_vec, implementation = "original")
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "M-MUSCLE, original")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(eval_StepFun_M(reg_mmuscle_original,2),type = "s",lwd = 2, col = "red")
lines(eval_StepFun_M(reg_mmuscle_original,1),type = "s",lwd = 2, col = "blue")
lines(eval_StepFun_M(reg_mmuscle_original,3),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )
