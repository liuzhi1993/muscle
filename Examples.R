library(muscle)
library(tictoc)

# generate data
n = 2048
alpha = 0.1
beta = 0.5
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

## MUSCLE
# simulate quantiles
tic()
q_muscle = simulQuantile_MUSCLE(n,alpha = alpha,beta = beta)
toc()

# segmentation with MUSCLE
tic()
reg_muscle = MUSCLE(Y,q_muscle,beta)
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(evalStepFun(reg_muscle),type = "s",lwd = 2, col = "red")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )

## MUSCLE-S
# segmentation with MUSCLE-S with m = 500
tic()
reg_muscle_s = MUSCLE(Y,q_muscle,beta,split = TRUE, m = 500)
toc()

# segmentation with MUSCLE-S with m = 500, print computation details
tic()
reg_muscle_s = MUSCLE(Y,q_muscle,beta,split = TRUE, m = 500, details = TRUE)
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "MUSCLE-S")
lines(1:n, signal, type = "s", lwd =2, col ="black")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )
lines(evalStepFun(reg_muscle_s),type = "s",lwd = 2, col = "blue")

## M-MUSCLE
beta_vec = c(0.25,0.5,0.75)
# simulate quantiles
q_mmuscle = simulQuantile_MMUSCLE(n,alpha = alpha, beta_vec = beta_vec)

# segmentation with M-MUSCLE
# this can be slow!
tic()
reg_mmuscle = MMUSCLE(Y,q_mmuscle,beta_vec)
toc()

plot(1:n,Y, type = "l", lwd = 2, col = "gray",ylim = c(-40,35)
     ,xlab = "",ylab = "", main = "M-MUSCLE")
lines(1:n, signal, type = "s", lwd =2, col ="black")
lines(eval_StepFun_M(reg_mmuscle,2),type = "s",lwd = 2, col = "red")
lines(eval_StepFun_M(reg_mmuscle,1),type = "s",lwd = 2, col = "orange")
lines(eval_StepFun_M(reg_mmuscle,3),type = "s",lwd = 2, col = "orange")
abline(v = c(390,666,1445), lty = 2, col= "green", lwd = 2 )
