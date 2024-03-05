## usethis namespace: start
#' @useDynLib MUSCLE, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' MUSCLE
#' @param Y Raw data
#' @param q quantiles, can be simulated by simulQuantile_MUSCLE
#' @param b target quantile, b = 0.5 for median
#' @param dyadic test on dyadic subintervals or full subintervals
#' @param split boolian variable indicates whether the data will be split by subsets with size m
#' @param m splitting size
#' @param deconv boolian variable, indicates whether deconvolution or not
#' @param lag shift lag, only for deconvolution
#' @return idealization of MUSCLE
#' @export
#'
MUSCLE <- function(Y, q, b = 0.5, dyadic = TRUE, split = FALSE, m = 2000, deconv = FALSE, lag = 192){
  n = length(Y)
  if(deconv == FALSE){
    lag = 0
    if(split == TRUE){
      print(paste("There are ",ceiling(n/m), "pieces:"))
      if(length(q)>=min(n,m)){
        print(paste("The 1 . piece is computing."))
        split_res = Split(n,m)
        num_split = length(split_res$start)
        start = split_res$start[1]
        end = split_res$end[1]
        size = end - start + 1
        if(dyadic == TRUE){
          res = .MUSCLE(Y[start:end],q[1:size],b,FALSE)
        }else{
          res = .MUSCLE_FULL(Y[start:end],q[1:size],b,FALSE)
        }
        num_cp_res = length(res$value)
        if(num_split>=2){
          for (i in 2:num_split) {
            print(paste("The", i,". piece is computing."))
            start = split_res$start[i]
            end = split_res$end[i]
            size = end - start + 1
            if(dyadic == TRUE){
              temp_res = .MUSCLE(Y[start:end],q[1:size],b,FALSE)
            }else{
              temp_res = .MUSCLE_FULL(Y[start:end],q[1:size],b,FALSE)
            }
            num_cp_temp_res = length(temp_res$value)
            c1 = res$last[1]
            c2 = res$last[2]
            d1 = temp_res$first[1]
            d2 = temp_res$first[2]
            if(c2>=d1 && c1<= d2){
              n = res$n + temp_res$n
              first = res$first
              if(num_cp_res>=2 && num_cp_temp_res>=2){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),b),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res==1 && num_cp_temp_res>=2){
                value = c(quantile(c(max(c1,d1),min(c2,d2)),b),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res>=2 && num_cp_temp_res==1){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),b))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }else{
                value = c(quantile(c(max(c1,d1),min(c2,d2)),b))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }
            }else{
              print(paste("There is a small jump betwenn ",i-1, ". and ", i, ". segments, alternative solution is computing."))
              n = res$n + temp_res$n
              first = res$first
              temp_start = res$left[num_cp_res]
              if(num_cp_temp_res == 1){
                temp_end = res$n + temp_res$n
              }else{
                temp_end = res$n + temp_res$left[2]-1
              }
              temp_size = temp_end - temp_start + 1
              print(paste("Problem segment start = ", temp_start, " end = ", temp_end))
              print(paste("Problem segment has length ", temp_size))
              if(temp_size <= length(q)){
                print(paste("Quantiles are sufficiently long."))
                if(dyadic == TRUE){
                  joint_res = .MUSCLE(Y[temp_start:temp_end],q[1:temp_size],b, TRUE)
                }else{
                  joint_res = .MUSCLE_FULL(Y[temp_start:temp_end],q[1:temp_size],b,TRUE)
                }
                if(length(joint_res$left)==1){
                  print(paste("Alternative solution is computed and jump deleted!"))
                  if(num_cp_res>=2 && num_cp_temp_res>=2){
                    value = c(res$value[1:(num_cp_res-1)],joint_res$value,temp_res$value[2:num_cp_temp_res])
                    left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                    last = temp_res$last
                  }else if(num_cp_res==1 && num_cp_temp_res>=2){
                    value = c(joint_res$value,temp_res$value[2:num_cp_temp_res])
                    left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                    last = temp_res$last
                  }else if(num_cp_res>=2 && num_cp_temp_res==1){
                    value = c(res$value[1:(num_cp_res-1)],joint_res$value)
                    left = c(res$left)
                    last = joint_res$last
                  }else{
                    value = c(joint_res$value)
                    left = c(res$left)
                    last = (joint_res$last)
                  }
                }else{
                  value = c(res$value, temp_res$value)
                  left = c(res$left, temp_res$left+res$n)
                  n = res$n + temp_res$n
                  first = res$first
                  last = temp_res$last
                }
              }else{
                value = c(res$value, temp_res$value)
                left = c(res$left, temp_res$left+res$n)
                n = res$n + temp_res$n
                first = res$first
                last = temp_res$last
              }
            }
            res = list(value = value, left = left, n = n, first = first, last = last)
            num_cp_res = length(res$value)
          }
        }
        return(res)
      }else{
        error("Quantile size too small!")
      }
    }else{
      if(n > length(q)){
        error("Sample size and quantile size do not match!")
      }else{
        q = q[1:n]
        if(dyadic == TRUE){
          return(.MUSCLE(Y, q, b, FALSE))
        }else{
          return(.MUSCLE_FULL(Y, q, b,FALSE))
        }
      }
    }
  }else{
    if(split == TRUE){
      print(paste("There are ",ceiling(n/m), "pieces:"))
      if(length(q)>=min(n,2000)){
        print(paste("The 1 . piece is computing."))
        split_res = Split(n,m)
        num_split = length(split_res$start)
        start = split_res$start[1]
        end = split_res$end[1]
        size = end - start + 1
        if(dyadic == TRUE){
          res = .DMUSCLE(Y[start:end],q[1:size],b,lag,FALSE)
        }else{
          res = .DMUSCLE_FULL(Y[start:end],q[1:size],b,lag,FALSE)
        }
        num_cp_res = length(res$value)
        if(num_split>=2){
          for (i in 2:num_split) {
            print(paste("The", i, ". piece is computing."))
            start = split_res$start[i]
            end = split_res$end[i]
            size = end - start + 1
            if(dyadic == TRUE){
              temp_res = .DMUSCLE(Y[start:end],q[1:size],b,lag,FALSE)
            }else{
              temp_res = .DMUSCLE_FULL(Y[start:end],q[1:size],b,lag,FALSE)
            }
            num_cp_temp_res = length(temp_res$value)
            c1 = res$last[1]
            c2 = res$last[2]
            d1 = temp_res$first[1]
            d2 = temp_res$first[2]
            if(c2>=d1 && c1<= d2){
              n = res$n + temp_res$n
              first = res$first
              if(num_cp_res>=2 && num_cp_temp_res>=2){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),b),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res==1 && num_cp_temp_res>=2){
                value = c(quantile(c(max(c1,d1),min(c2,d2)),b),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res>=2 && num_cp_temp_res==1){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),b))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }else{
                value = c(quantile(c(max(c1,d1),min(c2,d2)),b))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }
            }else{
              print(paste("There is a small jump betwenn ",i-1, ". and ", i, ". segments, alternative solution is computing."))
              n = res$n + temp_res$n
              first = res$first
              temp_start = res$left[num_cp_res]
              if(num_cp_temp_res == 1){
                temp_end = res$n + temp_res$n
              }else{
                temp_end = res$n + temp_res$left[2]-1
              }
              temp_size = temp_end - temp_start + 1
              print(paste("Problem segment start = ", temp_start, " end = ", temp_end))
              print(paste("Problem segment has length ", temp_size))
              if(temp_size <= length(q)){
                print(paste("Quantiles are sufficiently long."))
                if(dyadic == TRUE){
                  joint_res = .DMUSCLE(Y[temp_start:temp_end],q[1:temp_size],b, lag,TRUE)
                }else{
                  joint_res = .DMUSCLE_FULL(Y[temp_start:temp_end],q[1:temp_size],b,lag,TRUE)
                }
                if(length(joint_res$left)==1){
                  print(paste("Alternative solution is computed and jump deleted!"))
                  if(num_cp_res>=2 && num_cp_temp_res>=2){
                    value = c(res$value[1:(num_cp_res-1)],joint_res$value,temp_res$value[2:num_cp_temp_res])
                    left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                    last = temp_res$last
                  }else if(num_cp_res==1 && num_cp_temp_res>=2){
                    value = c(joint_res$value,temp_res$value[2:num_cp_temp_res])
                    left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                    last = temp_res$last
                  }else if(num_cp_res>=2 && num_cp_temp_res==1){
                    value = c(res$value[1:(num_cp_res-1)],joint_res$value)
                    left = c(res$left)
                    last = joint_res$last
                  }else{
                    value = c(joint_res$value)
                    left = c(res$left)
                    last = (joint_res$last)
                  }
                }else{
                  value = c(res$value, temp_res$value)
                  left = c(res$left, temp_res$left+res$n)
                  n = res$n + temp_res$n
                  first = res$first
                  last = temp_res$last
                }
              }else{
                value = c(res$value, temp_res$value)
                left = c(res$left, temp_res$left+res$n)
                n = res$n + temp_res$n
                first = res$first
                last = temp_res$last
              }
            }
            res = list(value = value, left = left, n = n, first = first, last = last)
            num_cp_res = length(res$value)
          }
        }
        return(res)
      }else{
        error("Quantile size too small!")
      }
    }else{
      if(n > length(q)){
        error("Sample size and quantile size do not match!")
      }else{
        q = q[1:n]
        if(dyadic == TRUE){
          return(.DMUSCLE(Y, q, b,lag))
        }else{
          return(.DMUSCLE_FULL(Y, q, b,lag))
        }
      }
    }
  }

}

#' MMUSCLE
#' @param Y Raw data
#' @param q_matrix quantiles, can be simulated by simulQuantile_MUSCLE
#' @param b_vec target quantile, b = 0.5 for median
#' @return idealization of MUSCLE
#' @export
#'
MMUSCLE <- function(Y,q_matrix, b_vec = c(0.5)){
  return(.MMUSCLE(Y,q_matrix,b_vec,FALSE))
}

#' Simulate MUSCLE quantiles
#' @param alpha Type I error of individual test
#' @param n sample size
#' @param b target quantile, b = 0.5 for median
#' @param deconv logical variable, TRUE for deconvolution (DMUSCLE)
#' @param Y data
#' @param lambda regularization parameter
#' @param m lag
#' @param sr sample rate
#' @param Kernel convolution kernel
#' @param E16 logical variable, TRUE for Element-16 measurement
#' @return simulated quantiles
#' @export
simulQuantile_MUSCLE = function(n, alpha = 0.1, b = 0.5, lambda = 1,deconv = FALSE,Kernel = NULL,Y = NULL, E16 = FALSE,sr = 20000,m=192){
  if(deconv == FALSE){
    if(alpha == 0.1 && b == 0.5){
      if(n == 1e3){
        return(q_1000)
      }else if(n == 2e3){
        return(q_2000)
      }else if(n == 5e3){
        return(q_5000)
      }else if(n == 1e4){
        return(q_10000)
      }else if(n == 1.5e4){
        return(q_15000)
      }else if(n == 2e4){
        return(q_20000)
      }else if(n == 2.5e4){
        return(q_25000)
      }else if(n == 3e4){
        return(q_30000)
      }else{
        return(.simulQuantile_MUSCLE(1-alpha, n, b))
      }
    }else{
      return(.simulQuantile_MUSCLE(1-alpha, n, b))
    }
  }else{
    res = rep(0,n)
    if(!is.null(Kernel)){
      m = length(Kernel)
      Kernelfft = fft(c(Kernel,rep(0,n-1)))/(n+m-1)
      r = 1000
      sd = 1 / sqrt(sum(Kernel^2))
      data = matrix(0,r,n)
      ACF = convolve(Kernel,Kernel, type = "open")[m:(2*m-1)]
      ACF = ACF/max(ACF)
      for (i in 1:r) {
        X = rnorm(n+m-1,0,sd)
        X_tilde = Re(fft(fft(X)*Kernelfft,inverse = TRUE))[m:(n+m-1)]
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n, b)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      return(res)
    }else if(E16 == TRUE){
      if(sr == 1250){
        m = 1024
      }else if(sr == 5000){
        m = 256
      }else if(sr == 10000){
        m = 128
      }else if(sr == 20000){
        m = 64
      }else if(sr == 50000){
        m = 200
      }else if(sr == 100000){
        m = 100
      }else if(sr == 200000){
        m = 50
      }else {
        m = 10
      }
      #m = 64
      #t = m/sr
      #x = (0:(3*m))/sr
      #Kernel = (x>=0)*(x<=t)*1/2*x^2+(x>t)*(x<2*t)*(-x^2+3*t*x-3/2*t^2)+(x>=2*t)*(x<=3*t)*(1/2*x^2-3*t*x+9/2*t^2)
      #Kernel = 1/t^3*Kernel
      #Kernel = Kernel/sum(Kernel)

      #Kernel = rep(1,11)
      #m = length(Kernel)
      ACF = rep(0,m)
      for (i in 0:(m-1)) {
        temp_Kernel = c(Kernel[(i+1):m],rep(0,i))
        ACF[i+1] = sum(Kernel*temp_Kernel)
      }
      ACF = ACF/max(ACF)

      ACF = dbacf(Y_conv, m, type = "correlation", plot = FALSE)$acf
      for (i in 1:m) {
        ACF[i] = max(ACF[i],0)
      }
      CovMatrix = matrix(0,n,n)
      for (i in 1:n) {
        CovMatrix[i,i:min(i+m-1,n)] = ACF[1:min(m,n-i+1)]
      }
      CovMatrix = CovMatrix + t(CovMatrix)
      for (i in 1:n) {
        CovMatrix[i,i] = CovMatrix[i,i]/2 + lambda
      }

      L = chol(CovMatrix)
      r = 1000
      data = matrix(0,r,n)
      for (i in 1:r) {
        X = rnorm(n,0,1)
        X_tilde = L%*%X;
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n, b)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      #res = .simulQuantile_DMUSCLE2(1-alpha,n,ACF,L)
      return(res)
    }else{
      ACF = dbacf(Y, m = m, type = "correlation", plot = F)$acf
      for (i in 1:(m+1)) {
        ACF[i] = max(ACF[i],0)
        print(ACF[i])
      }

      CovMatrix = matrix(0,n,n)
      for (i in 1:n) {
        CovMatrix[i,i:min(i+m-1,n)] = ACF[1:min(m,n-i+1)]
      }
      CovMatrix = CovMatrix + t(CovMatrix)
      for (i in 1:n) {
        CovMatrix[i,i] = CovMatrix[i,i]/2 + lambda
      }
      #return(CovMatrix)
      r = 1000
      data = matrix(n,r,n)
      L = t(chol(CovMatrix))
      print(norm(L%*%t(L)-CovMatrix))
      for (i in 1:r) {
        if (i %%100 == 0){
          print(paste(i/r*100,"% simulated."))
        }
        X = rnorm(n,0,1)
        X_tilde = L%*%X;
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n, b)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      #res = .simulQuantile_DMUSCLE2(1-alpha,n,ACF,L)
      return(res)
    }
  }
}

#' Split
#' @param n sample size
#' @param m splitting size
#' @export
Split <- function(n, m){
  if(n<=m){
    start_index = c(1)
    end_index = c(n)
  }else{
    q = n%/%m
    r = n%%m
    if(r==0){
      start_index = m*0:(q-1)+1
      end_index = m*(1:q)
    }else{
      if(q==1){
        if(r >= 0.75*m){
          start_index = c(1,m+1)
          end_index = c(m,n)
        }else{
          start_index = c(1,(r+m)%/%2+1)
          end_index = c((r+m)%/%2,n)
        }
      }else{
        if(r >= 0.75*m){
          start_index = seq(1,q*m+1,m)
          end_index = c(seq(m,q*m,m),n)
        }else{
          start_index = c(seq(1,(q-1)*m+1,m),(q-1)*m+1+(r+m)%/%2)
          end_index = c(seq(m,(q-1)*m,m),(q-1)*m+(r+m)%/%2,n)
        }
      }
    }
  }
  return(list(start = start_index, end = end_index))
}

#' logg
#' @param x numeric vector
#' @export
logg <- function(x)
  return(.logg(x))

#' maxmin of local errors
#' @param left left boundaries of true signals
#' @param left_hat left boundaries of etimation
#' @param n sample size
#' @export
D = function(left, left_hat, n){
  left = c(left,n)
  left_hat = c(left_hat,n)
  temp = rep(0,length(left))
  for (i in 1:length(left)) {
    temp[i] = min(abs(left[i]-left_hat))
  }
  return(max(temp)/n)
}


#' Over estimation rate
#' @param left left boundaries of true signals
#' @param left_hat left boundaries of etimation
#' @export
OER = function(left,left_hat){
  K = length(left) - 1
  K_hat = length(left_hat) - 1
  return(max(K_hat-K,0)/max(K_hat,1))
}

#' FDR
#' @param left left boundaries of true signals
#' @param left_hat left boundaries of etimation
#' @param n sample size
#' @export
FDR = function(left,left_hat,n){
  K = length(left) - 1
  K_hat = length(left_hat) - 1
  left_hat = c(left_hat,n)
  left = left[left!=1]
  FD = 0
  if(K_hat == 0){
    return(FD)
  }else{
    for (i in 2:K_hat) {
      lb = (left_hat[i]+left_hat[i-1])/2
      rb = (left_hat[i]+left_hat[i+1])/2
      if(sum(left>=lb & left<=rb)==0){
        FD = FD + 1
      }
    }
  }
  return(FD/(K_hat+1))
}

#' V-measure
#' @param left left boundaries of true signals
#' @param left_hat left boundaries of etimation
#' @param n sample size
#' @export
V = function(left,left_hat,n){
  S = length(left)
  S_hat = length(left_hat)
  A = matrix(0,length(left),length(left_hat))
  left = c(left,n+1)
  left_hat = c(left_hat,n+1)
  j = 1
  for (i in 1:S) {
    for (k in left[i]:(left[i+1]-1)) {
      if(k>=left_hat[j] && k<left_hat[j+1]){
        A[i,j] = A[i,j] + 1
      }else{
        j = j + 1
        A[i,j] = A[i,j] + 1
      }
    }
  }
  rs = rowSums(A)
  cs = colSums(A)
  H_CK = 0
  H_C = 0
  H_KC = 0
  H_K = 0
  for (j in 1:S_hat) {
    H_CK = H_CK - sum(A[,j]/n*logg(A[,j]/cs[j]))
  }
  H_C = -sum(rs/n*logg(rs/n))
  for (i in 1:S) {
    H_KC = H_KC - sum(A[i,]/n*logg(A[i,]/rs[i]))
  }
  H_K = - sum(cs/n*logg(cs/n))
  if(H_C == 0){
    h = 1
  }else{
    h = 1 - H_CK/H_C
  }
  if(H_K == 0){
    c = 1
  }else{
    c = 1 - H_KC/H_K
  }
  V = 2*h*c/(h+c)
  return(V)
}

# evaluate setp function
#' evaluate step function
#' @param stepF idealization result
#' @export
evalStepFun <- function(stepF)
{
  ret <- rep(0, stepF$n)
  left <- c(stepF$left, stepF$n)
  for (i in 1:length(stepF$left)) {
    ret[left[i]:left[i+1]] <- stepF$value[i]
  }
  ret
}

# create 'teeth' function
#' teeth function
#' @param n sample size
#' @param K number of change points
#' @param h higher level
#' @export
teethfun <- function(n, K, h=3)
{
  u <- rep(0, n)
  s0 <- n/(K+1)
  i <- 1
  s <- s0
  while (i <= n) {
    u[i:round(s)] <- h*(1 - (round(s/s0)%%2))
    i <- round(s) + 1
    s <- s + s0
  }
  u
}

