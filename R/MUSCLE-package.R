## usethis namespace: start
#' @useDynLib musclePST, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' Segmentation with MUSCLE.
#' @param Y observations.
#' @param q quantiles, can be simulated by simulQuantile_MUSCLE. It overwrites alpha if provided.
#' @param alpha type I error of each multiscale test.
#' @param beta target quantile, beta = 0.5 for median.
#' @param dyadic boolean variable, indicates test statistics obtained from dyadic subintervals or on full subintervals.
#' @param split boolean variable, indicates whether the data will be split by subsets with size m.
#' @param m splitting size.
#' @param details boolean variable, indicates whether computation details will be displayed.
#' @param deconv boolean variable, indicates whether the deconvolution results will be computed.
#' @param lag shift lag (the size of dependency), only for deconvolution (deconv = TRUE).
#' @return multiscale quantile segmentation results.
#' @export
#'
MUSCLE <- function(Y, q, alpha = 0.1, beta = 0.5, dyadic = TRUE, split = FALSE, m = 2000,
                   details = FALSE, deconv = FALSE, lag = 192){
  if(!is.numeric(Y) || length(Y) == 0L || any(!is.finite(Y))){
    stop("Y must be a non-empty finite numeric vector.")
  }
  if(length(beta) != 1L || !is.numeric(beta) || !is.finite(beta) ||
     beta <= 0 || beta >= 1){
    stop("beta must be a single number strictly between 0 and 1.")
  }
  if(missing(q)) {
    q = simulQuantile_MUSCLE(n = length(Y), alpha = alpha, beta = beta, 
                             exact = FALSE, lambda = 1, 
                             deconv = deconv, lag = lag)
  }
  if(!is.numeric(q) || any(!is.finite(q))){
    stop("q must be a finite numeric vector.")
  }
  if(length(dyadic) != 1L || length(split) != 1L || length(details) != 1L ||
     length(deconv) != 1L || !is.logical(dyadic) || !is.logical(split) ||
     !is.logical(details) || !is.logical(deconv)){
    stop("dyadic, split, details, and deconv must be single logical values.")
  }
  if(!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 1){
    stop("m must be a positive number.")
  }
  if(!is.numeric(lag) || length(lag) != 1L || !is.finite(lag) || lag < 0){
    stop("lag must be a non-negative number.")
  }
  n = length(Y)
  if(deconv == FALSE){
    if(split == TRUE){
      print(paste("There are ",ceiling(n/m), "pieces:"))
      if(length(q)>=min(n,m)){
        print(paste("The 1 . piece is computing."))
        split_res = split(n,m)
        num_split = length(split_res$start)
        start = split_res$start[1]
        end = split_res$end[1]
        size = end - start + 1
        if(dyadic == TRUE){
          res = .MUSCLE(Y[start:end],q[1:size],beta,FALSE,details)
        }else{
          res = .MUSCLE_FULL(Y[start:end],q[1:size],beta,FALSE,details)
        }
        num_cp_res = length(res$value)
        if(num_split>=2){
          for (i in 2:num_split) {
            print(paste("The", i,". piece is computing."))
            start = split_res$start[i]
            end = split_res$end[i]
            size = end - start + 1
            if(dyadic == TRUE){
              temp_res = .MUSCLE(Y[start:end],q[1:size],beta,FALSE,details)
            }else{
              temp_res = .MUSCLE_FULL(Y[start:end],q[1:size],beta,FALSE,details)
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
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),beta),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res==1 && num_cp_temp_res>=2){
                value = c(quantile(c(max(c1,d1),min(c2,d2)),beta),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res>=2 && num_cp_temp_res==1){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),beta))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }else{
                value = c(quantile(c(max(c1,d1),min(c2,d2)),beta))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }
            }else{
              if(details == TRUE){
                print(paste("There is a small jump between ",i-1, ". and ", i, ". segments, alternative solution is computing."))
              }
              n = res$n + temp_res$n
              first = res$first
              temp_start = res$left[num_cp_res]
              if(num_cp_temp_res == 1){
                temp_end = res$n + temp_res$n
              }else{
                temp_end = res$n + temp_res$left[2]-1
              }
              temp_size = temp_end - temp_start + 1
              if(details == TRUE){
                print(paste("Problem segment start = ", temp_start, " end = ", temp_end))
                print(paste("Problem segment has length ", temp_size))
              }
              if(temp_size <= length(q)){
                if(details == TRUE){
                  print(paste("Quantiles are sufficiently long."))
                }
                if(dyadic == TRUE){
                  joint_res = .MUSCLE(Y[temp_start:temp_end],q[1:temp_size],beta, TRUE,details)
                }else{
                  joint_res = .MUSCLE_FULL(Y[temp_start:temp_end],q[1:temp_size],beta,TRUE,details)
                }
                if(length(joint_res$left)==1){
                  if(details == TRUE){
                    print(paste("Alternative solution is computed and jump deleted!"))
                  }
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
        if(details == TRUE){
          stop("Quantile size too small!")
        }
      }
    }else{
      if(n > length(q)){
        stop("Sample size and quantile size do not match!")
      }else{
        q = q[1:n]
        if(dyadic == TRUE){
          return(.MUSCLE(Y, q, beta, FALSE,details))
        }else{
          return(.MUSCLE_FULL(Y, q, beta,FALSE,details))
        }
      }
    }
  }else{
    if(split == TRUE){
      print(paste("There are ",ceiling(n/m), "pieces:"))
      if(length(q)>=min(n,2000)){
        print(paste("The 1 . piece is computing."))
        split_res = split(n,m)
        num_split = length(split_res$start)
        start = split_res$start[1]
        end = split_res$end[1]
        size = end - start + 1
        if(dyadic == TRUE){
          res = .DMUSCLE(Y[start:end],q[1:size],beta,lag,FALSE,details)
        }else{
          res = .DMUSCLE_FULL(Y[start:end],q[1:size],beta,lag,FALSE,details)
        }
        num_cp_res = length(res$value)
        if(num_split>=2){
          for (i in 2:num_split) {
            print(paste("The", i, ". piece is computing."))
            start = split_res$start[i]
            end = split_res$end[i]
            size = end - start + 1
            if(dyadic == TRUE){
              temp_res = .DMUSCLE(Y[start:end],q[1:size],beta,lag,FALSE,details)
            }else{
              temp_res = .DMUSCLE_FULL(Y[start:end],q[1:size],beta,lag,FALSE,details)
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
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),beta),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res==1 && num_cp_temp_res>=2){
                value = c(quantile(c(max(c1,d1),min(c2,d2)),beta),
                          temp_res$value[2:num_cp_temp_res])
                left = c(res$left,temp_res$left[2:num_cp_temp_res]+res$n)
                last = temp_res$last
              }else if(num_cp_res>=2 && num_cp_temp_res==1){
                value = c(res$value[1:(num_cp_res-1)],quantile(c(max(c1,d1),min(c2,d2)),beta))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }else{
                value = c(quantile(c(max(c1,d1),min(c2,d2)),beta))
                left = c(res$left)
                last = c(max(c1,d1),min(c2,d2))
              }
            }else{
              if(details == TRUE){
                print(paste("There is a small jump between ",i-1, ". and ", i, ". segments, alternative solution is computing."))
              }
              n = res$n + temp_res$n
              first = res$first
              temp_start = res$left[num_cp_res]
              if(num_cp_temp_res == 1){
                temp_end = res$n + temp_res$n
              }else{
                temp_end = res$n + temp_res$left[2]-1
              }
              temp_size = temp_end - temp_start + 1
              if(details == TRUE){
                print(paste("Problem segment start = ", temp_start, " end = ", temp_end))
                print(paste("Problem segment has length ", temp_size))
              }
              if(temp_size <= length(q)){
                if(details == TRUE){
                  print(paste("Quantiles are sufficiently long."))
                }
                if(dyadic == TRUE){
                  joint_res = .DMUSCLE(Y[temp_start:temp_end],q[1:temp_size],beta, lag,TRUE,details)
                }else{
                  joint_res = .DMUSCLE_FULL(Y[temp_start:temp_end],q[1:temp_size],beta,lag,TRUE,details)
                }
                if(length(joint_res$left)==1){
                  if(details == TRUE){
                    print(paste("Alternative solution is computed and jump deleted!"))
                  }
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
        stop("Quantile size too small!")
      }
    }else{
      if(n > length(q)){
        stop("Sample size and quantile size do not match!")
      }else{
        q = q[1:n]
        if(dyadic == TRUE){
          return(.DMUSCLE(Y, q, beta,lag,FALSE,details))
        }else{
          return(.DMUSCLE_FULL(Y, q, beta,lag,FALSE,details))
        }
      }
    }
  }

}

#' Segmentation with MMUSCLE.
#' @param Y observations.
#' @param q_matrix quantile matrix, can be simulated by simulQuantile_MMUSCLE. It overwrites alpha if provided.
#' @param beta_vec target quantile vector, beta = 0.25, 0.5 and 0.75 stand for first quartile, median and third quartile respectively.
#' @return multiple multiscale quantiles segmentation results.
#' @export
#'
MMUSCLE <- function(Y, q_matrix, alpha = 0.1, beta_vec = c(0.25,0.5,0.75)){
  if(!is.numeric(Y) || length(Y) == 0L || any(!is.finite(Y))){
    stop("Y must be a non-empty finite numeric vector.")
  }
  if(!is.numeric(beta_vec) || length(beta_vec) == 0L ||
     any(!is.finite(beta_vec)) || any(beta_vec <= 0) || any(beta_vec >= 1)){
    stop("beta_vec must contain numbers strictly between 0 and 1.")
  }
  if(missing(q_matrix)) {
    q_matrix = simulQuantile_MMUSCLE(n = length(Y), alpha = alpha, 
                                     beta_vec = beta_vec)
  }
  if(!is.numeric(q_matrix) || length(q_matrix) == 0L ||
     any(!is.finite(q_matrix)) || ncol(as.matrix(q_matrix)) < length(Y)){
    stop("q_matrix must contain finite quantiles for every observation.")
  }
  return(.MMUSCLE(Y,q_matrix,beta_vec,FALSE,FALSE))
}

#' Simulate MUSCLE and D-MUSCLE quantiles
#' @param alpha type I error of each multiscale test.
#' @param n sample size.
#' @param beta target quantile, beta = 0.5 for median regression.
#' @param exact boolean variable, if exact = TRUE, quantiles will be simulated;
#'        otherwise, quantiles will be computed via appprximation
#' @param deconv boolean variable, deconv = TRUE for deconvolution (DMUSCLE).
#' @param Y observations, only required for deconvolution.
#' @param lambda regularization parameter, only required for deconvolution.
#' @param lag shift lag (the size of dependency), only for deconvolution (deconv = TRUE).
#' @param sr sampling rate, default value is 20k, only required for deconvolution.
#' @param Kernel convolution kernel, only required for deconvolution.
#' @param ACF autocorrelation function, only required for deconvolution.
#' @param E16 A logical variable, TRUE for Element-16 measurement, only required for deconvolution.
#' @return Simulated quantiles.
#' @export
simulQuantile_MUSCLE = function(n, alpha = 0.1, beta = 0.5, exact = FALSE, lambda = 1, deconv = FALSE,
                                Kernel = NULL, ACF = NULL, Y = NULL, E16 = FALSE,
                                sr = 20000, lag = 192){
  if(length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 1){
    stop("n must be a positive number.")
  }
  if(length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
     alpha <= 0 || alpha >= 1){
    stop("alpha must be a single number strictly between 0 and 1.")
  }
  if(length(beta) != 1L || !is.numeric(beta) || !is.finite(beta) ||
     beta <= 0 || beta >= 1){
    stop("beta must be a single number strictly between 0 and 1.")
  }
  if(deconv == FALSE){
    if(beta == 0.5 && n<=30000){
      if(alpha < 0.1){
        return(.simulQuantile_MUSCLE(1-alpha, n, beta))
      }else if(alpha == 0.1){
        return(q_30000_0_1_0_5[1:n])
      }else if(alpha < 0.2){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.1)*q_30000_0_2_0_5[1:n]+(0.2-alpha)*q_30000_0_1_0_5[1:n]))
        }
      }else if(alpha == 0.2){
        return(q_30000_0_2_0_5[1:n])
      }else if(alpha < 0.3){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.2)*q_30000_0_3_0_5[1:n]+(0.3-alpha)*q_30000_0_2_0_5[1:n]))
        }
      }else if(alpha == 0.3){
        return(q_30000_0_3_0_5[1:n])
      }else if(alpha < 0.4){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.3)*q_30000_0_4_0_5[1:n]+(0.4-alpha)*q_30000_0_3_0_5[1:n]))
        }
      }else if(alpha == 0.4){
        return(q_30000_0_4_0_5[1:n])
      }else if(alpha < 0.5){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.4)*q_30000_0_5_0_5[1:n]+(0.5-alpha)*q_30000_0_4_0_5[1:n]))
        }
      }else if(alpha == 0.5){
        return(q_30000_0_5_0_5[1:n])
      }else if(alpha < 0.6){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.5)*q_30000_0_6_0_5[1:n]+(0.6-alpha)*q_30000_0_5_0_5[1:n]))
        }
      }else if(alpha == 0.6){
        return(q_30000_0_6_0_5[1:n])
      }else if(alpha < 0.7){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.6)*q_30000_0_7_0_5[1:n]+(0.7-alpha)*q_30000_0_6_0_5[1:n]))
        }
      }else if(alpha == 0.7){
        return(q_30000_0_7_0_5[1:n])
      }else if(alpha < 0.8){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.7)*q_30000_0_8_0_5[1:n]+(0.8-alpha)*q_30000_0_7_0_5[1:n]))
        }
      }else if(alpha == 0.8){
        return(q_30000_0_8_0_5[1:n])
      }else if(alpha < 0.9){
        if(exact == TRUE){
          return(.simulQuantile_MUSCLE(1-alpha, n, beta))
        }else{
          return(10*((alpha-0.8)*q_30000_0_9_0_5[1:n]+(0.9-alpha)*q_30000_0_8_0_5[1:n]))
        }
      }else if(alpha == 0.9){
        return(q_30000_0_9_0_5[1:n])
      }else{
        return(.simulQuantile_MUSCLE(1-alpha, n, beta))
      }
    }else{
      return(.simulQuantile_MUSCLE(1-alpha, n, beta))
    }
  }else{
    res = rep(0,n)
    if(!is.null(Kernel)){
      lag = length(Kernel)
      Kernelfft = fft(c(Kernel,rep(0,n-1)))/(n+lag-1)
      r = 50/min(alpha,1-alpha)
      sd = 1 / sqrt(sum(Kernel^2))
      data = matrix(0,r,n)
      ACF = convolve(Kernel,Kernel, type = "open")[lag:(2*lag-1)]
      ACF = ACF/max(ACF)
      for (i in 1:r) {
        X = rnorm(n+lag-1,0,sd)
        X_tilde = Re(fft(fft(X)*Kernelfft,inverse = TRUE))[lag:(n+lag-1)]
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      return(res)
    }else if(!is.null(ACF)){
      for (i in 1:(lag+1)) {
        ACF[i] = max(ACF[i],0)
        print(ACF[i])
      }

      CovMatrix = matrix(0,n,n)
      for (i in 1:n) {
        CovMatrix[i,i:min(i+lag-1,n)] = ACF[1:min(lag,n-i+1)]
      }
      CovMatrix = CovMatrix + t(CovMatrix)

      for (i in 1:n) {
        CovMatrix[i,i] = CovMatrix[i,i]/2
      }
      CovMatrix.eigen = eigen(CovMatrix)
      D = CovMatrix.eigen$values
      V = CovMatrix.eigen$vectors
      D_sqrt = rep(0,length(D))
      for (i in 1:n) {
        if(D[i] < 0){
          print(paste("D[",i,"]=",D[i],"is deleted!"))
        }
        D[i] = max(D[i],0)
        D_sqrt[i] = sqrt(D[i])
      }
      D_sqrt = diag(D_sqrt)
      L = V%*%D_sqrt
      #print(norm(L%*%t(L)-CovMatrix))

      r = 50/min(alpha,1-alpha)
      data = matrix(0,r,n)
      for (i in 1:r) {
        if (i %%100 == 0){
          print(paste(i/r*100,"% simulated."))
        }
        X = rnorm(n,0,1)
        X_tilde = L%*%X;
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      #res = .simulQuantile_DMUSCLE2(1-alpha,n,ACF,L)
      return(res)
    }else if(E16 == TRUE){
      if(sr == 1250){
        lag = 1024
      }else if(sr == 5000){
        lag = 256
      }else if(sr == 10000){
        lag = 128
      }else if(sr == 20000){
        lag = 64
      }else if(sr == 50000){
        lag = 200
      }else if(sr == 100000){
        lag = 100
      }else if(sr == 200000){
        lag = 50
      }else {
        lag = 10
      }
      #lag = 64
      #t = lag/sr
      #x = (0:(3*lag))/sr
      #Kernel = (x>=0)*(x<=t)*1/2*x^2+(x>t)*(x<2*t)*(-x^2+3*t*x-3/2*t^2)+(x>=2*t)*(x<=3*t)*(1/2*x^2-3*t*x+9/2*t^2)
      #Kernel = 1/t^3*Kernel
      #Kernel = Kernel/sum(Kernel)

      #Kernel = rep(1,11)
      #lag = length(Kernel)
      ACF = rep(0,lag)
      for (i in 0:(lag-1)) {
        temp_Kernel = c(Kernel[(i+1):lag],rep(0,i))
        ACF[i+1] = sum(Kernel*temp_Kernel)
      }
      ACF = ACF/max(ACF)

      if(is.null(Y)){
        stop("Y is required for the E16 deconvolution branch.")
      }
      ACF = dbacf::dbacf(Y, lag, type = "correlation", plot = FALSE)$acf
      for (i in 1:lag) {
        ACF[i] = max(ACF[i],0)
      }
      CovMatrix = matrix(0,n,n)
      for (i in 1:n) {
        CovMatrix[i,i:min(i+lag-1,n)] = ACF[1:min(lag,n-i+1)]
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
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      #res = .simulQuantile_DMUSCLE2(1-alpha,n,ACF,L)
      return(res)
    }else{
      ACF = dbacf::dbacf(Y, m = lag, type = "correlation", plot = F)$acf
      for (i in 1:(lag+1)) {
        ACF[i] = max(ACF[i],0)
        #print(ACF[i])
      }

      CovMatrix = matrix(0,n,n)
      for (i in 1:n) {
        CovMatrix[i,i:min(i+lag-1,n)] = ACF[1:min(lag,n-i+1)]
      }
      CovMatrix = CovMatrix + t(CovMatrix)

      for (i in 1:n) {
        CovMatrix[i,i] = CovMatrix[i,i]/2
      }
      CovMatrix.eigen = eigen(CovMatrix)
      D = CovMatrix.eigen$values
      V = CovMatrix.eigen$vectors
      D_sqrt = rep(0,length(D))
      for (i in 1:n) {
        #if(D[i] < 0){
        #  print(paste("D[",i,"]=",D[i],"is deleted!"))
        #}
        D[i] = max(D[i],0)
        D_sqrt[i] = sqrt(D[i])
      }
      D_sqrt = diag(D_sqrt)
      L = V%*%D_sqrt

      #for (i in 1:n) {
      #  CovMatrix[i,i] = CovMatrix[i,i]/2 + lambda
      #}
      #L = t(chol(CovMatrix))

      #print(norm(L%*%t(L)-CovMatrix))

      r = 50/min(alpha,1-alpha)
      data = matrix(0,r,n)
      for (i in 1:r) {
        X = rnorm(n,0,1)
        X_tilde = L%*%X;
        data[i,] = .simulQuantile_DMUSCLE(X_tilde, ACF, n)
      }
      for (i in 1:n) {
        res[i] = quantile(data[,i],1-alpha)
      }
      return(res)
    }
  }
}

#' Simulate MMUSCLE quantiles.
#' @param alpha type I error of each multiscale test.
#' @param n sample size.
#' @param beta_vec target quantile vector, target quantile vector, beta = 0.25,
#' 0.5 and 0.75 stand for first quartile, median and third quartile respectively.
#' @export
simulQuantile_MMUSCLE = function(n, alpha = 0.1, beta_vec = c(0.25,0.5,0.75)){
  if(!is.numeric(beta_vec) || length(beta_vec) == 0L ||
     any(!is.finite(beta_vec)) || any(beta_vec <= 0) || any(beta_vec >= 1)){
    stop("beta_vec must contain numbers strictly between 0 and 1.")
  }
  num_beta = length(beta_vec)
  res = matrix(0,num_beta,n)
  for(i in 1:num_beta){
    res[i,] = .simulQuantile_MUSCLE(1-alpha, n, beta_vec[i])
  }
  return(res)
}


#' Split n observations into pieces with size m
#' @param n sample size
#' @param m splitting size
#' @export
split <- function(n, m){
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

#' Returns log(x) for positive x and zero for x = 0.
#' @param x numeric vector
#' @export
logg <- function(x)
  return(.logg(x))

#' Compute localization errors of estimated change points.
#' @param left true change points.
#' @param left_hat estimated change points.
#' @param n sample size.
#' @export
Local_errors = function(left, left_hat, n){
  left = c(left,n)
  left_hat = c(left_hat,n)
  temp = rep(0,length(left))
  for (i in 1:length(left)) {
    temp[i] = min(abs(left[i]-left_hat))
  }
  return(max(temp)/n)
}


#' Compute over estimation rate (OER) of change points estimation.
#' @param left true change points.
#' @param left_hat estimated change points.
#' @export
OER = function(left,left_hat){
  K = length(left) - 1
  K_hat = length(left_hat) - 1
  return(max(K_hat-K,0)/max(K_hat,1))
}

#' Compute false discovery rate (FDR) of change points estimation.
#' @param left true change points.
#' @param left_hat estimated change points.
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

#' Compute V-measure of change points estimation.
#' @param left true change points.
#' @param left_hat estimated change points.
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

#' Compute segmentation result.
#' @param stepF segmentation result.
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

#' Compute segmentation result (only for MMUSCLE).
#' @param stepF segmentation result.
#' @param m index of segmentation result.
#' @export
eval_StepFun_M <- function(stepF, m){
  n = stepF$n
  ret = c()
  num_seg = length(stepF$left)
  if(num_seg>1){
    for(i in 1:(num_seg-1)){
      ret = c(ret, rep(stepF$value[m,i],stepF$left[i+1]-stepF$left[i]))
    }
    ret = c(ret, rep(stepF$value[m,num_seg],n-stepF$left[num_seg]+1))
  }else{
    ret = rep(stepF$value[m,1],n)
  }
  return(ret)
}

#' create 'teeth' function
#' @param n sample size.
#' @param K number of change points.
#' @param h higher level.
#' @export
teethfun <- function(n, K, h=3){
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
