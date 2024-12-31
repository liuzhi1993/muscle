#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
using namespace Rcpp;
using namespace std;

double penfsC(int n, int len, NumericVector &q){
  double res = q[n] + sqrt(2*log(exp(1)*(n+1)/(len)));
  return res;
}

double penfsCC(int n, int len, double q){
  double res = q + sqrt(2*log(exp(1)*(n+1)/(len)));
  return res;
}

double f(double x,double b, double q){
  double res = (x*log(x/b)+(1-x)*log((1-x)/(1-b)))-q;
  return res;
}

double f_diff(double x, double b, double q){
  double res = log(x*(1-b)/((1-x)*b));
  return res;
}
double x1(double q, double b){
  double x = 1e-6;
  double x_new;
  double diff = 1;
  while(diff>=1e-6){
    x_new = x - f(x, b, q)/f_diff(x,b,q);
    // Rcout << "f(x,b,q) = " << f(x,b,q) << "\n";
    // Rcout << "f_diff(x,b,q) =" << f_diff(x,b,q) << "\n";
    // Rcout << "x_new = << " << x_new << "\n";
    diff = abs(x_new -x);
    x = x_new;
    if(x_new < 0){
      return(0);
    }
  }
  return x_new;
}
double x2(double q, double b){
  double x = 1 - 1e-6;
  double x_new;
  double diff = 1;
  while(diff>=1e-6){
    x_new = x - f(x, b, q)/f_diff(x,b,q);
    diff = abs(x_new -x);
    x = x_new;
    if(x_new >1){
      return(1);
    }
  }
  return x_new;
}

struct WT{
  vector<int> IDX;
  vector<vector<int>> RANK;
  vector<vector<int>> LR;

  int n;
  int level_max;
  WT(NumericVector &Y){
    n = Y.size();
    level_max = ceil(log2(n));
    int N = pow(2,level_max+1)-1;
    vector<int> LR_row (N);
    LR = vector<vector<int>> (2, LR_row);
    vector<int> RANK_row (n);
    RANK = vector<vector<int>> (level_max, RANK_row);
    IDX = vector<int> (n);
    //YY = vector<double> (n);
    left_and_right_bounds(LR,n,level_max);
    vector<int> V0(n);
    vector<int> V1(n);

    IDX = sort_indexes(Y);
    int old_Idx,Idx,i,k,lower_bound,upper_bound,middle_bound,Idx_left,Idx_right,Idx_middel;
    for(i = 0; i<n; i++){
      V0[IDX[i]] = i;
    }

    for(i = 0;i<level_max-1;i++){
      if(i%2 == 0){
        for(Idx = pow(2,i+1)-1;Idx<pow(2,i+2)-2;Idx=Idx+2){
          old_Idx = (Idx-1) >> 1;
          lower_bound = LR[0][old_Idx];
          upper_bound = LR[1][old_Idx];
          middle_bound = LR[1][Idx];
          Idx_left = 0;
          Idx_right = 0;
          Idx_middel = 0;

          k = lower_bound;
          if(V0[k] <= middle_bound){
            V1[lower_bound+Idx_left] = V0[k];
            ++Idx_left;
            RANK[i][k] = 1;
          }else{
            V1[middle_bound+1+Idx_right] = V0[k];
            ++Idx_right;
            RANK[i][k] = 0;
          }

          for(k = lower_bound+1;k<=upper_bound;k++){
            if(V0[k] <= middle_bound){
              V1[lower_bound+Idx_left] = V0[k];
              ++Idx_left;
              RANK[i][k] = RANK[i][k-1] + 1;
            }else{
              V1[middle_bound+1+Idx_right] = V0[k];
              ++Idx_right;
              RANK[i][k] = RANK[i][k-1];
            }
          }
        }
      }else{
        for(Idx = pow(2,i+1)-1;Idx<pow(2,i+2)-2;Idx=Idx+2){
          old_Idx = (Idx-1) >> 1;
          lower_bound = LR[0][old_Idx];
          upper_bound = LR[1][old_Idx];
          middle_bound = LR[1][Idx];
          Idx_left = 0;
          Idx_right = 0;
          Idx_middel = 0;

          k = lower_bound;
          if(V1[k] <= middle_bound){
            V0[lower_bound+Idx_left] = V1[k];
            ++Idx_left;
            RANK[i][k] = 1;
          }else{
            V0[middle_bound+1+Idx_right] = V1[k];
            ++Idx_right;
            RANK[i][k] = 0;
          }

          for(k = lower_bound+1;k<=upper_bound;k++){
            if(V1[k] <= middle_bound){
              V0[lower_bound+Idx_left] = V1[k];
              ++Idx_left;
              RANK[i][k] = RANK[i][k-1] + 1;
            }else{
              V0[middle_bound+1+Idx_right] = V1[k];
              ++Idx_right;
              RANK[i][k] = RANK[i][k-1];
            }
          }
        }
      }
    }

    i = level_max - 1;
    if(i%2 == 0){
      for(Idx = pow(2,i)-1; Idx<=pow(2,i+1)-2; Idx++){
        lower_bound = LR[0][Idx];
        upper_bound = LR[1][Idx];
        k = lower_bound;
        if(upper_bound == lower_bound){
          V1[k] = V0[k];
          RANK[i][k] = 1;
        }else{
          if(V0[k] == lower_bound){
            V1[k] = V0[k];
            V1[k+1] = V0[k+1];
            RANK[i][k] = 1;
            RANK[i][k+1] = 1;
          }else{
            V1[k] = V0[k+1];
            V1[k+1] = V0[k];
            RANK[i][k] = 0;
            RANK[i][k+1] = 1;
          }
        }
      }
    }else{
      for(Idx = pow(2,i)-1; Idx<=pow(2,i+1)-2; Idx++){
        lower_bound = LR[0][Idx];
        upper_bound = LR[1][Idx];
        k = lower_bound;
        if(upper_bound == lower_bound){
          V0[k] = V1[k];;
          RANK[i][k] = 1;
        }else{
          if(V1[k] == lower_bound){
            V0[k] = V1[k];
            V0[k+1] = V1[k+1];
            RANK[i][k] = 1;
            RANK[i][k+1] = 1;
          }else{
            V0[k] = V1[k+1];
            V0[k+1] = V1[k];
            RANK[i][k] = 0;
            RANK[i][k+1] = 1;
          }
        }
      }
    }
  }

  void left_and_right_bounds (vector<vector<int>> &LR, int n, int l) const {
    LR[0][0] = 0;
    LR[1][0] = n-1;

    for(int k = 0;k<=l-1;k++){
      for(int i = pow(2,k)-1; i<=pow(2,k+1)-2;i++){
        LR[0][2*i+1] = LR[0][i];
        LR[1][2*i+2] = LR[1][i];
        LR[1][2*i+1] = (LR[1][2*i+2]+LR[0][2*i+1])/2;
        LR[0][2*i+2] = fmin(LR[1][2*i+1]+1,LR[1][2*i+2]);
      }
    }
  }
  //
  vector<int> sort_indexes(const NumericVector &Y) const {

    // initialize original index locations
    int n = Y.size();
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&Y](size_t i1, size_t i2) {return Y[i1] < Y[i2];});

    return idx;
  }

  int rngQuantile(int k, int l, int r) const {
    int level = 0;
    int idx = 0;
    int nr_0_left, nr_0_right, old_bound;
    while (level < level_max) {
      if (l == LR[0][idx]) {
        nr_0_left = 0;
        nr_0_right = RANK[level][r];
      } else {
        nr_0_left = RANK[level][l-1];
        nr_0_right = RANK[level][r] - nr_0_left;
      }

      if (k <= nr_0_right) {
        idx = (idx << 1) + 1;
        l = LR[0][idx] + nr_0_left;
        r = LR[0][idx] + RANK[level][r] - 1;
      } else {
        old_bound = LR[0][idx];
        idx = (idx << 1) + 2;
        k -=  nr_0_right;
        l += LR[0][idx] - old_bound - nr_0_left;
        r += LR[0][idx]- old_bound - RANK[level][r];
      }
      ++level;
    }
    return l;
  }

  // int rngCounting(int l, int r, double v){
  //   int level = 0;
  //   int idx = 0;
  //   int res  = 0;
  //   int res_max = r-l+1;
  //   int nr_0_left, nr_0_right,old_bound;
  //
  //   while (level < level_max) {
  //     if (l == LR[0][idx]) {
  //       nr_0_left = 0;
  //       nr_0_right = RANK[level][r];
  //     } else {
  //       nr_0_left = RANK[level][l-1];
  //       nr_0_right = RANK[level][r] - nr_0_left;
  //     }
  //
  //     if (v <= YY[IDX[LR[1][idx*2+1]]]) {
  //       if(nr_0_right == 0){
  //         return res;
  //       }
  //       idx = (idx << 1) + 1;
  //       l = LR[0][idx] + nr_0_left;
  //       r = LR[0][idx] + RANK[level][r] - 1;
  //     } else {
  //       old_bound = LR[0][idx];
  //       idx = (idx << 1) + 2;
  //       res = res + nr_0_right;
  //       if(res >= res_max){
  //         return res;
  //       }
  //       l += LR[0][idx] - old_bound - nr_0_left;
  //       r += LR[0][idx] - old_bound - RANK[level][r];
  //     }
  //     ++level;
  //   }
  //   return res;
  // }
};

// [[Rcpp::export(.MUSCLE)]]
List MUSCLE(NumericVector &Y, NumericVector &q, double beta, bool test, bool details){
   //Rcout << "MUSCLE is running, please wait!\n";
   int n = Y.size();
   int imax,lq_idx, uq_idx,lk,uk, llen, ulen, llen_idx, ulen_idx,lq_max_idx, uq_max_idx,idx,power,l,lp,up;
   double lb_max, ub_max, lb, ub, q_tilde, aux_Y, aux_mu, aux_c,aux_lb, aux_ub,lq,uq, lq_max,uq_max;
   double lb_first, ub_first, lb_last, ub_last;
   NumericVector C(n);  //optimal cost
   NumericVector J(n);  //number of jumps
   NumericVector mu(n); //optimal value of step function
   NumericVector L(n);  //leftmost end points
   for (int i = 0; i < n; ++i){
     C[i] = INFINITY;
     J[i] = n;
     L[i] = 0;
   }
   double Q_max = fmax(-log(beta),-log(1-beta)); //maximum of f on [0,1]

   NumericVector Q_tilde(n);
   NumericVector Q_tilde_max(n);
   IntegerVector Q_tilde_max_idx(n);
   // bool get_p_min = false;
   bool get_p_min = false;
   int p_min = 1;
   for(int i = 0; i<n; i++){
     l = i + 1;
     for(int j = 0; j<n; j++){
       if(j==i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
         if(get_p_min == false && Q_tilde[j]<= Q_max -1e-6){
           get_p_min = true;
           p_min = floor(log2(l));
         }
       }else if(j>i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
       }else{
         Q_tilde[j] = R_NegInf;
       }
     }
     Q_tilde_max[i] = max(Q_tilde);
     Q_tilde_max_idx[i] = which_max(Q_tilde);
   }

   l = ceil(log2(n));
   WT Tree = WT(Y);

   //First search
   imax = 0;
   lb_max = R_NegInf;
   ub_max = R_PosInf;
   for (int i = 0; i< n; ++i) {
     lb = R_NegInf;
     ub = R_PosInf;
     power = fmax(floor(log2(i)),0);
     // Rcout << "power = " << power << "\n";
     // Rcout << "i = " << i << "\n";
     for (int p = power ; p>= p_min; p--) {
       l = pow(2,p);
       // Rcout << "l = " << l << "\n";
       for (int k = 0; k <= i-l+1; k = k+ 2) {
         // Rcout << "k = " << k << "\n";
         idx = fmin(fmax(1,(int)ceil(l*beta)),l);
         aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k,k+l-1)]];
         // Rcout << "aux_Y = " << aux_Y << "\n";
         //aux_Y = get_Quantile(idx,k,k+l-1,R,lr,RANK);
         q_tilde = pow(penfsC(i,l,q),2)/(2*l);
         // Rcout << "q_tilde = " << q_tilde << "\n";
         if(q_tilde<=Q_max-1e-6){
           lq = x1(q_tilde,beta);
           lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
           if(beta==0.5){
             uq_idx = l-lq_idx;
           }else{
             uq = x2(q_tilde,beta);
             uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
           }
           aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k,k+l-1)]];
           //aux_lb = get_Quantile(lq_idx,k,k+l-1,R,lr,RANK);
           if(aux_lb>lb){
             lb = aux_lb;
             lk = k;
             llen = l;
             lp = p;
             llen_idx = Q_tilde_max_idx[llen-1];
           }
           // Rcout << "lb = " << lb << "\n";
           aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k,k+l-1)]];
           //aux_ub = get_Quantile(uq_idx,k,k+l-1,R,lr,RANK);
           if(aux_ub<ub){
             ub = aux_ub;
             uk = k;
             ulen = l;
             up = p;
             ulen_idx = Q_tilde_max_idx[ulen-1];
           }
           // Rcout << "ub = " << ub << "\n";
         }
       }
     }
     if(lb <= ub){
       lb_first = lb;
       ub_first = ub;
       lb_last = lb;
       ub_last = ub;
       if(aux_Y > ub){
         mu[i] = ub;
       }else if(aux_Y < lb){
         mu[i] = lb;
       }else{
         mu[i] = aux_Y;
       }
       // Rcout << "mu[i] = " << mu[i] << "\n";
       // Rcout << "lb = " << lb << "\n";
       // Rcout << "ub = " << ub << "\n";
       C[i] = sum((Y[Range(0,i)]-mu[i])*(beta-ifelse(Y[Range(0,i)]<mu[i],1.0,0.0)));
       // Rcout << "C[i]" << i << "=" << C[i] << "\n";
       L[i] = 1;
       J[i] = 0;
     }else{
       q_tilde = Q_tilde_max[llen-1];
       // Rcout << "llen = " << llen << "\n";
       if(q_tilde<=Q_max-1e-6){
         lq_max = x1(q_tilde,beta);
         lq_max_idx = fmin(fmax(1,ceil(llen*lq_max)),llen);
         lb_max = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx,lk,lk+llen-1)]];
       }

       q_tilde = Q_tilde_max[ulen-1];
       // Rcout << "ulen = " << ulen << "\n";
       if(q_tilde<=Q_max-1e-6){
         uq_max = x2(q_tilde,beta);
         uq_max_idx = fmax(fmin(ulen,floor(ulen*uq_max)),1);
         ub_max = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx,uk,uk+ulen-1)]];
       }
       // Rcout << "lb_max =" << lb_max << " ub_max = "<< ub_max << "\n";
       if(lb_max>ub_max){
         // Rcout << "i = " << i << "\n";
         // Rcout << "No further solution!\n";
         break;
       }
     }
   }
   if(test == true && std::isinf(C[n-1])){
     if(details == true){
       Rcout << "Fast test is done! No alternative solution!\n";
     }
     IntegerVector left(2);
     left[0] = 1;
     left[1] = 1;
     return List::create(Named("left") = left);
   }
   // return List::create(Named("n") = n);
   // Later searchs
   int imin = 0;
   int njmp = 0;
   int jmin, aux_imin, aux_imax,i_llen_idx,j_llen_idx,i_ulen_idx,j_ulen_idx,lq_max_idx_i,uq_max_idx_i,lq_max_idx_j;
   bool isstop = FALSE;
   double lq_max_i,lb_max_i,uq_max_i, ub_max_i,lq_max_j,lb_max_j,uq_max_j,ub_max_j, uq_max_idx_j;

   while (std::isinf(C[n-1]) && njmp < n) {
     jmin = n;
     aux_imin = n;
     aux_imax = 0;
     lb_max_i = R_NegInf;
     ub_max_i = R_PosInf;
     for (int i=imin; i<n;i++) {
       // Rcout << "i = " << i << "\n";
       if(J[i] == njmp){
         jmin = fmin(jmin,i);
       }else if(std::isinf(C[i])){
         isstop = FALSE;
         lb_max_j = R_NegInf;
         ub_max_j = R_PosInf;
         for (int j=i-1; j>=jmin;j--) {
           if(J[j] == njmp){
             lb = R_NegInf;
             ub = R_PosInf;
             power = floor(log2(i-j));
             for (int p= power; p >= p_min;p--) {
               l = pow(2,p);
               for (int k=j; k<= i-l; k = k+2) {
                 idx = fmin(fmax(1,(int)ceil(l*beta)),l);
                 aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k+1,k+l)]];
                 //aux_Y = get_Quantile(idx,k+1,k+l,R,lr,RANK);
                 q_tilde = pow(penfsC(i-j-1,l,q),2)/(2*l);;
                 if(q_tilde<=Q_max-1e-6){
                   lq = x1(q_tilde,beta);
                   lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
                   if(beta==0.5){
                     uq_idx = l-lq_idx;
                   }else{
                     uq = x2(q_tilde,beta);
                     uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
                   }
                   aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k+1,k+l)]];
                   //aux_lb = get_Quantile(lq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_lb>lb){
                     lb = aux_lb;
                     lk = k;
                     llen = l;
                     lp = p;
                     i_llen_idx = n - imin - 1;
                     j_llen_idx = i - jmin - 1;
                   }
                   aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k+1,k+l)]];
                   //aux_ub = get_Quantile(uq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_ub<ub){
                     ub = aux_ub;
                     uk = k;
                     ulen = l;
                     up = p;
                     i_ulen_idx = n - imin - 1;
                     j_ulen_idx = i - jmin - 1;
                   }
                 }
               }
             }
             if(lb<=ub){
               lb_last = lb;
               ub_last = ub;
               if(aux_Y > ub){
                 aux_mu = ub;
               }else if(aux_Y < lb){
                 aux_mu = lb;
               }else{
                 aux_mu = aux_Y;
               }
               if(i < aux_imin){
                 aux_imin = i;
               }
               if(i > aux_imax){
                 aux_imax = i;
               }
               aux_c = sum((Y[Range(j+1,i)]-aux_mu)*(beta-ifelse(Y[Range(j+1,i)]<aux_mu,1.0,0.0)))+C[j];
               if(aux_c < C[i]){
                 C[i] = aux_c;
                 mu[i] = aux_mu;
                 L[i] = j+2;
                 J[i] = njmp + 1;
                 // Rcout << "C[i]" << i << "=" << C[i] << "\n";
               }
             }else{
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_i = x1(q_tilde,beta);
                 lq_max_idx_i = fmin(fmax(1,(int)ceil(llen*lq_max_i)),llen);
                 lb_max_i = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_i,lk+1,lk+llen)]];
                 //lb_max_i = get_Quantile(lq_max_idx_i,lk+1,lk+llen,R,lr,RANK);
               }

               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_i = x2(q_tilde,beta);
                 uq_max_idx_i = fmax(fmin(ulen,(int)floor(ulen*uq_max_i)),1);
                 ub_max_i = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_i,uk+1,uk+ulen)]];
                 //ub_max_i = get_Quantile(uq_max_idx_i,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_i>ub_max_i && j >= imax){
                 isstop = TRUE;
               }
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_j = x1(q_tilde,beta);
                 lq_max_idx_j = fmin(fmax(1,(int)ceil(llen*lq_max_j)),llen);
                 lb_max_j = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_j,lk+1,lk+llen)]];
                 //lb_max_j = get_Quantile(lq_max_idx_j,lk+1,lk+llen,R,lr,RANK);
               }
               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_j = x2(q_tilde,beta);
                 uq_max_idx_j = fmax(fmin(ulen,(int)floor(ulen*uq_max_j)),1);
                 ub_max_j = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_j,uk+1,uk+ulen)]];
                 //ub_max_j = get_Quantile(uq_max_idx_j,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_j>ub_max_j){
                 break;
               }
             }
           }
         }
         if(std::isinf(C[i])&&isstop){
           break;
         }
       }
     }
     imin = aux_imin;
     imax = aux_imax;
     njmp = njmp + 1;
     //Rcout << "Jump from" << njmp-1 << "to" << njmp << "\n";
   }

   NumericVector value(njmp+1);
   IntegerVector left(njmp+1);
   int cnt = n-1;
   for (int i = njmp; i >= 0; --i) {
     left[i] = L[cnt];
     value[i] = mu[cnt];
     cnt = L[cnt] - 2;
   }
   NumericVector first(2);
   NumericVector last(2);
   first[0] = lb_first;
   first[1] = ub_first;
   last[0] = lb_last;
   last[1] = ub_last;

   return List::create(Named("value") = value, Named("left") = left,
                       Named("n") = n, Named("first") = first,
                       Named("last") = last);
 }


// [[Rcpp::export(.DMUSCLE)]]

List DMUSCLE(NumericVector &Y, NumericVector &q, double beta, int lag, bool test, bool details){
   //Rcout << "MUSCLE is running, please wait!\n";
   int n = Y.size();
   if (n <= lag) {
    Rprintf("Error in DMUSCLE: the length of data should > the length of kernel! \n");
    return -1;
   }
   int imax,lq_idx, uq_idx,lk,uk, llen, ulen, llen_idx, ulen_idx,lq_max_idx, uq_max_idx,idx,power,l,lp,up;
   double lb_max, ub_max, lb, ub, q_tilde, aux_Y, aux_mu, aux_c,aux_lb, aux_ub,lq,uq, lq_max,uq_max;
   double lb_first, ub_first, lb_last, ub_last;
   NumericVector C(n);  //optimal cost
   NumericVector J(n);  //number of jumps
   NumericVector mu(n); //optimal value of step function
   NumericVector L(n);  //leftmost end points
   for (int i = 0; i < n; ++i){
     C[i] = INFINITY;
     J[i] = n;
     L[i] = 0;
   }
   double Q_max = fmax(-log(beta),-log(1-beta)); //maximum of f on [0,1]
   NumericVector Q_tilde(n);
   NumericVector Q_tilde_max(n);
   IntegerVector Q_tilde_max_idx(n);
   // bool get_p_min = false;
   int p_min = 1;
   for(int i = 0; i<n; i++){
     l = i + 1;
     // Rcout << "i = " << i << "\n";
     for(int j = 0; j<n; j++){
       if(j>=i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
       }else{
         Q_tilde[j] = R_NegInf;
       }
       //Rcout << Q_tilde;
     }
     //Rcout << "\n";
     Q_tilde_max[i] = max(Q_tilde);
     Q_tilde_max_idx[i] = which_max(Q_tilde);
   }

   //return List::create(Named("Q") = Q_tilde_max);
   l = ceil(log2(n));
   WT Tree = WT(Y);

   //First search
   imax = 0;
   lb_max = R_NegInf;
   ub_max = R_PosInf;
   for (int i = 0; i< n; ++i) {
     //Rcout << "i = " << i;
     lb = R_NegInf;
     ub = R_PosInf;
     power = floor(log2(i));
     for (int p = p_min; p<= power; p++) {
       l = pow(2,p);
       for (int k = lag; k <= i-l+1; k = k+ 2) {
         //Rcout << "k = " << k << "\n";
         idx = fmin(fmax(1,(int)ceil(l*beta)),l);
         aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k,k+l-1)]];
         //aux_Y = get_Quantile(idx,k,k+l-1,R,lr,RANK);
         q_tilde = pow(penfsC(i,l,q),2)/(2*l);;
         if(q_tilde<=Q_max-1e-6){
           lq = x1(q_tilde,beta);
           lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
           if(beta==0.5){
             uq_idx = l-lq_idx;
           }else{
             uq = x2(q_tilde,beta);
             uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
           }
           aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k,k+l-1)]];
           //aux_lb = get_Quantile(lq_idx,k,k+l-1,R,lr,RANK);
           if(aux_lb>lb){
             lb = aux_lb;
             lk = k;
             llen = l;
             lp = p;
             llen_idx = Q_tilde_max_idx[llen-1];
           }
           aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k,k+l-1)]];
           //aux_ub = get_Quantile(uq_idx,k,k+l-1,R,lr,RANK);
           if(aux_ub<ub){
             ub = aux_ub;
             uk = k;
             ulen = l;
             up = p;
             ulen_idx = Q_tilde_max_idx[ulen-1];
           }
         }
       }
     }
     if(lb <= ub){
       lb_first = lb;
       ub_first = ub;
       lb_last = lb;
       ub_last = ub;
       if(aux_Y > ub){
         mu[i] = ub;
       }else if(aux_Y < lb){
         mu[i] = lb;
       }else{
         mu[i] = aux_Y;
       }
       C[i] = sum((Y[Range(0,i)]-mu[i])*(beta-ifelse(Y[Range(0,i)]<mu[i],1.0,0.0)));
       // Rcout << "Cost at "<< i << "= "<<  C[i] << "\n";
       // Rcout << "C[i]" << i << "=" << C[i] << "\n";
       L[i] = 1;
       J[i] = 0;
     }else{
       q_tilde = Q_tilde_max[llen-1];
       // Rcout << "llen = " << llen << " q_tilde = " << q_tilde << "\n";
       if(q_tilde<=Q_max-1e-6){
         lq_max = x1(q_tilde,beta);
         lq_max_idx = fmin(fmax(1,ceil(llen*lq_max)),llen);
         lb_max = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx,lk,lk+llen-1)]];
       }
       //Rcout << "lq_max = " << lq_max << " ";
       //Rcout << "lq_max_idx = " << lq_max_idx << " ";
       //Rcout << "lb_max = " << lb_max << " ";

       q_tilde = Q_tilde_max[ulen-1];
       if(q_tilde<=Q_max-1e-6){
         uq_max = x2(q_tilde,beta);
         uq_max_idx = fmax(fmin(ulen,floor(ulen*uq_max)),1);
         ub_max = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx,uk,uk+ulen-1)]];
       }
       //Rcout << "ub_max = " << ub_max << "\n";
       if(lb_max>ub_max){
         break;
       }
     }
   }
   if(test == true && std::isinf(C[n-1])){
     if(details == true){
       Rcout << "Fast test is done! No alternative solution!\n";
     }
     IntegerVector left(2);
     left[0] = 1;
     left[1] = 1;
     return List::create(Named("left") = left);
   }

   // Later searchs
   int imin = 0;
   int njmp = 0;
   int jmin, aux_imin, aux_imax,i_llen_idx,j_llen_idx,i_ulen_idx,j_ulen_idx,lq_max_idx_i,uq_max_idx_i,lq_max_idx_j;
   bool isstop = FALSE;
   double lq_max_i,lb_max_i,uq_max_i, ub_max_i,lq_max_j,lb_max_j,uq_max_j,ub_max_j, uq_max_idx_j;

   while (std::isinf(C[n-1]) && njmp < n) {
     jmin = n;
     aux_imin = n;
     aux_imax = 0;
     lb_max_i = R_NegInf;
     ub_max_i = R_PosInf;
     for (int i=imin; i<n;i++) {
       // Rcout << "i = " << i << "\n";
       if(J[i] == njmp){
         jmin = fmin(jmin,i);
       }else if(std::isinf(C[i])){
         isstop = FALSE;
         lb_max_j = R_NegInf;
         ub_max_j = R_PosInf;
         for (int j=i-1; j>=jmin;j--) {
           if(J[j] == njmp){
             lb = R_NegInf;
             ub = R_PosInf;
             power = floor(log2(i-j));
             for (int p= p_min; p <= power; p++) {
               l = pow(2,p);
               for (int k=j+lag; k<= i-l; k = k+2) {
                 idx = fmin(fmax(1,(int)ceil(l*beta)),l);
                 aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k+1,k+l)]];
                 //aux_Y = get_Quantile(idx,k+1,k+l,R,lr,RANK);
                 q_tilde = pow(penfsC(i-j-1,l,q),2)/(2*l);;
                 if(q_tilde<=Q_max-1e-6){
                   lq = x1(q_tilde,beta);
                   lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
                   if(beta==0.5){
                     uq_idx = l-lq_idx;
                   }else{
                     uq = x2(q_tilde,beta);
                     uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
                   }
                   aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k+1,k+l)]];
                   //aux_lb = get_Quantile(lq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_lb>lb){
                     lb = aux_lb;
                     lk = k;
                     llen = l;
                     lp = p;
                     i_llen_idx = n - imin - 1;
                     j_llen_idx = i - jmin - 1;
                   }
                   aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k+1,k+l)]];
                   //aux_ub = get_Quantile(uq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_ub<ub){
                     ub = aux_ub;
                     uk = k;
                     ulen = l;
                     up = p;
                     i_ulen_idx = n - imin - 1;
                     j_ulen_idx = i - jmin - 1;
                   }
                 }
               }
             }
             if(lb<=ub){
               lb_last = lb;
               ub_last = ub;
               if(aux_Y > ub){
                 aux_mu = ub;
               }else if(aux_Y < lb){
                 aux_mu = lb;
               }else{
                 aux_mu = aux_Y;
               }
               if(i < aux_imin){
                 aux_imin = i;
               }
               if(i > aux_imax){
                 aux_imax = i;
               }
               aux_c = sum((Y[Range(j+1,i)]-aux_mu)*(beta-ifelse(Y[Range(j+1,i)]<aux_mu,1.0,0.0)))+C[j];
               if(aux_c < C[i]){
                 C[i] = aux_c;
                 mu[i] = aux_mu;
                 L[i] = j+2;
                 J[i] = njmp + 1;
                 // Rcout << "C[i]" << i << "=" << C[i] << "\n";
               }
             }else{
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_i = x1(q_tilde,beta);
                 lq_max_idx_i = fmin(fmax(1,(int)ceil(llen*lq_max_i)),llen);
                 lb_max_i = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_i,lk+1,lk+llen)]];
                 //lb_max_i = get_Quantile(lq_max_idx_i,lk+1,lk+llen,R,lr,RANK);
               }

               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_i = x2(q_tilde,beta);
                 uq_max_idx_i = fmax(fmin(ulen,(int)floor(ulen*uq_max_i)),1);
                 ub_max_i = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_i,uk+1,uk+ulen)]];
                 //ub_max_i = get_Quantile(uq_max_idx_i,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_i>ub_max_i && j >= imax){
                 isstop = TRUE;
               }
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_j = x1(q_tilde,beta);
                 lq_max_idx_j = fmin(fmax(1,(int)ceil(llen*lq_max_j)),llen);
                 lb_max_j = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_j,lk+1,lk+llen)]];
                 //lb_max_j = get_Quantile(lq_max_idx_j,lk+1,lk+llen,R,lr,RANK);
               }
               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_j = x2(q_tilde,beta);
                 uq_max_idx_j = fmax(fmin(ulen,(int)floor(ulen*uq_max_j)),1);
                 ub_max_j = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_j,uk+1,uk+ulen)]];
                 //ub_max_j = get_Quantile(uq_max_idx_j,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_j>ub_max_j){
                 break;
               }
             }
           }
         }
         if(std::isinf(C[i])&&isstop){
           break;
         }
       }
     }
     imin = aux_imin;
     imax = aux_imax;
     njmp = njmp + 1;
     // Rcout << "Jump from" << njmp-1 << "to" << njmp << "\n";
   }

   NumericVector value(njmp+1);
   IntegerVector left(njmp+1);
   int cnt = n-1;
   for (int i = njmp; i >= 0; --i) {
     left[i] = L[cnt];
     value[i] = mu[cnt];
     cnt = L[cnt] - 2;
   }
   NumericVector first(2);
   NumericVector last(2);
   first[0] = lb_first;
   first[1] = ub_first;
   last[0] = lb_last;
   last[1] = ub_last;

   return List::create(Named("value") = value, Named("left") = left,
                       Named("n") = n, Named("first") = first,
                       Named("last") = last);
 }


// [[Rcpp::export(.MMUSCLE)]]
List MMUSCLE(NumericVector &Y, NumericMatrix &q_matrix, NumericVector & beta_vec, bool test, bool details){
   //Rcout << "MUSCLE is running, please wait!\n";
   int num_beta = beta_vec.size();
   int n = Y.size();
   int imax,lq_max_idx, uq_max_idx,idx,power,l;
   double lb_max, ub_max, q_tilde,aux_lb, aux_ub, lq_max,uq_max;
   NumericVector aux_Y(num_beta);
   NumericVector lb_first(num_beta);
   NumericVector ub_first(num_beta);
   NumericVector lb_last(num_beta);
   NumericVector ub_last(num_beta);
   bool no_empty = true;
   // bool no_solution = false;
   bool no_solution_right = false;
   // bool no_solution_left = false;
   NumericVector lb(num_beta);
   NumericVector lp(num_beta);
   NumericVector lk(num_beta);
   NumericVector lq(num_beta);
   NumericVector llen(num_beta);
   NumericVector llen_idx(num_beta);
   NumericVector lq_idx(num_beta);

   NumericVector ub(num_beta);
   NumericVector up(num_beta);
   NumericVector uk(num_beta);
   NumericVector uq(num_beta);
   NumericVector ulen(num_beta);
   NumericVector ulen_idx(num_beta);
   NumericVector uq_idx(num_beta);

   NumericMatrix C(num_beta,n);  //optimal cost for each quantile
   NumericVector Cost(n); // optimal cost, equals to mean of C
   NumericVector J(n);  //number of jumps
   NumericMatrix mu(num_beta,n); //optimal value of step function
   NumericVector L(n);  //leftmost end points
   for (int i = 0; i < n; ++i){
     for(int j = 0; j < num_beta;++j){
       C(j,i) = INFINITY;
     }
     J[i] = n;
     L[i] = 0;
     Cost[i] = INFINITY;
   }
   NumericVector Q_max(num_beta);
   NumericVector l_min(num_beta);
   for(int i = 0;i < num_beta; ++i){
     Q_max[i] = fmax(-log(beta_vec[i]),-log(1-beta_vec[i]));
     l_min[i] = floor(1/fmin(beta_vec[i],1-beta_vec[i]));
     // Rcout << "Q_max[i] = " << Q_max[i] << "\n";
     // Rcout << "l_min[i] = " << l_min[i] << "\n";
   }
   NumericMatrix Q_tilde(num_beta,n);
   NumericMatrix Q_tilde_max(num_beta,n);
   IntegerMatrix Q_tilde_max_idx(num_beta,n);
   bool get_p_min = false;
   int p_min = 1;
   // Rcout << "p_min = " << p_min << "\n";
   for(int i = 0; i<n; i++){
     l = i + 1;
     for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
       for(int j = 0; j<n; j++){
         if(j==i){
           // Q_tilde(idx_beta,j) = pow(penfsC(j,l,q(idx_beta,_)),2)/(2*l);
           Q_tilde(idx_beta,j) = pow(penfsCC(j,l,q_matrix(idx_beta,j)),2)/(2*l);
           if(get_p_min == false && Q_tilde(idx_beta,j)< Q_max[i]-1e-6){
             get_p_min = true;
             p_min = floor(log2(l));
           }
         }else if(j>i){
           Q_tilde(idx_beta,j) = pow(penfsCC(j,l,q_matrix(idx_beta,j)),2)/(2*l);
         }else{
           Q_tilde(idx_beta,j) = R_NegInf;
         }
       }
       Q_tilde_max(idx_beta,i) = max(Q_tilde(idx_beta,_));
       Q_tilde_max_idx(idx_beta,i) = which_max(Q_tilde(idx_beta,_));
     }
   }

   l = ceil(log2(n));
   WT Tree = WT(Y);
   // return List::create(Named("n") = n);

   //First search
   imax = 0;
   lb_max = R_NegInf;
   ub_max = R_PosInf;
   for (int i = 0; i< n; ++i) {
     for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
       lb[idx_beta] = R_NegInf;
       ub[idx_beta] = R_PosInf;
     }
     power = fmax(floor(log2(i)),0);
     // Rcout << "i = " << i << "\n";
     // Rcout << "power = " << power << "\n";
     // Rcout << "p_min = " << p_min << "\n";
     for (int p = p_min; p <= power; p++) {
       l = pow(2,p);
       for (int k = 0; k <= i-l+1; k = k+ 2) {
         // Rcout << "k = " << k << "\n";
         for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
           // Rcout << "idx_beta = " << idx_beta << "\n";
           // Rcout << "b[idx_beta] = " << b[idx_beta] << "\n";
           idx = fmin(fmax(1,(int)ceil(l*beta_vec[idx_beta])),l);
           aux_Y[idx_beta] = Y[Tree.IDX[Tree.rngQuantile(idx,k,k+l-1)]];
           // Rcout << "aux_Y[idx_beta] = " << aux_Y[idx_beta] << "\n";
           if(l >= l_min[idx_beta]){
             // Rcout << "l = " << l << "\n";
             // q_tilde = pow(penfsC(i,l,q(idx_beta,_)),2)/(2*l);
             q_tilde = pow(penfsCC(i,l,q_matrix(idx_beta,i)),2)/(2*l);
             // Rcout << "q_tilde = " << q_tilde << "\n";
             // Rcout << "Q_max[idx_beta] =" << Q_max[idx_beta] << "\n";
             if(q_tilde<=Q_max[idx_beta]-1e-6){
               lq[idx_beta] = x1(q_tilde,beta_vec[idx_beta]);
               // Rcout << "lq[idx_beta] = " << lq[idx_beta] << "\n";
               lq_idx[idx_beta] = fmin(fmax(1,(int)ceil(l*lq[idx_beta])),l);
               // Rcout << "lq_idx[idx_beta] = " << lq_idx[idx_beta] << "\n";
               if(beta_vec[idx_beta]==0.5){
                 uq_idx[idx_beta] = l-lq_idx[idx_beta];
               }else{
                 uq[idx_beta] = x2(q_tilde,beta_vec[idx_beta]);
                 uq_idx[idx_beta] = fmax(fmin(l,(int)floor(l*uq[idx_beta])),1);
                 // Rcout << "uq[idx_beta] = " << uq[idx_beta] << "\n";
                 // Rcout << "uq_idx[idx_beta] = " << uq_idx[idx_beta] << "\n";
               }
               if(lq[idx_beta] > 0){
                 aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx[idx_beta],k,k+l-1)]];
               }else{
                 aux_lb = R_NegInf;
               }
               if(aux_lb>lb[idx_beta]){
                 lb[idx_beta] = aux_lb;
                 lk[idx_beta] = k;
                 llen[idx_beta] = l;
                 lp[idx_beta] = p;
                 llen_idx[idx_beta] = Q_tilde_max_idx(idx_beta,llen[idx_beta]-1);
               }

               // Rcout << "lb[idx_beta] =" << lb[idx_beta] << "\n";
               if(uq[idx_beta] < 1){
                 aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx[idx_beta],k,k+l-1)]];
               }else{
                aux_ub = R_PosInf;
               }
               //aux_ub = get_Quantile(uq_idx,k,k+l-1,R,lr,RANK);
               if(aux_ub<ub[idx_beta]){
                 ub[idx_beta] = aux_ub;
                 uk[idx_beta] = k;
                 ulen[idx_beta] = l;
                 up[idx_beta] = p;
                 ulen_idx[idx_beta] = Q_tilde_max_idx(idx_beta,ulen[idx_beta]-1);
               }
               // Rcout << "ub[idx_beta] =" << ub[idx_beta] << "\n";
             }
           }
         }
       }
     }
     for(int idx_beta = 0; idx_beta<num_beta; idx_beta++){
       if(ub[idx_beta]<lb[idx_beta]){
         no_empty = false;
         // Rcout << "Empty \n";
         // Rcout << "idx_beta = " << idx_beta << "\n";
         // Rcout << "ub[idx_beta]" << ub[idx_beta] << "\n";
         // Rcout << "lb[idx_beta]" << lb[idx_beta] << "\n";
         break;
       }
     }
     if(no_empty == true){
       for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
         lb_first[idx_beta] = lb[idx_beta];
         ub_first[idx_beta] = ub[idx_beta];
         lb_last[idx_beta] = lb[idx_beta];
         ub_last[idx_beta] = ub[idx_beta];
         // Rcout << "aux_Y[idx_beta] = " << aux_Y[idx_beta] <<"\n";
         // Rcout << "lb[idx_beta] = " << lb[idx_beta] <<"\n";
         // Rcout << "ub[idx_beta] = " << ub[idx_beta] <<"\n";
         if(aux_Y[idx_beta] > ub[idx_beta]){
           mu(idx_beta,i) = ub[idx_beta];
         }else if(aux_Y[idx_beta] < lb[idx_beta]){
           mu(idx_beta,i) = lb[idx_beta];
         }else{
           mu(idx_beta,i) = aux_Y[idx_beta];
         }
         // Rcout << "No empty \n";
         // Rcout << "idx_beta = " << idx_beta << "\n";
         // Rcout << "mu[idx_beta,i] = " << mu(idx_beta,i) << "\n";
         C(idx_beta,i) = sum((Y[Range(0,i)]-mu(idx_beta,i))*(beta_vec[idx_beta]-ifelse(Y[Range(0,i)]<mu(idx_beta,i),1.0,0.0)));
         // Rcout << "C[idx_beta,i] = " << C(idx_beta,i) << "\n";
       }
       Cost[i] = mean(C(_,i));
       L[i] = 1;
       J[i] = 0;
       // Rcout << "Cost[i]" << i << "=" << Cost[i] << "\n";
     }else{
       no_solution_right = false;
       for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
         q_tilde = Q_tilde_max(idx_beta, llen[idx_beta]-1);
         // Rcout << "llen[idx_beta] = " << llen[idx_beta] << "\n";
         if(q_tilde<=Q_max[idx_beta]-1e-6){
           lq_max = x1(q_tilde,beta_vec[idx_beta]);
           lq_max_idx = fmin(fmax(1,ceil(llen[idx_beta]*lq_max)),llen[idx_beta]);
           if(lq_max>0){
             lb_max = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx,lk[idx_beta],lk[idx_beta]+llen[idx_beta]-1)]];
           }
         }

         q_tilde = Q_tilde_max(idx_beta,ulen[idx_beta]-1);
         // Rcout << "ulen[idx_beta] = " << ulen[idx_beta] << "\n";

         if(q_tilde<=Q_max[idx_beta]-1e-6){
           uq_max = x2(q_tilde,beta_vec[idx_beta]);
           uq_max_idx = fmax(fmin(ulen[idx_beta],floor(ulen[idx_beta]*uq_max)),1);
           if(uq_max<1){
            ub_max = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx,uk[idx_beta],uk[idx_beta]+ulen[idx_beta]-1)]];
           }
         }
         // Rcout << "lb_max = " << lb_max << "\n";
         // Rcout << "ub_max = " << ub_max << "\n";
         if(lb_max>ub_max){
           no_solution_right = true;
           // Rcout << "i = " << i << "\n";
           // Rcout << "No solution right!\n";
           break;
         }
       }

       if(no_solution_right == true){
         break;
       }
     }
   }
   if(test == true && std::isinf(Cost[n-1])){
     if(details == true){
       Rcout << "Fast test is done! No alternative solution!\n";
     }
     IntegerVector left(2);
     left[0] = 1;
     left[1] = 1;
     return List::create(Named("left") = left);
   }
   // return List::create(Named("n") = n);

   // Rcout << "Later searchs!\n";
   // Later searchs
   int imin = 0;
   int njmp = 0;
   int jmin, aux_imin, aux_imax,lq_max_idx_i,uq_max_idx_i,lq_max_idx_j;
   bool isstop = false;
   bool isstop_j = false;
   double lq_max_i,lb_max_i,uq_max_i, ub_max_i,lq_max_j,lb_max_j,uq_max_j,ub_max_j, uq_max_idx_j;
   NumericVector i_llen_idx(num_beta);
   NumericVector j_llen_idx(num_beta);
   NumericVector i_ulen_idx(num_beta);
   NumericVector j_ulen_idx(num_beta);
   NumericVector aux_mu(num_beta);
   double aux_c = 0;
   // Rcout << "Cost[n-1] = " << Cost[n-1] << "\n";
   while (std::isinf(Cost[n-1]) && njmp < n) {
     jmin = n;
     aux_imin = n;
     aux_imax = 0;
     lb_max_i = R_NegInf;
     ub_max_i = R_PosInf;
     for (int i=imin; i<n;i++) {
       no_empty = true;
       // Rcout << "i = " << i << "\n";
       if(J[i] == njmp){
         jmin = fmin(jmin,i);
       }else if(std::isinf(Cost[i])){
         isstop = FALSE;
         lb_max_j = R_NegInf;
         ub_max_j = R_PosInf;
         for (int j=i-1; j>=jmin;j--) {
           // Rcout << "j = " << j << "\n";
           if(J[j] == njmp){
             for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
               lb[idx_beta] = R_NegInf;
               ub[idx_beta] = R_PosInf;
             }
             power = floor(log2(i-j));
             for (int p = power; p >= p_min;p--) {
               l = pow(2,p);
               // Rcout << "l = " << l << "\n";
               for (int k=j; k<= i-l; k = k+2) {
                 // Rcout << "k = " << k << "\n";
                 for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                   idx = fmin(fmax(1,(int)ceil(l*beta_vec[idx_beta])),l);
                   // Rcout << "idx = "<< idx << "\n";
                   aux_Y[idx_beta] = Y[Tree.IDX[Tree.rngQuantile(idx,k+1,k+l)]];
                   // Rcout << "aux_Y[idx_beta] = " << aux_Y[idx_beta] << "\n";
                   // aux_Y = get_Quantile(idx,k+1,k+l,R,lr,RANK);
                   if(l >= l_min[idx_beta]){
                     // q_tilde = pow(penfsC(i-j-1,l,q(idx_beta,_)),2)/(2*l);
                     q_tilde = pow(penfsCC(i-j-1,l,q_matrix(idx_beta,i-j-1)),2)/(2*l);
                     if(q_tilde<=Q_max[idx_beta]-1e-6){
                       lq[idx_beta] = x1(q_tilde,beta_vec[idx_beta]);
                       lq_idx[idx_beta] = fmin(fmax(1,(int)ceil(l*lq[idx_beta])),l);
                       if(beta_vec[idx_beta]==0.5){
                         uq_idx[idx_beta] = l-lq_idx[idx_beta];
                       }else{
                         uq[idx_beta] = x2(q_tilde,beta_vec[idx_beta]);
                         uq_idx[idx_beta] = fmax(fmin(l,(int)floor(l*uq[idx_beta])),1);
                       }
                       if(lq[idx_beta] >0){
                         aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx[idx_beta],k+1,k+l)]];
                       }else{
                         aux_lb = R_NegInf;
                       }
                       // aux_lb = get_Quantile(lq_idx,k+1,k+l,R,lr,RANK);
                       if(aux_lb>lb[idx_beta]){
                         lb[idx_beta] = aux_lb;
                         lk[idx_beta] = k;
                         llen[idx_beta] = l;
                         lp[idx_beta] = p;
                         i_llen_idx[idx_beta] = n - imin - 1;
                         j_llen_idx[idx_beta] = i - jmin - 1;
                       }
                       if(uq[idx_beta]<1){
                         aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx[idx_beta],k+1,k+l)]];
                       }else{
                         aux_ub = R_PosInf;
                       }
                       //aux_ub = get_Quantile(uq_idx,k+1,k+l,R,lr,RANK);
                       if(aux_ub<ub[idx_beta]){
                         ub[idx_beta] = aux_ub;
                         uk[idx_beta] = k;
                         ulen[idx_beta] = l;
                         up[idx_beta] = p;
                         i_ulen_idx[idx_beta] = n - imin - 1;
                         j_ulen_idx[idx_beta] = i - jmin - 1;
                       }
                     }
                   }
                 }
               }
             }
             for(int idx_beta = 0; idx_beta<num_beta; idx_beta++){
               if(ub[idx_beta]<lb[idx_beta]){
                 no_empty = false;
                 // Rcout << "Empty \n";
                 // Rcout << "idx_beta = " << idx_beta << "\n";
                 // Rcout << "ub[idx_beta]" << ub[idx_beta] << "\n";
                 // Rcout << "lb[idx_beta]" << lb[idx_beta] << "\n";
                 break;
               }
             }
             if(no_empty == true){
               for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                 lb_first[idx_beta] = lb[idx_beta];
                 ub_first[idx_beta] = ub[idx_beta];
                 lb_last[idx_beta] = lb[idx_beta];
                 ub_last[idx_beta] = ub[idx_beta];
                 // Rcout << "aux_Y[idx_beta] = " << aux_Y[idx_beta] <<"\n";
                 // Rcout << "lb[idx_beta] = " << lb[idx_beta] <<"\n";
                 // Rcout << "ub[idx_beta] = " << ub[idx_beta] <<"\n";
                 if(aux_Y[idx_beta] > ub[idx_beta]){
                   aux_mu[idx_beta] = ub[idx_beta];
                 }else if(aux_Y[idx_beta] < lb[idx_beta]){
                   aux_mu[idx_beta] = lb[idx_beta];
                 }else{
                   aux_mu[idx_beta] = aux_Y[idx_beta];
                 }
                 // Rcout << "No empty \n";
                 // Rcout << "idx_beta = " << idx_beta << "\n";
                 // Rcout << "aux_mu[idx_beta,i] = " << aux_mu[idx_beta] << "\n";
                 // C(idx_beta,i) = sum((Y[Range(j+1,i)]-aux_mu[idx_beta])*(b-ifelse(Y[Range(j+1,i)]<aux_mu[idx_beta],1.0,0.0)));
                 // Rcout << "C[idx_beta,i] = " << C(idx_beta,i) << "\n";
               }
               if(i < aux_imin){
                 aux_imin = i;
               }
               if(i > aux_imax){
                 aux_imax = i;
               }
               aux_c = 0;
               for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                  aux_c = aux_c + sum((Y[Range(j+1,i)]-aux_mu[idx_beta])*(beta_vec[idx_beta]-ifelse(Y[Range(j+1,i)]<aux_mu[idx_beta],1.0,0.0)));
               }
               // Rcout << "cost from j+1 to i = " << aux_c/num_beta << "\n";
               aux_c = aux_c/num_beta + Cost[j];
               // Rcout << "total cost from 1 to i = " << aux_c << "\n";
               if(aux_c < Cost[i]){
                 for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                   mu(idx_beta,i) = aux_mu[idx_beta];
                 }
                 Cost[i] = aux_c;
                 L[i] = j+2;
                 J[i] = njmp + 1;
                 // Rcout << "Cost[i]" << i << "=" << Cost[i] << "\n";
               }
             }else{
               for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                 q_tilde = Q_tilde_max[llen[idx_beta]-1];
                 if(q_tilde<=Q_max[idx_beta]-1e-6){
                   lq_max_i = x1(q_tilde,beta_vec[idx_beta]);
                   lq_max_idx_i = fmin(fmax(1,(int)ceil(llen[idx_beta]*lq_max_i)),llen[idx_beta]);
                   if(lq_max_i>0){
                     lb_max_i = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_i,lk[idx_beta]+1,lk[idx_beta]+llen[idx_beta])]];
                     // lb_max_i = get_Quantile(lq_max_idx_i,lk+1,lk+llen,R,lr,RANK);
                   }
                 }

                 q_tilde = Q_tilde_max[ulen[idx_beta]-1];
                 if(q_tilde<=Q_max[idx_beta]-1e-6){
                   uq_max_i = x2(q_tilde,beta_vec[idx_beta]);
                   uq_max_idx_i = fmax(fmin(ulen[idx_beta],(int)floor(ulen[idx_beta]*uq_max_i)),1);
                   if(uq_max_i<1){
                     ub_max_i = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_i,uk[idx_beta]+1,uk[idx_beta]+ulen[idx_beta])]];
                     //ub_max_i = get_Quantile(uq_max_idx_i,uk+1,uk+ulen,R,lr,RANK);
                   }
                 }
                 if(lb_max_i>ub_max_i && j >= imax){
                   isstop = TRUE;
                   break;
                 }
               }

               for(int idx_beta = 0; idx_beta < num_beta; idx_beta++){
                 q_tilde = Q_tilde_max[llen[idx_beta]-1];
                 if(q_tilde<=Q_max[idx_beta]-1e-6){
                   lq_max_j = x1(q_tilde,beta_vec[idx_beta]);
                   lq_max_idx_j = fmin(fmax(1,(int)ceil(llen[idx_beta]*lq_max_j)),llen[idx_beta]);
                   if(lq_max_j>0){
                     lb_max_j = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_j,lk[idx_beta]+1,lk[idx_beta]+llen[idx_beta])]];
                     //lb_max_j = get_Quantile(lq_max_idx_j,lk+1,lk+llen,R,lr,RANK);
                   }
                 }
                 q_tilde = Q_tilde_max[ulen[idx_beta]-1];
                 if(q_tilde<=Q_max[idx_beta]-1e-6){
                   uq_max_j = x2(q_tilde,beta_vec[idx_beta]);
                   uq_max_idx_j = fmax(fmin(ulen[idx_beta],(int)floor(ulen[idx_beta]*uq_max_j)),1);
                   if(uq_max_j<1){
                     ub_max_j = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_j,uk[idx_beta]+1,uk[idx_beta]+ulen[idx_beta])]];
                     //ub_max_j = get_Quantile(uq_max_idx_j,uk+1,uk+ulen,R,lr,RANK);
                   }
                 }
                 if(lb_max_j>ub_max_j){
                   isstop_j = true;
                   break;
                 }
                 if(isstop_j == true){
                   break;
                 }
               }
             }
           }
         }
         if(std::isinf(Cost[i])&&isstop){
           break;
         }
       }
     }
     imin = aux_imin;
     imax = aux_imax;
     njmp = njmp + 1;
     // Rcout << "Jump from" << njmp-1 << "to" << njmp << "\n";
   }

   NumericMatrix value(num_beta,njmp+1);
   IntegerVector left(njmp+1);
   int cnt = n-1;
   for (int i = njmp; i >= 0; --i) {
     left[i] = L[cnt];
     for(int idx_beta = 0;idx_beta<num_beta;idx_beta++){
       value(idx_beta,i) = mu(idx_beta,cnt);
     }
     cnt = L[cnt] - 2;
   }
   NumericMatrix first(num_beta,2);
   NumericMatrix last(num_beta,2);
   for(int idx_beta = 0;idx_beta<num_beta;idx_beta++){
     first(idx_beta,0) = lb_first[idx_beta];
     first(idx_beta,1) = ub_first[idx_beta];
     last(idx_beta,0) = lb_last[idx_beta];
     last(idx_beta,1) = ub_last[idx_beta];
   }
   return List::create(Named("value") = value, Named("left") = left,
                       Named("n") = n, Named("first") = first,
                       Named("last") = last, Named("Mu") = mu, Named("J")= J,Named("L") = L);
 }


// [[Rcpp::export(.MUSCLE_FULL)]]
List MUSCLE_FULL(NumericVector &Y, NumericVector &q, double beta, bool test, bool details){
    //Rcout << "MUSCLE FULL is running, please wait!\n";
    int n = Y.size();
    int imax,lq_idx, uq_idx,lk,uk, llen, ulen, llen_idx, ulen_idx,lq_max_idx, uq_max_idx,idx,l;
    double lb_max, ub_max, lb, ub, q_tilde, aux_Y, aux_mu, aux_c,aux_lb, aux_ub,lq,uq, lq_max,uq_max;
    double lb_first, ub_first, lb_last, ub_last;
    NumericVector C(n);  //optimal cost
    NumericVector J(n);  //number of jumps
    NumericVector mu(n); //optimal value of step function
    NumericVector L(n);  //leftmost end points
    for (int i = 0; i < n; ++i){
     C[i] = INFINITY;
     J[i] = n;
     L[i] = 0;
    }
    double Q_max = fmax(-log(beta),-log(1-beta)); //maximum of f on [0,1]
    NumericVector Q_tilde(n);
    NumericVector Q_tilde_max(n);
    IntegerVector Q_tilde_max_idx(n);
    bool get_l_min = false;
    int l_min;
    for(int i = 0; i<n; i++){
     l = i + 1;
     for(int j = 0; j<n; j++){
       if(j==i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
         if(get_l_min == false && Q_tilde[j] <= Q_max -1e-6){
           get_l_min = true;
           l_min = l;
         }
       }
       if(j>=i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
       }else{
         Q_tilde[j] = R_NegInf;
       }
     }
     Q_tilde_max[i] = max(Q_tilde);
     Q_tilde_max_idx[i] = which_max(Q_tilde);
    }

    //return List::create(Named("value") = Q_tilde_max);
    //int k_step = 1;
    int k_step;
    if(n<= 500){
     k_step = 2;
    }else if(n<=1000){
     k_step = floor(log(n)/1.5);
    }else{
     k_step = floor(log(n)*2);
    }
    //Rcout << "k_step = " << k_step << "\n";
    // for(int j = 0;j<n;j++){
    //   if(Q_tilde(j,j)<=Q_max-1e-6){
    //     l_min = j+1;
    //     break;
    //   }
    // }

    WT Tree = WT(Y);

    //First search
    imax = 0;
    lb_max = R_NegInf;
    ub_max = R_PosInf;
    for (int i = 0; i< n; ++i) {
     //Rcout << i << "\n";
     lb = R_NegInf;
     ub = R_PosInf;
     for (int l = i+1; l >=l_min; l--) {
       for (int k = 0; k <= i-l+1; k = k + fmax(floor(l/2),1)) {
         //fmax(floor(l/2),1)
         idx = fmin(fmax(1,(int)ceil(l*beta)),l);
         aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k,k+l-1)]];
         //aux_Y = get_Quantile(idx,k,k+l-1,R,lr,RANK);
         q_tilde = pow(penfsC(i,l,q),2)/(2*l);
         if(q_tilde<=Q_max-1e-6){
           lq = x1(q_tilde,beta);
           lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
           if(beta==0.5){
             uq_idx = l-lq_idx;
           }else{
             uq = x2(q_tilde,beta);
             uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
           }
           aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k,k+l-1)]];
           //aux_lb = get_Quantile(lq_idx,k,k+l-1,R,lr,RANK);
           if(aux_lb>lb){
             lb = aux_lb;
             lk = k;
             llen = l;
             llen_idx = Q_tilde_max_idx[llen-1];
           }
           aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k,k+l-1)]];
           //aux_ub = get_Quantile(uq_idx,k,k+l-1,R,lr,RANK);
           if(aux_ub<ub){
             ub = aux_ub;
             uk = k;
             ulen = l;
             ulen_idx = Q_tilde_max_idx[ulen-1];
           }
         }
       }
     }
     if(lb <= ub){
       lb_first = lb;
       ub_first = ub;
       lb_last = lb;
       ub_last = ub;
       if(aux_Y > ub){
         mu[i] = ub;
       }else if(aux_Y < lb){
         mu[i] = lb;
       }else{
         mu[i] = aux_Y;
       }
       C[i] = sum((Y[Range(0,i)]-mu[i])*(beta-ifelse(Y[Range(0,i)]<mu[i],1.0,0.0)));
       L[i] = 1;
       J[i] = 0;
     }else{
       q_tilde = Q_tilde_max[llen-1];
       if(q_tilde<=Q_max-1e-6){
         lq_max = x1(q_tilde,beta);
         lq_max_idx = fmin(fmax(1,ceil(llen*lq_max)),llen);
         lb_max = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx,lk,lk+llen-1)]];
       }

       q_tilde = Q_tilde_max[ulen-1];
       if(q_tilde<=Q_max-1e-6){
         uq_max = x2(q_tilde,beta);
         uq_max_idx = fmax(fmin(ulen,floor(ulen*uq_max)),1);
         ub_max = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx,uk,uk+ulen-1)]];
       }
       if(lb_max>ub_max){
         break;
       }
     }
    }
    if(test == true && std::isinf(C[n-1])){
      if(details == true){
        Rcout << "Fast test is done! No alternative solution!\n";
      }
      IntegerVector left(2);
      left[0] = 1;
      left[1] = 1;
      return List::create(Named("left") = left);
    }

    // Later searchs
    int imin = 0;
    int njmp = 0;
    int jmin, aux_imin, aux_imax,i_llen_idx,j_llen_idx,i_ulen_idx,j_ulen_idx,lq_max_idx_i,uq_max_idx_i,lq_max_idx_j;
    bool isstop = FALSE;
    double lq_max_i,lb_max_i,uq_max_i, ub_max_i,lq_max_j,lb_max_j,uq_max_j,ub_max_j, uq_max_idx_j;

    while (std::isinf(C[n-1]) && njmp < n) {
     jmin = n;
     aux_imin = n;
     aux_imax = 0;
     lb_max_i = R_NegInf;
     ub_max_i = R_PosInf;
     for (int i=imin; i<n;i++) {
       //Rcout << i << "\n";
       if(J[i] == njmp){
         jmin = fmin(jmin,i);
       }else if(std::isinf(C[i])){
         isstop = FALSE;
         lb_max_j = R_NegInf;
         ub_max_j = R_PosInf;
         for (int j=i-1; j>=jmin;j--) {
           if(J[j] == njmp){
             lb = R_NegInf;
             ub = R_PosInf;
             for (int l=i-j; l >=l_min;l--) {
               for (int k=j; k<= i-l; k = k + fmax(floor(l/2),1)) {
                 //fmax(floor(l/2),1)
                 idx = fmin(fmax(1,(int)ceil(l*beta)),l);
                 aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k+1,k+l)]];
                 //aux_Y = get_Quantile(idx,k+1,k+l,R,lr,RANK);
                 q_tilde = pow(penfsC(i-j-1,l,q),2)/(2*l);
                 if(q_tilde<=Q_max-1e-6){
                   lq = x1(q_tilde,beta);
                   lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
                   if(beta==0.5){
                     uq_idx = l-lq_idx;
                   }else{
                     uq = x2(q_tilde,beta);
                     uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
                   }
                   aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k+1,k+l)]];
                   //aux_lb = get_Quantile(lq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_lb>lb){
                     lb = aux_lb;
                     lk = k;
                     llen = l;
                     i_llen_idx = n - imin - 1;
                     j_llen_idx = i - jmin - 1;
                   }
                   aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k+1,k+l)]];
                   //aux_ub = get_Quantile(uq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_ub<ub){
                     ub = aux_ub;
                     uk = k;
                     ulen = l;
                     i_ulen_idx = n - imin - 1;
                     j_ulen_idx = i - jmin - 1;
                   }
                 }
               }
             }
             if(lb<=ub){
               lb_last = lb;
               ub_last = ub;
               if(aux_Y > ub){
                 aux_mu = ub;
               }else if(aux_Y < lb){
                 aux_mu = lb;
               }else{
                 aux_mu = aux_Y;
               }
               if(i < aux_imin){
                 aux_imin = i;
               }
               if(i > aux_imax){
                 aux_imax = i;
               }
               aux_c = sum((Y[Range(j+1,i)]-aux_mu)*(beta-ifelse(Y[Range(j+1,i)]<aux_mu,1.0,0.0)))+C[j];
               if(aux_c < C[i]){
                 C[i] = aux_c;
                 mu[i] = aux_mu;
                 L[i] = j+2;
                 J[i] = njmp + 1;
               }
             }else{
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_i = x1(q_tilde,beta);
                 lq_max_idx_i = fmin(fmax(1,(int)ceil(llen*lq_max_i)),llen);
                 lb_max_i = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_i,lk+1,lk+llen)]];
                 //lb_max_i = get_Quantile(lq_max_idx_i,lk+1,lk+llen,R,lr,RANK);
               }

               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_i = x2(q_tilde,beta);
                 uq_max_idx_i = fmax(fmin(ulen,(int)floor(ulen*uq_max_i)),1);
                 ub_max_i = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_i,uk+1,uk+ulen)]];
                 //ub_max_i = get_Quantile(uq_max_idx_i,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_i>ub_max_i && j >= imax){
                 isstop = TRUE;
               }
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_j = x1(q_tilde,beta);
                 lq_max_idx_j = fmin(fmax(1,(int)ceil(llen*lq_max_j)),llen);
                 lb_max_j = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_j,lk+1,lk+llen)]];
                 //lb_max_j = get_Quantile(lq_max_idx_j,lk+1,lk+llen,R,lr,RANK);
               }
               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_j = x2(q_tilde,beta);
                 uq_max_idx_j = fmax(fmin(ulen,(int)floor(ulen*uq_max_j)),1);
                 ub_max_j = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_j,uk+1,uk+ulen)]];
                 //ub_max_j = get_Quantile(uq_max_idx_j,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_j>ub_max_j){
                 break;
               }
             }
           }
         }
         if(std::isinf(C[i])&&isstop){
           break;
         }
       }
     }
     imin = aux_imin;
     imax = aux_imax;
     njmp = njmp + 1;
     // Rcout << "Jump from" << njmp-1 << "to" << njmp << "\n";
    }

    NumericVector value(njmp+1);
    IntegerVector left(njmp+1);
    int cnt = n-1;
    for (int i = njmp; i >= 0; --i) {
     left[i] = L[cnt];
     value[i] = mu[cnt];
     cnt = L[cnt] - 2;
    }
    NumericVector first(2);
    NumericVector last(2);
    first[0] = lb_first;
    first[1] = ub_first;
    last[0] = lb_last;
    last[1] = ub_last;
    return List::create(Named("value") = value, Named("left") = left,
                       Named("n") = n, Named("first") = first,
                       Named("last") = last);
}


// [[Rcpp::export(.DMUSCLE_FULL)]]
List DMUSCLE_FULL(NumericVector &Y, NumericVector &q, double beta, int lag, bool test, bool details){
   //Rcout << "MUSCLE FULL is running, please wait!\n";
   int n = Y.size();
   int imax,lq_idx, uq_idx,lk,uk, llen, ulen, llen_idx, ulen_idx,lq_max_idx, uq_max_idx,idx,l;
   double lb_max, ub_max, lb, ub, q_tilde, aux_Y, aux_mu, aux_c,aux_lb, aux_ub,lq,uq, lq_max,uq_max;
   double lb_first, ub_first, lb_last, ub_last;
   NumericVector C(n);  //optimal cost
   NumericVector J(n);  //number of jumps
   NumericVector mu(n); //optimal value of step function
   NumericVector L(n);  //leftmost end points
   for (int i = 0; i < n; ++i){
     C[i] = INFINITY;
     J[i] = n;
     L[i] = 0;
   }
   double Q_max = fmax(-log(beta),-log(1-beta)); //maximum of f on [0,1]
   NumericVector Q_tilde(n);
   NumericVector Q_tilde_max(n);
   IntegerVector Q_tilde_max_idx(n);
   bool get_l_min = false;
   int l_min;
   for(int i = 0; i<n; i++){
     l = i + 1;
     for(int j = 0; j<n; j++){
       if(j==i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
         if(get_l_min == false && Q_tilde[j] <= Q_max -1e-6){
           get_l_min = true;
           l_min = l;
         }
       }
       if(j>=i){
         Q_tilde[j] = pow(penfsC(j,l,q),2)/(2*l);
       }else{
         Q_tilde[j] = R_NegInf;
       }
     }
     Q_tilde_max[i] = max(Q_tilde);
     Q_tilde_max_idx[i] = which_max(Q_tilde);
   }

   //return List::create(Named("value") = Q_tilde_max);
   //int k_step = 1;
   int k_step;
   if(n<= 500){
     k_step = 2;
   }else if(n<=1000){
     k_step = floor(log(n)/1.5);
   }else{
     k_step = floor(log(n)*2);
   }
   //Rcout << "k_step = " << k_step << "\n";
   // for(int j = 0;j<n;j++){
   //   if(Q_tilde(j,j)<=Q_max-1e-6){
   //     l_min = j+1;
   //     break;
   //   }
   // }

   WT Tree = WT(Y);

   //First search
   imax = 0;
   lb_max = R_NegInf;
   ub_max = R_PosInf;
   for (int i = 0; i< n; ++i) {
     //Rcout << i << "\n";
     lb = R_NegInf;
     ub = R_PosInf;
     for (int l = i+1; l >=l_min; l--) {
       for (int k = lag; k <= i-l+1; k = k + fmax(floor(l/2),1)) {
         //fmax(floor(l/2),1)
         idx = fmin(fmax(1,(int)ceil(l*beta)),l);
         aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k,k+l-1)]];
         //aux_Y = get_Quantile(idx,k,k+l-1,R,lr,RANK);
         q_tilde = pow(penfsC(i,l,q),2)/(2*l);
         if(q_tilde<=Q_max-1e-6){
           lq = x1(q_tilde,beta);
           lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
           if(beta==0.5){
             uq_idx = l-lq_idx;
           }else{
             uq = x2(q_tilde,beta);
             uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
           }
           aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k,k+l-1)]];
           //aux_lb = get_Quantile(lq_idx,k,k+l-1,R,lr,RANK);
           if(aux_lb>lb){
             lb = aux_lb;
             lk = k;
             llen = l;
             llen_idx = Q_tilde_max_idx[llen-1];
           }
           aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k,k+l-1)]];
           //aux_ub = get_Quantile(uq_idx,k,k+l-1,R,lr,RANK);
           if(aux_ub<ub){
             ub = aux_ub;
             uk = k;
             ulen = l;
             ulen_idx = Q_tilde_max_idx[ulen-1];
           }
         }
       }
     }
     if(lb <= ub){
       lb_first = lb;
       ub_first = ub;
       lb_last = lb;
       ub_last = ub;
       if(aux_Y > ub){
         mu[i] = ub;
       }else if(aux_Y < lb){
         mu[i] = lb;
       }else{
         mu[i] = aux_Y;
       }
       C[i] = sum((Y[Range(0,i)]-mu[i])*(beta-ifelse(Y[Range(0,i)]<mu[i],1.0,0.0)));
       L[i] = 1;
       J[i] = 0;
     }else{
       q_tilde = Q_tilde_max[llen-1];
       if(q_tilde<=Q_max-1e-6){
         lq_max = x1(q_tilde,beta);
         lq_max_idx = fmin(fmax(1,ceil(llen*lq_max)),llen);
         lb_max = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx,lk,lk+llen-1)]];
       }

       q_tilde = Q_tilde_max[ulen-1];
       if(q_tilde<=Q_max-1e-6){
         uq_max = x2(q_tilde,beta);
         uq_max_idx = fmax(fmin(ulen,floor(ulen*uq_max)),1);
         ub_max = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx,uk,uk+ulen-1)]];
       }
       if(lb_max>ub_max){
         break;
       }
     }
   }
   if(test == true && std::isinf(C[n-1])){
     if(details == true){
       Rcout << "Fast test is done! No alternative solution!\n";
     }
     IntegerVector left(2);
     left[0] = 1;
     left[1] = 1;
     return List::create(Named("left") = left);
   }
   // Later searchs
   int imin = 0;
   int njmp = 0;
   int jmin, aux_imin, aux_imax,i_llen_idx,j_llen_idx,i_ulen_idx,j_ulen_idx,lq_max_idx_i,uq_max_idx_i,lq_max_idx_j;
   bool isstop = FALSE;
   double lq_max_i,lb_max_i,uq_max_i, ub_max_i,lq_max_j,lb_max_j,uq_max_j,ub_max_j, uq_max_idx_j;

   while (std::isinf(C[n-1]) && njmp < n) {
     jmin = n;
     aux_imin = n;
     aux_imax = 0;
     lb_max_i = R_NegInf;
     ub_max_i = R_PosInf;
     for (int i=imin; i<n;i++) {
       //Rcout << i << "\n";
       if(J[i] == njmp){
         jmin = fmin(jmin,i);
       }else if(std::isinf(C[i])){
         isstop = FALSE;
         lb_max_j = R_NegInf;
         ub_max_j = R_PosInf;
         for (int j=i-1; j>=jmin;j--) {
           if(J[j] == njmp){
             lb = R_NegInf;
             ub = R_PosInf;
             for (int l=i-j; l >=l_min;l--) {
               for (int k=j+lag; k<= i-l; k = k + fmax(floor(l/2),1)) {
                 //fmax(floor(l/2),1)
                 idx = fmin(fmax(1,(int)ceil(l*beta)),l);
                 aux_Y = Y[Tree.IDX[Tree.rngQuantile(idx,k+1,k+l)]];
                 //aux_Y = get_Quantile(idx,k+1,k+l,R,lr,RANK);
                 q_tilde = pow(penfsC(i-j-1,l,q),2)/(2*l);
                 if(q_tilde<=Q_max-1e-6){
                   lq = x1(q_tilde,beta);
                   lq_idx = fmin(fmax(1,(int)ceil(l*lq)),l);
                   if(beta==0.5){
                     uq_idx = l-lq_idx;
                   }else{
                     uq = x2(q_tilde,beta);
                     uq_idx = fmax(fmin(l,(int)floor(l*uq)),1);
                   }
                   aux_lb = Y[Tree.IDX[Tree.rngQuantile(lq_idx,k+1,k+l)]];
                   //aux_lb = get_Quantile(lq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_lb>lb){
                     lb = aux_lb;
                     lk = k;
                     llen = l;
                     i_llen_idx = n - imin - 1;
                     j_llen_idx = i - jmin - 1;
                   }
                   aux_ub = Y[Tree.IDX[Tree.rngQuantile(uq_idx,k+1,k+l)]];
                   //aux_ub = get_Quantile(uq_idx,k+1,k+l,R,lr,RANK);
                   if(aux_ub<ub){
                     ub = aux_ub;
                     uk = k;
                     ulen = l;
                     i_ulen_idx = n - imin - 1;
                     j_ulen_idx = i - jmin - 1;
                   }
                 }
               }
             }
             if(lb<=ub){
               lb_last = lb;
               ub_last = ub;
               if(aux_Y > ub){
                 aux_mu = ub;
               }else if(aux_Y < lb){
                 aux_mu = lb;
               }else{
                 aux_mu = aux_Y;
               }
               if(i < aux_imin){
                 aux_imin = i;
               }
               if(i > aux_imax){
                 aux_imax = i;
               }
               aux_c = sum((Y[Range(j+1,i)]-aux_mu)*(beta-ifelse(Y[Range(j+1,i)]<aux_mu,1.0,0.0)))+C[j];
               if(aux_c < C[i]){
                 C[i] = aux_c;
                 mu[i] = aux_mu;
                 L[i] = j+2;
                 J[i] = njmp + 1;
               }
             }else{
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_i = x1(q_tilde,beta);
                 lq_max_idx_i = fmin(fmax(1,(int)ceil(llen*lq_max_i)),llen);
                 lb_max_i = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_i,lk+1,lk+llen)]];
                 //lb_max_i = get_Quantile(lq_max_idx_i,lk+1,lk+llen,R,lr,RANK);
               }

               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_i = x2(q_tilde,beta);
                 uq_max_idx_i = fmax(fmin(ulen,(int)floor(ulen*uq_max_i)),1);
                 ub_max_i = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_i,uk+1,uk+ulen)]];
                 //ub_max_i = get_Quantile(uq_max_idx_i,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_i>ub_max_i && j >= imax){
                 isstop = TRUE;
               }
               q_tilde = Q_tilde_max[llen-1];
               if(q_tilde<=Q_max-1e-6){
                 lq_max_j = x1(q_tilde,beta);
                 lq_max_idx_j = fmin(fmax(1,(int)ceil(llen*lq_max_j)),llen);
                 lb_max_j = Y[Tree.IDX[Tree.rngQuantile(lq_max_idx_j,lk+1,lk+llen)]];
                 //lb_max_j = get_Quantile(lq_max_idx_j,lk+1,lk+llen,R,lr,RANK);
               }
               q_tilde = Q_tilde_max[ulen-1];
               if(q_tilde<=Q_max-1e-6){
                 uq_max_j = x2(q_tilde,beta);
                 uq_max_idx_j = fmax(fmin(ulen,(int)floor(ulen*uq_max_j)),1);
                 ub_max_j = Y[Tree.IDX[Tree.rngQuantile(uq_max_idx_j,uk+1,uk+ulen)]];
                 //ub_max_j = get_Quantile(uq_max_idx_j,uk+1,uk+ulen,R,lr,RANK);
               }
               if(lb_max_j>ub_max_j){
                 break;
               }
             }
           }
         }
         if(std::isinf(C[i])&&isstop){
           break;
         }
       }
     }
     imin = aux_imin;
     imax = aux_imax;
     njmp = njmp + 1;
     //Rcout << "Jump from" << njmp-1 << "to" << njmp << "\n";
   }

   NumericVector value(njmp+1);
   IntegerVector left(njmp+1);
   int cnt = n-1;
   for (int i = njmp; i >= 0; --i) {
     left[i] = L[cnt];
     value[i] = mu[cnt];
     cnt = L[cnt] - 2;
   }
   NumericVector first(2);
   NumericVector last(2);
   first[0] = lb_first;
   first[1] = ub_first;
   last[0] = lb_last;
   last[1] = ub_last;
   return List::create(Named("value") = value, Named("left") = left,
                       Named("n") = n, Named("first") = first,
                       Named("last") = last);
 }

double ff(double x,double b,double q){
  double res;
  if(x<1 && x>0){
    res = (x*log(x/b)+(1-x)*log((1-x)/(1-b)))-q;
  }else if(x == 0){
    res = -log(1-b)-q;
  }else if(x == 1){
    res = -log(b)-q;
  }else{
    res = -q;
  }
  return res;
}

NumericVector cumsum(NumericVector &x){
  // initialize the result vector
  NumericVector res(x.size());
  std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}

double quantileC(NumericVector &x, double q) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}

// [[Rcpp::export(.simulQuantile_MUSCLE)]]
NumericVector simulQuantile_MUSCLE(double p, int n, double beta) {
  Rcout << "Simulating quantiles and this can take long time...";
  int R = (int)(50/std::min(p,1-p));
  NumericMatrix data(R,n);
  NumericVector res(n);
  for (int r = 0; r<R;r++) {
    NumericVector X = rbinom(n,1,beta);
    NumericVector S = cumsum(X);
    S.push_front(0);
    NumericVector um(n);
    for(int i = 0; i < n; i++){
      um[i] = R_NegInf;
    }
    for (int I = 0; I < n; I++) {
      NumericVector pen(I+1);
      for(int i = 0; i< I+1; i++){
        pen[i] = sqrt(2*log(exp(1)*(I+1)/(1.0*(i+1))));
      }
      for (int l = 1; l<= (I+1); l++) {
        double X_bar = (S[I+1]-S[I-l+1])/l;
        um[l-1] = std::max(um[l-1],sqrt(2*l*ff(X_bar,beta,0)));
      }
      NumericVector tmp = (um[Range(0,I)]-pen);
      data(r,I) = max(tmp);
    }
    if(r % 20 == 0){
      Rcout << 1.0*r/R*100 << "% are simulated!\n";
    }
  }
  for (int I = 0; I < n;I++) {
    NumericVector temp = data( _ ,I);
    res[I] = quantileC(temp,p);
  }
  return(res);
}

// [[Rcpp::export(.logg)]]
NumericVector logg(NumericVector &x){
  int n = x.size();
  NumericVector res(n);
  for(int i=0;i<n;i++){
    if(x[i] == 0){
      res[i] = 0;
    }else{
      res[i] = log(x[i]);
    }
  }
  return(res);
}

NumericVector Multi(NumericMatrix &L, NumericVector &X){
  int n = X.size();
  NumericVector res(n);
  for(int i = 0;i<n;i++){
    for(int j = 0;j<n;j++){
      res[i] = res[i] + L(i,j)*X[j];
    }
  }
  return(res);
}

// [[Rcpp::export(.simulQuantile_DMUSCLE)]]
NumericVector simulQuantile_DMUSCLE(NumericVector &X, NumericVector &ACF, int n) {
  NumericVector S = cumsum(X);
  S.push_front(0);
  NumericVector um(n);
  NumericVector pen(n);
  NumericVector res(n);
  NumericVector sd(n);
  double sumACF = sum(ACF);
  int m = ACF.size();
  NumericVector ACF_adj(m);
  ACF_adj[0] = ACF[0];
  for(int i =1; i<m;i++){
    ACF_adj[i] = 2*ACF[i];
  }
  NumericVector SACF = cumsum(ACF_adj);
  SACF = cumsum(SACF);
  for(int i = 0;i<m;i++){
    sd[i] = SACF[i];
  }
  for(int i = m;i<n;i++){
    sd[i] = sd[m-1] + 2*(i-(m-1))*(sumACF - ACF[0]);
  }
  for(int i = 0;i<n;i++){
    sd[i] = sqrt(sd[i]);
  }
  //return sd;
  for(int i = 0; i < n; i++){
    um[i] = R_NegInf;
  }
  for (int I = 0; I < n; I++) {
    for(int i = 0; i< I+1; i++){
      pen[i] = sqrt(2*log(exp(1)*(I+1)/(1.0*(i+1))));
    }
    for (int l = 1; l<= (I+1); l++) {
      double X_diff = (S[I+1]-S[I-l+1]);
      um[l-1] = std::max(um[l-1],abs(X_diff)/sd[l-1]);
    }
    NumericVector tmp = um[Range(0,I)]-pen;
    res[I] = max(tmp);
  }
  return(res);
}

