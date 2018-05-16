/*
*Multi-session abundace and log-linear regression model* 
*for Rainbow Trout in the Lower Bow River, Alberta *
*Chris L. Cahill 16 May 2018*
*See also Hiroki Ito's translations of BUGS examples in*
*Bayesian Population Analysis (Kéry and Schaub 2012) to Stan code,* 
*available at: https://github.com/stan-dev/example-models/tree/master/BPA
*/
data {            
  int<lower=0> M; //number of rows in augmented datasets
  int<lower=0> J; //number of capture occassions 
  int<lower=0> nyear; //number of sessions
  int<lower=0, upper=1> y_array[M, J, nyear]; //3d capture history array
  real X[nyear]; //year  
  real km[nyear];  //length of stretch fished (unstandardized)  
}
transformed data {
  int<lower=0> s[nyear,M];  // totals in each row
  int<lower=0> Nobs[nyear]; //number of individuals observed each year

   for(i in 1:nyear){
     Nobs[i] = 0;
     for(j in 1:M) {
       s[i,j] = sum(y_array[j, ,i]);
       if (s[i,j] > 0){
        Nobs[i] = Nobs[i] + 1; 
       }
     }
   }
}
parameters {
  matrix<lower=0, upper=1>[nyear,J] mean_p; //detection probs
  real<lower=0> beta0; 
  real beta1; 
  vector<lower=0>[nyear] sigma; //weakly informative prior 
  vector [M] sigma_raw[nyear]; //non-centering terms
}
transformed parameters {
  matrix[nyear, J] mean_logit_p;
  real<lower=0, upper=1> p_array[nyear, M, J];
  vector<lower=0, upper=1> [nyear] psi; 
  real<lower=0> lambda[nyear];  
  vector<lower=-16, upper=16> [M] eps[nyear];

  for(i in 1:nyear){
    lambda[i] = exp( beta0 + beta1*( X[i]-X[4] )); //log-linear regression 
    psi[i] = lambda[i]/M;  //psi[i] is a derived variable
    eps[i,] = sigma[i]*sigma_raw[i,]; 
    for(j in 1:M){
      for(k in 1:J){
        mean_logit_p[i,k] = logit(mean_p[i,k]); //define the mean logit of cap prob per year
        p_array[i,j,k] = inv_logit(mean_logit_p[i,k] + eps[i,j]); //add the individual het terms
      }
    }
  } 
}
model {  
  //priors
  //capture probs[nyear, J] implicitly ~ uniform(0,1) 
  beta0 ~ normal(0,10); 
  beta1 ~ normal(0,1); 
  sigma ~ normal(0,1.25); //informative prior for individual heterogeneity

  for(i in 1:nyear){
     sigma_raw[i] ~ normal(0,1); 
    for (j in 1:M) { 
      if (s[i,j] > 0) {
        //z[i] == 1 == critter was observed in the field:
        target += bernoulli_lpmf(1 | psi[i]) + bernoulli_lpmf(y_array[j, ,i] | p_array[i,j,]);
      } else { //nind_obs[i]==0
        vector[2] lp;
        //z[i] == 1 == critter exists (1|psi), but was not observed:
        lp[1] = bernoulli_lpmf(1 | psi[i]) + bernoulli_lpmf(0 | p_array[i,j,]);
      
        //z[i] == 0 == critter doesn't exist (0|psi): 
        lp[2] = bernoulli_lpmf(0 | psi[i]); 
    
        target += log_sum_exp(lp);
      } 
    }
  }
}
generated quantities {
  //derived quantities:
  int G_N[nyear];  //Group or Session N
  vector[nyear] G_N_km; // N / Distance electro-fished
  int<lower=0, upper=1> z_matrix[nyear,M]; 

  for(i in 1:nyear){
    for(j in 1:M){
      if(s[i,j] > 0) {
       //z[i,j] == 1 == critter was observed in the field:
       z_matrix[i,j] = 1;  
      } else { //never detected
       real pr; //prob never detected given it was present
       pr = prod(rep_row_vector(1.0, J) - to_row_vector(p_array[i,j,])); 
       z_matrix[i,j] = bernoulli_rng(psi[i]*pr/(psi[i]*pr + (1-psi[i])));
      }
    }
  }
  for(i in 1:nyear){
    G_N[i] = sum(z_matrix[i,]);
    G_N_km[i] = G_N[i]/km[i]; 
  } 
}
