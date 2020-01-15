data {
  int<lower=0> N;
  int<lower=0> n;
  int<lower=0> y[N,n];
  real<lower=0> s_gmpr[n];
  real<lower=0> g[N];
  int<lower = 0, upper = 1> ip[n];
  real<lower=0> beta_mu;
  real<lower=0> alpha_mu;
  real<lower=0> beta_delta;
  real<lower=0> alpha_delta;
  real<lower=0> beta_phi;
  real<lower=0> alpha_phi;
}

transformed data {
  // this block basically splits the input and IPed 
  // samples into two different matrices for 
  // easier computation. 

  real<lower=0> s[n] = s_gmpr;
  int m = n/2; //number of samples per condition input/ip
  int k = 1; int l = 1;
  int<lower=0> y_input[N,m];
  int<lower=0> y_ip[N,m];
  for (j in 1:n){
    if (ip[j] == 0){
      for (i in 1:N){
        y_input[i,k] = y[i,j];
      }
      k+=1;
    } else {
      for (i in 1:N){
        y_ip[i,l] = y[i,j];
      }
      l+=1;
    }
  }
}


parameters {
  real <lower=0, upper=1> pi;
  real <lower = 0> mu[N];
  real <lower = 0> phi[N];
  real <lower = 0> delta[N];
}



model {
  vector[2] contributions;

  for (j in 1:m){
    for (i in 1:N) {
      //distribution of y_input
      target += neg_binomial_2_lpmf(y_input[i,j] | s[j] * g[i] * mu[i], phi[i]);

      //distribution of y_ip
      contributions[1] = log(pi) + neg_binomial_2_lpmf(y_ip[i,j] | s[j + m] * mu[i] * g[i] * delta[i], phi[i]); //mu[i] *
      contributions[2] = log(1 - pi) + normal_lpdf(y_ip[i,j] | 0, 0.00000000001); 
      target += log_sum_exp(contributions);
    }
  }
  

  phi ~ normal(0, 4);

  // Gamma: mean = alpha/beta, var = alpha / beta^2
  mu ~ gamma(alpha_mu, beta_mu);
  // just some reasonable values: large parameters are allowed, but the mean is around 0.4. So there is regularisation downards. 
  delta ~ gamma(0.2, 0.5);
  
}



/*comments: 

=======================================================================
one could think about having a prior on the reversal of phi. However, it seems to be working as it is. 

this would be the code: 
real <lower = 0> reciprocal_sqrt_phi[N];

transformed parameters {
  real <lower = 0> phi[N];
  for (i in 1:N){
    phi[i] = 1 / (reciprocal_sqrt_phi[i]^2);
  }
}
=======================================================================


=======================================================================
the priors for mu and delta are somewhat arbitrary. 
I found them to be working well, but maybe they should be adjusted. 
They definitely have an impact on the estimation. 
=======================================================================

*/

