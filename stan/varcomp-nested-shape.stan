data {
  int N;        // number of obs
  int N_J;      // number of higher groups
  int N_K;      // number of low groups
  vector[N] y;  // response (trait value)
  int J[N];     // labels for higher groups
  int K[N];     // labels for lower groups
  real sigma_within;
} 

parameters {
  real<lower=0> sigma_J;
  real<lower=0> sigma_K;
  
  real theta;
  vector[N_J] eta_J;
  vector[N_K] eta_K;
}

transformed parameters {
  vector[N] yhat;
  vector[N_J] mu_J;
  vector[N_K] dev_K;

  mu_J = theta + eta_J * sigma_J;
  dev_K = eta_K * sigma_K;
  
  for(i in 1:N) {
    yhat[i] = mu_J[J[i]] + dev_K[K[i]];
  }
}

model {
  eta_J ~ normal(0, 1);
  eta_K ~ normal(0, 1);
  
  theta ~ normal(0, 100);
  sigma_J ~ cauchy(0, 50);
  sigma_K ~ cauchy(0, 50);
  
  y ~ normal(yhat, sigma_within);
}

generated quantities {
  real var_within;
  real var_J;
  real var_K;
  real var_total;
  
  var_within = pow(sigma_within, 2);
  var_J = pow(sigma_J, 2);
  var_K = pow(sigma_K, 2);
  var_total = var_within + var_J + var_K;
}
