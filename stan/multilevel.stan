
data {
  int<lower=0> n;
  int<lower=0> k;
  int<lower=0> nsites;
  int<lower=0> site[n];
  vector[n] y;
  matrix[nsites,k] X;       // model matrix
}

transformed data {
  real mean_x1;
  real mean_x2;
  real mean_x3;
  
  real sd_x1;
  real sd_x2;
  real sd_x3;
  
  mean_x1 = mean(X[,2]);
  mean_x2 = mean(X[,3]);
  mean_x3 = mean(X[,4]);
  
  sd_x1 = sd(X[,2]);
  sd_x2 = sd(X[,3]);
  sd_x3 = sd(X[,4]);
}

parameters {
  vector[nsites] alpha;
  vector[k] beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_y;
}

transformed parameters {
  vector[nsites] alpha_hat;
  vector[n] y_hat;
  
  alpha_hat = X * beta;
  
  for (i in 1:n)
    y_hat[i] = alpha[site[i]];
}

model {
  sigma_alpha ~ cauchy(0, 10);
  sigma_y ~ cauchy(0, 10);
  beta ~ normal(0, 10);

  alpha ~ normal(alpha_hat, sigma_alpha);
  y ~ normal(y_hat, sigma_y);
}

generated quantities {
  vector[k-1] beta_std;

  beta_std[1] = beta[2] * sd_x1 / sd(alpha);
  beta_std[2] = beta[3] * sd_x2 / sd(alpha);
  beta_std[3] = beta[4] * sd_x3 / sd(alpha);
}
