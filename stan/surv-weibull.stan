functions {
  real weib_lpdf(real t, real l, real b) {
    return b * log(l) + log(b) + (b-1)*log(t) - (l*t)^b;
  }
  
  real weib_surv(real t, real l, real b) {
    return exp(-(l*t)^b);
  }
  
  real weib_haz(real t, real l, real b) {
    return l^b * b * t^(b-1);
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_pred;
  real t[N];
  vector[N_pred] t_new;
}

parameters {
  real<lower=0> l;
  real<lower=0> b;
}

model {
  l ~ cauchy(0, 30);
  b ~ cauchy(0, 30);
  
  for(i in 1:N) {
    t[i] ~ weib(l, b);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  vector[N_pred] haz_pred;
  
  for (i in 1:N) {
    log_lik[i] = weib_lpdf(t[i] | l, b);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] = weib_surv(t_new[j], l, b);
    haz_pred[j] = weib_haz(t_new[j], l, b);
  }
}
