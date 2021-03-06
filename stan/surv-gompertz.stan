functions {
  real gompertz_lpdf(real t, real l, real g) {
    return log(l) + g * t - (l/g) * (exp(g*t)-1);
  }
  
  real gompertz_surv(real t, real l, real g) {
    return exp(-(l/g) * (exp(g*t)-1));
  }
  
  real gompertz_haz(real t, real l, real g) {
    return l * exp(g*t);
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
  real<lower=0> g;
}

model {
  l ~ cauchy(0, 30);
  g ~ cauchy(0, 30);

  for(i in 1:N) {
    t[i] ~ gompertz(l, g);
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N_pred] surv_pred;
  vector[N_pred] haz_pred;
  
  for (i in 1:N){
    log_lik[i] = gompertz_lpdf(t[i] | l, g);
  }
  
  for(j in 1:N_pred) {
    surv_pred[j] = gompertz_surv(t_new[j], l, g);
    haz_pred[j] = gompertz_haz(t_new[j], l, g);
  }
}
