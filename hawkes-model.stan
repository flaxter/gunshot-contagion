functions {
  real gaussian(real x, real lengthscale) {
    return(1/(2*pi())*exp(-.5 * x^2/lengthscale^2)/lengthscale^2);
  }
}
data {
  int<lower=1> n;
  int hour[n];
  matrix[n,2] Space;
  vector[n] Time; // NB: these must be sorted from smallest to largest!
  real time_window;
  real space_window;

  vector[n] muST;
  real muSTintegral;
}
transformed data {
  matrix[n,n] timeD;
  matrix[n,n] spaceD;
  vector[n] timeD2;
  for(i in 1:n) {
    for(j in 1:n) {
      timeD[i,j] <- -(Time[i] - Time[j]);
      spaceD[i,j] <- distance(Space[i], Space[j]);
    }
    timeD2[i] <- -(time_window - Time[i]); 
  }
}
parameters {
  real<lower=0> lengthscaleS;
  real<lower=0> lengthscaleT;
  real<lower=0> a;
  real<lower=0> mu;
  real<lower=0> mu0;
}
transformed parameters {
  vector[n] ll;
  real lp;
  ll <- mu0 + muST * mu;
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      ll[i] <- ll[i] + (a * lengthscaleT * exp(timeD[i,j] * lengthscaleT) * gaussian(spaceD[i,j],lengthscaleS));
    }
  }
  lp <- sum(log(ll)) - muSTintegral * mu - mu0 * space_window * time_window +
          a * (sum(exp(timeD2*lengthscaleT)) - n);
}

model {
  increment_log_prob(lp);
    
  lengthscaleS ~ normal(0,10);
  lengthscaleT ~ normal(0,10);
  a ~  normal(0,10);
  mu0 ~ normal(0,1);
  mu ~ normal(0,1);
}

generated quantities {
  vector[n] background;
  real lengthscale_minutes;
  lengthscale_minutes <- 24*60/lengthscaleT;
  background <- (mu0 + muST *mu ) ./ ll;
}
