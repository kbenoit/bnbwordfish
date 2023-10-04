stanPoisson1d <- '
data {
int<lower=1> D; // number of docs
int<lower=1> T; // number of terms/cats
int<lower=1> N; // number of obs = number of docs * number of cats
int<lower=1> doc[N]; // doc id
int<lower=1> term[N]; // cat id
int<lower=0> y[N]; // counts as one long vector
}

parameters {
vector[D] alpha;

real mu_psi;
real<lower=0> sigma_psi;
vector[T-1] psi_rest; 

real<lower=0> sigma_beta;
vector[T] beta; 

vector[D] theta;
}

transformed parameters {
vector[T] psi;
vector[N] logmu;

for(t in 1:(T-1)){
psi[t] = psi_rest[t];
}
psi[T] = 0;

for(n in 1:N){
logmu[n] = alpha[doc[n]] + psi[term[n]] + theta[doc[n]] * beta[term[n]];
} 
}


model {

alpha ~ normal(0,10);

mu_psi ~ normal(0,5);
sigma_psi ~ cauchy(0,5);
psi_rest ~ normal(mu_psi,sigma_psi);

sigma_beta ~ cauchy(0,5);
beta ~ normal(0,sigma_beta);

theta ~ normal(0,1);

y ~ poisson_log(logmu);
}

'

stanNegbin1dWord <- '
data {
int<lower=1> D; // number of docs
int<lower=1> T; // number of terms/cats
int<lower=1> N; // number of obs = number of docs * number of cats
int<lower=1> doc[N]; // doc id
int<lower=1> term[N]; // cat id
int<lower=0> y[N]; // counts as one long vector
}


parameters {
vector[D] alpha;

real mu_psi;
real<lower=0> sigma_psi;
vector[T-1] psi_rest; 

real<lower=0> sigma_beta;
vector[T] beta; 

vector[D] theta;

vector<lower=0,upper=200>[T] phi_inv;
}

transformed parameters {
vector[T] psi;
vector[N] logmu;
vector[N] philong;

for(t in 1:(T-1)) {
    psi[t] = psi_rest[t];
}
psi[T] = 0;

for(n in 1:N){
logmu[n] = alpha[doc[n]] + psi[term[n]] + theta[doc[n]] * beta[term[n]];
philong[n] = 1/phi_inv[term[n]];
} 
}


model {

alpha ~ normal(0,10);

mu_psi ~ normal(0,5);
sigma_psi ~ cauchy(0,5);
psi_rest ~ normal(mu_psi,sigma_psi);

sigma_beta ~ cauchy(0,5);
beta ~ normal(0,sigma_beta);

theta ~ normal(0,1);


y ~ neg_binomial_2_log(logmu,philong);
}

'

# 2d Poisson - model, including some items loading on both dimensions (BothYes)
stanPoisson2dBothYes <- '
data {
int<lower=1> D; // number of docs
int<lower=1> T; // number of terms/cats
int<lower=1> N; // number of obs = number of docs * number of cats
int<lower=1> doc[N]; // doc id
int<lower=1> term[N]; // cat id
int<lower=0> y[N]; // counts as one long vector
int<lower=1> nbz; // number of terms for which both betas zero
int<lower=1> nd1; // number of terms for which only beta1 non-zero
int<lower=1> nd2; // number of terms for which only beta2 non-zero
}

parameters {
real mu_psi;
real<lower=0> sigma_psi;
real<lower=0> sigma_beta1;
real<lower=0> sigma_beta2;
real<lower=0> sigma_theta1;
real<lower=0> sigma_theta2;
real mu_beta1;
real mu_beta2;
vector[D] alpha;
vector[T-1] psi_rest; 
vector[T-2-nbz-nd2] beta1_rest; 
vector[T-2-nbz-nd1] beta2_rest; 
vector[D] theta1;
vector[D] theta2;
}

transformed parameters {
vector[T] psi;
vector[T] beta1;
vector[T] beta2;
vector[N] logmu;

psi[1] = 0;
for(t in 2:T){
psi[t] = psi_rest[t-1];
}


beta1[1] = 1;
beta2[1] = 0;
beta1[2] = 0;
beta2[2] = 1;

for(t in 3:(2+nbz)){
beta1[t] = 0;
beta2[t] = 0;
}
for(t in (3+nbz):(2+nbz+nd1)){
beta1[t] = beta1_rest[t-2-nbz];
beta2[t] = 0;
}
for(t in (3+nbz+nd1):(2+nbz+nd1+nd2)){
beta1[t] = 0;
beta2[t] = beta2_rest[t-2-nbz-nd1];
}
for(t in (3+nbz+nd1+nd2):T){
beta1[t] = beta1_rest[t-2-nbz-nd2];
beta2[t] = beta2_rest[t-2-nbz-nd1];
}

for(n in 1:N){
logmu[n] = alpha[doc[n]] + psi[term[n]] + theta1[doc[n]] * beta1[term[n]] + theta2[doc[n]] * beta2[term[n]];
} 
}


model {
// local variables

mu_psi ~ normal(0,5);
sigma_psi ~ cauchy(0,5);
mu_beta1 ~ normal(0,5);
mu_beta2 ~ normal(0,5);
sigma_beta1 ~ cauchy(0,5);
sigma_beta2 ~ cauchy(0,5);
sigma_theta1 ~ cauchy(0,5);
sigma_theta2 ~ cauchy(0,5);

alpha ~ normal(0,10);
psi_rest ~ normal(mu_psi,sigma_psi);
beta1_rest ~ normal(mu_beta1,sigma_beta1);
beta2_rest ~ normal(mu_beta2,sigma_beta2);
theta1 ~ normal(0,sigma_theta1);
theta2 ~ normal(0,sigma_theta2);

y ~ poisson_log(logmu);
}

'


# 2d Poisson model, no items loading on both dimensions (BothNo)
stanPoisson2dBothNo <- '
data {
int<lower=1> D; // number of docs
int<lower=1> T; // number of terms/cats
int<lower=1> N; // number of obs = number of docs * number of cats
int<lower=1> doc[N]; // doc id
int<lower=1> term[N]; // cat id
int<lower=0> y[N]; // counts as one long vector
int<lower=1> nbz; // number of terms for which both betas zero
int<lower=1> nd1; // number of terms for which only beta1 non-zero
int<lower=1> nd2; // number of terms for which only beta2 non-zero
}

parameters {
real mu_psi;
real<lower=0> sigma_psi;
real<lower=0> sigma_beta1;
real<lower=0> sigma_beta2;
real<lower=0> sigma_theta1;
real<lower=0> sigma_theta2;
real mu_beta1;
real mu_beta2;
vector[D] alpha;
vector[T-1] psi_rest; 
vector[T-2-nbz-nd2] beta1_rest; 
vector[T-2-nbz-nd1] beta2_rest; 
vector[D] theta1;
vector[D] theta2;
}

transformed parameters {
vector[T] psi;
vector[T] beta1;
vector[T] beta2;
vector[N] logmu;

psi[1] = 0;
for(t in 2:T){
psi[t] = psi_rest[t-1];
}


beta1[1] = 1;
beta2[1] = 0;
beta1[2] = 0;
beta2[2] = 1;

for(t in 3:(2+nbz)){
beta1[t] = 0;
beta2[t] = 0;
}
for(t in (3+nbz):(2+nbz+nd1)){
beta1[t] = beta1_rest[t-2-nbz];
beta2[t] = 0;
}
for(t in (3+nbz+nd1):(2+nbz+nd1+nd2)){
beta1[t] = 0;
beta2[t] = beta2_rest[t-2-nbz-nd1];
}
for(n in 1:N){
logmu[n] = alpha[doc[n]] + psi[term[n]] + theta1[doc[n]] * beta1[term[n]] + theta2[doc[n]] * beta2[term[n]];
}
}


model {
// local variables


mu_psi ~ normal(0,5);
sigma_psi ~ cauchy(0,5);
mu_beta1 ~ normal(0,5);
mu_beta2 ~ normal(0,5);
sigma_beta1 ~ cauchy(0,5);
sigma_beta2 ~ cauchy(0,5);
sigma_theta1 ~ cauchy(0,5);
sigma_theta2 ~ cauchy(0,5);

alpha ~ normal(0,10);
psi_rest ~ normal(mu_psi,sigma_psi);
beta1_rest ~ normal(mu_beta1,sigma_beta1);
beta2_rest ~ normal(mu_beta2,sigma_beta2);
theta1 ~ normal(0,sigma_theta1);
theta2 ~ normal(0,sigma_theta2);

y ~ poisson_log(logmu);
}
}

'
