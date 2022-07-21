data {
    int<lower=0> N;                         // quantidade de paises
    int<lower=0> D;                         // número de preditores
    int<lower=0> L;                         // quantidade de levels da hierarquia
    int<lower=0> y[N];                      // número de mortes em cada pais
    int<lower=0> pop[N];                    // população em cada pais
    int<lower=0> ll[N];                     // identificador da classe da hierarquia
    row_vector[D] x[N];                     // vetor de preditores
}
parameters {
    real mu[D];
    real<lower=0> sigma[D];
    vector[D] beta[L];
}
model {
    for (d in 1:D) {
        mu[d] ~ normal(0, 5);
        sigma[d] ~ exponential(5);
        for (l in 1:L) {
            beta[l,d] ~ normal(mu[d], sigma[d]);
        }
    }
    for (n in 1:N) {
        y[n] ~ binomial_logit(pop[n], x[n] * beta[ll[n]]);
    }
}
generated quantities {
    vector[N] y_rep;
    for (n in 1:N) {
        y_rep[n] = binomial_rng(pop[n],inv_logit(x[n] * beta[ll[n]]));
    }
}
