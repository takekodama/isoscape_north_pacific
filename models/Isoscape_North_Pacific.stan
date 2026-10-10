// =============================================================================
// Stan Model: Linear Environmental Drivers & Site RE for POM d15N
// =============================================================================

data {
  int<lower=1> N;                           // サンプル数
  int<lower=1> N_site;                      // サイト（海域）の数
  array[N] int<lower=1, upper=N_site> site; // 各サンプルのサイトID

  // 環境要因（全て標準化済みを推奨）
  vector[N] nitrate;                        // 硝酸塩濃度
  vector[N] temp;                           // 水温
  vector[N] chla;                           // log(Chl-a)
  vector[N] sal;                            // 塩分

  vector[N] d15N_pom;                       // POMのd15N観測値 (per mil)
}

parameters {
  real alpha0;                              // 切片 (ベースライン d15N)

  // 各環境要因の回帰係数
  real gamma_nitrate;
  real gamma_temp;
  real gamma_chla;
  real gamma_sal;

  // サイト間のランダム効果（非中心化）
  vector[N_site] site_raw;
  real<lower=0> sigma_site;

  real<lower=0> sigma;                      // 観測残差の標準偏差
}

transformed parameters {
  vector[N_site] site_effect = site_raw * sigma_site;
  vector[N] d15N_expected;

  for (i in 1:N) {
    // 直線回帰モデル + サイトランダム効果
    d15N_expected[i] = alpha0
                       + gamma_nitrate * nitrate[i]
                       + gamma_temp * temp[i]
                       + gamma_chla * chla[i]
                       + gamma_sal * sal[i]
                       + site_effect[site[i]];
  }
}

model {
  // --- 事前分布 (Priors) ---
  alpha0        ~ normal(5.0, 3.0);
  gamma_nitrate ~ normal(0, 1);
  gamma_temp    ~ normal(0, 1);
  gamma_chla    ~ normal(0, 1);
  gamma_sal     ~ normal(0, 1);

  site_raw      ~ normal(0, 1);
  sigma_site    ~ exponential(1);

  sigma         ~ cauchy(0, 2.5);

  // --- 尤度 (Likelihood) ---
  d15N_pom ~ normal(d15N_expected, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(d15N_pom[n] | d15N_expected[n], sigma);
  }
}

