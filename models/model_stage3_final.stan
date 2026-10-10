data {
  int<lower=0> N;
  int<lower=0> N_site;

  vector[N] d15N;
  vector[N] POM_obs; // Model 1 から渡された予測値

  vector[N] BL;

  // サイトごとの環境共変量
  vector[N_site] temp_site;
  vector[N_site] nitrate_site;
  vector[N_site] chla_site;
  vector[N_site] sal_site;

  array[N] int<lower=1> site;
  array[N] int<lower=1> species;
  int<lower=1> N_species;
}

parameters {
  // 各サイトの真のPOMベースライン（POM_obs を中心とする潜在変数）
  vector[N_site] POM_true;
  real<lower=0> sigma_pom;

  // ---- Fish Model ----
  real beta0;
  real beta_pom;      // POM d15N に対する応答係数（通常 1 前後）
  real beta_BL;       // 体長効果

  // 魚の栄養段階/同位体濃縮に与える環境の直接効果
  real beta_temp;
  real beta_nitrate;
  real beta_chla;
  real beta_sal;

  vector[N_species] species_raw;
  real<lower=0> sigma_species;

  real<lower=0> sigma; // 魚モデルの観測残差
}

transformed parameters {
  vector[N_species] species_effect = species_raw * sigma_species;
}

model {
  // ====================
  // Priors
  // ====================
  sigma_pom ~ normal(0.5, 0.2); // 観測誤差の広がりを適度に制約

  beta0 ~ normal(0, 5);
  beta_pom ~ normal(1, 0.2);   // ベースラインの反映係数（事前分布で1付近に誘導）

  beta_BL ~ normal(0, 1);
  beta_temp ~ normal(0, 1);
  beta_nitrate ~ normal(0, 1);
  beta_chla ~ normal(0, 1);
  beta_sal ~ normal(0, 1);

  species_raw ~ normal(0, 1);
  sigma_species ~ exponential(1);
  sigma ~ exponential(1);

  // ====================
  // 1. POM Observation / Latent Model
  // ====================
  // Model 1の予測値 (POM_obs) を平均とする潜在変数 POM_true の定義
  for (n in 1:N) {
    POM_obs[n] ~ normal(POM_true[site[n]], sigma_pom);
  }

  // ====================
  // 2. Fish Model
  // ====================
  for (n in 1:N) {
    d15N[n] ~ normal(
      beta0
      + beta_pom * POM_true[site[n]]
      + beta_BL * BL[n]
      + beta_temp * temp_site[site[n]]
      + beta_nitrate * nitrate_site[site[n]]
      + beta_chla * chla_site[site[n]]
      + beta_sal * sal_site[site[n]]
      + species_effect[species[n]],
      sigma
    );
  }
}
