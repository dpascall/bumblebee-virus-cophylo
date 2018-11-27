//a comparative model looking at the associations between hosts and viruses
//per Hadfield et al - A Tale of Two Phylogenies marginalising
//phylogenetic uncertainty see that paper for descriptions of all terms
//and formulas for calculation of covariance matrices
//code based on and inspired by that of Diogo Melo - originally licensed
//under the MIT Licence 
//github.com/diogro/stanAnimal
//marginalisation of phylogenetic uncertainty implemented per solution
//on Stan forums
//discourse.mc-stan.org/t/sampling-over-an-uncertain-variance-covariance-matrix
//pooled binomial likelihood implemented per Thompson 1971
data {
  int<lower=1> datapoints;
  int<lower=1> virusspecies;
  int<lower=1> hostspecies;
  int<lower=1> virushostcombinations;
  int<lower=1> pools;
  int<lower=1> spatialcompositions;
  int<lower=1> matrices;
  int<lower=0, upper=1> run_estimation;

  int<lower=0, upper=1> y[datapoints];
  int<lower=1, upper=hostspecies> host[datapoints];
  int<lower=1, upper=virusspecies> virus[datapoints];
  int<lower=1, upper=virushostcombinations> combination[datapoints];
  int<lower=1, upper=pools> pool[datapoints];
  int<lower=1, upper=spatialcompositions> spatialcomposition[datapoints];
  int<lower=1, upper=datapoints> ID[datapoints];
  int<lower=1> samples[datapoints];
  
  //for posterior prediction of expected prevalence in each host virus combination
  int<lower=1, upper=hostspecies> host_pred[virushostcombinations];
  int<lower=1, upper=virusspecies> virus_pred[virushostcombinations];
  int<lower=1, upper=virushostcombinations> combination_pred[virushostcombinations];
  
  matrix[hostspecies, hostspecies] HostPhy[matrices];
  matrix[virusspecies, virusspecies] VirusPhy[matrices];
  matrix[virushostcombinations, virushostcombinations] HostInter[matrices];
  matrix[virushostcombinations, virushostcombinations] VirusInter[matrices];
  matrix[virushostcombinations, virushostcombinations] CoevoInter[matrices];
}
transformed data{
  real log_unif;
  matrix[hostspecies, hostspecies] CD_HostPhy[matrices];
  matrix[virusspecies, virusspecies] CD_VirusPhy[matrices];
  matrix[virushostcombinations, virushostcombinations] CD_HostInter[matrices];
  matrix[virushostcombinations, virushostcombinations] CD_VirusInter[matrices];
  matrix[virushostcombinations, virushostcombinations] CD_CoevoInter[matrices];
  
  //for constant correction to posterior density
  log_unif = -log(matrices);
  
  //for speed cholesky decompose all covariance matrices
  for (i in 1:matrices) {
    CD_HostPhy[i] = cholesky_decompose(HostPhy[i]);
    CD_VirusPhy[i] = cholesky_decompose(VirusPhy[i]);
    CD_HostInter[i] = cholesky_decompose(HostInter[i]);
    CD_VirusInter[i] = cholesky_decompose(VirusInter[i]);
    CD_CoevoInter[i] = cholesky_decompose(CoevoInter[i]);
  }
}
parameters{
  real<lower=0> sigma_virus;
  real<lower=0> sigma_virusphy;
  real<lower=0> sigma_host;
  real<lower=0> sigma_hostphy;
  real<lower=0> sigma_virusinter;
  real<lower=0> sigma_hostinter;
  real<lower=0> sigma_inter;
  real<lower=0> sigma_coevointer;
  real<lower=0> sigma_pool;
  real<lower=0> sigma_spatial;
  real<lower=0> sigma_residual;
  
  //non-centering
  vector[virusspecies] virus_tilde;
  vector[hostspecies] host_tilde;
  vector[virushostcombinations] inter_tilde;
  vector[pools] pool_tilde;
  vector[spatialcompositions] spatial_tilde;
  vector[datapoints] residual_tilde;
  vector[virusspecies] vphy_tilde;
  vector[hostspecies] hphy_tilde;
  vector[virushostcombinations] vinter_tilde;
  vector[virushostcombinations] hinter_tilde;
  vector[virushostcombinations] coevointer_tilde;
  
  //latent mean
  real alpha; 
}
transformed parameters{
  //predicted probability for each datapoint for each matrix
  vector<lower=0, upper=1>[datapoints] predprob[matrices];
  vector<lower=0, upper=1>[datapoints] theta[matrices];
  
  //predicted means for each partially pooled effect
  vector[virusspecies] mean_virus;
  vector[hostspecies] mean_host;
  vector[virushostcombinations] mean_inter;
  vector[pools] mean_pool;
  vector[spatialcompositions] mean_spatial;
  vector[datapoints] mean_residual;
  vector[virusspecies] mean_virusphy[matrices];
  vector[hostspecies] mean_hostphy[matrices];
  vector[virushostcombinations] mean_virusinter[matrices];
  vector[virushostcombinations] mean_hostinter[matrices];
  vector[virushostcombinations] mean_coevointer[matrices];
  
  //vector for calculating marginalised log likelihoods
  vector[matrices] lp;
  
  //add exact correction for marginalisation
  lp = rep_vector(log_unif, matrices);
  
  //non-centred calculation of means for non-phylogenetically associated effects
  mean_virus = sqrt(sigma_virus) * virus_tilde;
  mean_host = sqrt(sigma_host) * host_tilde;
  mean_inter = sqrt(sigma_inter) * inter_tilde;
  mean_pool = sqrt(sigma_pool) * pool_tilde;
  mean_spatial = sqrt(sigma_spatial) * spatial_tilde;
  mean_residual = sqrt(sigma_residual) * residual_tilde;
  
  //non-centered calculation of means for phylogenetically associated effects
  for (i in 1:matrices) {
    mean_virusphy[i] = sqrt(sigma_virusphy) * (CD_VirusPhy[i] * vphy_tilde);
    mean_hostphy[i] = sqrt(sigma_hostphy) * (CD_HostPhy[i] * hphy_tilde);
    mean_virusinter[i] = sqrt(sigma_virusinter) * (CD_VirusInter[i] * vinter_tilde);
    mean_hostinter[i] = sqrt(sigma_hostinter) * (CD_HostInter[i] * hinter_tilde);
    mean_coevointer[i] = sqrt(sigma_coevointer) * (CD_CoevoInter[i] * coevointer_tilde);
  }
  
  //precalculate expected probability for each data point
  for (i in 1:datapoints) {
    for (j in 1:matrices) {
      //predprob saved for later posterior predictive simulations
      predprob[j,i] = inv_logit(alpha + mean_virus[virus[i]] + mean_host[host[i]] + mean_inter[combination[i]] + mean_virusphy[j,virus[i]] + mean_hostphy[j,host[i]] + mean_virusinter[j,combination[i]] + mean_hostinter[j,combination[i]] + mean_coevointer[j,combination[i]] + mean_pool[pool[i]] + mean_spatial[spatialcomposition[i]] + mean_residual[ID[i]]);
      theta[j,i] = 1 - (1 - predprob[j,i])^samples[i];
    }
  }
  
  //calculate model likelihood marginalising over the phylogenetic uncertainty
  if (run_estimation==1) {
    for (i in 1:matrices) {
      //likelihood of at least 1 success in n trials with underlying
      //probability of success theta is modified Bernoulli with success
      //probability = 1-(1-p)^n -- quantity precalculated above in theta
      //see Thompson 1971 Biometrics
      lp[i] = lp[i] + bernoulli_lpmf(y | theta[i]);
    }
  }
}
model {
  alpha ~ logistic(0,1);
  
  vphy_tilde ~ normal(0,1);
  hphy_tilde ~ normal(0,1);
  vinter_tilde ~ normal(0,1);
  hinter_tilde ~ normal(0,1);
  coevointer_tilde ~ normal(0,1);
  virus_tilde ~ normal(0,1);
  host_tilde ~ normal(0,1);
  inter_tilde ~ normal(0,1);
  pool_tilde ~ normal(0,1);
  spatial_tilde ~ normal(0,1);
  residual_tilde ~ normal(0,1);
  
  sigma_virus ~ exponential(1);
  sigma_virusphy ~ exponential(1);
  sigma_host ~ exponential(1);
  sigma_hostphy ~ exponential(1);
  sigma_virusinter ~ exponential(1);
  sigma_hostinter ~ exponential(1);
  sigma_coevointer ~ exponential(1);
  sigma_inter ~ exponential(1);
  sigma_pool ~ exponential(1);
  sigma_spatial ~ exponential(1);
  sigma_residual ~ exponential(1);

  //increment the log posterior by the summed likelihood
  if (run_estimation==1) {
    target += log_sum_exp(lp);
  }
}
generated quantities {  
  //define ICC parameters
  real<lower=0> denominator;
  real<lower=0, upper=1> ICC_virus;
  real<lower=0, upper=1> ICC_virusphy;
  real<lower=0, upper=1> ICC_host;
  real<lower=0, upper=1> ICC_hostphy;
  real<lower=0, upper=1> ICC_virusinter;
  real<lower=0, upper=1> ICC_hostinter;
  real<lower=0, upper=1> ICC_inter;
  real<lower=0, upper=1> ICC_coevointer;
  real<lower=0, upper=1> ICC_pool;
  real<lower=0, upper=1> ICC_spatial;
  real<lower=0, upper=1> ICC_residual;
  real<lower=0, upper=1> ICC_nonphylogenetic;
  real<lower=0, upper=1> ICC_phylogenetic;
  
  //define parameters for data simulation
  int<lower=0, upper=1> y_sim[datapoints];
  
  //define parameters for predicted means
  real<lower=0, upper=1> y_pred[virushostcombinations];

  //define parameters for model comparision 
  vector[matrices] uncertainty_log_lik;
  vector[datapoints] log_lik;
  real calc;
  
  //define for the case where matrices == 1 and calc 
  //would otherwise be undefined
  calc = 0;
  
  //generate denominator for ICC calculations
  denominator = sigma_virus + sigma_virusphy + sigma_host + sigma_hostphy + sigma_virusinter + sigma_hostinter + sigma_coevointer + sigma_inter + sigma_pool + sigma_spatial + sigma_residual + (pi()^2)/3;
  
  //generate proportion of variance explained by each model component, as well as all phylogenetically and all non-phylogenetically associated effects
  ICC_virus = sigma_virus/denominator;
  ICC_virusphy = sigma_virusphy/denominator;
  ICC_host = sigma_host/denominator;
  ICC_hostphy = sigma_hostphy/denominator;
  ICC_virusinter = sigma_virusinter/denominator;
  ICC_hostinter = sigma_hostinter/denominator;
  ICC_inter = sigma_inter/denominator;
  ICC_coevointer = sigma_coevointer/denominator;
  ICC_pool = sigma_pool/denominator;
  ICC_spatial = sigma_spatial/denominator;
  ICC_residual = sigma_residual/denominator;
  ICC_nonphylogenetic = (sigma_virus + sigma_host + sigma_inter + sigma_spatial + sigma_pool)/denominator;
  ICC_phylogenetic = (sigma_virusphy + sigma_hostphy + sigma_virusinter + sigma_hostinter + sigma_coevointer)/denominator;

  //generate pointwise log-likelihoods for PSIS-LOO model comparison - see Vehtari, Gelman and Gabry 2017
  for (i in 1:datapoints) {
    uncertainty_log_lik = rep_vector(log_unif, matrices);
    for (j in 1:matrices) {
      //likelihood of at least 1 success in n trials with underlying probability of success theta is modified Bernoulli with p = 1-(1-theta)^n -- see Thompson 1971 Biometrics
      uncertainty_log_lik[j] = uncertainty_log_lik[j] + bernoulli_lpmf(y[i] | theta[j,i]);
	}
	log_lik[i] = log_sum_exp(uncertainty_log_lik);
  }

  //generate data simulations, using infection probability averaged over
  //each tree (VCV) in the phylogenetic uncertainty case
  if (matrices == 1) {
    for(i in 1:datapoints) {
      y_sim[i] = bernoulli_rng(theta[1,i]);
    }
  } else {
    for(i in 1:datapoints) {
      calc = 0;
      for (j in 1:matrices) {
        calc = calc + predprob[j,i]/matrices;
      }
      y_sim[i] = bernoulli_rng(1 - (1 - calc)^samples[i]);
    }
  }
  
  //generate posterior predicted means, using infection probability averaged over
  //each tree (VCV) in the phylogenetic uncertainty case
  if (matrices == 1) {
    for(i in 1:virushostcombinations) {
      y_pred[i] = inv_logit(alpha + mean_virus[virus_pred[i]] + mean_host[host_pred[i]] + mean_inter[combination_pred[i]] + mean_virusphy[1,virus_pred[i]] + mean_hostphy[1,host[i]] + mean_virusinter[1,combination_pred[i]] + mean_hostinter[1,combination_pred[i]] + mean_coevointer[1,combination_pred[i]] + normal_rng(0,sqrt(sigma_pool)) + normal_rng(0,sqrt(sigma_spatial)) + normal_rng(0,sqrt(sigma_residual)));
    }
  } else {
    for(i in 1:virushostcombinations) {
      calc = 0;
      for (j in 1:matrices) {
        calc = calc + (mean_virusphy[j,virus_pred[i]] + mean_hostphy[j,host[i]] + mean_virusinter[j,combination_pred[i]] + mean_hostinter[j,combination_pred[i]] + mean_coevointer[j,combination_pred[i]])/matrices;
      }
      calc = calc + alpha + mean_virus[virus_pred[i]] + mean_host[host_pred[i]] + mean_inter[combination_pred[i]] + normal_rng(0,sqrt(sigma_pool)) + normal_rng(0,sqrt(sigma_spatial)) + normal_rng(0,sqrt(sigma_residual));
      y_pred[i] = inv_logit(calc);
    }
  }
}
