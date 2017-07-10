library(rstan)
library(foreach)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() + 2)

mapCat <- function(f, coll) {
  return(unlist(lapply(coll, f)))
}

get_stan_draws <- function(gen_model, gen_data) {
  gendata <- stan(gen_model,
                  data=gen_data,
                  control=list(max_treedepth = 12),
                  algorithm='Fixed_param')
  gdf <- as.data.frame(gendata)
}

gen_percentiles <- function(gen_model_file, gen_data, data_vars, data_rows, outcome_vars, outcome_rows, 
                            target_model_file, num_replicates=1000) {
  get_stan_fit_quantiles <- function(draw, param_names, target_model, target_data,
                                     num_chains=1, num_iters=1000) {
    library(rstan)
    target_fit = sampling(target_model, data=target_data,
                          chains=num_chains, iter=num_iters,
                          control=list(adapt_delta=0.98, max_treedepth=14));
    target_matrix = as.matrix(target_fit);
    
    rankings = rep(NA, length(param_names));
    mode(rankings) = "numeric";
    for (i in 1:length(param_names)) {
      theta0 = draw[,param_names[i]];
      samples_theta0 = target_matrix[,param_names[i]];
      quantile = sum(samples_theta0 <= theta0)/length(samples_theta0);
      rankings[i] = quantile;
    }
    names(rankings) = param_names;
    return(rankings);
  }
  
  gdf <- get_stan_draws(gen_model_file, gen_data)
  getVarNames <-function(name) {
    return(names(gdf)[grep(paste("^", name, sep=""), names(gdf))])
  }
  
  data_names <- mapCat(getVarNames, data_vars)
  outcome_names <- mapCat(getVarNames, outcome_vars)
  nonParamNames <- union(union(data_names, outcome_names), c("lp__"))
  paramNames <- setdiff(names(gdf), nonParamNames)
  target_model <- stan_model("hierarchical_logistic_regression.stan")
  
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)
  
  replicates_matrix <- foreach(i=1:num_replicates, .combine=rbind) %do% {
    names <- c(data_vars, outcome_vars)
    rows <- c(data_rows, outcome_rows)
    for (j in seq_along(names)) {
      dv <- names[j]
      r <- rows[j]
      if (r > 1) {
        gen_data[[dv]] <-  matrix(unlist(gdf[i,][getVarNames(dv)], use.names=F), nrow=r)
      } else {
        gen_data[[dv]] <-  unlist(gdf[i,][getVarNames(dv)], use.names=F)
      }
    }
    get_stan_fit_quantiles(gdf[i,], paramNames, target_model, gen_data)
  }
  
  colnames(replicates_matrix) = paramNames
  dput(replicates_matrix, file=paste(target_model_file, "matrix", sep="."))
  stopCluster(cl)
}
#===================================================================


setwd("~/scm/mygmo")
d <- list(N=8, K=2, J=2, L=3, jj=c(1, 1, 2, 2, 1, 1, 2, 2))
replicates_matrix <- gen_percentiles("models/hls_gendata.stan", d,
                                     c("x", "u"), c(d$N, d$J),
                                     c("y"), c(1),
                                     "models/hierarchical_logistic_regression.stan")
hist(replicates_matrix[,-c(1, 4)])
