library(rstan)
library(foreach)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() + 2)

set.seed(42)

mapCat <- function(f, coll) {
  return(unlist(lapply(coll, f)))
}

get_stan_draws <- function(gen_model, gen_data) {
  gendata <- stan(gen_model,
                  data=gen_data,
                  seed=1234,
                  control=list(max_treedepth = 12),
                  algorithm='Fixed_param')
  return(as.data.frame(gendata))
}

add_draws_to_constants <- function(draws, data_vars, data_rows, constants) {
  getVarNames <-function(name) {
    return(names(draws)[grep(paste("^", name, sep=""), names(draws))])
  }
  
  for (j in seq_along(data_vars)) {
    dv <- data_vars[j]
    r <- data_rows[j]
    if (r > 1) {
      constants[[dv]] <-  matrix(unlist(draws[getVarNames(dv)], use.names=F), nrow=r)
    } else {
      constants[[dv]] <-  unlist(draws[getVarNames(dv)], use.names=F)
    }
  }
  return(constants)
}

get_stan_fit_quantiles <- function(draw, param_names, target_model, target_data,
                                   num_chains=1, num_iters=1000) {
  library(rstan)
  target_fit = sampling(target_model, data=target_data,
                        chains=num_chains, iter=num_iters,
                        seed=123433,
                        control=list(adapt_delta=0.99, max_treedepth=14));
  target_matrix = as.matrix(target_fit);
  
  rankings = rep(NA, length(param_names));
  mode(rankings) = "numeric";
  for (i in 1:length(param_names)) {
    theta0 = draw[,param_names[i]];
    samples_theta0 = target_matrix[,param_names[i]];
    quantile = sum(samples_theta0 > theta0)/length(samples_theta0);
    rankings[i] = quantile;
  }
  names(rankings) = param_names;
  return(rankings);
}

gen_percentiles <- function(gen_model_file, gen_data, data_vars, data_rows, 
                            target_model_file, num_replicates=1000) {
  gdf <- get_stan_draws(gen_model_file, gen_data)
  getVarNames <-function(name) {
    return(names(gdf)[grep(paste("^", name, sep=""), names(gdf))])
  }
  data_names <- mapCat(getVarNames, data_vars)
  nonParamNames <- union(data_names, c("lp__"))
  paramNames <- setdiff(names(gdf), nonParamNames)
  target_model <- stan_model(target_model_file)
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)

  replicates_matrix <- foreach(i=1:num_replicates, .combine=rbind,
                               .export = c("get_stan_fit_quantiles",
                                           "add_draws_to_constants")) %dopar% {
    target_data <- add_draws_to_constants(gdf[i,], data_vars, data_rows, gen_data)
    get_stan_fit_quantiles(gdf[i,], paramNames, target_model, target_data)
  }
  
  colnames(replicates_matrix) = paramNames
  dput(replicates_matrix, file=paste(target_model_file, "matrix", sep="."))
  stopCluster(cl)
  return(replicates_matrix)
}

#+================3===


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
#===================================================================

setwd("~/scm/mygmo") #XXX probably need to change this for your machine


#========== HLR ============
#d <- list(N=8, K=2, J=2, L=3, jj=c(1, 1, 2, 2, 1, 1, 2, 2))
#d$x <- matrix(rnorm(d$N * d$K, 0, 5), nrow=d$N)
#d$u <- matrix(rnorm(d$J * d$L, 0, 5), nrow=d$J)
#data_vars <- c("y")
#data_rows <- c(1)
#rm <- gen_percentiles("models/hls_gendata.stan", d,
#                      data_vars, data_rows,
#                      "models/hierarchical_logistic_regression.stan",
#                      num_replicates = 256)
#hist(rm[,-c(5, 7, 9, 12)], breaks=400)

#========= 8 schools NCP ============
#d <- list(J=8, K=2, sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#rm2 <- gen_percentiles("models/gen_8schools.stan", d, c("y"), c(1), "models/8schools.stan")
#hist(rm2, breaks=400)

#========= Lin Regr ======
d <- list(N=25, X=rnorm(25, 0, 5))
lin_regr_rm <- gen_percentiles("models/gen_lin_regr.stan", d, c("y"), c(1), "models/lin_regr.stan")
hist(lin_regr_rm, breaks=400)

#========= 8 Schools CP - see divergences, spike at 1 for tau
#d <- list(J=8, K=2, sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#rm <- gen_percentiles("models/gen_8schools.stan", d,
#                      c("y"), c(1),
#                      "models/8schools_centered.stan")

# See how many divergences centered and non-centered (beta); if both -
#   0. Plot divergent and non-divergent separately on marginal plots - use pairs plots
#   1. maybe tighten priors on gamma 

#draws <- get_stan_draws("models/hls_gendata.stan", d)
#target_data <- add_draws_to_constants(draws[1,], data_vars, data_rows, d)
#fit <- stan(file="models/hierarchical_logistic_regression.stan", data=target_data,
#            control=list(adapt_delta=0.95, max_treedepth=15),
#            iter = 2000, warmup = 1000, chains = 4, seed=12345)
#params <- c("z", "L_Omega", "tau", "gamma", "sigma")
#for (i in 1:length(params)) {
#  for (j in (i+1):length(params)) {
#    if (j > i && j <= length(params)) {
#      print(paste0(params[i], params[j]))
#      png(paste0(params[i], params[j], ".png"), width = 2000, height = 1500)
#      pairs(fit, pars=c(params[i], params[j]))
#      dev.off()
#    }
#  }
#}

#c1 <- read.csv2("diagnostic_6.csv", comment.char = "#", sep = ",", dec = ".")
#c1 <- c1[1000:nrow(c1),]
#c1 <- c1[order(c1$divergent__),]
#colors <- c(add.alpha(c("maroon"), 0.02),
#            add.alpha(c("lime green"), 0.4))[unclass(c1$divergent__ + 1)]
#syms <- 21
#pairs(~ sigma + tau.1 + tau.2 + gamma.1.1 + L_Omega.1 + z.1.1, data=c1, pch = syms,
#      col = colors, bg = colors)
#sum(c1$divergent)
#pairs(~ tau.2 + gamma.1.1 + gamma1.2 + gamma1.3 + gamma2.1 + gamma2.2 + L_Omega.1 + z.1.1, data=c1, pch = syms,
#      col = colors, bg = colors)


# ============== fabs vs while loop test
#m <- as.matrix(stan("models/fabs.stan", algorithm = 'Fixed_param', chains = 4, iter = 1e8/2))

#fhist <- hist(m[,'fabs_normal'], plot=F, breaks = (0:1000) * 30 / 1000.0)
#whist <- hist(m[,'while_normal'], plot=F, breaks = (0:1000) * 30 / 1000.0)
#plot(fhist, col=add.alpha('blue', 0.2), main="", xlab="theta.1", yaxt='n', ann=FALSE, xlim = c(0, 1))
#plot(whist, col=add.alpha('red', 0.2), add=T)

#minus <- m[,'fabs_normal'] - m[,'while_normal']
#hist(minus)
#quantile(minus, probs=seq(5, 95, 5)/100)
