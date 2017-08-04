library(rstan)
library(foreach)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() + 2)

set.seed(42)

mapCat <- function(f, coll) {unlist(lapply(coll, f))}

ensure_model = function(model) { 
  if (typeof(model) == "character")
    model = stan_model(model)
  model
}

get_stan_draws <- function(gen_model, gen_data, num_draws = 4000, num_chains = 4) {
  gendata <- sampling(ensure_model(gen_model),
                  data=gen_data,
                  #seed=1234,
                  #seed=23489432,
                  chains = num_chains,
                  warmup = 0,
                  iter = num_draws / num_chains,
                  control=list(max_treedepth = 12),
                  algorithm='Fixed_param')
  as.data.frame(gendata)
}

add_draws_to_constants <- function(draws, data_vars, data_rows, constants) {
  getVarNames <-function(name) {
    names(draws)[grep(paste("^", name, sep=""), names(draws))]
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
  constants
}

get_posterior_thetas <- function(target_model, target_data,
                                   num_chains=1, num_iters=2000) {
  library(rstan)
  target_model = ensure_model(target_model)
  target_fit = sampling(target_model, data=target_data,
                        chains=num_chains, iter=num_iters,
                        #seed=123433,
                        #seed=23489432,
                        control=list(adapt_delta=0.99, max_treedepth=14));
  sampler_params <- get_sampler_params(target_fit, inc_warmup = FALSE)
  div <- sum(sapply(sampler_params, function(s) {return(sum(s[,"divergent__"]))}))
  s <- summary(target_fit)$summary
  model_file = target_model@model_name
  write(paste(s[, "n_eff"], collapse=", "), file=paste0(model_file, ".n_effs"), append=T)
  if (div > 10)
    stop(paste0("Too many divergences! ", div))
  list(as.data.frame(target_fit), div)
}

get_quantiles <- function(param_names, op) {
  function(cell) {
    theta0s <- cell[[1]]
    posterior_thetas <- cell[[2]]
    rankings = rep(NA, length(param_names));
    mode(rankings) = "numeric";
    for (i in 1:length(param_names)) {
      theta0 = theta0s[,param_names[i]];
      samples_theta0 = posterior_thetas[,param_names[i]];
      quantile = mean(op(samples_theta0, theta0))
      rankings[i] = quantile;
    }
    names(rankings) = param_names;
    rankings
  }
}

get_param_names <- function(gdf, data_vars) {
  getVarNames <-function(name) {
    return(names(gdf)[grep(paste("^", name, sep=""), names(gdf))])
  }
  data_names <- mapCat(getVarNames, data_vars)
  nonParamNames <- union(data_names, c("lp__"))
  setdiff(names(gdf), nonParamNames)
}

gen_replications <- function(gen_model_file, gen_data, data_vars, data_rows, 
                            target_model_file, num_replicates=1000) {
  gdf <- get_stan_draws(gen_model_file, gen_data, num_draws=num_replicates)

  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)
  if (num_replicates > nrow(gdf)) {
    stop(paste0("Only ", nrow(gdf), " draws for ", num_replicates, " replicates!"))
  }
  target_model = ensure_model(target_model_file)
  replications <- foreach(i=1:num_replicates, #.combine=rbind,
                               .export = c("get_posterior_thetas", "get_param_names",
                                           "mapCat", "ensure_model",
                                           "add_draws_to_constants")) %dopar% {
    target_data <- add_draws_to_constants(gdf[i,], data_vars, data_rows, gen_data)
    results <- get_posterior_thetas(target_model, target_data)
    list(theta0=gdf[i,get_param_names(gdf, data_vars)], posterior_thetas=results[[1]])
  }
  
  dput(replications, file=paste(target_model@model_name, "replications", sep="."))
  stopCluster(cl)
  replications
}

gen_percentiles <- function(gen_model_file, gen_data, data_vars, data_rows, 
                            target_model_file, num_replicates=1000) {
  replications <- gen_replications(gen_model_file, gen_data, data_vars, data_rows,
                                   target_model_file, num_replicates)
  paramNames <- get_param_names(replications, data_vars)
  sapply(replications, get_quantiles(paramNames))
}

#+================3===

add.alpha <- function(col, alpha=1) {
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdfreplications <- function(replications, data_vars, name) {
  paramNames <- get_param_names(replications[[1]][[2]], data_vars)
  quants_lte <- sapply(replications, get_quantiles(paramNames, function(x, y) {return(x <= y)}))
  quants_gt <- sapply(replications, get_quantiles(paramNames, function(x, y) {return(x > y)}))
  plots <- function(quants_lte, quants_gt, breaks) {
    if (length(paramNames) > 4)
      paramNames <- paramNames[0:4]
    par(mfrow=c(length(paramNames) + 1, 2))
    par(mar=c(4,2,2,1)+0.1)
    if (length(paramNames) > 1) {
      for (pn in paramNames) {
        hist(quants_lte[pn,], breaks = breaks, main=NULL, xlab=pn)
        hist(quants_gt[pn,], breaks = breaks, main=NULL, xlab=pn)
      }
    }

    hist(quants_lte, breaks = breaks, main=NULL)
    hist(quants_gt, breaks = breaks, main=NULL)
  }
  
  pdf(paste(name, "pdf", sep="."))
  widths <- c(0.05, 0.01, 0.005)
  for (binwidth in widths) {
    plots(quants_lte, quants_gt, seq(0, 1, binwidth))
  }
  dev.off()
}

inlapdfreplications <- function(replications, data_vars, name) {
  paramNames <- get_param_names(replications[[1]][[2]], data_vars)
  quants_lte <- sapply(replications, get_quantiles(paramNames, function(x, y) {return(x <= y)}))
  quants_gt <- sapply(replications, get_quantiles(paramNames, function(x, y) {return(x > y)}))
  plots <- function(quants_lte, quants_gt) {
    if (length(paramNames) > 4)
      paramNames <- paramNames[0:4]
    par(mfrow=c(length(paramNames) + 1, 2))
    par(mar=c(4,2,2,1)+0.1)
    for (pn in paramNames) {
      INLA::inla.ks.plot(quants_lte[pn,], punif)
      INLA::inla.ks.plot(quants_gt[pn,], punif)
    }
    INLA::inla.ks.plot(quants_lte, punif)
    INLA::inla.ks.plot(quants_gt, punif)
  }

  pdf(paste("inla", name, "pdf", sep="."))
  plots(quants_lte, quants_gt)
  dev.off()
}

#===================================================================

setwd("~/scm/mygmo") #XXX probably need to change this for your machine

# ============= Dan unknown mean and variance
# N_reps= 10000;
# N_samp = 1000;
# 
# prior_shape = 1;
# prior_rate = 1e-3;
# n_data=100;
# prior_mean=0;
# prior_prec=1;
# 
# do_posterior_prop_unknown_mean_unknown_prec = function(n_data, N_samp, true_mean , true_prec, prior_mean, prior_prec,prior_shape,prior_rate) {
#   data = rnorm(n_data, mean = true_mean, sd = 1/sqrt(true_prec))
#   x_bar = mean(data)
#   ss_x = sum((data - x_bar)^2)
#   
#   post_mean = (prior_prec*prior_mean + n_data*x_bar)/(prior_prec+n_data)
#   post_prec = prior_prec + n_data
#   post_shape = prior_shape + n_data/2
#   post_rate = prior_rate +ss_x/2 +0.5*n_data*prior_prec/(n_data + prior_prec)*(x_bar - prior_mean)^2
#   
#   prec_samps = rgamma(N_samp,shape=post_shape,rate=post_rate)
#   mean_samps = rnorm(N_samp,mean=post_mean, sd = 1/sqrt(post_prec*prec_samps))
#   
#   prec_quant = mean(prec_samps <= true_prec)
#   mean_quant = mean(mean_samps <= true_mean)
#   
#   return( list(mean= mean_quant, prec=prec_quant))
# }
# 
# output_mean = rep(NA,N_reps)
# output_prec = rep(NA,N_reps)
# for(i in 1:N_reps) { 
#   true_prec = rgamma(1,shape=prior_shape,rate=prior_rate)
#   true_mean = rnorm(1, mean=prior_mean, sd = 1/sqrt(true_prec))
#   tmp = do_posterior_prop_unknown_mean_unknown_prec(n_data, N_samp, true_mean , true_prec, prior_mean, prior_prec,prior_shape,prior_rate)
#   output_mean[i] =  tmp$mean
#   output_prec[i] = tmp$prec
# }
# 
# INLA::inla.ks.plot(output_mean,punif)
# INLA::inla.ks.plot(output_prec,punif)

#========== HLR ============
# d <- list(N=8, K=2, J=2, L=3, jj=c(1, 1, 2, 2, 1, 1, 2, 2))
# d$x <- matrix(rnorm(d$N * d$K, 0, 5), nrow=d$N)
# d$u <- matrix(rnorm(d$J * d$L, 0, 5), nrow=d$J)
# hlr_reps <- gen_replications("models/hls_gendata.stan", d, c("y"), c(1), "models/hierarchical_logistic_regression.stan", num_replicates = 256)
# hlr_reps_minus <- lapply(hlr_reps, function(s) {return(list(s[[1]], s[[2]][,-c(5, 7, 8, 22, 25)]))})
# pdfreplications(hlr_reps_minus, c("y"), "hlrm")
# inlapdfreplications(hlr_reps_minus, c("y"), "hlrm")

#========= 8 schools NCP ============
# d8 <- list(J=8, K=2, sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
# school_ncp_reps <- gen_replications("models/gen_8schools.stan", d8, c("y"), c(1), "models/8schools.stan")
# pdfreplications(school_ncp_reps, c("y"), "8schools_ncp")
# inlapdfreplications(school_ncp_reps, c("y"), "8schools_ncp")

#========= 8 schools NCP off prior============
#d8 <- list(J=8, K=2, sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
#school_ncp_reps <- gen_replications("models/gen_8schools_p.stan", d8, c("y"), c(1), "models/8schools.stan")
#school_ncp_off_reps <- school_ncp_reps
#pdfreplications(school_ncp_off_reps, c("y"), "8schools_ncp_off")
#inlapdfreplications(school_ncp_off_reps, c("y"), "8schools_ncp_off")

#========= Lin Regr ====================
#d <- list(N=25, X=rnorm(25, 0, 5))
#lin_regr_rm <- gen_replications("models/gen_lin_regr.stan", d, c("y"), c(1), "models/lin_regr.stan", num_replicates = 1000)
#pdfreplications(lin_regr_rm, c("y"), "lin_regr")
#inlapdfreplications(lin_regr_rm, c("y"), "lin_regr")

#========= Lin Regr constant sigma=======
d <- list(N=25, X=rnorm(25, 0, 5))
gen_lin_regr_c_model <- stan_model("models/gen_lin_regr_c.stan")
lin_regr_c_model <- stan_model("models/lin_regr_c.stan")
for (i in 11:100) {
  lin_regr_rm_c2 <- gen_replications(gen_lin_regr_c_model, d, c("y"), c(1), lin_regr_c_model, num_replicates = 1000)
  pdfreplications(lin_regr_rm_c2, c("y"), paste0("lin_regr_c", i))
}
x <- 1
lin_regr_rm_c20 <- gen_replications(gen_lin_regr_c_model, d, c("y"), c(1), lin_regr_c_model, num_replicates = 1000*x)
pdfreplications(lin_regr_rm_c2, c("y"), paste0("lin_regr_c", x))
inlapdfreplications(lin_regr_rm_c2, c("y"), paste0("lin_regr_c", x))

#========= Lin Regr fabs sigma======
#d <- list(N=25, X=rnorm(25, 0, 5))
#lin_regr_rm_fabs <- gen_replications("models/gen_lin_regr_fabs.stan", d, c("y"), c(1), "models/lin_regr.stan", num_replicates = 1000)
#pdfreplications(lin_regr_rm_fabs, c("y"), "lin_regr_fabs")

#========= 8 Schools CP - see divergences, spike at 1 for tau
# d8 <- list(J=8, K=2, sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
# school_cp_reps <- gen_replications("models/gen_8schools.stan", d8, c("y"), c(1), "models/8schools_centered.stan")
#pdfreplications(school_cp_reps, c("y"), "8schools_cp")
#inlapdfreplications(school_cp_reps, c("y"), "8schools_cp")

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
