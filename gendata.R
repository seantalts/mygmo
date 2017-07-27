library(rstan)
library(foreach)
library(doParallel)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() + 2)

set.seed(42)

mapCat <- function(f, coll) {
  return(unlist(lapply(coll, f)))
}

get_stan_draws <- function(gen_model, gen_data, num_draws = 4000, num_chains = 4) {
  gendata <- stan(gen_model,
                  data=gen_data,
                  #seed=1234,
                  seed=23489432,
                  chains = num_chains,
                  warmup = 0,
                  iter = num_draws / num_chains,
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

get_posterior_thetas <- function(theta0_draws, target_model, target_data,
                                   num_chains=1, num_iters=1000) {
  library(rstan)
  target_fit = sampling(target_model, data=target_data,
                        chains=num_chains, iter=num_iters,
                        #seed=123433,
                        seed=23489432,
                        control=list(adapt_delta=0.99, max_treedepth=14));
  return(as.data.frame(target_fit));
}

get_quantiles <- function(param_names) {
  return(function(cell) {
    theta0s <- cell[[1]]
    posterior_thetas <- cell[[2]]
    rankings = rep(NA, length(param_names));
    mode(rankings) = "numeric";
    for (i in 1:length(param_names)) {
      theta0 = theta0s[,param_names[i]];
      samples_theta0 = posterior_thetas[,param_names[i]];
      quantile = sum(samples_theta0 > theta0)/length(samples_theta0);
      rankings[i] = quantile;
    }
    names(rankings) = param_names;
    return(rankings);
  })
}

get_param_names <- function(gdf, data_vars) {
  getVarNames <-function(name) {
    return(names(gdf)[grep(paste("^", name, sep=""), names(gdf))])
  }
  data_names <- mapCat(getVarNames, data_vars)
  nonParamNames <- union(data_names, c("lp__"))
  paramNames <- setdiff(names(gdf), nonParamNames)
  return(paramNames)
}

gen_replications <- function(gen_model_file, gen_data, data_vars, data_rows, 
                            target_model_file, num_replicates=1000) {
  gdf <- get_stan_draws(gen_model_file, gen_data, num_draws=num_replicates)
  target_model <- stan_model(target_model_file)
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)
  if (num_replicates > nrow(gdf)) {
    print(paste0("Only ", nrow(gdf), " draws for ", num_replicates, " replicates!"))
    return(NULL)
  }
  replications <- foreach(i=1:num_replicates, #.combine=rbind,
                               .export = c("get_posterior_thetas", "get_param_names",
                                           "mapCat",
                                           "add_draws_to_constants")) %dopar% {
    target_data <- add_draws_to_constants(gdf[i,], data_vars, data_rows, gen_data)
    list(gdf[i,get_param_names(gdf, data_vars)],
         get_posterior_thetas(gdf[i,], target_model, target_data))
  }
  
  dput(replications, file=paste(target_model_file, "replications", sep="."))
  stopCluster(cl)
  return(replications)
}

gen_percentiles <- function(gen_model_file, gen_data, data_vars, data_rows, 
                            target_model_file, num_replicates=1000) {
  replications <- gen_replications(gen_model_file, gen_data, data_vars, data_rows,
                                   target_model_file, num_replicates)
  paramNames <- get_param_names(replications, data_vars)
  quantile_matrix <- vapply(replications, get_quantiles(paramNames))
  return(quantile_matrix)
}

#+================3===

add.alpha <- function(col, alpha=1) {
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
#d <- list(N=25, X=rnorm(25, 0, 5))
#lin_regr_rm <- gen_percentiles("models/gen_lin_regr.stan", d, c("y"), c(1), "models/lin_regr.stan", num_replicates = 10000)
#hist(lin_regr_rm, breaks=500)

#========= Lin Regr constant sigma======
d <- list(N=25, X=rnorm(25, 0, 5))
lin_regr_rm_c2 <- gen_replications("models/gen_lin_regr_c.stan", d, c("y"), c(1), "models/lin_regr_c.stan", num_replicates = 10000)
quants_gt2 <- sapply(lin_regr_rm_c2, get_quantiles(c("alpha", "beta")))
plots <- function(quants_lte, quants_gt, breaks) {
  par(mfrow=c(3,2))
  hist(quants_lte["alpha",], breaks = breaks)
  hist(quants_gt["alpha",], breaks = breaks)
  hist(quants_lte["beta",], breaks = breaks)
  hist(quants_gt["beta",], breaks = breaks)
  hist(quants_lte, breaks = breaks)
  hist(quants_gt, breaks = breaks)
}

pdf("sassy.pdf")
widths <- c(0.05, 0.01, 0.005)
for (binwidth in widths) {
  plots(quants_lte, quants_gt, seq(0, 1, binwidth))
}
for (binwidth in widths) {
  plots(quants_lte2, quants_gt2, seq(0, 1, binwidth))
}
dev.off()


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
