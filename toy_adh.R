rm(list=ls())
suppressPackageStartupMessages({
  library(stringr)
  library(spdep)
  library(dirmult)
  library(adespatial)
  library(mvtnorm)
  library(vegan)
  library(MEMseq)
})

# ---- parameter setup ----
ad = T # alpha or beta analysis
p_mix = 1 # prob of tumor in tissue area, and vice versa
ad_type = "h" # entropy
bd_type = "bray" # Bray-Curtis
method = "seq" # mem selection method
spa_pat = simulate_spa_pat() # generate simulation spatial pattern
coords = spa_pat[[1]] # spatial coordinates
region = spa_pat[[2]] # tumor/tissue indicator

# Specify autocorrelation, Gaussian process with squared-exponential kernel
dmat = as.matrix(stats::dist(coords))
sigma = exp(- (dmat^2) / (2 * (2^2))) # autocorrelation variance

nsim = 500 # number of Monte Carlo replication
seed_start = 0 # seed
p = 50 # number of taxa
gamma = 0 # region effect in alpha diversity analysis
gamma_b = 0 # region effect in beta diversity analysis
alpha0 = 1 # baseline taxa intensity
L = 2000 # baseline library size
nperm = 999 # number of permutation
set.seed(0)
w = rnorm(p, mean = 0, sd = 1) # factor of spatial autocorrelation on each taxa
w2 = rnorm(p, mean = 0, sd = 1) # factor of region effect on each taxa

# Specify MEM
nb = dnearneigh(coords, 0, 1)
lw = nb2listw(nb)
mem = mem(lw, MEM.autocor = "positive")

res_df = NULL

# ---- Monte Carlo simulation ----

cat("Null case:\n")

timestart = Sys.time()
res_df = rep()

for(i in 1:nsim + seed_start){
  seed = i
  suppressMessages({
    out = simulate_one(region = region, sigma = sigma, ad = ad,
                       ad_type = bd_type, bd_type = bd_type,
                       L = L, p_mix = p_mix,
                       method = method,
                       p = p, alpha0 = alpha0,
                       gamma = gamma, gamma_b = gamma_b,
                       w = w, w2 = w2,
                       lw = lw, mem = mem, coords = coords,
                       nperm = nperm, seed = seed)
  })
  df_i = data.frame(p_emp = out$p_emp,
                    num_mem = out$num_mem,
                    select_time = out$select_time,
                    analysis_time = out$analysis_time,
                    total_time = out$total_time)
  res_df = rbind(res_df, df_i)
}

cat("Monte Carlo simulation time:\n")
Sys.time() - timestart # 5 mins

cat("P-value:\n")
summary(c(res_df$p_emp))
cat("Type I error:\n")
mean(res_df$p_emp < 0.05) # 0.03

cat("Number of selected MEM:\n")
summary(c(res_df$num_mem))
cat("Selection time:\n")
summary(c(res_df$select_time))
cat("Analysis time:\n")
summary(c(res_df$analysis_time))
cat("Total (selection + ANOVA) time:\n")
summary(c(res_df$total_time))

