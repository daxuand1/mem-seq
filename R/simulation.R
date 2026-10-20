#' Simulation data spatial pattern. Refer to the paper.
#' @param a1 The width of the large rectangle.
#' @param b1 The height of the large rectangle.
#' @param a2 The width of the rectangles on the foremost left and right.
#' @param a3 The width of the tumor rectangles.
#' @param b2 The height of the rectangle at the top/bottom.
#'
#' @return A list containing the spatial pattern results.
#' @export
simulate_spa_pat = function(a1 = 10, b1 = 16, a2 = 2,
                            a3 = 2, b2 = 3){
  b3 = b1 - 2 * b2
  a4 = a1 - 2 * a2 - 2 * a3

  dt = rbind(cbind(x = rep(1:a1, b2),
                   y = rep(1:b2, each = a1),
                   fac = 0), # bottom
             cbind(x = rep(1:a2, b3),
                   y = rep(1:b3 + b2, each = a2),
                   fac = 0), # mid left normal
             cbind(x = rep(1:a3 + a2, b3),
                   y = rep(1:b3 + b2, each = a3),
                   fac = 1), # mid left tumor
             cbind(x = rep(1:a3 + a2 + a3 + a4, b3),
                   y = rep(1:b3 + b2, each = a3),
                   fac = 1), # mid right tumor
             cbind(x = rep(1:a2 + a2 + a3 + a4 + a3, b3),
                   y = rep(1:b3 + b2, each = a2),
                   fac = 0), # mid right normal
             cbind(x = rep(1:a1, b2),
                   y = rep(1:b2 + b2 + b3, each = a1),
                   fac = 0) # upper
  )
  coords = dt[, 1:2]
  region = dt[, 3]

  return(list(coords = coords,
              region = region))
}

#' Single simulation.
#' @param region Tumor/tissue indicator.
#' @param sigma Autocorrelation variance.
#' @param ad Indicator for alpha or beta analysis. True or False.
#' @param ad_type Type of alpha diversity metric. It can be "ls" (library size), "h" (entropy), or "s" (Simpson diversity).
#' @param bd_type Type of beta diversity metric used in vegdist in the vegan package.
#' @param L Baseline library size.
#' @param p_mix Probability of tumor in tissue area, and vice versa.
#' @param method MEM selection method. It can be "seq" (sequential), "fwd" (forward selection), and "glb" (global).
#' @param p Number of taxa.
#' @param alpha0 Baseline taxa intensity.
#' @param gamma Region effect in alpha diversity analysis.
#' @param gamma_b Region effect in beta diversity analysis.
#' @param w Factor of spatial autocorrelation on each taxa.
#' @param w2 Factor of region effect on each taxa.
#' @param lw A listw object from the spdep package.
#' @param mem A numeric matrix of candidate Moran's eigenvectors.
#' @param coords Spatial coordinates of spots.
#' @param nperm Number of permutation.
#' @param seed Random seed.
#'
#' @return A list containing the simulation results.
#' @export
simulate_one = function(region, sigma, ad = T,
                        ad_type = "h", bd_type = "bray",
                        L = 2000, p_mix = 1,
                        method = "seq",
                        p = 50, alpha0 = 1,
                        gamma = 0, gamma_b = 0,
                        w = 1, w2 = 1,
                        lw = NULL, mem = NULL, coords = NULL,
                        nperm = 999, seed = 1) {

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop(
      "Package 'mvtnorm' is required for simulate_one(). ",
      "Install it with install.packages('mvtnorm').",
      call. = FALSE
    )
  }

  if (!requireNamespace("dirmult", quietly = TRUE)) {
    stop(
      "Package 'dirmult' is required for simulate_one(). ",
      "Install it with install.packages('dirmult').",
      call. = FALSE
    )
  }

  # seed = 15
  set.seed(seed)

  # number of sample
  n = length(region)

  # base concentration parameters for taxa
  if(length(alpha0) == 1) alpha_vec = rep(alpha0, p) else alpha_vec = alpha0

  # mixture of tumor and normal tissue
  region_final = region
  if(p_mix != 1){
    region_final[region == 0] = stats::rbinom(sum(region == 0), 1, 1 - p_mix)
    region_final[region == 1] = stats::rbinom(sum(region == 1), 1, p_mix)
  }

  # latent sample effects b ~ MVN(0, sigma)
  b = as.numeric(mvtnorm::rmvnorm(1, sigma = sigma))

  # library size
  lib_size = round(L + region_final * 100 * gamma +
                     b * 100 + stats::rnorm(n, 0, 100))

  counts = matrix(NA, nrow = n, ncol = p)
  for(i in seq_len(n)){
    alpha_i = alpha_vec * exp(w * b[i] + gamma_b * w2 * region_final[i])
    p_i = dirmult::rdirichlet(1, alpha_i)
    counts[i, ] = as.numeric(stats::rmultinom(1, size = lib_size[i], prob = p_i))
  }
  rownames(counts) = paste0("s", seq_len(n))
  colnames(counts) = paste0("t", seq_len(p))

  if(ad){
    y = apply(counts, 1, function(x){

      if(ad_type == "ls") return(sum(x))

      if(ad_type == "r") return(sum(x > 0))

      x = x[x > 0]
      prop = x / sum(x)
      if(ad_type == "s"){
        return(1 - sum(prop^2))
      }else{
        return(- sum(prop * log(prop)))
      }
    })
    y = y - mean(y)
    outcome = y
  }else{

    if(bd_type == "jaccard"){
      outcome = vegan::vegdist(counts, "canberra", T)
    }else{
      outcome = vegan::vegdist(counts, bd_type)
    }

    y = outcome
  }

  timestart = Sys.time()

  if(method == "null"){
    fm = stats::as.formula("outcome ~ region_final")
    num_mem = 0
    select_time = 0
  }else{
    if(method == "fwd"){
      D = as.matrix(outcome)
      centerM = diag(n) - 1/n
      K = -0.5 * centerM %*% (D * D) %*% centerM
      eK = eigen(K, symmetric = TRUE)
      ord = order(abs(eK$values), decreasing = T)
      ind = which(cumsum(abs(eK$values)[ord]) >=
                    sum(abs(eK$values)) * 0.9)[1]
      val = sqrt(abs(eK$values[ord[1:ind]]))
      vec = eK$vectors[, ord[1:ind]]
      y = t(t(vec) * val)

      if(is.null(lw)){
        nb = spdep::dnearneigh(coords, 0, 1)
        lw = spdep::nb2listw(nb)
      }
      ms = mem_fwd(y, region_final, lw, nperm)
    }else if(method == "glb"){
      if(is.null(mem)){
        nb = spdep::dnearneigh(coords, 0, 1)
        lw = spdep::nb2listw(nb)
        mem = mem(lw, MEM.autocor = "positive")
      }
      ms = mem_glb(y, region_final, mem, nperm)
    }else{
      if(is.null(mem)){
        nb = spdep::dnearneigh(coords, 0, 1)
        lw = spdep::nb2listw(nb)
        mem = mem(lw, MEM.autocor = "positive")
      }
      ms = mem_seq(y, region_final, mem, nperm)
    }

    num_mem = ms$num_mem
    if(num_mem == 0){
      fm = stats::as.formula("outcome ~ region_final")
    }else{
      cond_terms = paste(ms$mem_names, collapse = " + ")
      fm = stats::as.formula(paste("outcome ~ region_final + Condition(",
                                   cond_terms, ")", sep = ""))
    }
    select_time = difftime(Sys.time(), timestart, units = "mins")
  }

  timestart = Sys.time()
  if(ad){
    mod_obs = vegan::rda(fm, data = mem)
  }else{
    mod_obs = vegan::dbrda(fm, data = mem)
  }
  an_obs = stats::anova(mod_obs, permutations = nperm)
  p_emp = an_obs$`Pr(>F)`[1]
  p_emp = ifelse(is.na(p_emp), 1, p_emp)

  analysis_time = difftime(Sys.time(), timestart, units = "mins")

  total_time = select_time + analysis_time

  list(p_emp = p_emp,
       num_mem = num_mem,
       select_time = select_time,
       analysis_time = analysis_time,
       total_time = total_time)
}
