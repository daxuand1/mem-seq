#' Sequential MEM Selection.
#' @param y A numeric matrix or dissimilarity matrix of outcome variables.
#' @param x A vector, matrix or data frame of covariates.
#' @param mem A numeric matrix of candidate Moran's eigenvectors.
#' @param nperm The number of permutations.
#' @param ncores Number of cores for permutation tests.
#'
#' @return A list containing the selection results.
#' @export
mem_seq = function(y, x, mem, nperm = 999, ncores = 1){

  mem = as.matrix(mem)
  mem = apply(mem, 2, function(x){x / sqrt(sum(x^2))})

  if(is.data.frame(x)){
    x = as.matrix(x)
  }else if(is.vector(x)){
    x = as.matrix(x, ncol = 1)
  }
  x = apply(x, 2, function(x){x - mean(x)})
  n = nrow(x)

  if("dist" %in% class(y)){
    glb_fit = vegan::dbrda(y ~ mem)
  }else{
    if(is.data.frame(y)){
      y = as.matrix(y)
    }else if(is.vector(y)){
      y = as.matrix(y, ncol = 1)
    }
    y = apply(y, 2, function(x){x - mean(x)})
    glb_fit = vegan::rda(y ~ mem)
  }

  an_glb = stats::anova(glb_fit, permutations = permute::how(nperm = nperm),
                 parallel = ncores)
  if(an_glb$`Pr(>F)`[1] <= 0.05){

    if("dist" %in% class(y)){
      d2 = -0.5 * as.matrix(y)^2
      dm = d2 %*% mem
      a_pool = colSums(mem * dm)
      tss = -sum(d2) / n
    }else{
      ym = crossprod(y, mem)   # d × p
      a_pool = colSums(ym^2)
      tss = sum(y^2)
    }

    ord = order(a_pool, decreasing = T)
    a_pool_ord = a_pool[ord]

    r2_pool = cumsum(a_pool_ord) / tss
    r2adj_pool = 1 - (1 - r2_pool) * (n - 1) / (n - 1:ncol(mem) - 1)
    r2adj_full = r2adj_pool[ncol(mem)]

    mem_num_y_r2adjfull = which.max(r2adj_pool >= r2adj_full)
    mem_num_y_r2adjmax = which.max(r2adj_pool)

    mem_names_1 = colnames(mem)[ord[1:mem_num_y_r2adjfull]]

    mem_names_2 =
      colnames(mem)[ord[(mem_num_y_r2adjfull + 1):mem_num_y_r2adjmax]]

    # check MEM for x

    mem_remain = mem[, !(colnames(mem) %in% mem_names_1)]

    xm = crossprod(x, mem_remain)   # d × p
    a_pool = colSums(xm^2)
    tss = sum(x^2)

    ord = order(a_pool, decreasing = T)
    a_pool_ord = a_pool[ord]

    r2_pool = cumsum(a_pool_ord) / tss
    r2adj_pool = 1 - (1 - r2_pool) * (n - 1) / (n - 1:ncol(mem_remain) - 1)
    r2adj_full = r2adj_pool[ncol(mem_remain)]

    mem_num_x_r2adjfull = which.max(r2adj_pool >= r2adj_full)
    mem_names_3 = colnames(mem_remain)[ord[1:mem_num_x_r2adjfull]]
    ind = mem_names_3 %in% mem_names_2

    mem_names = unique(c(mem_names_1, mem_names_3[ind]))
    num_mem = length(mem_names)
  }else{
    num_mem = 0
    mem_names = NULL
  }

  return(list(num_mem = num_mem,
              mem_names = mem_names))
}

#' Forward MEM Selection. A wrapper of mem.select from the adespatial package.
#' @param y A numeric matrix of outcome variables.
#' @param x A vector, matrix or data frame of covariates.
#' @param lw A listw object from the spdep package.
#' @param nperm The number of permutations.
#'
#' @return A list containing the selection results.
#' @export
mem_fwd = function(y, x, lw, nperm = 999){

  ms = adespatial::mem.select(y, lw, "positive", "FWD",
                  nperm.global = nperm,
                  nperm = nperm)
  selected_mem = ms$MEM.select

  if(!is.null(selected_mem)){
    num_mem = ncol(selected_mem)
    mem_names = colnames(selected_mem)
  }else{
    num_mem = 0
    mem_names = NULL
  }

  return(list(num_mem = num_mem,
              mem_names = mem_names))
}

#' Global MEM Selection.
#' @param y A numeric matrix or dissimilarity matrix of outcome variables.
#' @param x A vector, matrix or data frame of covariates.
#' @param mem A numeric matrix of candidate Moran's eigenvectors.
#' @param nperm The number of permutations.
#' @param ncores Number of cores for permutation tests.
#'
#' @return A list containing the selection results.
#' @export
mem_glb = function(y, x, mem, nperm = 999, ncores = 1){

  mem = as.matrix(mem)

  if("dist" %in% class(y)){
    glb_fit = vegan::dbrda(y ~ mem)
  }else{
    if(is.data.frame(y)){
      y = as.matrix(y)
    }else if(is.vector(y)){
      y = as.matrix(y, ncol = 1)
    }
    y = apply(y, 2, function(x){x - mean(x)})
    glb_fit = vegan::rda(y ~ mem)
  }

  an_glb = stats::anova(glb_fit, permutations = permute::how(nperm = nperm),
                 parallel = ncores)
  if(an_glb$`Pr(>F)`[1] <= 0.05){
    num_mem = ncol(mem)
    mem_names = colnames(mem)
  }else{
    num_mem = 0
    mem_names = NULL
  }

  return(list(num_mem = num_mem,
              mem_names = mem_names))
}
