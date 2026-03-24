#' Generate a valid probability table for attribute Xa
#'
#' @keywords internal
basGen <- function(x){
  aux <- runif(x[1])
  aux <- aux / sum(aux)
  aux
}

#' Generate a valid conditional probability table for
#' attribute Xa given its parent Xb
#' This is a two entry table
#'
#' @keywords internal
parGen <- function(x){
  p <- list()
  for (i in 1:x[2]){
    aux <- runif(x[1])
    aux <- aux / sum(aux)
    p[[i]] <- aux
  }
  unlist(p)
}


#' Generate a valid conditional probability table for
#' attribute Xa given its parents Xb and Xc
#' This a three entry table
#'
#' @keywords internal
PARGEN <- function(x){
  lon <- length(x)
  p <- list()
  if (lon == 1){
    p <- basGen(x)
  }else if(lon == 2){
    p <- parGen(x)
  }else{
    p <- replicate(x[3], parGen(x), simplify = FALSE)
  }
  unlist(p)
}


#' Generates an individual (potential solution)
#' @param W Number of nodes in BN structure
#' @param X Cardinality of CPTs
#' @param Y Number of parameters per node
#' @param Z BN Structure with CPT slots
#'
#' @keywords internal
individual <- function(W, X, Y, Z){
  ind <- as.vector(unlist(lapply(X, PARGEN)))
  unname(unlist(lapply(X, PARGEN)))
  k <- 1
  id_params <- grep(".params", names(Z))
  for (i in 1:W){
    for (j in 1:Y[[i]]){
      Z[[id_params]][[i]][j] <- ind[k]
      k <- k + 1
    }
  }
  Z
}

#' Replace a vector of parameters in a CPT
#'
#' @keywords internal
vec2net <- function(vec, W, X, Y, Z){
  ind <- vec
  k <- 1
  id_params <- grep(".params", names(Z))
  for (i in 1:W){
    for (j in 1:Y[[i]]){
      Z[[id_params]][[i]][j] <- ind[k]
      k <- k + 1
    }
  }
  Z
}

#' Fast version of vec2net using precomputed offsets
#'
#' @param vec parameter vector
#' @param Z template bnc_bn object
#' @param id_params index of .params in Z
#' @param offsets precomputed list with start/end indices per node
#' @keywords internal
vec2netFast <- function(vec, Z, id_params, offsets){
  for (i in seq_along(offsets)) {
    Z[[id_params]][[i]][] <- vec[offsets[[i]]]
  }
  Z
}

#' Precompute offsets for vec2netFast
#'
#' @keywords internal
precompute_offsets <- function(Y) {
  offsets <- vector("list", length(Y))
  k <- 1L
  for (i in seq_along(Y)) {
    offsets[[i]] <- k:(k + Y[[i]] - 1L)
    k <- k + Y[[i]]
  }
  offsets
}

#' Retrieve Families for a Given Node in Bayesian Network
#'
#' This function identifies and returns the family for a given node in a Bayesian network.
#' The family is a set of nodes, which includes the node itself, its parents, and the class.
#'
#' @param node A node in the bayesian network.
#' @param df An adjacency list representing the network, formatted as a data frame.
#' @param class.str The name of the class node.
#'
#' @keywords internal
#' @return A vector containing the node and its family members.
get_family <- function(node, df, class.str) {
  # Identify parents of the given node.
  parents <- df$from[df$to == node]

  # Exclude the class node from parents, and ensure it's the last family member.
  parents <- setdiff(parents, class.str)
  if (node != class.str) {
    parents <- c(parents, class.str)
  }

  c(node, parents)
}



#' Generates a population of potential solutions
#' @param pobSize Number of desired individuals at initial population
#'
#' @keywords internal
population <- function(pobSize, W, X, Y, Z){
  pob <- list()
  pob <- replicate(pobSize, individual(W, X, Y, Z), simplify = FALSE)
  pob
}



#' Transform a matrix into a vector
#'
#' @keywords internal
vec <- function(mt){unname(unlist(mt))}


#' Reparation to satisfy the constraints for (\eqn{\theta_{a, b}}) in [0 , 1]
#'
#' @keywords internal
reflect <- function(vector){
  l <- which(vector < 0)
  u <- which(vector > 1)
  if (length(c(l, u)) <= 0){
    vector
  }else{
    vector[l] <- abs(vector[l]) %% 1
    vector[u] <- 1 - vector[u] %% 1
  }
  vector
}


#' Repair to keep the sum of row vectors equal to 1
#' In a two entry way table
#'
#' @keywords internal
lev2 <- function(x, index){
  lista <- c()
  for (i in 1:x[2]){
    lista <- c(lista, rep(index, x[1]))
    index <- index + 1
  }
  lista
}

#' Repair to keep the sum of row vectors equal to 1
#' In a three entry way table
#'
#' @keywords internal
lev3 <- function(x, index){
  lista <- c()
  for (i in 1:x[3]){
    for (j in 1:x[2]){
      lista <- c(lista, rep(index, x[1]))
      index <- index + 1
    }
  }
  lista
}

#' Repair to keep the sum of row vectors equal to 1
#' Gives BN structure
#'
#' @keywords internal
strRep <- function(x){
  index <- 1
  lista <- c()
  for (i in 1:length(x)){
    depth <- length(x[[i]])
    if (depth == 2){
      lista <- c(lista, lev2(x[[i]], index))
      index <- lista[length(lista)] + 1
    } else if (depth == 3){
      lista <- c(lista, lev3(x[[i]], index))
      index <- lista[length(lista)] + 1
    } else if (depth == 1){
      lista <- c(lista, rep(index, x[[i]]))
      index <- lista[length(lista)] + 1
    } else{
      cat('There is a not valid CPT \n')
    }
  }
  lista
}

#' Repair to keep the sum of row vectors equal to 1
#' over a mutant vector
#' @param x is the mutant vector
#' @param s is the netwotk structure given by \code{strRep}
#' @keywords internal
keepSum <- function(x, s){
  for (i in 1:s[length(s)]){
    x[which(s == i)] <- x[which(s == i)] / sum(x[which(s == i)])
  }
  x
}

#' Fast version of keepSum using precomputed group indices
#'
#' @param x the mutant vector
#' @param group_indices list of precomputed index vectors for each group
#' @keywords internal
keepSumFast <- function(x, group_indices){
  for (idx in group_indices){
    x[idx] <- x[idx] / sum(x[idx])
  }
  x
}



#' Parameter adaptation as in JADE
#' Normal distribution
#'
#' @keywords internal
randN <- function(Q, mCR){
  CR <- rnorm(Q, mCR, 0.1)
  if (CR < 0 ){
    CR <- 0
  } else if( CR > 1){
    CR <- 1
  }
  CR
}

#' Parameter adaptation as in JADE
#' Cauchy distribution
#'
#' @keywords internal
randC <- function(Q, mF){
  F <- rcauchy(Q, mF, 0.05)
  if (F > 0 & F <= 1){
    F
  } else if(F > 1){
    F <- 1
    F
  } else if(F <= 0){
    randC(Q, mF)
  }
}

#' Parameter adaptation as in JADE
#' Lehmer mean
#'
#' @keywords internal
meanL <- function(sF){
  sum(sF^2) / sum(sF)
}

#' Parameter adaptation as in JADE
#'
#' @keywords internal
meanWL <- function(S, diff){
  wk <- diff / sum(diff)
  sum(wk * S^2) /sum(wk * S)
}


#' JADE uses a mutation strategy called ”DE/current-to- p best”
#' Vectors are selected from the best p vectors in the current
#' population
#'
#' @keywords internal
getPBest <- function(x, n = 30) {
  which(x >= -sort(-x, partial = n)[n])
}

#' Apply a custom edgelist to a bnclassify BN structure
#'
#' @param bn A bnc_dag object from bnclassify
#' @param edgelist A matrix with 2 columns (from, to) representing directed edges
#' @param class.name The name of the class variable
#'
#' @keywords internal
apply_edgelist <- function(bn, edgelist, class.name) {
  if (!is.matrix(edgelist) && !is.data.frame(edgelist)) {
    stop("'edgelist' must be a matrix or data.frame with 2 columns.")
  }
  if (ncol(edgelist) != 2) {
    stop("'edgelist' must have exactly 2 columns (from, to).")
  }

  df <- as.data.frame(edgelist, stringsAsFactors = FALSE)
  colnames(df) <- c("from", "to")
  all_nodes <- unique(c(df$from, df$to))
  if (!class.name %in% all_nodes) {
    stop("'class.name' must appear in the edgelist nodes.")
  }

  # Replace edges in the DAG using named access
  bn$.dag$edges <- edgelist

  # Set call_struct to tan_cl so lp() processes it correctly
  bn$.call_struct[[1]] <- "tan_cl"

  # Rebuild families from the edgelist
  unique_nodes <- unique(c(df$from, df$to))
  unique_nodes <- unique_nodes[order(unique_nodes != class.name, decreasing = TRUE)]
  families <- lapply(unique_nodes, function(node) get_family(node, df, class.name))
  names(families) <- unique_nodes
  bn$.families <- families

  bn
}


#' Precompute data indices for fast CLL evaluation
#'
#' @param data The dataset
#' @param Z A bnc_bn object
#'
#' @keywords internal
precompute_cll_data <- function(data, Z) {
  params <- Z$.params
  class_var <- Z$.class
  node_names <- names(params)
  n <- nrow(data)
  class_levels <- dimnames(params[[class_var]])[[1]]
  n_classes <- length(class_levels)
  observed_class_idx <- match(as.character(data[[class_var]]), class_levels)

  feature_nodes <- setdiff(node_names, class_var)
  feature_info <- vector("list", length(feature_nodes))
  names(feature_info) <- feature_nodes

  for (nm in feature_nodes) {
    cpt <- params[[nm]]
    dn <- dimnames(cpt)
    var_names <- names(dn)
    dims <- dim(cpt)
    n_dims <- length(dims)
    class_dim <- which(var_names == class_var)
    non_class_dims <- setdiff(seq_len(n_dims), class_dim)

    fixed_idx <- matrix(0L, nrow = n, ncol = n_dims)
    for (d in non_class_dims) {
      fixed_idx[, d] <- match(as.character(data[[var_names[d]]]), dn[[d]])
    }

    strides <- cumprod(c(1L, dims[-n_dims]))

    partial <- rep(1L, n)
    for (d in non_class_dims) {
      partial <- partial + (fixed_idx[, d] - 1L) * strides[d]
    }

    feature_info[[nm]] <- list(
      partial = partial,
      class_stride = strides[class_dim]
    )
  }

  list(
    feature_info = feature_info,
    observed_class_idx = observed_class_idx,
    n = n,
    n_classes = n_classes,
    feature_nodes = feature_nodes,
    class_var = class_var
  )
}


#' Fast Conditional Log-Likelihood without DAG verification
#'
#' @param net A bnc_bn object
#' @param cll_data Precomputed data from precompute_cll_data
#'
#' @keywords internal
fastCLL <- function(net, cll_data) {
  params <- net$.params
  n <- cll_data$n
  n_classes <- cll_data$n_classes

  log_prior <- log(as.numeric(params[[cll_data$class_var]]))
  log_joint <- matrix(rep(log_prior, each = n), nrow = n, ncol = n_classes)

  for (nm in cll_data$feature_nodes) {
    cpt_vec <- as.numeric(params[[nm]])
    info <- cll_data$feature_info[[nm]]
    partial <- info$partial
    cs <- info$class_stride

    for (ci in seq_len(n_classes)) {
      log_joint[, ci] <- log_joint[, ci] + log(cpt_vec[partial + (ci - 1L) * cs])
    }
  }

  log_norm <- log_joint - matrixStats::rowLogSumExps(log_joint)
  sum(log_norm[cbind(seq_len(n), cll_data$observed_class_idx)])
}


#' Compute predictive accuracy.
#'
#' @param x A vector of predicted labels.
#' @param y A vector of true labels.
#' @export
#'
#' @examples
#' data(car)
#' dpl.lshade <- lshade(NP = 20, G = 25, data = car, class.name = names(car)[7], c = 0.1,
#' structure = "tan", pB = 0.05, edgelist = NULL, verbose = 5)
#' p <- predict(dpl.lshade$Best, car)
#' accuracy(p, car$class)
accuracy <- function(x, y) {
  if (length(x) != length(y)) {
    stop("The 'predicted' and 'actual' vectors must have the same length.")
  }

  # Check if the inputs are character or factor vectors
  if (!(is.character(x) || is.factor(x)) ||
      !(is.character(y) || is.factor(y))) {
    stop("Both arguments, 'predicted' and 'actual', must be character or factor vectors.")
  }

  # Convert character inputs to factor if necessary
  if (is.character(x)) {
    x <- factor(x)
  }

  if (is.character(y)) {
    y <- factor(y)
  }

  # Calculate the total number of observations
  total_observations <- length(x)

  # Calculate the number of correctly classified observations
  correct_predictions <- sum(x == y)

  # Calculate and return the accuracy
  accuracy <- correct_predictions / total_observations
  return(accuracy)
}




