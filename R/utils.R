#' Generate a valid probability table for attribute Xa
#'
#' @keywords internal
basGen = function(x){
  aux = runif(x[1])
  aux = aux/sum(aux)
  aux
}

#' Generate a valid conditional probability table for
#' attribute Xa given its parent Xb
#' This is a two entry table
#'
#' @keywords internal
parGen = function(x){
  p = list()
  for (i in 1:x[2]){
    aux = runif(x[1])
    aux = aux/sum(aux)
    p[[i]] = aux
  }
  unlist(p)
}


#' Generate a valid conditional probability table for
#' attribute Xa given its parents Xb and Xc
#' This a three entry table
#'
#' @keywords internal
PARGEN = function(x){
  lon = length(x)
  p = list()
  if (lon == 1){
    p = basGen(x)
  }else if(lon == 2){
    p = parGen(x)
  }else{
    p = replicate(x[3], parGen(x), simplify = FALSE)
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
individual = function(W,X,Y,Z){
  ind = as.vector(unlist(lapply(X, PARGEN)))
  unname(unlist(lapply(X, PARGEN)))
  k = 1
  id_params = grep(".params", names(Z))
  for (i in 1:W){
    for (j in 1:Y[[i]]){
      Z[[id_params]][[i]][j] = ind[k]
      k = k + 1
    }
  }
  Z
}

#' Replace a vector of parameters in a CPT
#'
#' @keywords internal
vec2net = function(vec,W,X,Y,Z){
  ind = vec
  k = 1
  id_params = grep(".params", names(Z))
  for (i in 1:W){
    for (j in 1:Y[[i]]){
      Z[[id_params]][[i]][j] = ind[k]
      k = k + 1
    }
  }
  Z
}


#' Generates a population of potential solutions
#' @param pobSize Number of desired individuals at initial population
#'
#' @keywords internal
population = function(pobSize, W, X, Y, Z){
  pob = list()
  pob = replicate(pobSize, individual(W,X,Y,Z), simplify=FALSE)
  pob
}



#' Transform a matrix into a vector
#'
#' @keywords internal
vec = function(mt){unname(unlist(mt))}


#' Reparation to satisfy the constraints for (\eqn{\theta_{a,b}}) in [0 , 1]
#'
#' @keywords internal
reflect = function(vector){
  l = which(vector < 0)
  u = which(vector > 1)
  if (length(c(l,u)) <= 0){
    vector
  }else{
    vector[l] = abs(vector[l]) %% 1
    vector[u] = 1 - vector[u] %% 1
  }
  vector
}


#' Repair to keep the sum of row vectors equal to 1
#' In a two entry way table
#'
#' @keywords internal
lev2 = function(x, index){
  lista = c()
  for (i in 1:x[2]){
    lista = c(lista, rep(index, x[1]))
    index = index + 1
  }
  lista
}

#' Repair to keep the sum of row vectors equal to 1
#' In a three entry way table
#'
#' @keywords internal
lev3 = function(x, index){
  lista = c()
  for (i in 1:x[3]){
    for (j in 1:x[2]){
      lista = c(lista, rep(index, x[1]))
      index = index + 1
    }
  }
  lista
}

#' Repair to keep the sum of row vectors equal to 1
#' Gives BN structure
#'
#' @keywords internal
strRep = function(x){
  index = 1
  lista = c()
  for (i in 1:length(x)){
    depth = length(x[[i]])
    if (depth == 2){
      lista = c(lista, lev2(x[[i]], index))
      index = lista[length(lista)] + 1
    } else if (depth == 3){
      lista = c(lista, lev3(x[[i]], index))
      index = lista[length(lista)] + 1
    } else if (depth == 1){
      lista = c(lista, rep(index, x[[i]]))
      index = lista[length(lista)] + 1
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
keepSum = function(x, s){
  for (i in 1:s[length(s)]){
    x[which(s == i)] = x[which(s == i)] / sum(x[which(s == i)])
  }
  x
}



#' Parameter adaptation as in JADE
#' Normal distribution
#'
#' @keywords internal
randN = function(Q, mCR){
  CR = rnorm(Q, mCR, 0.1)
  if (CR < 0 ){
    CR = 0
  } else if( CR > 1){
    CR = 1
  }
  CR
}

#' Parameter adaptation as in JADE
#' Cauchy distribution
#'
#' @keywords internal
randC = function(Q, mF){
  F = rcauchy(Q, mF, 0.05)
  if (F > 0 & F <= 1){
    F
  } else if(F > 1){
    F = 1
    F
  } else if(F <= 0){
    randC(Q, mF)
  }
}

#' Parameter adaptation as in JADE
#' Lehmer mean
#'
#' @keywords internal
meanL = function(sF){
  sum(sF^2) / sum(sF)
}

#' Parameter adaptation as in JADE
#'
#' @keywords internal
meanWL = function(S, diff){
  wk = diff / sum(diff)
  sum(wk * S^2) /sum(wk * S)
}


#' JADE uses a mutation strategy called ”DE/current-to- p best”
#' Vectors are selected from the best p vectors in the current
#' population
#'
#' @keywords internal
getPBest = function(x, n=30) {
  which(x >= -sort(-x, partial=n)[n])
}

#' Compute predictive accuracy.
#'
#' @param x A vector of predicted labels.
#' @param y A vector of true labels.
#' @export
#'
#' @examples
#' data(car)
#' dpl.lshade <- lshade(NP=20, G=25, data = car, class.name = names(car)[7], c = 0.1,
#' structure = "tan", pB=0.05, edgelist = NULL, verbose = 5)
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
