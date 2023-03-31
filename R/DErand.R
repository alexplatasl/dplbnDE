#' Discriminative parameter learning of BN by DE variant rand/k/.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by Differential Evolution with variant rand/k
#'
#' @name DErand
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations that may be performed before the algorithm is halted.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param F A numeric. Mutation factor. Default is 0.5.
#' @param CR A numeric. Cross over factor. Default is 0.7.
#' @param mutation.pairs A positive integer giving the number of pairs (1 or 2) used in the mutation step.
#' @param crossover A character. Crossover type among binomial (bin) or exponential (exp).
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure.
#' that will replace de learned structure.
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return An object of class \code{DE}, which is a list with the following components:
#' \item{Best}{A \code{bnc_bn} object with the best individual in the final population, i.e., the bayesian network with the best fitness at the end of evolution.}
#' \item{BestCLL}{A numeric specifying the Conditional Log-Likelihood of the best individual.}
#' \item{pobFinal}{A list of \code{bnc_bn} objects with the final population, i.e., a set of bayesian networks with optimized parameters at the end of evolution.}
#' \item{CLLPobFinal}{A numeric vector specifying the Conditional Log-Likelihood of the final population.}
#' \item{N.evals}{An integer giving the total number of evaluations.}
#' \item{convergence}{A numeric vector giving the maximum Conditional Log-Likelihood at each generation.}
#' \item{evaluations}{An integer vector giving the total number of evaluations at each generation.}
#'
#' @examples
#' # Load data
#' data(car)
#' # Parameter learning with "rand/2/bin" variant
#' dpl.rand2bin <- DErand(NP=25, G=40, data = car, class.name = names(car)[7], F = 0.5,
#' CR = 0.5, mutation.pairs = 2, crossover = "bin", structure = "tan",edgelist = NULL,
#' verbose = 10)
#' # Print results
#' print(dpl.rand2bin)
#' \dontrun{plot(dpl.rand2bin)}



DErand <- function(NP=40, G=100, data, class.name, F = 0.5, CR = 0.7, mutation.pairs = c(1,2),
                 crossover = c("bin","exp"), structure = c("nb","tancl","hc"),
                 edgelist=NULL, verbose=25, ...){

  if (NP <=5){
    warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE)
    NP <- 40
  }

  if (G <=1){
    warning("'G' <= 1; set to default value 100\n", immediate. = TRUE)
    G <- 100
  }

  if (F < 0 | F > 2) {
    warning("'F' not in [0,2]; set to default value 0.5\n", immediate. = TRUE)
    F <- 0.5
  }

  if (CR < 0 | CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.7\n", immediate. = TRUE)
    CR <- 0.7
  }

  if (mutation.pairs < 1 | mutation.pairs > 2) {
    warning("'mutation.pairs' not in {1,2}; set to default value 1\n",
            immediate. = TRUE)
    mutation.pairs <- 1
  }

  neval = 0
  record_CLL = c()
  record_evals = c()

  # Algorithm to learn structure
  if(!is.null(edgelist)){
    structure <- "nb"
  }
  structure <- match.arg(structure)
  if (structure == "tancl"){
    bn = bnclassify::tan_cl(class.name, data, ...)
  }else if(structure == "hc"){
    bn = bnclassify::tan_hc(class.name, data, ...)
  } else {
    bn = bnclassify::nb(class.name, data)
  }


  # To replace BN structure from adjacency list (if provided)
  if (!is.null(edgelist)){
    bn[[2]][[2]] = edgelist
  }

  # Start CPTs
  Z = bnclassify::lp(bn, data, smooth = 0.01)
  W = length(Z$.params)
  X = lapply(Z$.params, dim)
  Y = sapply(Z$.params, length)
  dim = sum(unlist(Y))
  COL =  strRep(X)

  # Differential Evolution
  # Starts population
  pop = population(NP, W, X,  Y, Z)
  # Fitness evaluation
  CLL = function(net) bnclassify::cLogLik(net, data)
  fitness = unlist(lapply(pop, CLL))
  neval = neval + NP
  record_evals = c(neval)
  best_idx = which.max(fitness)
  best = pop[[best_idx]]
  record_CLL = c(record_CLL, fitness[best_idx])

  for(i in seq_len(G)){
    for(j in seq_len(NP)){

      # Number of pairs used in mutation strategy
      if (mutation.pairs == 1){
        ### best/1
        idxs = seq_len(NP)[-j]
        R = pop[sample(idxs,3)]
        r0 = R[[1]]; r1 = R[[2]]; r2 = R[[3]]

        # mutant vector
        mutant = vec(r0$.params) + F * (vec(r1$.params) - vec(r2$.params))
      }else{
        ### best/2
        idxs = seq_len(NP)[-j]
        R = pop[sample(idxs,5)]
        r0 = R[[1]]; r1 = R[[2]]; r2 = R[[3]]; r3 = R[[4]]; r4 = R[[5]]

        # mutant vector
        mutant = vec(r0$.params) + F * (vec(r1$.params) - vec(r2$.params) + vec(r3$.params) - vec(r4$.params))
      }

      # Repair 1
      mutant = reflect(mutant)


      # Crossover type
      crossover <- match.arg(crossover)
      if (crossover == "bin"){
        ### /bin
        cross_points = runif(dim) > CR
        if(!all(cross_points)){
          cross_points[sample(1:dim,1)] = TRUE
        }
        trial = vec(pop[[j]]$.params);
        trial[cross_points] = mutant[cross_points] # Cross Over
      }else {
        ### /exp
        trial = vec(pop[[j]]$.params);
        k = sample(1:dim,1)
        trial[k] = mutant[k] # Cross Over
        L = 1
        beta = runif(1, 0, 1)
        while (beta <= CR && L < (dim - 1)) {
          L = L + 1
          beta = runif(1, 0, 1)
          k = 1 + (k %% dim)
          trial[k] = mutant[k] # Cross Over
        }
      }

      # Repair 2
      trial = keepSum(trial, COL)
      trial = vec2net(trial, W, X, Y, Z)
      f = CLL(trial)
      neval = neval + 1
      if(f > fitness[j]){
        fitness[j] = f
        pop[[j]] = trial
        if(f > fitness[best_idx]){
          best_idx = j
          best = trial
        }  # end if
      } # end if
    } # end j loop (NP)
    record_CLL = c(record_CLL, fitness[best_idx])
    record_evals = c(record_evals, neval)

    # Print best fitness each verbose generations
    if (verbose > 0 & (i %% verbose) == 0){
      cat("Gen: ", i,"\t CLL= ",fitness[best_idx], "\t NP= ", NP, "\n")
    }
  } # end i loop (Generations)

  outDE = list(Best=best,
               BestCLL=fitness[best_idx],
               pobFinal= pop,
               CLLPobFinal = fitness,
               N.evals = neval,
               convergence = record_CLL,
               evaluations = record_evals)
  attr(outDE, "class") <- "DE"
  outDE
}
