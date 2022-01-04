#' Discriminative parameter learning of BN by JADE.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by Adaptive Differential Evolution with optional external Archive
#'
#' @name jade
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations that may be performed before the algorithm is halted.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param c A numeric. An adaptation parameter. Default is 0.1.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param pB A numeric. JADE mutation strategy.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure
#' that will replace de learned structure.
#' @param archive A logical. If TRUE, trial vector r2 is randomly selected from the union of
#' the current population and the external archive.
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return A list with the following components.
#' \item{Best}{A \code{bnc_bn} object with the best individual in the final population.}
#' \item{BestCLL}{A numeric specifying the Conditional Log-Likelihood of the best individual.}
#' \item{pobFinal}{A list of \code{bnc_bn} objects with the final population.}
#' \item{CLLPobFinal}{A numeric vector specifying the Conditional Log-Likelihood of the final population.}
#' \item{N.evals}{An integer giving the total number of evaluations.}
#' \item{convergence}{A numeric vector giving the maximum Conditional Log-Likelihood at each generation.}
#' \item{evaluations}{An integer vector giving the total number of evaluations at each generation.}


jade <- function(NP=40, G=100, data, class.name, c = 0.1, structure = c("nb","tancl","hc"),
                   pB=0.05, edgelist=NULL, archive = TRUE, verbose=25, ...){

  if (NP <=5){
    warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE)
    NP <- 40
  }

  if (G <=1){
    warning("'G' <= 1; set to default value 100\n", immediate. = TRUE)
    G <- 100
  }

  if (c < 0 | c > 0.2){
    warning("'c' not in [0,0.2]; set to default value 0.1\n", immediate. = TRUE)
    c <- 0.1
  }

  if (pB <= 0 | pB > 1){
    warning("'pB' not in (0,1]; set to default value 0.05\n", immediate. = TRUE)
    pB <- 0.1
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
  # Default muCR y muF
  mCR = 0.68; mF = 0.5
  # Starts Archive
  Archive = list()
  # p-best  individuals to choose
  pBest = NP * pB

  for(i in seq_len(G)){
    SF = c(); SCR = c()
    for(j in seq_len(NP)){
      CR = randN(1,mCR);
      F = randC(1, mF)
      ### rand/1
      idxs = seq_len(NP)[-j]
      idbs = getPBest(fitness, pBest)
      idrs = seq_len(NP+length(Archive))[-j]

      xi = pop[[j]]
      xp = pop[[sample(idbs,1)]]
      r1 = pop[[sample(idxs,1)]]
      r2 = c(pop, Archive)[[sample(idrs,1)]]
      # mutant vector
      mutant = vec(xi$.params) + F * (vec(xp$.params) - vec(xi$.params)) + F * (vec(r1$.params) - vec(r2$.params))
      # Repair 1
      mutant = reflect(mutant)
      ### /bin
      cross_points = runif(dim) > CR
      if(!all(cross_points)){
        cross_points[sample(1:dim,1)] = TRUE
      }
      trial = vec(pop[[j]]$.params);
      trial[cross_points] = mutant[cross_points] # Cross Over
      # Repair 2
      trial = keepSum(trial, COL)
      trial = vec2net(trial, W, X, Y, Z)
      f = CLL(trial)
      neval = neval + 1
      if(f > fitness[j]){
        fitness[j] = f
        if (archive){
          Archive = c(Archive, pop[j])
        }
        pop[[j]] = trial
        SCR = c(SCR, CR); SF = c(SF, F)
        if(f > fitness[best_idx]){
          best_idx = j
          best = trial
        }  # end if
      } # end if
    } # end j loop (NP)
    record_CLL = c(record_CLL, fitness[best_idx])
    record_evals = c(record_evals, neval)

    # At random removes solutions from A, such that |A| <= NP
    if (length(Archive) > 0){
      Archive = Archive[sample(seq_len(length(Archive)), min(length(Archive),NP))]
    }

    # Updates mCR and mF values
    if (length(SCR) > 0){
      mCR = ((1 - c)* mCR) + (c * mean(SCR))
    }else{
      mCR = 0.9
    }
    if (length(SF)> 0 ){
      mF = ((1 - c)* mF) + (c * meanL(SF))
    }else{
      mF = 0.9
    }

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
