#' Discriminative parameter learning of BN by SHADEILS.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by SHADE with Iterative Local Search. Combines SHADE (adaptive DE)
#' as global explorer with Solis-Wets local search for intensification.
#' Includes partial restart mechanism when stagnation is detected.
#'
#' @name shadeils
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations that may be performed before the algorithm is halted.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param c Numeric. An adaptation parameter for SHADE memory. Default is 0.1.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param pB Numeric. Percentage of best individuals for mutation strategy. Default is 0.05.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure
#' that will replace de learned structure.
#' @param ls_budget_ratio Numeric. Fraction of NP used as budget for each local search call. Default is 0.3.
#' @param ls_freq Positive integer. Apply local search every ls_freq generations. Default is 5.
#' @param threshold_improvement Numeric. Minimum relative improvement before stagnation counter increments. Default is 0.01.
#' @param max_iters_no_improve Positive integer. Consecutive stagnation cycles before partial restart. Default is 3.
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return An object of class \code{DE}, which is a list with the following components:
#' \item{Best}{A \code{bnc_bn} object with the best individual in the final population.}
#' \item{BestCLL}{A numeric specifying the Conditional Log-Likelihood of the best individual.}
#' \item{pobFinal}{A list of \code{bnc_bn} objects with the final population.}
#' \item{CLLPobFinal}{A numeric vector specifying the Conditional Log-Likelihood of the final population.}
#' \item{N.evals}{An integer giving the total number of evaluations.}
#' \item{convergence}{A numeric vector giving the maximum Conditional Log-Likelihood at each generation.}
#' \item{evaluations}{An integer vector giving the total number of evaluations at each generation.}
#'
#' @references Molina, D., LaTorre, A., Herrera, F. (2018). SHADE with Iterative Local Search
#' for Large-Scale Global Optimization. IEEE CEC 2018.
#'
#' @examples
#' data(car)
#' dpl.shadeils <- shadeils(NP = 20, G = 25, data = car, class.name = names(car)[7],
#' c = 0.1, pB = 0.05, ls_freq = 5, verbose = 5)
#' print(dpl.shadeils)
#' \dontrun{plot(dpl.shadeils)}


shadeils <- function(NP = 40, G = 100, data, class.name, c = 0.1,
                     structure = c("nb", "tancl", "tan", "hc"),
                     pB = 0.05, edgelist = NULL, ls_budget_ratio = 0.3,
                     ls_freq = 5, threshold_improvement = 0.01,
                     max_iters_no_improve = 3, verbose = 25, ...) {

  if (NP <= 5) {
    warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE)
    NP <- 40
  }
  if (G <= 1) {
    warning("'G' <= 1; set to default value 100\n", immediate. = TRUE)
    G <- 100
  }
  if (c < 0 || c > 0.2) {
    warning("'c' not in [0, 0.2]; set to default value 0.1\n", immediate. = TRUE)
    c <- 0.1
  }
  if (pB <= 0 || pB > 1) {
    warning("'pB' not in (0, 1]; set to default value 0.05\n", immediate. = TRUE)
    pB <- 0.05
  }

  neval <- 0
  record_CLL <- c()
  record_evals <- c()

  # Structure learning
  if (is.null(structure)) structure <- "nb"
  structure <- match.arg(structure)
  if (structure == "tan") structure <- "tancl"
  if (structure == "tancl") {
    bn <- bnclassify::tan_cl(class.name, data, ...)
  } else if (structure == "hc") {
    bn <- bnclassify::tan_hc(class.name, data, ...)
  } else {
    bn <- bnclassify::nb(class.name, data)
  }

  if (!is.null(edgelist)) {
    bn <- apply_edgelist(bn, edgelist, class.name)
  }

  # CPTs and precomputed structures
  Z <- bnclassify::lp(bn, data, smooth = 0.01)
  W <- length(Z$.params)
  X <- lapply(Z$.params, dim)
  Y <- sapply(Z$.params, length)
  dim <- sum(unlist(Y))
  COL <- strRep(X)
  group_indices <- lapply(seq_len(max(COL)), function(i) which(COL == i))
  id_params <- grep(".params", names(Z))
  offsets <- precompute_offsets(Y)
  cll_data <- precompute_cll_data(data, Z)

  # SHADE setup
  H <- 6
  ki <- 1
  terminal <- 0
  pop <- population(NP, W, X, Y, Z)
  CLL_fn <- function(net) fastCLL(net, cll_data)
  fitness <- unlist(lapply(pop, CLL_fn))
  neval <- neval + NP
  record_evals <- c(neval)
  best_idx <- which.max(fitness)
  best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])
  mCR <- rep(0.5, H); mF <- mCR
  Archive <- list()
  pBest <- max(1, round(NP * pB))

  # Solis-Wets state for best individual
  sw_rho <- 0.5
  sw_bias <- rep(0, dim)
  sw_ns <- 0L; sw_nf <- 0L

  # Stagnation tracking
  iters_no_improve <- 0

  for (i in seq_len(G)) {
    best_before <- fitness[best_idx]

    # ============================
    # PHASE 1: SHADE generation
    # ============================
    SF <- c(); SCR <- c(); improvement <- c()
    for (j in seq_len(NP)) {
      ri <- sample.int(H, 1)
      CR <- if (mCR[ri] == terminal) 0 else randN(1, mCR[ri])
      F <- randC(1, mF[ri])

      idxs <- seq_len(NP)[-j]
      idrs <- seq_len(NP + length(Archive))[-j]
      idbs <- getPBest(fitness, pBest)

      xi <- pop[[j]]
      xp <- pop[[sample(idbs, 1)]]
      r1 <- pop[[sample(idxs, 1)]]
      r2 <- c(pop, Archive)[[sample(idrs, 1)]]

      mutant <- vec(xi$.params) + F * (vec(xp$.params) - vec(xi$.params)) +
        F * (vec(r1$.params) - vec(r2$.params))
      mutant <- reflect(mutant)

      cross_points <- runif(dim) > CR
      if (!all(cross_points)) cross_points[sample(1:dim, 1)] <- TRUE
      trial <- vec(pop[[j]]$.params)
      trial[cross_points] <- mutant[cross_points]

      trial <- keepSumFast(trial, group_indices)
      trial <- vec2netFast(trial, Z, id_params, offsets)
      f <- CLL_fn(trial)
      neval <- neval + 1

      if (f > fitness[j]) {
        improvement <- c(improvement, abs(f - fitness[j]))
        fitness[j] <- f
        Archive <- c(Archive, pop[j])
        pop[[j]] <- trial
        SCR <- c(SCR, CR); SF <- c(SF, F)
        if (f > fitness[best_idx]) {
          best_idx <- j
          best <- trial
        }
      }
    }

    # Archive pruning
    if (length(Archive) > NP) {
      Archive <- Archive[sample(seq_len(length(Archive)), NP)]
    }

    # Update SHADE memory
    if (length(SCR) > 0 && length(SF) > 0) {
      if (mCR[ki] == terminal || max(SCR) == 0) {
        mCR[ki] <- terminal
      } else {
        mCR[ki] <- meanWL(SCR, improvement)
      }
      mF[ki] <- meanWL(SF, improvement)
      ki <- (ki %% H) + 1
    }

    # ============================
    # PHASE 2: Solis-Wets LS
    # ============================
    if (i %% ls_freq == 0) {
      ls_budget <- max(10, round(ls_budget_ratio * NP))
      best_vec <- vec(best$.params)
      sw_result <- solisWetsSearch(best_vec, fitness[best_idx], sw_rho, sw_bias,
                                   sw_ns, sw_nf, ls_budget, group_indices,
                                   Z, id_params, offsets, cll_data)
      neval <- neval + sw_result$evals
      sw_rho <- sw_result$rho
      sw_bias <- sw_result$bias
      sw_ns <- sw_result$n_success
      sw_nf <- sw_result$n_fail

      if (sw_result$f_vec > fitness[best_idx]) {
        new_net <- vec2netFast(sw_result$vec, Z, id_params, offsets)
        pop[[best_idx]] <- new_net
        fitness[best_idx] <- sw_result$f_vec
        best <- new_net
      }
    }

    record_CLL <- c(record_CLL, fitness[best_idx])
    record_evals <- c(record_evals, neval)

    # ============================
    # PHASE 3: Stagnation restart
    # ============================
    if (best_before != 0) {
      rel_improve <- (fitness[best_idx] - best_before) / abs(best_before)
    } else {
      rel_improve <- fitness[best_idx] - best_before
    }

    if (rel_improve < threshold_improvement) {
      iters_no_improve <- iters_no_improve + 1
    } else {
      iters_no_improve <- 0
    }

    if (iters_no_improve >= max_iters_no_improve) {
      saved_best <- best
      saved_fitness <- fitness[best_idx]
      for (k in seq_len(NP)) {
        if (k != best_idx) {
          pop[[k]] <- individual(W, X, Y, Z)
          fitness[k] <- CLL_fn(pop[[k]])
          neval <- neval + 1
        }
      }
      mCR <- rep(0.5, H); mF <- rep(0.5, H)
      sw_rho <- 0.5; sw_bias <- rep(0, dim); sw_ns <- 0L; sw_nf <- 0L
      iters_no_improve <- 0
    }

    if (verbose > 0 && (i %% verbose) == 0) {
      cat("Gen: ", i, "\t CLL= ", fitness[best_idx], "\t NP= ", NP, "\n")
    }
  }

  outDE <- list(Best = best, BestCLL = fitness[best_idx],
                pobFinal = pop, CLLPobFinal = fitness,
                N.evals = neval, convergence = record_CLL,
                evaluations = record_evals)
  attr(outDE, "class") <- "DE"
  outDE
}
