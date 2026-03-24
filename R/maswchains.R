#' Discriminative parameter learning of BN by MA-SW-Chains.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by a Memetic Algorithm with Solis-Wets Chains. Each individual maintains
#' its own persistent Solis-Wets state. After each DE generation, the most
#' promising individuals receive longer local search chains.
#'
#' @name maswchains
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations that may be performed before the algorithm is halted.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param F A numeric. Mutation factor for DE/rand/1. Default is 0.5.
#' @param CR A numeric. Cross over factor. Default is 0.9.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure
#' that will replace de learned structure.
#' @param ls_evals_ratio Numeric. Fraction of NP used as total LS budget per generation. Default is 0.5.
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
#' @references Molina, D., Lozano, M., Herrera, F. (2010). MA-SW-Chains: Memetic Algorithm
#' based on Local Search Chains. Evolutionary Computation 18(1).
#'
#' @examples
#' data(car)
#' dpl.masw <- maswchains(NP = 20, G = 25, data = car, class.name = names(car)[7],
#' verbose = 5)
#' print(dpl.masw)
#' \dontrun{plot(dpl.masw)}


maswchains <- function(NP = 40, G = 100, data, class.name, F = 0.5, CR = 0.9,
                       structure = c("nb", "tancl", "tan", "hc"),
                       edgelist = NULL, ls_evals_ratio = 0.5,
                       verbose = 25, ...) {

  if (NP <= 5) {
    warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE)
    NP <- 40
  }
  if (G <= 1) {
    warning("'G' <= 1; set to default value 100\n", immediate. = TRUE)
    G <- 100
  }
  if (F < 0 || F > 2) {
    warning("'F' not in [0, 2]; set to default value 0.5\n", immediate. = TRUE)
    F <- 0.5
  }
  if (CR < 0 || CR > 1) {
    warning("'CR' not in [0, 1]; set to default value 0.9\n", immediate. = TRUE)
    CR <- 0.9
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

  # Initialize population
  pop <- population(NP, W, X, Y, Z)
  CLL_fn <- function(net) fastCLL(net, cll_data)
  fitness <- unlist(lapply(pop, CLL_fn))
  neval <- neval + NP
  record_evals <- c(neval)
  best_idx <- which.max(fitness)
  best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])

  # Solis-Wets state per individual
  sw_states <- lapply(seq_len(NP), function(k) {
    list(rho = 0.5, bias = rep(0, dim), n_success = 0L, n_fail = 0L,
         improvement = 0, n_applications = 0L)
  })

  for (i in seq_len(G)) {
    # ============================
    # PHASE 1: DE/rand/1/bin
    # ============================
    for (j in seq_len(NP)) {
      idxs <- seq_len(NP)[-j]
      R <- pop[sample(idxs, 3)]
      r0 <- R[[1]]; r1 <- R[[2]]; r2 <- R[[3]]
      mutant <- vec(r0$.params) + F * (vec(r1$.params) - vec(r2$.params))
      mutant <- reflect(mutant)

      cross_points <- runif(dim) > CR
      if (!all(cross_points)) cross_points[sample(1:dim, 1)] <- TRUE
      trial <- vec(pop[[j]]$.params)
      trial[cross_points] <- mutant[cross_points]

      trial <- keepSumFast(trial, group_indices)
      trial_net <- vec2netFast(trial, Z, id_params, offsets)
      f <- CLL_fn(trial_net)
      neval <- neval + 1

      if (f > fitness[j]) {
        # Reset SW state if individual changed substantially
        old_vec <- vec(pop[[j]]$.params)
        if (sqrt(sum((trial - old_vec)^2)) > 0.1 * sqrt(dim)) {
          sw_states[[j]] <- list(rho = 0.5, bias = rep(0, dim),
                                 n_success = 0L, n_fail = 0L,
                                 improvement = sw_states[[j]]$improvement,
                                 n_applications = sw_states[[j]]$n_applications)
        }
        fitness[j] <- f
        pop[[j]] <- trial_net
        if (f > fitness[best_idx]) {
          best_idx <- j
          best <- trial_net
        }
      }
    }

    # ============================
    # PHASE 2: SW-Chains
    # ============================
    ls_budget_total <- max(10, round(ls_evals_ratio * NP))

    # Compute intensity per individual
    ranks <- rank(-fitness)  # 1 = best
    rank_score <- (NP - ranks) / NP
    hist_scores <- sapply(sw_states, function(s) s$improvement)
    max_hist <- max(abs(hist_scores), 1e-30)
    hist_norm <- hist_scores / max_hist
    intensity <- 0.6 * rank_score + 0.4 * pmax(hist_norm, 0)
    intensity <- pmax(intensity, 1e-10)
    prob_ls <- intensity / sum(intensity)

    ls_remaining <- ls_budget_total
    while (ls_remaining > 0) {
      # Roulette selection
      sel <- sample.int(NP, 1, prob = prob_ls)

      chain_budget <- min(max(5, round(ls_remaining * prob_ls[sel])), ls_remaining)

      sel_vec <- vec(pop[[sel]]$.params)
      sw <- sw_states[[sel]]
      f_before <- fitness[sel]

      sw_result <- solisWetsSearch(sel_vec, f_before, sw$rho, sw$bias,
                                   sw$n_success, sw$n_fail, chain_budget,
                                   group_indices, Z, id_params, offsets, cll_data)

      neval <- neval + sw_result$evals
      ls_remaining <- ls_remaining - sw_result$evals

      mejora <- sw_result$f_vec - f_before
      sw_states[[sel]]$rho <- sw_result$rho
      sw_states[[sel]]$bias <- sw_result$bias
      sw_states[[sel]]$n_success <- sw_result$n_success
      sw_states[[sel]]$n_fail <- sw_result$n_fail
      sw_states[[sel]]$improvement <- sw_states[[sel]]$improvement + mejora
      sw_states[[sel]]$n_applications <- sw_states[[sel]]$n_applications + 1L

      if (sw_result$f_vec > fitness[sel]) {
        new_net <- vec2netFast(sw_result$vec, Z, id_params, offsets)
        pop[[sel]] <- new_net
        fitness[sel] <- sw_result$f_vec
        if (sw_result$f_vec > fitness[best_idx]) {
          best_idx <- sel
          best <- new_net
        }
      }

      # Reset if rho too small
      if (sw_states[[sel]]$rho < 1e-10) {
        sw_states[[sel]]$rho <- 0.5
        sw_states[[sel]]$bias <- rep(0, dim)
        sw_states[[sel]]$n_success <- 0L
        sw_states[[sel]]$n_fail <- 0L
      }
    }

    record_CLL <- c(record_CLL, fitness[best_idx])
    record_evals <- c(record_evals, neval)

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
