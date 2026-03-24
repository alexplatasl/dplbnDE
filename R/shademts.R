#' Discriminative parameter learning of BN by SHADE-MTS.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by SHADE hybridized with MTS-LS1 (Multiple Trajectory Search Local Search 1).
#' SHADE handles global exploration while MTS-LS1 performs coordinate-wise
#' local search on the best individual, exploiting the near-separable
#' structure of BN parameter spaces.
#'
#' @name shademts
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param c Numeric. An adaptation parameter. Default is 0.1.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param pB Numeric. Percentage of best individuals for mutation strategy. Default is 0.05.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure.
#' @param ls_freq Positive integer. Apply MTS-LS1 every ls_freq generations. Default is 3.
#' @param ls_max_dims Positive integer. Maximum dimensions to perturb per LS call. Default is 200.
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return An object of class \code{DE}.
#'
#' @references Tseng, L.Y., Chen, C. (2008). Multiple Trajectory Search for
#' Large Scale Global Optimization. IEEE CEC 2008.
#'
#' @examples
#' data(car)
#' dpl.smts <- shademts(NP = 20, G = 25, data = car, class.name = names(car)[7],
#' c = 0.1, pB = 0.05, ls_freq = 3, verbose = 5)
#' print(dpl.smts)
#' \dontrun{plot(dpl.smts)}


shademts <- function(NP = 40, G = 100, data, class.name, c = 0.1,
                     structure = c("nb", "tancl", "tan", "hc"),
                     pB = 0.05, edgelist = NULL, ls_freq = 3,
                     ls_max_dims = 200, verbose = 25, ...) {

  if (NP <= 5) { warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE); NP <- 40 }
  if (G <= 1) { warning("'G' <= 1; set to default value 100\n", immediate. = TRUE); G <- 100 }
  if (c < 0 || c > 0.2) { warning("'c' not in [0, 0.2]; set to default value 0.1\n", immediate. = TRUE); c <- 0.1 }
  if (pB <= 0 || pB > 1) { warning("'pB' not in (0, 1]; set to default value 0.05\n", immediate. = TRUE); pB <- 0.05 }

  neval <- 0; record_CLL <- c(); record_evals <- c()

  # Structure learning
  if (is.null(structure)) structure <- "nb"
  structure <- match.arg(structure)
  if (structure == "tan") structure <- "tancl"
  if (structure == "tancl") { bn <- bnclassify::tan_cl(class.name, data, ...)
  } else if (structure == "hc") { bn <- bnclassify::tan_hc(class.name, data, ...)
  } else { bn <- bnclassify::nb(class.name, data) }
  if (!is.null(edgelist)) bn <- apply_edgelist(bn, edgelist, class.name)

  Z <- bnclassify::lp(bn, data, smooth = 0.01)
  W <- length(Z$.params); X <- lapply(Z$.params, dim)
  Y <- sapply(Z$.params, length); dim <- sum(unlist(Y))
  COL <- strRep(X)
  group_indices <- lapply(seq_len(max(COL)), function(i) which(COL == i))
  id_params <- grep(".params", names(Z))
  offsets <- precompute_offsets(Y)
  cll_data <- precompute_cll_data(data, Z)

  # SHADE setup (no LPSR, fixed NP)
  H <- 6; ki <- 1; terminal <- 0
  mCR <- rep(0.5, H); mF <- rep(0.5, H)
  pBest_n <- max(1, round(NP * pB))

  pop <- population(NP, W, X, Y, Z)
  CLL_fn <- function(net) fastCLL(net, cll_data)
  fitness <- unlist(lapply(pop, CLL_fn))
  neval <- NP
  record_evals <- c(neval)
  best_idx <- which.max(fitness); best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])
  Archive <- list()

  # MTS-LS1 step sizes (per dimension)
  SR <- rep(0.4, dim)

  for (i in seq_len(G)) {
    # ============================
    # SHADE generation
    # ============================
    SF <- c(); SCR <- c(); improvement <- c()
    for (j in seq_len(NP)) {
      ri <- sample.int(H, 1)
      CR <- if (mCR[ri] == terminal) 0 else randN(1, mCR[ri])
      F <- randC(1, mF[ri])

      idxs <- seq_len(NP)[-j]
      idrs <- seq_len(NP + length(Archive))[-j]
      idbs <- getPBest(fitness, pBest_n)

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
        if (f > fitness[best_idx]) { best_idx <- j; best <- trial }
      }
    }

    if (length(Archive) > NP)
      Archive <- Archive[sample(seq_len(length(Archive)), NP)]

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
    # MTS-LS1 on best individual
    # ============================
    if (i %% ls_freq == 0) {
      ls_budget <- min(dim, ls_max_dims) * 2  # ~2 evals per dim attempted
      best_vec <- vec(best$.params)
      mts_result <- mtsLS1(best_vec, fitness[best_idx], SR,
                           ls_max_dims, ls_budget,
                           group_indices, Z, id_params, offsets, cll_data)
      neval <- neval + mts_result$evals
      SR <- mts_result$SR

      if (mts_result$f_vec > fitness[best_idx]) {
        new_net <- vec2netFast(mts_result$vec, Z, id_params, offsets)
        pop[[best_idx]] <- new_net
        fitness[best_idx] <- mts_result$f_vec
        best <- new_net
      }
    }

    record_CLL <- c(record_CLL, fitness[best_idx])
    record_evals <- c(record_evals, neval)

    if (verbose > 0 && (i %% verbose) == 0)
      cat("Gen: ", i, "\t CLL= ", fitness[best_idx], "\t NP= ", NP, "\n")
  }

  outDE <- list(Best = best, BestCLL = fitness[best_idx],
                pobFinal = pop, CLLPobFinal = fitness,
                N.evals = neval, convergence = record_CLL,
                evaluations = record_evals)
  attr(outDE, "class") <- "DE"
  outDE
}
