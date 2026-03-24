#' Discriminative parameter learning of BN by NL-SHADE-RSP.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by NL-SHADE-RSP (Non-Linear population size reduction SHADE with
#' Ranked-based Selective Pressure). Uses exponential population reduction
#' and fitness-ranked selection of donor vectors.
#'
#' @name nlshadersp
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
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return An object of class \code{DE}.
#'
#' @references Stanovov, V., Akhmedova, S., Semenkin, E. (2021). NL-SHADE-RSP
#' Algorithm with Adaptive Archive and Selective Pressure for CEC 2021. IEEE CEC 2021.
#'
#' @examples
#' data(car)
#' dpl.nlsr <- nlshadersp(NP = 20, G = 25, data = car, class.name = names(car)[7],
#' c = 0.1, pB = 0.05, verbose = 5)
#' print(dpl.nlsr)
#' \dontrun{plot(dpl.nlsr)}


nlshadersp <- function(NP = 40, G = 100, data, class.name, c = 0.1,
                       structure = c("nb", "tancl", "tan", "hc"),
                       pB = 0.05, edgelist = NULL, verbose = 25, ...) {

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

  # NL-SHADE-RSP parameters
  H <- 20; ki <- 1; terminal <- 0
  mF <- rep(0.5, H); mCR <- rep(0.5, H)
  N.init <- NP; N.min <- 4
  MAX_NFE <- G * NP; NFE <- 0

  pop <- population(NP, W, X, Y, Z)
  CLL_fn <- function(net) fastCLL(net, cll_data)
  fitness <- unlist(lapply(pop, CLL_fn))
  neval <- NP; NFE <- NP
  record_evals <- c(neval)
  best_idx <- which.max(fitness); best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])
  Archive <- list()
  pBest_n <- max(1, round(NP * pB))

  for (i in seq_len(G)) {
    SF <- c(); SCR <- c(); improvement <- c()
    ratio <- NFE / MAX_NFE

    # Ranked selective pressure weights for r1 selection
    ranks <- rank(-fitness)  # 1 = best
    rsp_weights <- 3 * (NP - ranks)
    rsp_weights <- pmax(rsp_weights, 1e-10)
    rsp_probs <- rsp_weights / sum(rsp_weights)

    for (j in seq_len(NP)) {
      ri <- sample.int(H, 1)
      CR <- if (mCR[ri] == terminal) 0 else randN(1, mCR[ri])
      F <- randC(1, mF[ri])

      # Adaptive archive: prob of selecting from archive
      archive_prob <- length(Archive) / (length(Archive) + NP + 1e-30)

      # Select r1 with ranked selective pressure
      r1_candidates <- seq_len(NP)[-j]
      r1_probs_adj <- rsp_probs[r1_candidates]
      r1_probs_adj <- r1_probs_adj / sum(r1_probs_adj)
      r1_idx <- r1_candidates[sample.int(length(r1_candidates), 1, prob = r1_probs_adj)]

      # Select r2 from pop or archive based on archive_prob
      if (runif(1) < archive_prob && length(Archive) > 0) {
        r2 <- Archive[[sample.int(length(Archive), 1)]]
      } else {
        r2_candidates <- seq_len(NP)[-c(j, r1_idx)]
        if (length(r2_candidates) == 0) r2_candidates <- seq_len(NP)[-j]
        r2 <- pop[[sample(r2_candidates, 1)]]
      }

      # p-best selection
      idbs <- getPBest(fitness, pBest_n)
      xp <- pop[[sample(idbs, 1)]]

      # current-to-pbest/1 mutation
      mutant <- vec(pop[[j]]$.params) +
        F * (vec(xp$.params) - vec(pop[[j]]$.params)) +
        F * (vec(pop[[r1_idx]]$.params) - vec(r2$.params))
      mutant <- reflect(mutant)

      # Binomial crossover
      cross_points <- runif(dim) > CR
      if (!all(cross_points)) cross_points[sample(1:dim, 1)] <- TRUE
      trial <- vec(pop[[j]]$.params)
      trial[cross_points] <- mutant[cross_points]

      trial <- keepSumFast(trial, group_indices)
      trial <- vec2netFast(trial, Z, id_params, offsets)
      f <- CLL_fn(trial)
      neval <- neval + 1; NFE <- NFE + 1

      if (f > fitness[j]) {
        improvement <- c(improvement, abs(f - fitness[j]))
        fitness[j] <- f
        Archive <- c(Archive, pop[j])
        pop[[j]] <- trial
        SCR <- c(SCR, CR); SF <- c(SF, F)
        if (f > fitness[best_idx]) { best_idx <- j; best <- trial }
      }
    }

    # Archive pruning
    if (length(Archive) > NP)
      Archive <- Archive[sample(seq_len(length(Archive)), NP)]

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

    # Non-linear population size reduction (NLPSR)
    nl_ratio <- ratio^(1 - ratio)
    New.NP <- round(N.init + (N.min - N.init) * nl_ratio)
    New.NP <- max(New.NP, N.min)
    if (NP > New.NP) {
      ord <- order(fitness, decreasing = TRUE)
      pop <- pop[ord]; fitness <- fitness[ord]
      pop <- pop[1:New.NP]; fitness <- fitness[1:New.NP]
      best_idx <- 1; best <- pop[[1]]
      if (length(Archive) > New.NP)
        Archive <- Archive[sample(seq_len(length(Archive)), New.NP)]
      NP <- New.NP
      pBest_n <- max(1, round(NP * pB))
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
