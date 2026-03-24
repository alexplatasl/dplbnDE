#' Discriminative parameter learning of BN by jSO.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by Algorithm jSO, an improved variant of iL-SHADE. Uses a weighted
#' mutation strategy (current-to-pbest-w/1) with phase-based parameter
#' constraints that favor exploration early and exploitation late.
#'
#' @name jso
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure.
#' @param verbose positive integer indicating the number of generations until the iteration progress should be printed.
#' @param ... other structure learning options from \link[bnclassify]{tan_cl} or \link[bnclassify]{tan_hc}.
#'
#' @export
#' @return An object of class \code{DE}.
#'
#' @references Brest, J., Maucec, M.S., Boskovic, B. (2017). Single Objective
#' Real-Parameter Optimization: Algorithm jSO. IEEE CEC 2017.
#'
#' @examples
#' data(car)
#' dpl.jso <- jso(NP = 20, G = 25, data = car, class.name = names(car)[7], verbose = 5)
#' print(dpl.jso)
#' \dontrun{plot(dpl.jso)}


jso <- function(NP = 40, G = 100, data, class.name,
                structure = c("nb", "tancl", "tan", "hc"),
                edgelist = NULL, verbose = 25, ...) {

  if (NP <= 5) { warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE); NP <- 40 }
  if (G <= 1) { warning("'G' <= 1; set to default value 100\n", immediate. = TRUE); G <- 100 }

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

  # jSO parameters
  H <- 5
  mF <- rep(0.3, H); mCR <- rep(0.8, H)
  mF[H] <- 0.9; mCR[H] <- 0.9  # last slot fixed
  ki <- 1
  p_min <- 0.125; p_max <- 0.25
  N.init <- NP; N.min <- 4
  MAX_NFE <- G * NP
  LPSR.1 <- (N.min - N.init) / MAX_NFE
  NFE <- 0

  pop <- population(NP, W, X, Y, Z)
  CLL_fn <- function(net) fastCLL(net, cll_data)
  fitness <- unlist(lapply(pop, CLL_fn))
  neval <- NP; NFE <- NP
  record_evals <- c(neval)
  best_idx <- which.max(fitness); best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])
  Archive <- list()

  for (i in seq_len(G)) {
    SF <- c(); SCR <- c(); improvement <- c()
    ratio <- NFE / MAX_NFE

    # p linearly increasing
    p_cur <- p_min + (p_max - p_min) * ratio
    pBest <- max(1, round(NP * p_cur))

    for (j in seq_len(NP)) {
      ri <- sample.int(H, 1)

      # Generate F with Cauchy(M_F[r], 0.1)
      F <- randC_jso(1, mF[ri])
      # Phase constraint on F
      if (ratio < 0.6 && F > 0.7) F <- 0.7

      # Generate CR with Normal(M_CR[r], 0.1)
      CR <- rnorm(1, mCR[ri], 0.1)
      CR <- max(0, min(1, CR))
      # Phase constraints on CR
      if (ratio < 0.25 && CR < 0.7) CR <- 0.7
      else if (ratio < 0.5 && CR < 0.6) CR <- 0.6

      # Weighted mutation: current-to-pbest-w/1
      if (ratio < 0.2) { Fw <- 0.7 * F
      } else if (ratio < 0.4) { Fw <- 0.8 * F
      } else { Fw <- 1.2 * F }

      idxs <- seq_len(NP)[-j]
      idrs <- seq_len(NP + length(Archive))[-j]
      idbs <- getPBest(fitness, pBest)

      xi <- pop[[j]]
      xp <- pop[[sample(idbs, 1)]]
      r1 <- pop[[sample(idxs, 1)]]
      r2 <- c(pop, Archive)[[sample(idrs, 1)]]

      mutant <- vec(xi$.params) + Fw * (vec(xp$.params) - vec(xi$.params)) +
        F * (vec(r1$.params) - vec(r2$.params))
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

    # Update memory (skip last slot which is fixed)
    if (length(SCR) > 0 && length(SF) > 0 && ki < H) {
      mF[ki] <- meanWL(SF, improvement)
      wk <- improvement / sum(improvement)
      mCR[ki] <- sum(wk * SCR)  # weighted arithmetic mean
      ki <- ki + 1
      if (ki >= H) ki <- 1  # skip H-th slot
    }

    # LPSR
    New.NP <- round(LPSR.1 * NFE + N.init, 0)
    New.NP <- max(New.NP, N.min)
    if (NP > New.NP) {
      ord <- order(fitness, decreasing = TRUE)
      pop <- pop[ord]; fitness <- fitness[ord]
      pop <- pop[1:New.NP]; fitness <- fitness[1:New.NP]
      best_idx <- 1; best <- pop[[1]]
      if (length(Archive) > New.NP)
        Archive <- Archive[sample(seq_len(length(Archive)), New.NP)]
      NP <- New.NP
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
