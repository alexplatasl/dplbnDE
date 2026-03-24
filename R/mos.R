#' Discriminative parameter learning of BN by MOS.
#'
#' Learn parameters of a Bayesian Network in a discriminative way
#' by Multiple Offspring Sampling. Dynamically combines SHADE (adaptive DE)
#' and Solis-Wets local search, adjusting participation of each technique
#' based on recent performance.
#'
#' @name mos
#'
#' @param NP positive integer giving the number of candidate solutions in the initial population.
#' @param G positive integer specifying the maximum number of generations that may be performed before the algorithm is halted.
#' @param data The data frame from which to learn the classifier.
#' @param class.name A character. Name of the class variable.
#' @param structure A character. Name of the structure learning function. "tan" uses Tree Augmented Network.
#'  "nb" uses Naive Bayes. "hc" uses Hill Climbing.
#' @param pB Numeric. Percentage of best individuals for SHADE mutation strategy. Default is 0.05.
#' @param edgelist A matrix. An optional edge list to use a custom BN structure
#' that will replace de learned structure.
#' @param n_splits Positive integer. Number of budget partitions. Default is 10.
#' @param min_participation Numeric. Minimum participation per technique. Default is 0.1.
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
#' @references LaTorre, A., Muelas, S., Pena, J.M. (2013). Large Scale Global Optimization:
#' Experimental Results with MOS-based Hybrid Algorithms. IEEE CEC 2013.
#'
#' @examples
#' data(car)
#' dpl.mos <- mos(NP = 20, G = 25, data = car, class.name = names(car)[7],
#' n_splits = 5, verbose = 5)
#' print(dpl.mos)
#' \dontrun{plot(dpl.mos)}


mos <- function(NP = 40, G = 100, data, class.name,
                structure = c("nb", "tancl", "tan", "hc"),
                pB = 0.05, edgelist = NULL, n_splits = 10,
                min_participation = 0.1, verbose = 25, ...) {

  if (NP <= 5) {
    warning("'NP' <= 5; set to default value 40\n", immediate. = TRUE)
    NP <- 40
  }
  if (G <= 1) {
    warning("'G' <= 1; set to default value 100\n", immediate. = TRUE)
    G <- 100
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
  best_idx <- which.max(fitness)
  best <- pop[[best_idx]]
  record_CLL <- c(fitness[best_idx])
  record_evals <- c(neval)

  # SHADE memory
  H <- 6; ki <- 1; terminal <- 0
  mCR <- rep(0.5, H); mF <- mCR
  Archive <- list()
  pBest_n <- max(1, round(NP * pB))

  # Solis-Wets state
  sw_rho <- 0.5; sw_bias <- rep(0, dim); sw_ns <- 0L; sw_nf <- 0L

  # Technique participation
  participation <- c(shade = 0.5, sw = 0.5)
  G_per_split <- max(1, floor(G / n_splits))

  for (s in seq_len(n_splits)) {
    G_shade <- max(1, round(participation["shade"] * G_per_split))
    sw_budget <- max(5, round(participation["sw"] * G_per_split * NP))

    # ============================
    # Technique 1: SHADE generations
    # ============================
    f_before_shade <- fitness[best_idx]

    for (g in seq_len(G_shade)) {
      SF <- c(); SCR <- c(); improvement <- c()
      for (j in seq_len(NP)) {
        ri <- sample.int(H, 1)
        CR <- if (mCR[ri] == terminal) 0 else randN(1, mCR[ri])
        Fval <- randC(1, mF[ri])

        idxs <- seq_len(NP)[-j]
        idrs <- seq_len(NP + length(Archive))[-j]
        idbs <- getPBest(fitness, pBest_n)

        xi <- pop[[j]]
        xp <- pop[[sample(idbs, 1)]]
        r1 <- pop[[sample(idxs, 1)]]
        r2 <- c(pop, Archive)[[sample(idrs, 1)]]

        mutant <- vec(xi$.params) + Fval * (vec(xp$.params) - vec(xi$.params)) +
          Fval * (vec(r1$.params) - vec(r2$.params))
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
          SCR <- c(SCR, CR); SF <- c(SF, Fval)
          if (f > fitness[best_idx]) {
            best_idx <- j
            best <- trial
          }
        }
      }

      if (length(Archive) > NP) {
        Archive <- Archive[sample(seq_len(length(Archive)), NP)]
      }
      if (length(SCR) > 0 && length(SF) > 0) {
        if (mCR[ki] == terminal || max(SCR) == 0) {
          mCR[ki] <- terminal
        } else {
          mCR[ki] <- meanWL(SCR, improvement)
        }
        mF[ki] <- meanWL(SF, improvement)
        ki <- (ki %% H) + 1
      }

      record_CLL <- c(record_CLL, fitness[best_idx])
      record_evals <- c(record_evals, neval)
    }

    quality_shade <- max(0, fitness[best_idx] - f_before_shade)

    # ============================
    # Technique 2: Solis-Wets
    # ============================
    f_before_sw <- fitness[best_idx]
    best_vec <- vec(best$.params)

    sw_result <- solisWetsSearch(best_vec, fitness[best_idx], sw_rho, sw_bias,
                                 sw_ns, sw_nf, sw_budget, group_indices,
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

    record_CLL <- c(record_CLL, fitness[best_idx])
    record_evals <- c(record_evals, neval)

    quality_sw <- max(0, fitness[best_idx] - f_before_sw)

    # ============================
    # Adjust participation
    # ============================
    total_quality <- quality_shade + quality_sw
    if (total_quality > 0) {
      participation["shade"] <- quality_shade / total_quality
      participation["sw"] <- quality_sw / total_quality
      participation <- pmax(participation, min_participation)
      participation <- participation / sum(participation)
    }

    if (verbose > 0) {
      gen_approx <- s * G_per_split
      if (gen_approx %% verbose == 0 || s == n_splits) {
        cat("Split: ", s, "/", n_splits, "\t CLL= ", fitness[best_idx],
            "\t SHADE=", round(participation["shade"], 2),
            "\t SW=", round(participation["sw"], 2), "\n")
      }
    }
  }

  outDE <- list(Best = best, BestCLL = fitness[best_idx],
                pobFinal = pop, CLLPobFinal = fitness,
                N.evals = neval, convergence = record_CLL,
                evaluations = record_evals)
  attr(outDE, "class") <- "DE"
  outDE
}
