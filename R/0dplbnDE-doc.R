#' Discriminative Parameter Learning of Bayesian Networks by Differential Evolution
#'
#' Implements Differential Evolution (DE) to train parameters of Bayesian Networks
#' (BN) for optimizing the Conditional Log-Likelihood (Discriminative Learning)
#' instead of the log-likelihood (Generative Learning). Any given BN structure
#' encodes assumptions about conditional independencies among the attributes and
#' will result  in error if they do not hold in the data. Such an error includes
#' the classification dimension. The main goal of Discriminative learning is minimize
#' this type of error.
#'
#'
#' DE variants:
#' Based on different strategies followed by the operators of DE, there are
#' different variants, which define the way in which the mutant and trial
#' vectors are generated. The most popular variant is called \code{DE/rand/1/bin}, where ``DE''
#' means Differential Evolution, the word ``rand'' indicates that the so-called base vector
#' is randomly chosen, ``1'' is the number of vector pairs (i.e., vector differences to be
#' calculated) chosen, and finally ``bin'' means that a binomial recombination is chosen.
#' The following is a list with the available variants:
#'
#' \itemize{
#' \item \code{\link{DErand}}: Implements \code{DE/rand/} variant with 1 or 2 pairs
#' of vector differences, and binomial or exponential recombination. (Price and Storn, 1996)
#' \item \code{\link{DEbest}}: Implements \code{DE/best/} variant with 1 or 2 pairs
#' of vector differences and binomial or exponential recombination. (Price and Storn, 1996)
#' \item \code{\link{jade}}: A variant that includes some mechanisms to decrease the
#' dependence to its parameter values such as F and CR. JADE uses a mutation strategy
#' called \code{DE/current-to-pbest}, where \eqn{p \in (0 , 1]}. Base vectors are selected
#' from the best \eqn{100p}% vectors in the current population. JADE, with the aim to
#' maintaining diversity, uses an optional external archive. (Zhang and Sanderson, 2009)
#' \item \code{\link{lshade}}: An improved version of JADE, LSHADE incorporates a
#' Linear Population Size Reduction (LPSR) in order to improve the performance. (Tanabe and Fukunaga, 2014)
#' }
#'
#' @name dplbnDE-package
#' @docType package
"_PACKAGE"
#' @importFrom stats median sd runif rnorm rcauchy
#' @importFrom graphics hist
#' @importFrom bnclassify tan_cl tan_hc nb lp cLogLik
#'
#' @references Price K and Storn R (1996), Minimizing the real functions of the
#' icec’96 contest by differential evolution. In \emph{Proc. of IEEE C. Evol.
#' Computat.}, pp. 842--844.
#'
#'   Zhang J and Sanderson A (2009). Jade: adaptive differential evolution with
#'   optional external archive. \emph{IEEE Trans. Evol. Comput.}, pp. 945--958.
#'
#'   Tanabe R and Fukunaga A (2014). Improving the search performance of shade
#'   using linear population size reduction. In \emph{Proc. of IEEE C. Evol.
#'   Computat.}, pp. 1658--1665.
NULL

#' Discriminative parameter learning of bayesian networks by differential evolution
#'
#' A list with Bayesian Networks, structure and parameters learned in a discriminative way by
#' diferential evolution returned by \code{\link{lshade}}, \code{\link{jade}},
#' \code{\link{DEbest}}, or \code{\link{DErand}} functions.
#'
#'@name DE
#'
#'
#'
#' @examples
#' data(car)
#' dpl.lshade <- lshade(NP = 40, G = 50, data = car, class.name = names(car)[7], c = 0.1,
#' structure = "tan", pB = 0.05, edgelist = NULL, verbose = 5)
#' dpl.lshade
NULL

#' Car Evaluation Data Set.
#'
#' Data set from the UCI repository:
#' \url{https://archive.ics.uci.edu/ml/datasets/Car+Evaluation}.
#'
#' @format A \code{data.frame} with 7 columns and 1728 rows.
#' @docType data
#' @name car
NULL
