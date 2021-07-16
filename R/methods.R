#' Print
#' @export
#'
#'
#' @param object A \code{\link{DE}} object.
#'
print.DE <- function(x,...){
  cat("Number of evaluations: \t", x$N.evals, "\n")
  cat("Final population size: \t", length(x$pobFinal), "\n\n")
  cat("Summary results of fitness in final population: \n\n")
  cat("Best CLL: \t", x$BestCLL, "\n")
  cat("Worst CLL: \t", min(x$CLLPobFinal), "\n")
  cat("Median: \t", median(x$CLLPobFinal), "\n")
  cat("Std. Dev.: \t", sd(x$CLLPobFinal), "\n")
}

#' Plot main results of evolution
#' @export
#'
#' @param object A \code{\link{DE}} object.
#'
plot.DE <- function(x,...){
  par(mfrow=c(1,2))
  hist(x$CLLPobFinal, xlab = "CLL",  main = "CLL of final population", breaks = sqrt(length(x$pobFinal)))
  plot(x$evaluations, x$convergence, type = "l", xlab = "Evaluations", ylab = "CLL",
       main = "Convergence plot")
}
