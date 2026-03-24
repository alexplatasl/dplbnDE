# Quick smoke test: 1 dataset, 1 structure, 1 seed, all algorithms
devtools::load_all("dplbnDE")
library(bnclassify)

data_dir <- "data"
df <- read.csv(file.path(data_dir, "heartDisc.csv"), stringsAsFactors = TRUE)
for (col in names(df)) if (!is.factor(df[[col]])) df[[col]] <- as.factor(df[[col]])
cn <- names(df)[ncol(df)]
n_atts <- ncol(df) - 1
G <- 15 * n_atts
NP <- 30  # smaller for smoke test

cat("Dataset: heart, Atts:", n_atts, ", G:", G, ", NP:", NP, "\n\n")

algorithms <- list(
  jade_archive = function() jade(NP=NP, G=G, data=df, class.name=cn, c=0.1, pB=0.05, archive=TRUE, verbose=0),
  lshade       = function() lshade(NP=NP, G=G, data=df, class.name=cn, c=0.1, pB=0.05, verbose=0),
  shadeils     = function() shadeils(NP=NP, G=G, data=df, class.name=cn, c=0.1, pB=0.05, ls_freq=5, verbose=0),
  mos          = function() mos(NP=NP, G=G, data=df, class.name=cn, pB=0.05, n_splits=10, verbose=0),
  jso          = function() jso(NP=NP, G=G, data=df, class.name=cn, verbose=0),
  nlshadersp   = function() nlshadersp(NP=NP, G=G, data=df, class.name=cn, c=0.1, pB=0.05, verbose=0),
  shademts     = function() shademts(NP=NP, G=G, data=df, class.name=cn, c=0.1, pB=0.05, ls_freq=3, verbose=0)
)

for (alg_name in names(algorithms)) {
  set.seed(42)
  t0 <- proc.time()[3]
  r <- algorithms[[alg_name]]()
  elapsed <- proc.time()[3] - t0
  acc <- accuracy(predict(r$Best, df), df[[cn]])
  cat(sprintf("%-15s CLL=%10.2f  Acc=%.4f  Evals=%5d  Time=%.1fs\n",
              alg_name, r$BestCLL, acc, r$N.evals, elapsed))
}

# WANBIA baseline
bn <- nb(cn, df)
fitted <- lp(bn, df, smooth = 1)
w_cll <- cLogLik(fitted, df)
w_acc <- accuracy(predict(fitted, df), df[[cn]])
cat(sprintf("%-15s CLL=%10.2f  Acc=%.4f\n", "wanbia", w_cll, w_acc))

cat("\nSmoke test PASSED!\n")
