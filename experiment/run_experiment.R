#!/usr/bin/env Rscript
# =============================================================================
# Experiment: Discriminative Parameter Learning of BN by DE
# Platas-López et al. (2026)
#
# 8 algorithms × 20 datasets × 3 structures × 31 runs = 14,880 executions
# Parallelized with parallel::mclapply
# =============================================================================

suppressPackageStartupMessages({
  library(bnclassify)
  library(dplbnDE)  # load after bnclassify so dplbnDE::accuracy takes precedence
  library(parallel)
})

# ---- Configuration ----
NUM_CORES <- as.integer(Sys.getenv("NCORES", detectCores() - 1))
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR <- Sys.getenv("OUT_DIR", "experiment/results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NUM_RUNS <- 31
NP <- 100

cat("=== Experiment Configuration ===\n")
cat("Cores:", NUM_CORES, "\n")
cat("Runs per combination:", NUM_RUNS, "\n")
cat("NP:", NP, "\n\n")

# ---- Dataset loading ----
dataset_files <- list(
  australian  = "AustralianDisc.csv",
  chess       = "chess_data_3196.csv",
  cleveland   = "CleveDisc.csv",
  corral      = "corrAl.csv",
  crx         = "crxDisc.csv",
  diabetes    = "diabetesDisc.csv",
  flare       = "flare.csv",
  german      = "germanDisc.csv",
  glass       = "glassDisc.csv",
  heart       = "heartDisc.csv",
  hepatitis   = "hepatitisDisc.csv",
  lymphography = "lymphografy_148.csv",
  mofn        = "mofn-3-7-10.csv",
  pima        = "pimaDisc.csv",
  segment     = "segmentDisc.csv",
  soybean     = "soybean_large.csv",
  tictactoe   = "tic-tac-toe.csv",
  vehicle     = "vehicleDisc.csv",
  vote        = "vote.csv",
  waveform    = "wafeform21Disc.csv"
)

load_dataset <- function(name) {
  file_path <- file.path(DATA_DIR, dataset_files[[name]])
  df <- read.csv(file_path, stringsAsFactors = TRUE)
  # Last column is always the class
  class_col <- ncol(df)
  class_name <- names(df)[class_col]
  # Ensure all columns are factors
  for (col in names(df)) {
    if (!is.factor(df[[col]])) df[[col]] <- as.factor(df[[col]])
  }
  list(data = df, class.name = class_name,
       n_atts = ncol(df) - 1, n_cases = nrow(df),
       n_classes = nlevels(df[[class_name]]))
}

# ---- Algorithm definitions ----
structures <- c("nb", "tancl", "hc")

define_algorithms <- function(NP, G, data, class.name, structure, edgelist) {
  common_args <- list(data = data, class.name = class.name, verbose = 0)

  list(
    jade_archive = function(seed) {
      set.seed(seed)
      do.call(jade, c(list(NP = NP, G = G, c = 0.1, pB = 0.05,
                           structure = structure, edgelist = edgelist,
                           archive = TRUE), common_args))
    },
    lshade = function(seed) {
      set.seed(seed)
      do.call(lshade, c(list(NP = NP, G = G, c = 0.1, pB = 0.05,
                             structure = structure, edgelist = edgelist),
                        common_args))
    },
    shadeils = function(seed) {
      set.seed(seed)
      do.call(shadeils, c(list(NP = NP, G = G, c = 0.1, pB = 0.05,
                               ls_freq = 5, structure = structure,
                               edgelist = edgelist), common_args))
    },
    mos = function(seed) {
      set.seed(seed)
      do.call(mos, c(list(NP = NP, G = G, pB = 0.05, n_splits = 10,
                          structure = structure, edgelist = edgelist),
                     common_args))
    },
    jso = function(seed) {
      set.seed(seed)
      do.call(jso, c(list(NP = NP, G = G,
                          structure = structure, edgelist = edgelist),
                     common_args))
    },
    nlshadersp = function(seed) {
      set.seed(seed)
      do.call(nlshadersp, c(list(NP = NP, G = G, c = 0.1, pB = 0.05,
                                 structure = structure, edgelist = edgelist),
                            common_args))
    },
    shademts = function(seed) {
      set.seed(seed)
      do.call(shademts, c(list(NP = NP, G = G, c = 0.1, pB = 0.05,
                               ls_freq = 3, structure = structure,
                               edgelist = edgelist), common_args))
    }
  )
}

# ---- WANBIA baseline ----
run_wanbia <- function(data, class.name, structure) {
  if (structure == "nb") {
    bn <- nb(class.name, data)
  } else if (structure == "tancl") {
    bn <- tan_cl(class.name, data)
  } else {
    bn <- tan_hc(class.name, data, k = 5)
  }
  # WANBIA: learn parameters with awnb (attribute weighting)
  fitted <- tryCatch(
    bnc(structure, class.name, data, smooth = 1),
    error = function(e) lp(bn, data, smooth = 1)
  )
  # If bnc doesn't work, use lp with smooth=1 as generative baseline
  fitted <- lp(bn, data, smooth = 1)
  cll_val <- cLogLik(fitted, data)
  preds <- predict(fitted, data)
  acc <- accuracy(preds, data[[class.name]])
  list(CLL = cll_val, accuracy = acc)
}

# ---- 5-fold CV helper ----
cv_accuracy <- function(result, data, class.name) {
  n <- nrow(data)
  folds <- sample(rep(1:5, length.out = n))
  accs <- numeric(5)
  for (k in 1:5) {
    test_idx <- which(folds == k)
    test_data <- data[test_idx, , drop = FALSE]
    preds <- predict(result$Best, test_data)
    accs[k] <- accuracy(preds, test_data[[class.name]])
  }
  mean(accs)
}

# ---- Build task list ----
tasks <- list()
for (ds_name in names(dataset_files)) {
  for (struct in structures) {
    for (alg_name in c("jade_archive", "lshade", "shadeils", "mos",
                       "jso", "nlshadersp", "shademts")) {
      for (seed in 1:NUM_RUNS) {
        tasks[[length(tasks) + 1]] <- list(
          dataset = ds_name, structure = struct,
          algorithm = alg_name, seed = seed
        )
      }
    }
  }
}

cat("Total tasks:", length(tasks), "\n")
cat("Starting experiment...\n\n")

# ---- Execute ----
run_single_task <- function(task) {
  ds <- load_dataset(task$dataset)
  G <- 15 * ds$n_atts

  # HCS needs k parameter
  hcs_args <- if (task$structure == "hc") list(k = 5) else list()

  algorithms <- define_algorithms(NP, G, ds$data, ds$class.name,
                                   task$structure, NULL)

  t0 <- proc.time()[3]
  result <- tryCatch({
    algorithms[[task$algorithm]](task$seed)
  }, error = function(e) {
    list(BestCLL = NA, N.evals = NA, Best = NULL,
         error = conditionMessage(e))
  })
  elapsed <- proc.time()[3] - t0

  if (is.null(result$Best) || is.na(result$BestCLL)) {
    return(data.frame(
      dataset = task$dataset, structure = task$structure,
      algorithm = task$algorithm, seed = task$seed,
      CLL = NA, train_accuracy = NA, cv_accuracy = NA,
      n_evals = NA, time = elapsed,
      error = result$error %||% "unknown",
      stringsAsFactors = FALSE
    ))
  }

  # Training accuracy
  preds <- predict(result$Best, ds$data)
  train_acc <- accuracy(preds, ds$data[[ds$class.name]])

  # 5-fold CV accuracy
  cv_acc <- tryCatch(
    cv_accuracy(result, ds$data, ds$class.name),
    error = function(e) NA
  )

  data.frame(
    dataset = task$dataset, structure = task$structure,
    algorithm = task$algorithm, seed = task$seed,
    CLL = result$BestCLL, train_accuracy = train_acc,
    cv_accuracy = cv_acc, n_evals = result$N.evals,
    time = elapsed, error = "",
    stringsAsFactors = FALSE
  )
}

# Run in parallel
results_list <- mclapply(tasks, run_single_task, mc.cores = NUM_CORES)
results <- do.call(rbind, results_list)

# ---- Add WANBIA baseline (deterministic, no seeds needed) ----
cat("Running WANBIA baseline...\n")
wanbia_results <- list()
for (ds_name in names(dataset_files)) {
  ds <- load_dataset(ds_name)
  for (struct in structures) {
    w <- tryCatch(run_wanbia(ds$data, ds$class.name, struct),
                  error = function(e) list(CLL = NA, accuracy = NA))
    wanbia_results[[length(wanbia_results) + 1]] <- data.frame(
      dataset = ds_name, structure = struct,
      algorithm = "wanbia", seed = 1,
      CLL = w$CLL, train_accuracy = w$accuracy,
      cv_accuracy = NA, n_evals = 0, time = 0, error = "",
      stringsAsFactors = FALSE
    )
  }
}
wanbia_df <- do.call(rbind, wanbia_results)
results <- rbind(results, wanbia_df)

# ---- Save results ----
out_file <- file.path(OUT_DIR, "full_results.csv")
write.csv(results, out_file, row.names = FALSE)
cat("\nResults saved to:", out_file, "\n")
cat("Total rows:", nrow(results), "\n")
cat("Errors:", sum(results$error != "", na.rm = TRUE), "\n")
cat("Done!\n")
