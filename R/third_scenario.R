#' Improved k-progressive simulation function 
#' 
#' @param count_matrices_dir Path to directory containing count matrices (.RDS files) OR a list of count matrices
#' @param adjm_file Path to adjacency matrix file (.txt format)
#' @param max_k Maximum k to test (default: NULL, uses all available matrices)
#' @param n_runs Number of simulation runs to perform for each k (default: 10)
#' @param n_cores Number of CPU cores to use (default: 15)
#' @param genes_on_rows Logical. TRUE if genes are on rows, FALSE if on columns (default: TRUE)
#' @param quantile_threshold Quantile threshold for cutoff_adjacency (default: 0.99)
#' @param output_file Optional path to save results
#' @param verbose Print progress messages (default: TRUE)
#' 
#' @return List with summary and detailed results across all k values
run_k_progressive_study <- function(count_matrices_dir,
                                    adjm_file, 
                                    max_k = NULL,
                                    n_runs = 10,
                                    n_cores = 15,
                                    genes_on_rows = TRUE,
                                    quantile_threshold = 0.99,
                                    output_file = NULL,
                                    verbose = TRUE) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scGraphVerse)
    library(INetTool)
    library(BiocParallel)
    library(glue)
    library(tools)
  })
  
  modules <- init_py(python_path = "/usr/bin/python3", required = TRUE)
  
  # Validate and load data
  if (!file.exists(adjm_file)) {
    stop("Adjacency matrix file not found: ", adjm_file)
  }
  
  # Load count matrices
  if (is.character(count_matrices_dir) && length(count_matrices_dir) == 1) {
    if (!dir.exists(count_matrices_dir)) {
      stop("Count matrices directory not found: ", count_matrices_dir)
    }
    rds_files <- list.files(count_matrices_dir, pattern = "\\.RDS$|\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) {
      stop("No RDS files found in directory: ", count_matrices_dir)
    }
    count_matrices <- lapply(rds_files, readRDS)
    names(count_matrices) <- tools::file_path_sans_ext(basename(rds_files))
  } else if (is.list(count_matrices_dir)) {
    count_matrices <- count_matrices_dir
  } else {
    stop("count_matrices_dir must be either a directory path or a list of matrices")
  }
  
  if (!genes_on_rows) {
    count_matrices <- lapply(count_matrices, t)
  }
  
  total_matrices <- length(count_matrices)
  if (is.null(max_k)) {
    max_k <- total_matrices
  } else {
    max_k <- min(max_k, total_matrices)
  }
  count_matrices <- count_matrices[1:max_k]
  
  # Load adjacency matrix
  adjm <- as.matrix(read.table(adjm_file))
  colnames(adjm) <- rownames(adjm)
  
  n_genes <- nrow(adjm)
  n_cells <- ncol(count_matrices[[1]])
  ratio <- n_genes / n_cells
  
  if (verbose) {
    cat("=== K-Progressive Study ===\n")
    cat("Testing k values: 1 to", max_k, "\n")
    cat("Genes:", n_genes, "| Cells:", n_cells, "| Runs per k:", n_runs, "\n")
    cat("Note: JRF will only be tested for k >= 2 (requires multiple datasets)\n\n")
  }
  
  # Ground truth community
  adj_comm <- community_path(adjm, plot = FALSE, verbose = FALSE, genes_path = Inf)
  
  all_results <- list()
  
  # Helper function for single analysis
  run_single_method <- function(matrices, method, method_modules = NULL) {
    # Convert list of matrices to MultiAssayExperiment
    mae <- create_mae(matrices)

    timing <- system.time({
      networks <- infer_networks(mae, method = method,
                                 grnboost_modules = method_modules,
                                 nCores = n_cores)
    })

    # generate_adjacency and cutoff_adjacency return SummarizedExperiment
    wadj_se <- generate_adjacency(networks)
    sym_wadj_se <- symmetrize(wadj_se, weight_function = "mean")
    binary_adj_se <- cutoff_adjacency(mae, sym_wadj_se, n = 2, method = method,
                                      grnboost_modules = method_modules,
                                      nCores = n_cores,
                                      quantile_threshold = quantile_threshold,
                                      debug = FALSE)

    # Extract assays for pscores
    binary_adj_list <- SummarizedExperiment::assays(binary_adj_se)

    scores <- pscores(adjm, binary_adj_list)$Statistics %>%
      mutate(Method = method,
             Time_seconds = timing["elapsed"],
             Time_minutes = timing["elapsed"] / 60)

    scores <- add_community_metrics(scores, binary_adj_list)

    return(list(scores = scores, binary_adj = binary_adj_se, sym_wadj = sym_wadj_se))
  }
  
  # other helpers
  create_consensus_scores <- function(binary_adj_se, sym_wadj_se, method_name, timing) {
    # create_consensus accepts SummarizedExperiment and returns SE, extract the matrix
    consensus_vote_se <- create_consensus(binary_adj_se, method = "vote")
    consensus_vote <- SummarizedExperiment::assays(consensus_vote_se)[[1]]

    consensus_union_se <- create_consensus(binary_adj_se, method = "union")
    consensus_union <- SummarizedExperiment::assays(consensus_union_se)[[1]]

    consensus_inet_se <- create_consensus(binary_adj_se, sym_wadj_se,
                                          method = "INet", threshold = 0.05, ncores = n_cores)
    consensus_inet <- SummarizedExperiment::assays(consensus_inet_se)[[1]]

    scores <- pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))$Statistics %>%
      mutate(Method = method_name,
             Predicted_Matrix = c("vote", "union", "inet"),
             Integration = "late",
             Time_seconds = timing["elapsed"],
             Time_minutes = timing["elapsed"] / 60)

    scores <- add_community_metrics(scores, list(consensus_vote, consensus_union, consensus_inet))

    return(scores)
  }
  
  add_community_metrics <- function(scores_df, consensus_matrices) {
    tryCatch({
      communities <- lapply(consensus_matrices, function(mat) {
        community_path(mat, plot = FALSE, verbose = FALSE, genes_path = Inf)
      })
      topscore <- community_similarity(adj_comm, communities)
      
      required_cols <- c("VI", "NMI", "ARI", "Modularity", "Communities", "Density", "Transitivity")
      for (col in required_cols) {
        if (!col %in% colnames(scores_df)) {
          scores_df[[col]] <- NA_real_
        }
      }
      
      if (!is.null(topscore$topology_measures) && nrow(topscore$topology_measures) > 0) {
        for (i in 1:min(nrow(topscore$topology_measures), nrow(scores_df))) {
          scores_df$Modularity[i] <- topscore$topology_measures$Modularity[i]
          scores_df$Communities[i] <- topscore$topology_measures$Communities[i]
          scores_df$Density[i] <- topscore$topology_measures$Density[i]
          scores_df$Transitivity[i] <- topscore$topology_measures$Transitivity[i]
        }
      }
      
      if (!is.null(topscore$community_metrics) && nrow(topscore$community_metrics) > 0) {
        for (i in 1:min(nrow(topscore$community_metrics), nrow(scores_df))) {
          scores_df$VI[i] <- topscore$community_metrics$VI[i]
          scores_df$NMI[i] <- topscore$community_metrics$NMI[i]
          scores_df$ARI[i] <- topscore$community_metrics$ARI[i]
        }
      }
    }, error = function(e) {
      if (verbose) cat("Community analysis failed:", e$message, "\n")
    })
    
    return(scores_df)
  }
  
  # Main loop over k values
  for (k in 1:max_k) {
    if (verbose) cat(glue("=== Testing k = {k} ===\n"))
    
    current_matrices <- count_matrices[1:k]
    k_results <- list()
    
    for (run_id in 1:n_runs) {
      if (verbose) cat(glue("K={k}, Run={run_id}/{n_runs}\n"))
      
      all_scores <- tibble()
      
      if (k == 1) {
        # === K=1: Only single-dataset methods (GENIE3, GRNBoost2) ===
        
        # GENIE3
        genie3_result <- run_single_method(current_matrices, "GENIE3")
        genie3_scores <- genie3_result$scores %>%
          mutate(Integration = "single", Predicted_Matrix = "single")
        all_scores <- bind_rows(all_scores, genie3_scores)
        
        # GRNBoost2
        grnboost_result <- run_single_method(current_matrices, "GRNBoost2", modules)
        grnboost_scores <- grnboost_result$scores %>%
          mutate(Integration = "single", Predicted_Matrix = "single")
        all_scores <- bind_rows(all_scores, grnboost_scores)
        
        # JRF is NOT run for k=1 
        if (verbose) cat("JRF skipped for k=1 (requires multiple datasets)\n")
        
      } else {
        # === K>1: All methods including multi-dataset strategies ===
        
        # Convert current_matrices to MultiAssayExperiment for API compatibility
        current_mae <- create_mae(current_matrices)

        # Prepare early matrix (earlyj takes MAE as input and returns MAE with merged data)
        # Extract the merged experiment as a list for run_single_method
        early_mae <- earlyj(current_mae, rowg = TRUE)
        early_matrix <- as.list(MultiAssayExperiment::experiments(early_mae))

        # GENIE3 - Late integration
        timing_genie3_late <- system.time({
          networks <- infer_networks(current_mae, method = "GENIE3", nCores = n_cores)
        })
        wadj_se <- generate_adjacency(networks)
        sym_wadj_se <- symmetrize(wadj_se, weight_function = "mean")
        binary_adj_se <- cutoff_adjacency(current_mae, sym_wadj_se, n = 3, method = "GENIE3",
                                          nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        genie3_late_scores <- create_consensus_scores(binary_adj_se, sym_wadj_se, "GENIE3", timing_genie3_late)
        all_scores <- bind_rows(all_scores, genie3_late_scores)

        # GENIE3 - Early integration
        genie3_early_result <- run_single_method(early_matrix, "GENIE3")
        genie3_early_scores <- genie3_early_result$scores %>%
          mutate(Integration = "early", Predicted_Matrix = "early")
        all_scores <- bind_rows(all_scores, genie3_early_scores)

        # GRNBoost2 - Late integration
        timing_grnboost_late <- system.time({
          networks <- infer_networks(current_mae, method = "GRNBoost2",
                                     grnboost_modules = modules, nCores = n_cores)
        })
        wadj_se <- generate_adjacency(networks)
        sym_wadj_se <- symmetrize(wadj_se, weight_function = "mean")
        binary_adj_se <- cutoff_adjacency(current_mae, sym_wadj_se, n = 3, method = "GRNBoost2",
                                          grnboost_modules = modules, nCores = n_cores,
                                          quantile_threshold = quantile_threshold, debug = FALSE)

        grnboost_late_scores <- create_consensus_scores(binary_adj_se, sym_wadj_se, "GRNBoost2", timing_grnboost_late)
        all_scores <- bind_rows(all_scores, grnboost_late_scores)
        
        # GRNBoost2 - Early integration
        grnboost_early_result <- run_single_method(early_matrix, "GRNBoost2", modules)
        grnboost_early_scores <- grnboost_early_result$scores %>%
          mutate(Integration = "early", Predicted_Matrix = "early")
        all_scores <- bind_rows(all_scores, grnboost_early_scores)
        
        # JRF - Joint modeling (only for k >= 2)
        jrf_scores <- tibble(
          Predicted_Matrix = c("vote", "union", "inet"),
          Method = "JRF", Integration = "joint",
          TP = 0, TN = 0, FP = 0, FN = 0,
          TPR = 0, FPR = 0, Precision = 0, F1 = 0, MCC = 0,
          VI = NA_real_, NMI = NA_real_, ARI = NA_real_,
          Modularity = NA_real_, Communities = NA_real_,
          Density = NA_real_, Transitivity = NA_real_,
          Time_seconds = NA_real_, Time_minutes = NA_real_
        )
        
        timing_jrf <- system.time({
          tryCatch({
            networks <- infer_networks(current_mae, method = "JRF", nCores = n_cores)
            if (!is.null(networks) && length(networks) > 0 && !is.null(networks[[1]])) {
              wadj_se <- generate_adjacency(networks)
              sym_wadj_se <- symmetrize(wadj_se, weight_function = "mean")
              binary_adj_se <- cutoff_adjacency(current_mae, sym_wadj_se, n = 3, method = "JRF",
                                                nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)

              jrf_scores <- create_consensus_scores(binary_adj_se, sym_wadj_se, "JRF", timing_jrf) %>%
                mutate(Integration = "joint")
            }
          }, error = function(e) {
            if (verbose) cat("JRF failed:", e$message, "\n")
          })
        })
        
        all_scores <- bind_rows(all_scores, jrf_scores)
      }
      
      # Add metadata to all scores
      all_scores <- all_scores %>%
        mutate(
          k = k,
          Run = run_id,
          Ratio = ratio,
          p = n_genes,
          n_cells = n_cells
        )
      
      k_results[[run_id]] <- all_scores
      
      BiocParallel::bpstop(BiocParallel::bpparam())
      if (run_id < n_runs) Sys.sleep(1)
    }
    
    all_results[[as.character(k)]] <- bind_rows(k_results, .id = "Run")
    if (verbose) cat(glue("Completed k = {k}\n\n"))
  }
  
  # Combine and summarize results
  final_df <- bind_rows(all_results, .id = "k")
  final_df$k <- as.numeric(final_df$k)
  
  summary_df <- final_df %>%
    group_by(k, Method, Integration, Predicted_Matrix) %>%
    summarise(
      across(where(is.numeric) & !matches("k|Run"), 
             list(mean = ~mean(.x, na.rm = TRUE), 
                  sd = ~sd(.x, na.rm = TRUE)), 
             .names = "{.col}_{.fn}"), 
      .groups = "drop"
    ) %>%
    arrange(k, Method, Integration, Predicted_Matrix)
  
  if (!is.null(output_file)) {
    write.table(final_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) cat("Results saved to:", output_file, "\n")
  }
  
  if (verbose) {
    cat("\n=== Study Completed ===\n")
    cat("Methods and integrations tested:\n")
    method_summary <- final_df %>%
      count(k, Method, Integration) %>%
      arrange(k, Method, Integration)
    print(method_summary)
  }
  
  return(list(
    summary = summary_df,
    detailed = final_df,
    parameters = list(
      max_k = max_k,
      n_runs = n_runs,
      n_genes = n_genes,
      n_cells = n_cells,
      methods_by_k = final_df %>% 
        distinct(k, Method, Integration) %>% 
        arrange(k, Method, Integration)
    )
  ))
}

# Example
# count_mat <- readRDS("path/sim_n100p500k5.RDS")
# results <- run_k_progressive_study(
#  count_matrices_dir = count_mat,
#  adjm_file = "path/adjm_p500.txt",
#  max_k = 5, 
#  n_runs = 10,
#  genes_on_rows = TRUE,
#  n_cores = 15,
#  verbose = TRUE,
#  output_file = "path/sim_n100p500k5_summary.txt"
# )

