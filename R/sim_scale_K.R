#' Simplified k-progressive simulation function 
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
  
  # Load required libraries
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scGraphVerse)
    library(INetTool)
    library(BiocParallel)
    library(glue)
    library(tools)
  })
  
  # Initialize Python environment for GRNBoost2
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
  
  # Handle gene orientation
  if (!genes_on_rows) {
    count_matrices <- lapply(count_matrices, t)
  }
  
  # Determine max_k
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
  
  # Calculate metrics
  n_genes <- nrow(adjm)
  n_cells <- ncol(count_matrices[[1]])
  ratio <- n_genes / n_cells
  
  if (verbose) {
    cat("=== K-Progressive Study ===\n")
    cat("Testing k values: 1 to", max_k, "\n")
    cat("Genes:", n_genes, "| Cells:", n_cells, "| Runs per k:", n_runs, "\n\n")
  }
  
  # Ground truth community
  adj_comm <- community_path(adjm, plot = FALSE, verbose = FALSE, genes_path = Inf)
  
  # Storage for all results
  all_results <- list()
  
  # Helper function to add community metrics safely
  add_community_metrics <- function(scores_df, consensus_matrices) {
    tryCatch({
      communities <- lapply(consensus_matrices, function(mat) {
        community_path(mat, plot = FALSE, verbose = FALSE, genes_path = Inf)
      })
      topscore <- community_similarity(adj_comm, communities)
      
      # Add missing columns
      required_cols <- c("VI", "NMI", "ARI", "Modularity", "Communities", "Density", "Transitivity")
      for (col in required_cols) {
        if (!col %in% colnames(scores_df)) {
          scores_df[[col]] <- NA_real_
        }
      }
      
      # Add topology measures
      if (!is.null(topscore$topology_measures) && nrow(topscore$topology_measures) > 0) {
        for (i in 1:min(nrow(topscore$topology_measures), nrow(scores_df))) {
          scores_df$Modularity[i] <- topscore$topology_measures$Modularity[i]
          scores_df$Communities[i] <- topscore$topology_measures$Communities[i]
          scores_df$Density[i] <- topscore$topology_measures$Density[i]
          scores_df$Transitivity[i] <- topscore$topology_measures$Transitivity[i]
        }
      }
      
      # Add community metrics
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
  
  # Main loop
  for (k in 1:max_k) {
    if (verbose) cat(glue("=== Testing k = {k} ===\n"))
    
    current_matrices <- count_matrices[1:k]
    k_results <- list()
    
    for (run_id in 1:n_runs) {
      if (verbose) cat(glue("K={k}, Run={run_id}/{n_runs}\n"))
      
      timing <- list()
      all_scores <- tibble()
      
      # === k=1: Special case ===
      if (k == 1) {
        
        # GENIE3
        timing[["GENIE3"]] <- system.time({
          networks <- infer_networks(current_matrices, method = "GENIE3", nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(current_matrices, sym_wadj, n = 1, method = "GENIE3", nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        scores <- pscores(adjm, binary_adj)$Statistics %>%
          mutate(Method = "GENIE3", Predicted_Matrix = "single")
        scores <- add_community_metrics(scores, binary_adj)
        
        # Duplicate as early and late
        genie3_scores <- bind_rows(
          scores %>% mutate(Predicted_Matrix = "early"),
          scores %>% mutate(Predicted_Matrix = "late")
        )
        all_scores <- bind_rows(all_scores, genie3_scores)
        
        # GRNBoost2
        timing[["GRNBoost2"]] <- system.time({
          networks <- infer_networks(current_matrices, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(current_matrices, sym_wadj, n = 1, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        scores <- pscores(adjm, binary_adj)$Statistics %>%
          mutate(Method = "GRNBoost2", Predicted_Matrix = "single")
        scores <- add_community_metrics(scores, binary_adj)
        
        # Duplicate as early and late
        grnboost_scores <- bind_rows(
          scores %>% mutate(Predicted_Matrix = "early"),
          scores %>% mutate(Predicted_Matrix = "late")
        )
        all_scores <- bind_rows(all_scores, grnboost_scores)
        
      } else {
        # === k>1: Full analysis ===
        
        early_matrix <- list(earlyj(current_matrices, rowg = TRUE))
        
        # GENIE3
        # Late integration
        timing[["GENIE3_late"]] <- system.time({
          networks <- infer_networks(current_matrices, method = "GENIE3", nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(current_matrices, sym_wadj, n = 3, method = "GENIE3", nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        consensus_vote <- create_consensus(binary_adj, method = "vote")
        consensus_union <- create_consensus(binary_adj, method = "union")
        consensus_inet <- create_consensus(binary_adj, sym_wadj, method = "INet", threshold = 0.05, ncores = n_cores)
        
        scores <- pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))$Statistics %>%
          mutate(Method = "GENIE3", Predicted_Matrix = c("vote", "union", "inet"))
        scores <- add_community_metrics(scores, list(consensus_vote, consensus_union, consensus_inet))
        all_scores <- bind_rows(all_scores, scores)
        
        # Early integration
        timing[["GENIE3_early"]] <- system.time({
          networks <- infer_networks(early_matrix, method = "GENIE3", nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(early_matrix, sym_wadj, n = 2, method = "GENIE3", nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        scores <- pscores(adjm, binary_adj)$Statistics %>%
          mutate(Method = "GENIE3", Predicted_Matrix = "early")
        scores <- add_community_metrics(scores, binary_adj)
        all_scores <- bind_rows(all_scores, scores)
        
        # GRNBoost2
        # Late integration
        timing[["GRNBoost2_late"]] <- system.time({
          networks <- infer_networks(current_matrices, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(current_matrices, sym_wadj, n = 3, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        consensus_vote <- create_consensus(binary_adj, method = "vote")
        consensus_union <- create_consensus(binary_adj, method = "union")
        consensus_inet <- create_consensus(binary_adj, sym_wadj, method = "INet", threshold = 0.05, ncores = n_cores)
        
        scores <- pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))$Statistics %>%
          mutate(Method = "GRNBoost2", Predicted_Matrix = c("vote", "union", "inet"))
        scores <- add_community_metrics(scores, list(consensus_vote, consensus_union, consensus_inet))
        all_scores <- bind_rows(all_scores, scores)
        
        # Early integration
        timing[["GRNBoost2_early"]] <- system.time({
          networks <- infer_networks(early_matrix, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores)
        })
        wadj <- generate_adjacency(networks)
        sym_wadj <- symmetrize(wadj, weight_function = "mean")
        binary_adj <- cutoff_adjacency(early_matrix, sym_wadj, n = 2, method = "GRNBoost2", grnboost_modules = modules, nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
        
        scores <- pscores(adjm, binary_adj)$Statistics %>%
          mutate(Method = "GRNBoost2", Predicted_Matrix = "early")
        scores <- add_community_metrics(scores, binary_adj)
        all_scores <- bind_rows(all_scores, scores)
        
        # JRF
        jrf_scores <- tibble(
          Predicted_Matrix = c("vote", "union", "inet"),
          Method = "JRF",
          TP = 0, TN = 0, FP = 0, FN = 0,
          TPR = 0, FPR = 0, Precision = 0, F1 = 0, MCC = 0,
          VI = NA_real_, NMI = NA_real_, ARI = NA_real_,
          Modularity = NA_real_, Communities = NA_real_,
          Density = NA_real_, Transitivity = NA_real_
        )
        
        timing[["JRF"]] <- system.time({
          tryCatch({
            networks <- infer_networks(current_matrices, method = "JRF", nCores = n_cores)
            
            # Process JRF
            if (!is.null(networks) && length(networks) > 0 && !is.null(networks[[1]])) {
              importance_cols <- grep("importance", names(networks[[1]]), value = TRUE)
              
              if (length(importance_cols) >= k) {
                # Use first k importance columns
                jrf_list <- lapply(importance_cols[1:k], function(col) {
                  df <- networks[[1]][, c("gene1", "gene2", col)]
                  colnames(df)[3] <- col
                  df
                })
                
                wadj <- generate_adjacency(jrf_list)
                sym_wadj <- symmetrize(wadj, weight_function = "mean")
                binary_adj <- cutoff_adjacency(current_matrices, sym_wadj, n = 3, method = "JRF", nCores = n_cores, quantile_threshold = quantile_threshold, debug = FALSE)
                
                consensus_vote <- create_consensus(binary_adj, method = "vote")
                consensus_union <- create_consensus(binary_adj, method = "union")
                consensus_inet <- create_consensus(binary_adj, sym_wadj, method = "INet", threshold = 0.1, ncores = n_cores)
                
                scores <- pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))$Statistics %>%
                  mutate(Method = "JRF", Predicted_Matrix = c("vote", "union", "inet"))
                scores <- add_community_metrics(scores, list(consensus_vote, consensus_union, consensus_inet))
                
                jrf_scores <- scores
              }
            }
          }, error = function(e) {
            if (verbose) cat("JRF failed:", e$message, "\n")
          })
        })
        
        # ALWAYS add JRF (even if failed)
        all_scores <- bind_rows(all_scores, jrf_scores)
      }
      
      # Add timing and metadata
      timing_df <- tibble(
        Method_Integration = names(timing),
        Time_seconds = sapply(timing, function(x) x["elapsed"])
      ) %>%
        separate(Method_Integration, into = c("Method", "Integration"), sep = "_", fill = "right") %>%
        mutate(Integration = case_when(
          is.na(Integration) ~ ifelse(k == 1, "single", "joint"),
          TRUE ~ Integration
        ))
      
      all_scores <- all_scores %>%
        mutate(
          k = k,
          Run = run_id,
          Ratio = ratio,
          p = n_genes,
          Integration = case_when(
            k == 1 ~ "single",
            Method == "JRF" ~ "joint",
            Predicted_Matrix == "early" ~ "early",
            Predicted_Matrix %in% c("vote", "union", "inet") ~ "late",
            TRUE ~ "unknown"
          )
        ) %>%
        left_join(timing_df, by = c("Method", "Integration")) %>%
        mutate(Time_minutes = Time_seconds / 60)
      
      k_results[[run_id]] <- all_scores
      
      # Cleanup
      BiocParallel::bpstop(BiocParallel::bpparam())
      if (run_id < n_runs) Sys.sleep(1)
    }
    
    all_results[[as.character(k)]] <- bind_rows(k_results, .id = "Run")
    if (verbose) cat(glue("Completed k = {k}\n\n"))
  }
  
  # Combine and summarize
  final_df <- bind_rows(all_results, .id = "k")
  final_df$k <- as.numeric(final_df$k)
  
  summary_df <- final_df %>%
    group_by(k, Method, Predicted_Matrix) %>%
    summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), .groups = "drop") %>%
    arrange(k, Method, Predicted_Matrix)
  
  # Save results
  if (!is.null(output_file)) {
    write.table(final_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) cat("Results saved to:", output_file, "\n")
  }
  
  if (verbose) {
    cat("\n=== Study Completed ===\n")
    cat("Methods found in results:\n")
    print(table(final_df$Method, final_df$k))
  }
  
  return(list(
    summary = summary_df,
    detailed = final_df,
    parameters = list(
      max_k = max_k,
      n_runs = n_runs,
      n_genes = n_genes,
      n_cells = n_cells
    )
  ))
}

count_mat <- readRDS("../data/simdata/sim_n100p500k5.RDS")
results <- run_k_progressive_study(
   count_matrices_dir = count_mat,
   adjm_file = "../data/adjacency/adjm_p500.txt",
   max_k = 5, 
   n_runs = 10,
   genes_on_rows = TRUE,
   n_cores = 15,
   verbose = TRUE,
   output_file = "../data/results/sim_n100p500k5_summary.txt"
 )
