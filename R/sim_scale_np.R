#' Enhanced simulation function for network inference methods comparison
#' 
#' @param adjm_file Path to adjacency matrix file (.txt format)
#' @param count_matrices_file Path to count matrices file (.RDS format)
#' @param n_runs Number of simulation runs to perform (default: 10)
#' @param n_cores Number of CPU cores to use (default: 15) - NOTE: JRF doesn't use this parameter
#' @param quantile_threshold Quantile threshold for cutoff_adjacency (default: 0.99)
#' @param inet_threshold Threshold parameter for INet consensus (default: 0.05 for late, 0.1 for JRF)
#' @param verbose Print progress messages (default: TRUE)
#' 
#' @return Data frame with aggregated results across all runs
run_simulation_study <- function(adjm_file,
                                 count_matrices_file, 
                                 n_runs = 10,
                                 n_cores = 15,
                                 quantile_threshold = 0.99,
                                 inet_threshold_late = 0.05,
                                 inet_threshold_jrf = 0.1,
                                 verbose = TRUE) {
  
  # Helper function for timing with messages
  time_step <- function(step_name, expr, verbose = TRUE) {
    if (verbose) cat("  → Starting:", step_name, "...")
    start_time <- Sys.time()
    result <- expr
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    if (verbose) {
      if (elapsed < 60) {
        cat(" ✓ Completed in", round(elapsed, 2), "seconds\n")
      } else if (elapsed < 3600) {
        cat(" ✓ Completed in", round(elapsed/60, 2), "minutes\n")
      } else {
        cat(" ✓ Completed in", round(elapsed/3600, 2), "hours\n")
      }
    }
    
    return(list(result = result, time = elapsed))
  }
  
  # Load required libraries
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scGraphVerse)
  })
  
  # Initialize Python environment for GRNBoost2
  modules <- init_py(python_path = "/usr/bin/python3", required = TRUE)
  
  # Validate input files
  if (!file.exists(adjm_file)) {
    stop("Adjacency matrix file not found: ", adjm_file)
  }
  if (!file.exists(count_matrices_file)) {
    stop("Count matrices file not found: ", count_matrices_file)
  }
  
  # Load data once
  if (verbose) cat("Loading data files...\n")
  adjm_step <- time_step("Loading adjacency matrix", {
    adjm <- as.matrix(read.table(adjm_file))
    colnames(adjm) <- rownames(adjm)
    adjm
  }, verbose)
  adjm <- adjm_step$result
  
  count_step <- time_step("Loading count matrices", {
    readRDS(count_matrices_file)
  }, verbose)
  count_matrices <- count_step$result
  
  # Prepare early integration matrix once
  early_step <- time_step("Preparing early integration matrix", {
    list(earlyj(count_matrices, rowg = TRUE))
  }, verbose)
  early_matrix <- early_step$result
  
  # Pre-calculate some metrics
  n_genes <- nrow(adjm)
  n_cells <- ncol(count_matrices[[1]])
  ratio <- n_genes / n_cells
  
  if (verbose) {
    cat("Simulation parameters:\n")
    cat("  - Genes (p):", n_genes, "\n")
    cat("  - Cells (n):", n_cells, "\n")
    cat("  - Ratio (p/n):", round(ratio, 3), "\n")
    cat("  - Runs:", n_runs, "\n")
    cat("  - Cores:", n_cores, "(Note: JRF uses C implementation and doesn't use nCores)\n\n")
  }
  
  # Storage for all runs
  all_results <- list()
  
  # Run simulations
  for (run_id in 1:n_runs) {
    if (verbose) cat("=== Running simulation", run_id, "of", n_runs, "===\n")
    
    # Storage for timing and results
    timing_results <- list()
    method_scores <- list()
    
    # Helper function to add community metrics
    add_community_metrics <- function(scores_df, adj_comm, pred_communities, method_name) {
      if (verbose) cat("  → Starting: Computing community similarity metrics...")
      start_time <- Sys.time()
      topscore <- community_similarity(adj_comm, pred_communities)
      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
      if (verbose) {
        if (elapsed < 60) {
          cat(" ✓ Completed in", round(elapsed, 2), "seconds\n")
        } else {
          cat(" ✓ Completed in", round(elapsed/60, 2), "minutes\n")
        }
      }
      
      # Add topology measures
      topology_step <- time_step("Adding topology measures", {
        topscore$topology_measures %>%
          rownames_to_column("Predicted_Matrix") %>%
          mutate(
            Predicted_Matrix = case_when(
              Predicted_Matrix == "Predicted_1" ~ scores_df$Predicted_Matrix[1],
              Predicted_Matrix == "Predicted_2" ~ scores_df$Predicted_Matrix[2], 
              Predicted_Matrix == "Predicted_3" ~ scores_df$Predicted_Matrix[3],
              TRUE ~ Predicted_Matrix
            ),
            Method = method_name,
            Ratio = ratio,
            p = n_genes
          ) %>%
          full_join(scores_df, by = c("Predicted_Matrix", "Method", "Ratio", "p"))
      }, verbose)
      scores_df <- topology_step$result
      
      # Add community metrics
      community_step <- time_step("Adding community metrics", {
        topscore$community_metrics %>%
          rownames_to_column("Predicted_Matrix") %>%
          mutate(
            Predicted_Matrix = case_when(
              Predicted_Matrix == "Predicted_1" ~ scores_df$Predicted_Matrix[1],
              Predicted_Matrix == "Predicted_2" ~ scores_df$Predicted_Matrix[2],
              Predicted_Matrix == "Predicted_3" ~ scores_df$Predicted_Matrix[3], 
              TRUE ~ Predicted_Matrix
            ),
            Method = method_name,
            Ratio = ratio,
            p = n_genes
          ) %>%
          full_join(scores_df, by = c("Predicted_Matrix", "Method", "Ratio", "p"))
      }, verbose)
      scores_df <- community_step$result
      
      return(scores_df)
    }
    
    # Helper function for processing late integration (GENIE3/GRNBoost2)
    process_late_integration <- function(networks, method_name, inet_thresh = inet_threshold_late) {
      # Generate and symmetrize adjacency matrices
      wadj_step <- time_step("Generating adjacency matrices", {
        generate_adjacency(networks)
      }, verbose)
      wadj_list <- wadj_step$result
      
      sym_step <- time_step("Symmetrizing adjacency matrices", {
        symmetrize(wadj_list, weight_function = "mean")
      }, verbose)
      sym_wadj_list <- sym_step$result
      
      # Apply cutoff - WITH nCores for GENIE3/GRNBoost2
      cutoff_step <- time_step("Applying cutoff with shuffled networks", {
        cutoff_adjacency(
          count_matrices = count_matrices,
          weighted_adjm_list = sym_wadj_list,
          n = 3,
          method = method_name,
          nCores = n_cores,  # OK for GENIE3/GRNBoost2
          quantile_threshold = quantile_threshold,
          grnboost_modules = if(method_name == "GRNBoost2") modules else NULL,
          debug = FALSE
        )
      }, verbose)
      binary_adj_list <- cutoff_step$result
      
      # Create consensus matrices
      vote_step <- time_step("Creating vote consensus", {
        create_consensus(binary_adj_list, method = "vote")
      }, verbose)
      consensus_vote <- vote_step$result
      
      union_step <- time_step("Creating union consensus", {
        create_consensus(binary_adj_list, method = "union")
      }, verbose)
      consensus_union <- union_step$result
      
      inet_step <- time_step("Creating INet consensus", {
        create_consensus(
          adj_matrix_list = binary_adj_list,
          weighted_list = sym_wadj_list,
          method = "INet",
          threshold = inet_thresh,
          ncores = n_cores
        )
      }, verbose)
      consensus_inet <- inet_step$result
      
      # Calculate scores
      pscores_step <- time_step("Calculating performance scores", {
        pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))
      }, verbose)
      scores <- pscores_step$result
      
      scores_df <- scores$Statistics %>%
        mutate(
          Predicted_Matrix = c("vote", "union", "inet"),
          Method = method_name,
          Ratio = ratio,
          p = n_genes
        )
      
      # Community detection
      comm_step <- time_step("Detecting communities for consensus matrices", {
        list(
          community_path(consensus_vote, plot = FALSE, verbose = FALSE),
          community_path(consensus_union, plot = FALSE, verbose = FALSE), 
          community_path(consensus_inet, plot = FALSE, verbose = FALSE)
        )
      }, verbose)
      communities <- comm_step$result
      
      return(list(
        scores = scores_df,
        communities = communities
      ))
    }
    
    # Helper function for processing early integration (GENIE3/GRNBoost2)
    process_early_integration <- function(networks, method_name) {
      # Generate and symmetrize adjacency matrices
      wadj_step <- time_step("Generating adjacency matrices", {
        generate_adjacency(networks)
      }, verbose)
      wadj_list <- wadj_step$result
      
      sym_step <- time_step("Symmetrizing adjacency matrices", {
        symmetrize(wadj_list, weight_function = "mean")
      }, verbose)
      sym_wadj_list <- sym_step$result
      
      # Apply cutoff - WITH nCores for GENIE3/GRNBoost2
      cutoff_step <- time_step("Applying cutoff with shuffled networks", {
        cutoff_adjacency(
          count_matrices = early_matrix,
          weighted_adjm_list = sym_wadj_list,
          n = 2,
          method = method_name,
          nCores = n_cores,  # OK for GENIE3/GRNBoost2
          quantile_threshold = quantile_threshold,
          grnboost_modules = if(method_name == "GRNBoost2") modules else NULL,
          debug = FALSE
        )
      }, verbose)
      binary_adj_list <- cutoff_step$result
      
      # Calculate scores
      pscores_step <- time_step("Calculating performance scores", {
        pscores(adjm, binary_adj_list)
      }, verbose)
      scores <- pscores_step$result
      
      scores_df <- scores$Statistics %>%
        mutate(
          Predicted_Matrix = "early",
          Method = method_name,
          Ratio = ratio,
          p = n_genes
        )
      
      # Community detection
      comm_step <- time_step("Detecting communities for early integration", {
        community_path(binary_adj_list[[1]], plot = FALSE, verbose = FALSE)
      }, verbose)
      communities <- comm_step$result
      
      return(list(
        scores = scores_df,
        communities = list(communities)
      ))
    }
    
    # Ground truth community for comparison
    adj_comm <- community_path(adjm, plot = FALSE, verbose = FALSE)
    
    # === GENIE3 ===
    if (verbose) cat("Running GENIE3...\n")
    
    # Late integration - WITH nCores
    genie3_late_step <- time_step("GENIE3 network inference (late integration)", {
      infer_networks(count_matrices, method = "GENIE3", nCores = n_cores)
    }, verbose)
    genie3_late_networks <- genie3_late_step$result
    timing_results[["GENIE3_late"]] <- list(elapsed = genie3_late_step$time)
    
    if (verbose) cat("  Processing GENIE3 late integration results...\n")
    genie3_late_results <- process_late_integration(genie3_late_networks, "GENIE3")
    genie3_scores <- add_community_metrics(
      genie3_late_results$scores, 
      adj_comm, 
      genie3_late_results$communities,
      "GENIE3"
    )
    
    # Early integration - WITH nCores
    genie3_early_step <- time_step("GENIE3 network inference (early integration)", {
      infer_networks(early_matrix, method = "GENIE3", nCores = n_cores)
    }, verbose)
    genie3_early_networks <- genie3_early_step$result
    timing_results[["GENIE3_early"]] <- list(elapsed = genie3_early_step$time)
    
    if (verbose) cat("  Processing GENIE3 early integration results...\n")
    genie3_early_results <- process_early_integration(genie3_early_networks, "GENIE3")
    genie3_early_scores <- add_community_metrics(
      genie3_early_results$scores,
      adj_comm,
      genie3_early_results$communities,
      "GENIE3"
    )
    
    genie3_scores <- rbind(genie3_scores, genie3_early_scores)
    
    # === GRNBoost2 ===
    if (verbose) cat("Running GRNBoost2...\n")
    
    # Late integration - WITH nCores
    grnboost_late_step <- time_step("GRNBoost2 network inference (late integration)", {
      infer_networks(
        count_matrices, 
        method = "GRNBoost2",
        grnboost_modules = modules,
        nCores = n_cores
      )
    }, verbose)
    grnboost_late_networks <- grnboost_late_step$result
    timing_results[["GRNBoost2_late"]] <- list(elapsed = grnboost_late_step$time)
    
    if (verbose) cat("  Processing GRNBoost2 late integration results...\n")
    grnboost_late_results <- process_late_integration(grnboost_late_networks, "GRNBoost2")
    grnboost_scores <- add_community_metrics(
      grnboost_late_results$scores,
      adj_comm,
      grnboost_late_results$communities,
      "GRNBoost2"
    )
    
    # Early integration - WITH nCores
    grnboost_early_step <- time_step("GRNBoost2 network inference (early integration)", {
      infer_networks(
        early_matrix,
        method = "GRNBoost2", 
        grnboost_modules = modules,
        nCores = n_cores
      )
    }, verbose)
    grnboost_early_networks <- grnboost_early_step$result
    timing_results[["GRNBoost2_early"]] <- list(elapsed = grnboost_early_step$time)
    
    if (verbose) cat("  Processing GRNBoost2 early integration results...\n")
    grnboost_early_results <- process_early_integration(grnboost_early_networks, "GRNBoost2")
    grnboost_early_scores <- add_community_metrics(
      grnboost_early_results$scores,
      adj_comm,
      grnboost_early_results$communities,
      "GRNBoost2"
    )
    
    grnboost_scores <- rbind(grnboost_scores, grnboost_early_scores)
    
    # === JRF (Joint Random Forest) ===
    if (verbose) cat("Running JRF (using C implementation, no nCores parameter)...\n")
    
    # JRF network inference - WITHOUT nCores parameter
    jrf_inference_step <- time_step("JRF network inference (joint integration)", {
      infer_networks(count_matrices, method = "JRF")  # No nCores parameter!
    }, verbose)
    jrf_networks <- jrf_inference_step$result
    timing_results[["JRF_joint"]] <- list(elapsed = jrf_inference_step$time)
    
    if (verbose) cat("  Processing JRF joint integration results...\n")
    
    # Process JRF results (joint integration)
    wadj_step <- time_step("Generating adjacency matrices", {
      generate_adjacency(jrf_networks)
    }, verbose)
    wadj_list <- wadj_step$result
    
    sym_step <- time_step("Symmetrizing adjacency matrices", {
      symmetrize(wadj_list, weight_function = "mean")
    }, verbose)
    sym_wadj_list <- sym_step$result
    
    # Apply cutoff - WITHOUT nCores for JRF
    cutoff_step <- time_step("Applying cutoff with shuffled networks", {
      cutoff_adjacency(
        count_matrices = count_matrices,
        weighted_adjm_list = sym_wadj_list,
        n = 3,
        method = "JRF",
        # nCores = n_cores,  # REMOVED - JRF doesn't use this parameter
        quantile_threshold = quantile_threshold,
        debug = FALSE
      )
    }, verbose)
    binary_adj_list <- cutoff_step$result
    
    # Create consensus matrices
    vote_step <- time_step("Creating vote consensus", {
      create_consensus(binary_adj_list, method = "vote")
    }, verbose)
    consensus_vote <- vote_step$result
    
    union_step <- time_step("Creating union consensus", {
      create_consensus(binary_adj_list, method = "union")
    }, verbose)
    consensus_union <- union_step$result
    
    inet_step <- time_step("Creating INet consensus", {
      create_consensus(
        adj_matrix_list = binary_adj_list,
        weighted_list = sym_wadj_list,
        method = "INet", 
        threshold = inet_threshold_jrf,
        ncores = n_cores  # INet consensus can still use cores
      )
    }, verbose)
    consensus_inet <- inet_step$result
    
    # Calculate scores
    pscores_step <- time_step("Calculating performance scores", {
      pscores(adjm, list(consensus_vote, consensus_union, consensus_inet))
    }, verbose)
    scores <- pscores_step$result
    
    jrf_scores <- scores$Statistics %>%
      mutate(
        Predicted_Matrix = c("vote", "union", "inet"),
        Method = "JRF",
        Ratio = ratio,
        p = n_genes
      )
    
    # Community detection
    comm_step <- time_step("Detecting communities for consensus matrices", {
      list(
        community_path(consensus_vote, plot = FALSE, verbose = FALSE),
        community_path(consensus_union, plot = FALSE, verbose = FALSE),
        community_path(consensus_inet, plot = FALSE, verbose = FALSE)
      )
    }, verbose)
    jrf_communities <- comm_step$result
    
    jrf_scores <- add_community_metrics(jrf_scores, adj_comm, jrf_communities, "JRF")
    
    # === Combine timing data ===
    timing_df <- data.frame(
      Method = names(timing_results),
      Time_in_Seconds = sapply(timing_results, function(x) x[["elapsed"]])
    ) %>%
      mutate(
        Time_in_Minutes = Time_in_Seconds / 60,
        Time_in_Hours = Time_in_Seconds / 3600,
        Ratio = ratio,
        p = n_genes
      ) %>%
      separate(Method, into = c("Method", "Integration"), sep = "_", fill = "right") %>%
      mutate(Integration = case_when(
        Method == "JRF" ~ "joint",
        TRUE ~ Integration
      ))
    
    # === Combine all scores ===
    all_scores <- rbind(genie3_scores, grnboost_scores, jrf_scores)
    
    # Add timing information
    all_scores <- all_scores %>%
      mutate(
        Integration = case_when(
          Method == "JRF" ~ "joint",
          Predicted_Matrix == "early" ~ "early", 
          Predicted_Matrix %in% c("vote", "union", "inet") ~ "late",
          TRUE ~ Predicted_Matrix
        )
      ) %>%
      left_join(
        timing_df %>% select(Method, Integration, Time_in_Seconds, Time_in_Minutes, Time_in_Hours),
        by = c("Method", "Integration")
      ) %>%
      mutate(Run = run_id) %>%
      select(-Integration)
    
    all_results[[run_id]] <- all_scores
    
    # Clean up BiocParallel workers
    BiocParallel::bpstop(BiocParallel::bpparam())
    if (run_id < n_runs) Sys.sleep(2)
    
    if (verbose) cat("Completed run", run_id, "\n\n")
  }
  
  # === Aggregate results across all runs ===
  if (verbose) cat("Aggregating results across", n_runs, "runs...\n")
  
  combined_results <- bind_rows(all_results)
  
  # Calculate summary statistics
  summary_results <- combined_results %>%
    group_by(Method, Predicted_Matrix) %>%
    summarise(
      across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(Method, Predicted_Matrix)
  
  if (verbose) {
    cat("Simulation completed!\n")
    cat("Summary of key metrics (mean ± sd):\n")
    
    key_metrics <- summary_results %>%
      select(Method, Predicted_Matrix, 
             TPR_mean, TPR_sd, Precision_mean, Precision_sd, 
             F1_mean, F1_sd, Time_in_Minutes_mean) %>%
      mutate(
        TPR = paste0(round(TPR_mean, 3), " ± ", round(TPR_sd, 3)),
        Precision = paste0(round(Precision_mean, 3), " ± ", round(Precision_sd, 3)),
        F1 = paste0(round(F1_mean, 3), " ± ", round(F1_sd, 3)),
        Time_min = round(Time_in_Minutes_mean, 1)
      ) %>%
      select(Method, Predicted_Matrix, TPR, Precision, F1, Time_min)
    
    print(key_metrics)
  }
  
  return(list(
    summary = summary_results,
    detailed = combined_results,
    parameters = list(
      n_genes = n_genes,
      n_cells = n_cells,
      n_runs = n_runs,
      n_cores = n_cores,
      quantile_threshold = quantile_threshold,
      note = "JRF uses C implementation and doesn't use nCores parameter"
    )
  ))
}

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p100.txt",
  count_matrices_file = "../data/simdata/sim_n100p100.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n100p100_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n100p100_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)


results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p200.txt",
  count_matrices_file = "../data/simdata/sim_n100p200.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n100p200_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n100p200_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p500.txt",
  count_matrices_file = "../data/simdata/sim_n100p500.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n100p500_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n100p500_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p700.txt",
  count_matrices_file = "../data/simdata/sim_n100p700.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n100p700_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n100p700_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 500

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p100.txt",
  count_matrices_file = "../data/simdata/sim_n500p100.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n500p100_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n500p100_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)



results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p200.txt",
  count_matrices_file = "../data/simdata/sim_n500p200.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n500p200_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n500p200_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p500.txt",
  count_matrices_file = "../data/simdata/sim_n500p500.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n500p500_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n500p500_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_simulation_study(
  adjm_file = "../data/adjacency/adjm_p700.txt",
  count_matrices_file = "../data/simdata/sim_n500p700.RDS", 
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE
)

# Access summary results
summary_df <- results$summary

# Access detailed results
detailed_df <- results$detailed

# Save results
write.table(summary_df, "../data/results/sim_n500p700_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "../data/results/sim_n500p700_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

