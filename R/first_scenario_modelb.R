#' Enhanced simulation function for binary network inference methods comparison
#'
#' @param adjm_file Path to adjacency matrix file (.txt format)
#' @param count_matrices_file Path to count matrices file (.RDS format)
#' @param n_runs Number of simulation runs to perform (default: 10)
#' @param n_cores Number of CPU cores to use (default: 15)
#' @param inet_threshold Threshold parameter for INet consensus (default: 0.05)
#' @param methods Vector of methods to run. Options: "PCzinb", "ZILGM" (default: c("PCzinb", "ZILGM"))
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Data frame with aggregated results across all runs
run_binary_simulation_study <- function(adjm_file,
                                       count_matrices_file,
                                       n_runs = 10,
                                       n_cores = 15,
                                       inet_threshold = 0.05,
                                       methods = c("PCzinb", "ZILGM"),
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

  # Validate input files
  if (!file.exists(adjm_file)) {
    stop("Adjacency matrix file not found: ", adjm_file)
  }
  if (!file.exists(count_matrices_file)) {
    stop("Count matrices file not found: ", count_matrices_file)
  }

  # Validate methods parameter
  available_methods <- c("PCzinb", "ZILGM")
  if (!all(methods %in% available_methods)) {
    invalid_methods <- methods[!methods %in% available_methods]
    stop("Invalid methods specified: ", paste(invalid_methods, collapse = ", "),
         ". Available methods are: ", paste(available_methods, collapse = ", "))
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
  count_matrices_list <- count_step$result

  # Convert to MultiAssayExperiment for compatibility with new API
  mae_step <- time_step("Converting to MultiAssayExperiment", {
    create_mae(count_matrices_list)
  }, verbose)
  count_matrices <- mae_step$result

  # Prepare early integration matrix once (earlyj returns a MAE with merged data)
  early_step <- time_step("Preparing early integration matrix", {
    earlyj(count_matrices, rowg = TRUE)
  }, verbose)
  early_matrix <- early_step$result

  # Pre-calculate some metrics
  n_genes <- nrow(adjm)
  # Extract first experiment from MultiAssayExperiment to get dimensions
  n_cells <- ncol(MultiAssayExperiment::experiments(count_matrices)[[1]])
  ratio <- n_genes / n_cells

  if (verbose) {
    cat("Simulation parameters:\n")
    cat("  - Genes (p):", n_genes, "\n")
    cat("  - Cells (n):", n_cells, "\n")
    cat("  - Ratio (p/n):", round(ratio, 3), "\n")
    cat("  - Runs:", n_runs, "\n")
    cat("  - Cores:", n_cores, "\n\n")
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

      # Try to compute community similarity with error handling
      topscore <- tryCatch({
        community_similarity(adj_comm, pred_communities, plot = FALSE)
      }, error = function(e) {
        if (verbose) cat(" ⚠ Error in community comparison:", e$message, "\n")
        # Return minimal structure to continue processing
        list(
          community_metrics = data.frame(
            VI = rep(NA, length(pred_communities)),
            NMI = rep(NA, length(pred_communities)),
            ARI = rep(NA, length(pred_communities)),
            row.names = paste0("Predicted_", 1:length(pred_communities))
          ),
          topology_measures = data.frame(
            Modularity = rep(NA, length(pred_communities)),
            n_communities = rep(NA, length(pred_communities)),
            Density = rep(NA, length(pred_communities)),
            Transitivity = rep(NA, length(pred_communities)),
            row.names = paste0("Predicted_", 1:length(pred_communities))
          )
        )
      })

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

    # Helper function for processing binary methods (late integration)
    process_binary_late_integration <- function(binary_se, method_name, inet_thresh = inet_threshold) {
      # Binary methods (PCzinb, ZILGM) return SummarizedExperiment with adjacency matrices as assays
      # No need to call generate_adjacency()

      # Create consensus matrices directly from binary networks (extract matrix from SE)
      vote_step <- time_step("Creating vote consensus", {
        consensus_se <- create_consensus(binary_se, method = "vote")
        SummarizedExperiment::assays(consensus_se)[[1]]
      }, verbose)
      consensus_vote <- vote_step$result

      union_step <- time_step("Creating union consensus", {
        consensus_se <- create_consensus(binary_se, method = "union")
        SummarizedExperiment::assays(consensus_se)[[1]]
      }, verbose)
      consensus_union <- union_step$result

      # Skip INet consensus for binary methods since they don't provide weighted matrices
      if (verbose) cat("  → Skipping: INet consensus (requires weighted matrices, not available for binary methods)\n")

      # Calculate scores
      pscores_step <- time_step("Calculating performance scores", {
        pscores(adjm, list(consensus_vote, consensus_union))
      }, verbose)
      scores <- pscores_step$result

      scores_df <- scores$Statistics %>%
        mutate(
          Predicted_Matrix = c("vote", "union"),
          Method = method_name,
          Ratio = ratio,
          p = n_genes
        )

      # Community detection
      comm_step <- time_step("Detecting communities for consensus matrices", {
        list(
          community_path(consensus_vote, plot = FALSE, verbose = FALSE),
          community_path(consensus_union, plot = FALSE, verbose = FALSE)
        )
      }, verbose)
      communities <- comm_step$result

      return(list(
        scores = scores_df,
        communities = communities
      ))
    }

    # Helper function for processing binary methods (early integration)
    process_binary_early_integration <- function(binary_se, method_name) {
      # Binary methods (PCzinb, ZILGM) return SummarizedExperiment with adjacency matrices as assays
      # No need to call generate_adjacency()

      # Extract first assay (early integration produces single matrix)
      binary_matrix <- SummarizedExperiment::assays(binary_se)[[1]]

      # Calculate scores directly
      pscores_step <- time_step("Calculating performance scores", {
        pscores(adjm, list(binary_matrix))
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
        community_path(binary_matrix, plot = FALSE, verbose = FALSE)
      }, verbose)
      communities <- comm_step$result

      return(list(
        scores = scores_df,
        communities = list(communities)
      ))
    }

    # Ground truth community for comparison
    adj_comm <- community_path(adjm, plot = FALSE, verbose = FALSE)

    # Initialize scores storage
    method_scores_list <- list()

    # === PCzinb ===
    if ("PCzinb" %in% methods) {
      if (verbose) cat("Running PCzinb...\n")

      # Late integration
      pczinb_late_step <- time_step("PCzinb network inference (late integration)", {
        infer_networks(count_matrices, method = "PCzinb", nCores = n_cores)
      }, verbose)
      pczinb_late_networks <- pczinb_late_step$result
      timing_results[["PCzinb_late"]] <- list(elapsed = pczinb_late_step$time)

      if (verbose) cat("  Processing PCzinb late integration results...\n")
      pczinb_late_results <- process_binary_late_integration(pczinb_late_networks, "PCzinb")
      pczinb_scores <- add_community_metrics(
        pczinb_late_results$scores,
        adj_comm,
        pczinb_late_results$communities,
        "PCzinb"
      )

      # Early integration
      pczinb_early_step <- time_step("PCzinb network inference (early integration)", {
        infer_networks(early_matrix, method = "PCzinb", nCores = n_cores)
      }, verbose)
      pczinb_early_networks <- pczinb_early_step$result
      timing_results[["PCzinb_early"]] <- list(elapsed = pczinb_early_step$time)

      if (verbose) cat("  Processing PCzinb early integration results...\n")
      pczinb_early_results <- process_binary_early_integration(pczinb_early_networks, "PCzinb")
      pczinb_early_scores <- add_community_metrics(
        pczinb_early_results$scores,
        adj_comm,
        pczinb_early_results$communities,
        "PCzinb"
      )

      pczinb_scores <- rbind(pczinb_scores, pczinb_early_scores)
      method_scores_list[["PCzinb"]] <- pczinb_scores
    }

    # === ZILGM ===
    if ("ZILGM" %in% methods) {
      if (verbose) cat("Running ZILGM...\n")

      # Late integration
      zilgm_late_step <- time_step("ZILGM network inference (late integration)", {
        infer_networks(count_matrices, method = "ZILGM", nCores = n_cores,
                       zilgm_params = list(family = "NBII", update_type = "IRLS",
                                           do_boot = TRUE, boot_num = 5, sym = "OR", nCores = 15, nlambda = 15))
      }, verbose)
      zilgm_late_networks <- zilgm_late_step$result
      timing_results[["ZILGM_late"]] <- list(elapsed = zilgm_late_step$time)

      if (verbose) cat("  Processing ZILGM late integration results...\n")
      zilgm_late_results <- process_binary_late_integration(zilgm_late_networks, "ZILGM")
      zilgm_scores <- add_community_metrics(
        zilgm_late_results$scores,
        adj_comm,
        zilgm_late_results$communities,
        "ZILGM"
      )

      # Early integration
      zilgm_early_step <- time_step("ZILGM network inference (early integration)", {
        infer_networks(early_matrix, method = "ZILGM", nCores = n_cores,
                       zilgm_params = list(family = "NBII", update_type = "IRLS",
                                           do_boot = TRUE, boot_num = 5, sym = "OR", nCores = 15, nlambda = 15))
      }, verbose)
      zilgm_early_networks <- zilgm_early_step$result
      timing_results[["ZILGM_early"]] <- list(elapsed = zilgm_early_step$time)

      if (verbose) cat("  Processing ZILGM early integration results...\n")
      zilgm_early_results <- process_binary_early_integration(zilgm_early_networks, "ZILGM")
      zilgm_early_scores <- add_community_metrics(
        zilgm_early_results$scores,
        adj_comm,
        zilgm_early_results$communities,
        "ZILGM"
      )

      zilgm_scores <- rbind(zilgm_scores, zilgm_early_scores)
      method_scores_list[["ZILGM"]] <- zilgm_scores
    }

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
      separate(Method, into = c("Method", "Integration"), sep = "_", fill = "right")

    # === Combine all scores ===
    all_scores <- do.call(rbind, method_scores_list)

    # Add timing information
    all_scores <- all_scores %>%
      mutate(
        Integration = case_when(
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
      inet_threshold = inet_threshold,
      note = "PCzinb and ZILGM return binary matrices directly"
    )
  ))
}

# Run simulations for different scenarios

# n=100 scenarios
results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p100.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p100.RDS",
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE,
  methods = "ZILGM"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p100_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p100_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p100.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p100.RDS",
  n_runs = 2,
  n_cores = 15,
  verbose = TRUE,
  methods = "PCzinb"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p100_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p100_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)





results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p200.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p200.RDS",
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE,
  methods = "ZILGM"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p200_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p200_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p200.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p200.RDS",
  n_runs = 10,
  n_cores = 15,
  verbose = TRUE,
  methods = "PCzinb"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p200_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p200_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)



results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p500.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p500.RDS",
  n_runs = 5,
  n_cores = 15,
  verbose = TRUE,
  methods = "ZILGM"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p500_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_zilgm_n100p500_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

results <- run_binary_simulation_study(
  adjm_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/adjacency/adjm_p500.txt",
  count_matrices_file = "/home/francescoc/Desktop/network/graphwork/analysis/simulation/simdata/sim_n100p500.RDS",
  n_runs = 5,
  n_cores = 15,
  verbose = TRUE,
  methods = "PCzinb"
)

# Save results
summary_df <- results$summary
detailed_df <- results$detailed
write.table(summary_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p500_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detailed_df, "/home/francescoc/Desktop/network/graphwork/analysis/simulation/results/sim_pczinb_n100p500_detailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)