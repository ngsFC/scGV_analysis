# scGraphVerse Benchmarking and Case Study Analysis

This repository contains a comprehensive analysis using the scGraphVerse package. The analysis includes data generation, simulation studies, and case study using PBMC data.

## Repository Structure

```
scGV_analysis/
├── R/                          # R scripts and analysis notebooks
│   ├── PBMC.gtruth.Rmd        # Data generation and ground truth creation
│   ├── case_study_JRF.Rmd     # Real data case study analysis
│   ├── sim_scale_np.R         # Simulation study: scaling by nodes (n) and genes (p)
│   └── sim_scale_K.R          # Simulation study: scaling by datasets (K)
├── data/                       # Data storage
│   ├── adjacency/             # Ground truth adjacency matrices
│   ├── simdata/               # Simulated count matrices
│   ├── results/               # Simulation results
│   └── PBMC.top600.RDS        # Processed PBMC data
└── output/                     # Analysis outputs
    ├── plots/                 # Generated plots and figures
    └── tables/                # Results tables and summaries
```

## Workflow Overview

The analysis follows a three-stage workflow:

### 1. 📊 Data Generation (`PBMC.gtruth.Rmd`)
**Purpose**: Generate ground truth networks and simulated data from real PBMC single-cell RNA-seq data.

**What it does**:
- Downloads PBMC data from Single Cell Atlas
- Creates ground truth adjacency matrices using STRING database at different scales:
  - **p100**: ~100 genes network
  - **p200**: ~200 genes network  
  - **p500**: ~500 genes network
  - **p700**: ~700 genes network
- Generates simulated count matrices using ZINB model with different parameters:
  - Various sample sizes (n=100, n=500)
  - Different cluster numbers (K=3, K=5)
- Saves processed PBMC data for case studies

**Key outputs**:
- Ground truth adjacency matrices (`data/adjacency/adjm_p*.txt`)
- Simulated count matrices (`data/simdata/sim_n*p*.RDS`)
- Network visualization plots (`output/gtruth_p*.png`)
- Processed PBMC data (`data/PBMC.top600.RDS`)

### 2. 🧪 Simulation Studies

#### 2.1 Node/Gene Scaling Study (`sim_scale_np.R`)
**Purpose**: Evaluate network inference methods across different network sizes and sample sizes.

**What it does**:
- Tests multiple network inference methods:
  - **JRF (Joint Random Forest)**: Main method of interest
  - **GENIE3**: Random Forest
  - **GRNBoost2**: Gradient boosting method
- Evaluates performance across:
  - Different network sizes (p=100, 200, 500, 700 genes)
  - Different sample sizes (n=100, 500 cells)
- Calculates performance metrics:
  - Precision, Recall, F1-score, MCC
- Runs multiple iterations (default: 10) for statistical robustness

**Key outputs**:
- Summary results (`data/results/sim_n*p*_summary.txt`)
- Detailed results (`data/results/sim_n*p*_detailed.txt`)

#### 2.2 Multiple dataset Scaling Study (`sim_scale_K.R`)
**Purpose**: Evaluate how network inference performance increasing the number of datasets.

**What it does**:
- Uses the K=5 simulated data
- Progressively tests inference on subsets: K=1,2,3,4,5
- Uses the same performance metrics as the scaling study

**Key outputs**:
- K-progressive results (`data/results/sim_n*p*k*_summary.txt`)

### 3. 🔬 Case Study Analysis (`case_study_JRF.Rmd`)
**Purpose**: Apply the best method (JRF) to real PBMC data.

**What it does**:
- Loads processed PBMC data (B cells from multiple donors)
- Applies JRF method to infer gene regulatory networks
- Performs community detection to identify functional modules
- Conducts pathway enrichment analysis using KEGG database
- Validates inferred edges against:
  - STRING protein-protein interaction database
  - Literature mining via PubMed
  
**Key outputs**:
- Network community plots (`output/plots/network_communities.png`)
- Pathway enrichment visualizations (`output/plots/pathway_enrichment_*.png`)
- Summary tables (`output/tables/*.csv`)
- Final results summary (`output/results_cases.txt`)

### Customization

**Adjusting computational resources**:
- Modify `n_cores` parameter in scripts (default: 15)
- Reduce `n_runs` for faster execution (default: 10)

**Testing different parameters**:
- Network inference thresholds in simulation functions
- Community detection parameters in case study
- Pathway enrichment significance cutoffs

**Support**:
- Check scGraphVerse documentation
- File issues on the package GitHub repository

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
