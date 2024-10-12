#  Random Initialize

## Overview

These functions are for initializing parameters for a
Mixed Markov Model. The functions include:

-   **InitializeClusters**: Initializes cluster assignments randomly.
-   **filterSequencesByCluster**: Extracts sequences for each cluster.
-   **InitialStateProb**: Calculates initial state probabilities.
-   **InitialTransitionMatrix**: Calculates initial transition matrices.
-   **InitializeParameters**: Initializes all parameters for the HMM.

## Functions

### InitializeClusters

**Description:** Initializes cluster assignments randomly.

**Parameters:** - `sequences`: List of sequences. - `nClusters`: Number
of clusters.

**Returns:** `IntegerVector` of cluster assignments.

### filterSequencesByCluster

**Description:** Extracts sequences for each cluster.

**Parameters:** - `sequences`: List of sequences. - `clusters`:
`IntegerVector` of cluster assignments. - `k`: Cluster index.

**Returns:** List of sequences for the specified cluster.

### InitialStateProb

**Description:** Calculates initial state probabilities.

**Parameters:** - `sequences`: List of sequences. - `categories`: Number
of categories.

**Returns:** `NumericVector` of initial state probabilities.

### InitialTransitionMatrix

**Description:** Calculates initial transition matrices.

**Parameters:** - `sequences`: List of sequences. - `nCategories`:
Number of categories.

**Returns:** `NumericMatrix` of initial transition probabilities.

### InitializeParameters

**Description:** Initializes all parameters for the HMM.

**Parameters:** - `sequences`: List of sequences.

## Usage

To use this package, follow these steps:

**Source the Rcpp code**:
`r    library(Rcpp)    sourceCpp("02_Initialize.cpp")`

### Example:

    # Example data
    sequences <- list(...)  # Replace with actual sequences
    nClusters <- 3  # Replace with actual number of clusters
    nCategories <- 4  # Replace with actual number of categories

    # Initialize parameters
    params <- InitializeParameters(sequences, nClusters, nCategories)

    # Access elements from the returned list
    assigned_clusters <- params$assigned_clusters
    mixing_prop <- params$mixing_prop
    initial_values_alpha <- params$initial_values_alpha
    transition_matrices <- params$transition_matrices
    grouped_sequences_by_clusters <- params$grouped_sequences_by_clusters
