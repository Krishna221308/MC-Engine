# MC-Engine: Monte Carlo Correlated Path Generator

This repository contains the foundational C++ code for simulating **correlated Geometric Brownian Motion (GBM)** paths.

## Overview

In multi-asset options pricing (like basket options), assets do not move independently; they are correlated. To simulate this correctly in a Monte Carlo engine, we must:
1. Define a target correlation matrix ($\Sigma$).
2. Decompose this matrix using the **Cholesky Decomposition** to get a lower triangular matrix ($L$) such that $L L^T = \Sigma$.
3. Generate independent standard normal random variables ($Z$).
4. Multiply them by the Cholesky factor to get correlated shocks ($X = L \cdot Z$).
5. Feed these correlated shocks into the discrete **Geometric Brownian Motion (GBM)** formula to generate price paths over time.

This project implements that exact pipeline from scratch in C++.

## Files Architecture

*   **`cholesky.hpp`**: Implements the Cholesky decomposition algorithm. It takes a symmetric, positive-definite correlation matrix and returns the lower triangular factor $L$.
*   **`correlated_normal.hpp`**: A random number generator class that uses the standard `<random>` library to draw independent normals, and then multiplies them by the Cholesky matrix $L$ to produce correlated random variables.
*   **`gbm.hpp`**: The simulation engine. It takes the initial asset parameters (Price, Drift $\mu$, Volatility $\sigma$) and generates complete simulated price paths across multiple time steps using the correlated shocks.
*   **`main.cpp`**: The orchestrator. It sets up the parameters, runs 20,000 Monte Carlo paths, extracts the first-step returns, calculates the **empirical correlation matrix**, and compares it to the target correlation matrix to verify mathematical correctness.

## Compilation and Execution

You can compile this code using any modern C++ compiler (C++11 or later is required for `<random>`):

```bash
# Compile with GCC (Optimized for speed)
g++ -O3 main.cpp -o mc_engine

# Run the executable
./mc_engine
```

## Expected Output

When running the executable, it will print:
1. The target correlation matrix you explicitly defined.
2. The Cholesky lower triangular matrix.
3. An empirical correlation matrix derived purely from the statistics of 20,000 independent simulated paths.

The empirical matrix should closely match the target matrix (typically within ~0.01 to 0.02 due to Monte Carlo variance), proving the underlying math works perfectly.