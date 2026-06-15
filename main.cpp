#include "cholesky.hpp"
#include "correlated_normal.hpp"
#include "gbm.hpp"
#include <iostream>
#include <vector>

Matrix empiricalCorrelation(const std::vector<std::vector<double>>& returns) {
    size_t n = returns.size();
    size_t m = returns[0].size();

    // Means
    std::vector<double> mean(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (double r : returns[i]) mean[i]+=r;
        mean[i] /= static_cast<double>(m);
    }

    // Std
    std::vector<double> stdev(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        double sumSq = 0.0;
        for (double r : returns[i]) sumSq += (r - mean[i])*(r - mean[i]);
        stdev[i] = std::sqrt(sumSq / static_cast<double>(m - 1));
    }

    // Correlation Matrix
    Matrix corr(n, std::vector<double> (n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double cov = 0.0;
            for (size_t k = 0; k < m; ++k) {
                cov += (returns[i][k] - mean[i]) * (returns[j][k] - mean[j]);
            }
            cov /= static_cast<double> (m-1);
            corr[i][j] = cov / (stdev[i] * stdev[j]);
        }
    }
    return corr; 
}

void printMatrix(const Matrix& M, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : M) {
        for (double v : row) {
            printf("%8.4f ", v);
        }
        std::cout << '\n';
    }
    std::cout<< '\n';
}

int main() {
    Matrix sigma = {
        {1.00, 0.50, 0.30},
        {0.50, 1.00, 0.20},
        {0.30, 0.20, 1.00},
    };
    printMatrix(sigma,"Target correlation matrix (Sigma)");

    // Cholesky
    Matrix L = choleskyDecompose(sigma);
    printMatrix(L, "Cholesky factor (L)");

    // Correlated normal generator
    CorrelatedNormalGenerator gen(L, /*seed=*/42);

    // Defining asset params
    std::vector<AssetParams> assets = {
        {100.0, 0.05, 0.20},
        {50.0,  0.03, 0.25},
        {200.0, 0.07, 0.15},
    };

    int numSteps = 50;
    double dt = 1.0 / 252.0;
    int numPaths = 20000;

    // Monte Carlo paths, collect log returns of the first step
    size_t n = assets.size();
    std::vector<std::vector<double>> firstStepReturns(n);
    for (size_t i = 0; i < n; ++i) firstStepReturns[i].reserve(numPaths);

    for (int p = 0; p < numPaths; ++p) {
        Matrix path = simulatePath(assets, gen, numSteps, dt);

        for (size_t i = 0; i < n; ++i) {
            double logreturn = std::log(path[1][i] / path[0][i]);
            firstStepReturns[i].push_back(logreturn);
        }
    }

    // Verification
    Matrix empCorr = empiricalCorrelation(firstStepReturns);
    printMatrix(empCorr, "Empirical Correlation Matrix (from simulation)");

    std::cout << "Compare this to the target Sigma above.\n";
    std::cout << "With numPaths = " << numPaths << ", values should be close (within ~0.01-0.02).\n";
 
    return 0;
}