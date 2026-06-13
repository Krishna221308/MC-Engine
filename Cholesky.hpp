#pragma once
#include <vector>
#include <cmath>
#include <stdexcept> 

using Matrix = std::vector<std::vector<double>>;

// Computes L using cholesky decomposition
// Sigma must be a symmetric positive-definite NxN matrix.

inline Matrix choleskyDecompose(const Matrix& sigma) {
    size_t n = sigma.size();
    Matrix L(n, std::vector<double> (n,0.0));
    for (size_t i=0;i<n;++i) {
        for (size_t j=0;j<=i;++j) {
            double sum=0.0;
            if (j==i) {
                for (size_t k=0;k<j;++k) {
                    sum += L[j][k]*L[j][k]
                }
                double val = sigma[i][i]-sum;
                if (val < 0.0) {
                    throw std::runtime_error(
                        "Matrix is not positive-definite (non positiive pivot encountered).")
                }
                L[i][j]=sqrt(val);
            } else {
                for (size_t k = 0; k < j; ++k) {
                    sum+=(L[i][k] * L[j][k]);
                }
                L[i][j]=(sigma[i][j]-sum)/L[j][j];
            }
        }
    }
    return L;
}