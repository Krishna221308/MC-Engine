#pragma once
#include "cholesky.hpp"
#include <random>

// Generates correlated standard-normal draws using a precomputed Cholesky factor L.
class CorrelatedNormalGenerator{
public:
    CorrelatedNormalGenerator(const Matrix& L, unsigned int seed) 
        : L_(L), n_(L.size()), rng_(seed), dist_(0.0,1.0) {}
    
    // Return L*Z where Z is a vector of independent N(0,1) draws.
    std::vector<double> next() {
        std::vector<double> z(n_);
        for (size_t i = 0; i < n_; ++i) {
            z[i] = dist_(rng_);
        }

        std::vector<double> x(n_,0.00);
        for (size_t i = 0; i < n_; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                x[i]+=L_[i][j]*z[j];
            }
        }
        return x;
    }           

private:
    Matrix L_;
    size_t n_;
    std::mt19937 rng_;
    std::normal_distribution<double> dist_;
};