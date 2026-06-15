#pragma once
#include "correlated_normal.hpp"
#include <vector>

// dS = mu*S*dt + sigma*S*dW
struct AssetParams {
    double S0;
    double mu;
    double sigma;
};

// Returns a matrix of shape [numSteps + 1][numAssets], including the inital price at t=0
inline Matrix simulatePath(const std::vector<AssetParams>& assets, CorrelatedNormalGenerator& gen, int numSteps, double dt) {
    size_t n=assets.size();
    Matrix path(numSteps+1,std::vector<double> (n));

    for (size_t i = 0; i < n; ++i) path[0][i]=assets[i].S0;

    for (int t = 1; t <= numSteps; ++t) {
        std::vector<double> shock = gen.next();

        for (size_t i = 0; i < n; ++i) {
            double mu = assets[i].mu;
            double sigma = assets[i].sigma;
            double prevS = path[t - 1][i];

            double drift = (mu - 0.5*sigma*sigma)*dt;
            double diffusion = sigma * std::sqrt(dt) * shock[i];

            path[t][i] = prevS * std::exp(drift + diffusion);
        }
    }
    return path;
}