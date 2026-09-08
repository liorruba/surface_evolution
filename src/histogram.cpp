// Log-binned histogram with a fixed number of bins per decade.
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/histogram.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"

// Bins from minBin to at least maxBin, binsPerDecade log-uniform bins per factor of ten.
Histogram::Histogram(double minBin, double maxBin, int binsPerDecade) {
        binsPerDecade = std::max(1, binsPerDecade);
        minBin = std::max(minBin, 1e-12);
        maxBin = std::max(maxBin, minBin * 1.0001);
        const int numOfBins = (int) std::ceil(std::log10(maxBin / minBin) * binsPerDecade) + 1;
        logMin = std::log10(minBin);
        logStep = 1.0 / binsPerDecade;
        bins.resize(numOfBins);
        for (int i = 0; i < numOfBins; i++)
                bins[i] = std::pow(10.0, logMin + i * logStep);
        counts.assign(numOfBins, 0);

        char logEntry[160];
        snprintf(logEntry, sizeof(logEntry), "Creating histogram from %g to %g with %d bins per decade (%d bins).", minBin, bins.back(), binsPerDecade, numOfBins - 1);
        addLogEntry(logEntry, false);
}

// Add a value (values outside the range are not counted)
void Histogram::add(double val) {
        if (!(val >= bins.front()) || val >= bins.back())
                return;
        const int k = (int) std::floor((std::log10(val) - logMin) / logStep);
        if (k >= 0 && k + 1 < (int) bins.size())
                counts[k]++;
}

// Print: one line per bin edge with the count of the bin starting at that edge (the last count is 0)
void Histogram::print(const char *path) {
        std::ofstream histogramFile;
        histogramFile.open(path, std::ios_base::out);

        if (histogramFile) {
                char logEntry[160];
                snprintf(logEntry, sizeof(logEntry), "Histogram file successfully created in path %s.", path);
                addLogEntry(logEntry, false);

                for (size_t i = 0; i < bins.size(); i++) {
                        histogramFile << bins[i] << "\t" << counts[i] << std::endl;
                }
        }
        else {
                char logEntry[160];
                snprintf(logEntry, sizeof(logEntry), "ERROR: could not create histogram file %s.", path);
                addLogEntry(logEntry, true);
        }
}
