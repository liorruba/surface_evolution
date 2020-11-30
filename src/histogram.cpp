// This file creates and manipulates a histogram "object".
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include "../include/histogram.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"


// Class constructor
Histogram::Histogram(double minBin, double maxBin, int numOfBins){
        char logEntry[100];
        sprintf(logEntry, "Creating histogram with minimum bin value %f and maximum bin value %f and size %d", minBin, maxBin, numOfBins);
        addLogEntry(logEntry);

        // calculate the exponent for the logvec function via change of base:
        double minBinExponent = log(minBin)/log(sqrt(2.0));
        double maxBinExponent = log(maxBin)/log(sqrt(2.0));

        // Create the bins and vals arrays:
        bins = logspace(minBinExponent, maxBinExponent, numOfBins, sqrt(2.0));

        // Initialize the counts vector to zero:
        std::vector<int> _counts(numOfBins, 0);
        counts = _counts;
}

// Add value to histogram
void Histogram::add(double val){
        if ((val < bins.front()) || (val > bins.back())) {
                // TODO: INCREASE HISTOGRAM SIZE.
        }
        else {
                // Iterate over histogram:
                for (size_t i = 1; i < bins.size(); i++) {
                        // If in bin:
                        if (val < bins[i]) {
                                counts[i-1]++;
                                break;
                        }
                }
        }
}

// Print histogram
void Histogram::print(const char *path){
        std::ofstream histogramFile;

        histogramFile.open(path, std::ios_base::out);

        if(histogramFile) {
                char logEntry[100];
                sprintf(logEntry, "Histogram file successfully created in path %s.",path);
                addLogEntry(logEntry);

                for (size_t i = 0; i < bins.size(); i++) {
                        histogramFile << bins[i] << "\t" << counts[i] << std::endl;
                }
        }

        else {
                char logEntry[100];
                sprintf(logEntry, "Cannot create histogram file in path %s.",path);
                addLogEntry(logEntry);
        }
}
