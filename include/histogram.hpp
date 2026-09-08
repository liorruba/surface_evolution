#pragma once
#include <vector>

// Log-binned histogram: bins[k] <= value < bins[k+1] is counted in counts[k]; a fixed number of bins per decade.
class Histogram {
  std::vector<double> bins;
  std::vector<long> counts;
  double logMin;
  double logStep;

  public:
    Histogram(double minBin, double maxBin, int binsPerDecade);
    void add(double val);
    void print(const char *path);
};
