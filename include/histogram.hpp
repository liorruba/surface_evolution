// This file creates and manipulates a histogram "object".
class Histogram {
  std::vector<double> bins;
  std::vector<int> counts;

  public:
    Histogram(double minBin, double maxBin, int numOfBins);
    void add(double val);
    void print(const char *path);
};
