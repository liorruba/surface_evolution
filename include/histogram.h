// This file creates and manipulates a histogram "object".
// A histogram is a 2-D array in which the first column is bins and the second column is counts.
typedef struct{
  double * bins;
  int * counts;
  int length;
} histogram;

// Create a histogram "object" with some range [minBin, maxBim] and size numOfBins.
histogram createHistogram(double minBin, double maxBin, int numOfBins);

// This function receives a value and an existing histogram and returns a histogram.
histogram addToHistogram(double *val, histogram hist);

// This function prints a histogram to file:
void printHistogram(histogram hist, char *path);
