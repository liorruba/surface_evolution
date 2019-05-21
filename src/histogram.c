// This file creates and manipulates a histogram "object".
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../include/histogram.h"
#include "../include/grid.h"
#include "../include/log.h"

// Create a histogram "object" with some range [minBin, maxBim] and size numOfBins.
histogram createHistogram(double minBin, double maxBin, int numOfBins){
  histogram hist;

  char logEntry[100];
  sprintf(logEntry, "Creating histogram with minimum bin value %f and maximum bin balue %f and size %d", minBin, maxBin, numOfBins);
  addLogEntry(logEntry);

  // Calulcate the exponent for the logvec function via change of base:
  double minBinExponent = log(minBin)/log(sqrt(2));
  double maxBinExponent = log(maxBin)/log(sqrt(2));

  // Create the bins and vals arrays:
  double * bins = logspace(minBinExponent, maxBinExponent, numOfBins, sqrt(2));

  // Set the histogram length:
  hist.length = numOfBins;
  // Allocate memory to the data array in histogram:
  hist.bins = (double *)malloc(sizeof(double) * (numOfBins));
  hist.counts = (int *)malloc(sizeof(int) * (numOfBins));

  for (int i = 0; i < numOfBins; i++){
    hist.bins[i] = bins[i];
    hist.counts[i] = 0;
  }

  return hist;
}

// This function recieves a value and an existing histogram and returns a histogram.
histogram addToHistogram(double *val, histogram hist){
  if ((*val < hist.bins[0]) || (*val > hist.bins[hist.length - 1])){
    // TODO: INCREASE HISTOGRAM SIZE.
  }
  else {
    // Iterate over histogram:
    for (int i = 1; i < hist.length; i++){
      // If in bin:
      if (*val < hist.bins[i]){
        hist.counts[i-1]++;
        break;
      }
    }
  }
  return hist;
}

// This function prints a histogram to file:
void printHistogram(histogram hist, char *path){
  FILE *histogramFile;

  histogramFile = fopen(path, "w+");
  if (histogramFile == NULL){
    char logEntry[100];
    sprintf(logEntry, "Cannot create histogram file in path %s.",path);
    addLogEntry(logEntry);
  }
  else {
    char logEntry[100];
    sprintf(logEntry, "Histogram file successfully created in path %s.",path);
    addLogEntry(logEntry);
  }

  for (int i = 0; i < hist.length; i++){
    fprintf(histogramFile, "%f\t",hist.bins[i]);
    fprintf(histogramFile, "%d\n",hist.counts[i]);
  }
  fclose(histogramFile);
}
