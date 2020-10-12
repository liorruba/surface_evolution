// This file contains various functions to help manage the code.
#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cctype>
#include<ctime>
#include<cinttypes>
#include<vector>
#include<fstream>
#include "../include/utility.hpp"

// Get current time
uint64_t get_time()
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e9 + t.tv_nsec;
}

// Generate a uniformly distributed random number:
double randU(double low, double high)
{
	return (high - low) * drand48() + low;
}

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
std::vector<double> linspace(double low, double high, int numberOfElements) {
	std::vector<double> linVec(numberOfElements);
	double res = (high - low) / (double) (numberOfElements - 1);

	for (int i = 0; i < numberOfElements; i++) {
		linVec[i] = low + i * res;
	}

	return linVec;
}

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
std::vector<double> logspace(double low, double high, int numberOfElements, double base) {
  std::vector<double> logVec(numberOfElements);
  std::vector<double> exponents = linspace(low, high, numberOfElements);

  for (int i = 0; i < numberOfElements; i++){
    logVec[i] = pow(base, exponents[i]);
  }

  return logVec;
}

// Linear interpolation. Evaluates y(reqX) by interpolation the vector y(x).
double linearInterp(std::vector<double> x, std::vector<double> y, double reqX){
  int i = 0;
  double xPrev = x[0], xNext = x[1], yPrev = y[0], yNext = y[1];

  // Check if in range:
  if (reqX < x.front() || reqX > x.back()){
    return 0;
  }

  while (!(reqX >= xPrev & reqX <= xNext)){
    xPrev = x[i]; xNext = x[i+1];
    yPrev = y[i]; yNext = y[i+1];
    i++;
  }

  // Linearly interpolate between points:
  double t = (reqX - xPrev)/(xNext - xPrev);

  return (1 - t) * yPrev + t * yNext;
}

// Surface normal
std::vector<double> surfNormal(std::vector<double> ptA, std::vector<double> ptB, std::vector<double> ptC){
  std::vector<double> v1(3);
  std::vector<double> v2(3);
  std::vector<double> res(3);

  v1[0] = ptB[0] - ptA[0];
  v1[1] = ptB[1] - ptA[1];
  v1[2] = ptB[2] - ptA[2];
  v2[0] = ptC[0] - ptA[0];
  v2[1] = ptC[1] - ptA[1];
  v2[2] = ptC[2] - ptA[2];

  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];

  double norm = vecNorm(res);
  res[0] = res[0] / norm;
  res[1] = res[1] / norm;
  res[2] = -res[2] / norm;
  return res;
}

// Calculate the norm of a vector
double vecNorm(std::vector<double> vec){
  return sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));
}

// Angle between xy plane and vector
double xyPlaneVecAngle(std::vector<double> vec){
  return asin(vec[2]/vecNorm(vec));
}
