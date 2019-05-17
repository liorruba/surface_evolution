// This file contains functions that create a 2D grid.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
double * linspace(double low, double high, int numberOfElements) {
	double *linVec = (double *)malloc(sizeof(double) * (numberOfElements));
	double res = (high - low) / (double) (numberOfElements - 1);

	for (int i = 0; i < numberOfElements; i++) {
		linVec[i] = low + i * res;
	}

	return linVec;
}

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
double * logspace(double low, double high, int numberOfElements, double base) {
  double * logVec = (double *)malloc(sizeof(double) * (numberOfElements));
  double * exponents = linspace(low, high, numberOfElements);

  for (int i = 0; i < numberOfElements; i++){
    logVec[i] = pow(base, exponents[i]);
  }

  return logVec;
}

// This function creates an array xmat to serve as the x coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridX(double *x, double *y, int lengthx, int lengthy){
	// Allocate memory to the array:
	double **xmat = (double **)malloc(sizeof(double) * (lengthy));
	for (int i = 0; i < lengthy; i++)
		xmat[i] = (double *)malloc(sizeof(double) * (lengthx));

	// Populate the array:
	for (int i = 0; i < lengthx; i++) {
		for (int j = 0; j < lengthy; j++) {
			xmat[i][j] = x[j];
		}
	}

	return xmat;
}

// This function creates an array ymat to serve as the y coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridY(double *x, double *y, int lengthx, int lengthy){
	// Allocate memory to the array:
	double **ymat = (double **)malloc(sizeof(double) * (lengthx));
	for (int i = 0; i < lengthx; i++)
		ymat[i] = (double *)malloc(sizeof(double) * (lengthy));

	// Populate the array:
	for (int i = 0; i < lengthx; i++) {
		 for (int j = 0; j < lengthy; j++) {
			ymat[j][i] = y[j];
		}
	}

	return ymat;
}
