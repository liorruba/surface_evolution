// This file contains functions that create a 2D grid.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
double * linspace(double low, double high, int numberOfElements);

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
double * logspace(double low, double high, int numberOfElements, double base);

// This function creates an array xmat to serve as the x coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridX(double *x, double *y, int lengthx, int lengthy);

// This function creates an array ymat to serve as the y coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridY(double *x, double *y, int lengthx, int lengthy);
