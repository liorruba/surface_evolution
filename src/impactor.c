#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../include/regolit_main.h"

double impactorMass(double *impactorDiameter, double *impactorDensity){

  double impactorRadius = *impactorDiameter;

  return 4/3 * pi * pow(impactorRadius ,3) * impactorDensity;
}
