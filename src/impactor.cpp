// Class for impactor object
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/crater.hpp"
#include "../include/utility.hpp"

// Initialize an impactor from a cumulative distribution
Impactor::Impactor() {
  double quantile = randU(0,1);
  radius = 0.5 * minimumImpactorDiameter * pow(quantile, -1/slope_b); // Set impactor radius from the cumulative distribution, meters
  mass = calcMass(radius);
  velocity = meanImpactVelocity;
}

// Initialize an impactor with given parameters
Impactor::Impactor(double _radius) {
  radius = _radius; // Set impactor radius from the cumulative distribution, meters
  mass = calcMass(_radius);
  velocity = meanImpactVelocity;
}

// Calculate impactor mass, kg 
double Impactor::calcMass(double radius){
  return 4.0 / 3.0 * M_PI * pow(radius ,3) * impactorDensity;
}
