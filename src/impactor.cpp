// Class for impactor object
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include "../include/regolit_main.hpp"
#include "../include/layer.hpp"
#include "../include/impactor.hpp"
#include "../include/crater.hpp"
#include "../include/utility.hpp"

// Initialize an impactor from a cumulative distribution
Impactor::Impactor() {
  double quantile = randU(0,1);
  radius = 0.5 * minimumImpactorDiameter * pow(quantile, -1/slope_b); // Set impactor radius from the cumulative distribution, meters
  velocity = meanImpactVelocity;
  density = impactorDensity;
  mass = calcMass(radius);
}

// Initialize an impactor with given radius
Impactor::Impactor(double _radius) {
  radius = _radius;
  velocity = meanImpactVelocity;
  density = impactorDensity;
  mass = calcMass(_radius);
}

// Initialize a fully specified impactor (ejecta fragments: target material at the landing speed)
Impactor::Impactor(double _radius, double _velocity, double _density) {
  radius = _radius;
  velocity = _velocity;
  density = _density;
  mass = calcMass(_radius);
}

// Calculate impactor mass, kg
double Impactor::calcMass(double _radius) const {
  return 4.0 / 3.0 * M_PI * pow(_radius, 3) * density;
}
