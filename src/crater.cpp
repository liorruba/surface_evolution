// Class for crater object
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include "../include/layer.hpp"
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/utility.hpp"
#include "../include/crater.hpp"

// First constructor: randomize impact location
Crater::Crater(Impactor impactor) : ejectaDistance(), ejectaThickness() {
  xPosition = randU(-regionWidth/2, regionWidth/2);
  yPosition = randU(-regionWidth/2, regionWidth/2);
  transientRadius = calcTransientCraterRadius(impactor.radius, impactor.mass, impactor.velocity);
  transientRadiusGravity = calcTransientCraterRadiusGravity(impactor.radius, impactor.mass, impactor.velocity);
  finalRadius = calcFinalCraterRadius();
  calcEjectaThickness(impactor);
  numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
  ejectedMass = Layer(0,0,0,0); // Mass ejected from the crater stored in a single layer. This will be used to redistribute material in ejecta.
}

// Second constructor: predetermined impact location
Crater::Crater(Impactor impactor, double _xPosition, double _yPosition){
  xPosition = _xPosition;
  yPosition = _yPosition;
  transientRadius = calcTransientCraterRadius(impactor.radius, impactor.mass, impactor.velocity);
  transientRadiusGravity = calcTransientCraterRadiusGravity(impactor.radius, impactor.mass, impactor.velocity);
  finalRadius = calcFinalCraterRadius();
  calcEjectaThickness(impactor);
  numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
}

// Third constructor: predetermined crater radius
Crater::Crater(double _xPosition, double _yPosition, double _finalRadius) : xPosition(_xPosition), yPosition(_yPosition), finalRadius(_finalRadius) {}

///////////////////
// Crater physical parameters
///////////////////
// Transient crater volume:
double Crater::calcTransientVolume(double impactorRadius, double impactorMass, double impactVelocity){
  double buff1 = (g * impactorRadius / pow(impactVelocity,2.0)) * pow(targetDensity/impactorDensity, -1.0/3.0);
  double buff2 = pow(Ybar/targetDensity/pow(impactVelocity, 2.0), (2.0 + mu)/2.0);

  return k1 * (impactorMass/targetDensity) * pow(buff1 + buff2, -3 * mu / (2 + mu));
}

// Transient crater volume (gravity regime):
double Crater::calcTransientVolumeGravity(double impactorRadius, double impactorMass, double impactVelocity){
  double buff1 = pow(g * impactorRadius / pow(impactVelocity, 2), (-3 * mu / (2 + mu)));
  double buff2 = pow(targetDensity / impactorDensity, mu / (2 + mu));

  return k1 * (impactorMass/targetDensity) * buff1 * buff2;
}

// Transient crater radius:
double Crater::calcTransientCraterRadius(double impactorRadius, double impactorMass, double impactVelocity){
  return pow(3 * calcTransientVolume(impactorRadius, impactorMass, impactVelocity) / M_PI, 1.0/3.0);
}

// Transient crater radius (gravity regime):
double Crater::calcTransientCraterRadiusGravity(double impactorRadius, double impactorMass, double impactVelocity){
  return pow(3 * calcTransientVolumeGravity(impactorRadius, impactorMass, impactVelocity) / M_PI, 1.0/3.0);
}

// Final crater radius:
double Crater::calcFinalCraterRadius(){
  return 1.18 * transientRadius;
}

//////////////////
// Ejecta profile:
//////////////////
void Crater::calcEjectaThickness(Impactor impactor){
  std::vector<double> z_model_shell_radius = linspace(0.1, transientRadius, numberOfZModelShells);
  double Ctg = 0.85; // Richardson 2009, Eq. 20
  double Cvpg = sqrt(2) / Ctg * (mu/(mu + 1));
  double transitionStr = targetDensity * pow(impactor.velocity, 2) * pow((g * impactor.radius/pow(impactor.velocity, 2)) * pow(impactorDensity/targetDensity,1.0/3.0), 2/(2+mu)); // Richardson 2007, Eq. 18
  double Cvps = Cvpg * sqrt(targetDensity * g * transientRadiusGravity / (Ybar + transitionStr)) * pow(transientRadiusGravity/transientRadius, 1/mu);
  double Kg = pow(Cvpg, 2);

  // Ejecta volume (Richardson 2009, Eq. 22):
  std::vector<double> ejectaVolume;
  for (size_t i = 0; i < (z_model_shell_radius.size()); i++){
    ejectaVolume.insert(ejectaVolume.begin(), Kg * M_PI * (pow(z_model_shell_radius[i+1], 3) - pow(z_model_shell_radius[i], 3)) / 0.8);
  }

  // Ejecta velocity, launch distance and area (Richardson 2007, Eq. 29, Richardson 2009, Eq. 17 and 23-25):
  double ejectaVelocityGravity;
  double finalEjectaVelocity;
  double horizontalVelocity;
  double verticalVelocity;
  for (size_t i = 0; i < (z_model_shell_radius.size() - 1); i++){
    ejectaVelocityGravity = Cvpg * sqrt(g * transientRadiusGravity) * pow(z_model_shell_radius[i] / transientRadiusGravity, -1/mu);
    finalEjectaVelocity = sqrt(pow(ejectaVelocityGravity, 2) - pow(Cvpg, 2) * g * z_model_shell_radius[i] - pow(Cvps, 2) * Ybar / targetDensity);
    horizontalVelocity = finalEjectaVelocity * cos(M_PI/180 * (55 - (20 * z_model_shell_radius[i]/transientRadius)));
    verticalVelocity = finalEjectaVelocity * sin(M_PI/180 * (55 - (20 * z_model_shell_radius[i]/transientRadius)));
    ejectaDistance.insert(ejectaDistance.begin(), z_model_shell_radius[i] + 2 * horizontalVelocity * verticalVelocity / g);
  }

  // Ejecta area and thickness 25-27:
  double ejectaArea;
  for (size_t i = (z_model_shell_radius.size() - 1); i > 0 ; i--){
    ejectaArea = M_PI * (pow(ejectaDistance[i], 2) - pow(ejectaDistance[i - 1], 2));
    ejectaThickness.insert(ejectaThickness.begin(), ejectaVolume[i]/ejectaArea);
  }
}
