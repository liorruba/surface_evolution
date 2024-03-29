// Class for crater object
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <complex>
#include "../include/layer.hpp"
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/utility.hpp"
#include "../include/crater.hpp"

// First constructor: randomize impact location
Crater::Crater(Impactor impactor) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness() {
        xLocation = randU(-regionWidth/2, regionWidth/2);
        yLocation = randU(-regionWidth/2, regionWidth/2);
        transientRadius = calcTransientCraterRadius(impactor);
        transientRadiusGravity = calcTransientCraterRadiusGravity(impactor);
        finalRadius = calcFinalCraterRadius();
        finalDepth = depthToDiameter * 2 * finalRadius;
        finalDepth_init = finalDepth;
        rimHeight = rimToDiameter * 2 * finalRadius;
        calcEjectaThickness(impactor);

        // Number of secondaries
        numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
}

// Second constructor: predetermined impact location
Crater::Crater(Impactor impactor, double _xLocation, double _yLocation) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness(){
        xLocation = _xLocation;
        yLocation = _yLocation;
        transientRadius = calcTransientCraterRadius(impactor);
        transientRadiusGravity = calcTransientCraterRadiusGravity(impactor);
        finalRadius = calcFinalCraterRadius();
        finalDepth = depthToDiameter * 2 * finalRadius;
        finalDepth_init = finalDepth;
        rimHeight = rimToDiameter * 2 * finalRadius;
        calcEjectaThickness(impactor);

        // Number of secondaries
        numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
}

// Third constructor: predetermined impact location and ejected mas
Crater::Crater(Impactor impactor, double _xLocation, double _yLocation, Layer _ejectedMass) : ejectedMass(_ejectedMass), ejectaDistance(), ejectaThickness(){
        xLocation = _xLocation;
        yLocation = _yLocation;
        transientRadius = calcTransientCraterRadius(impactor);
        transientRadiusGravity = calcTransientCraterRadiusGravity(impactor);
        finalRadius = calcFinalCraterRadius();
        finalDepth = depthToDiameter * 2 * finalRadius;
        finalDepth_init = finalDepth;
        rimHeight = rimToDiameter * 2 * finalRadius;
        calcEjectaThickness(impactor);

        // Number of secondaries
        numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
}

// Fourth constructor: predetermined crater radius
// This type of crater has no ejecta
Crater::Crater(double _xLocation, double _yLocation, double _finalRadius) : xLocation(_xLocation), yLocation(_yLocation), finalRadius(_finalRadius), ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness() {
        finalDepth = depthToDiameter * 2 * finalRadius;
        rimHeight = rimToDiameter * 2 * finalRadius;
        numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
        finalDepth = depthToDiameter * 2 * finalRadius;
        finalDepth_init = finalDepth;

        // Number of secondaries
        numberOfSecondaries = 0;
}

///////////////////
// Crater physical parameters
///////////////////
// Transient crater volume:
double Crater::calcTransientVolume(Impactor impactor){
        double buff1 = (g * impactor.radius / pow(impactor.velocity,2.0)) * pow(targetDensity/impactorDensity, -1.0/3.0);
        double buff2 = pow(Ybar/targetDensity/pow(impactor.velocity, 2.0), (2.0 + mu)/2.0);

        return k1 * (impactor.mass/targetDensity) * pow(buff1 + buff2, -3 * mu / (2 + mu));
}

// Transient crater volume (gravity regime):
double Crater::calcTransientVolumeGravity(Impactor impactor){
        double buff1 = pow(g * impactor.radius / pow(impactor.velocity, 2), (-3 * mu / (2 + mu)));
        double buff2 = pow(targetDensity / impactorDensity, mu / (2 + mu));

        return k1 * (impactor.mass/targetDensity) * buff1 * buff2;
}

// Transient crater radius:
double Crater::calcTransientCraterRadius(Impactor impactor){
        return pow(3 * calcTransientVolume(impactor) / M_PI, 1.0/3.0);
}

// Transient crater radius (gravity regime):
double Crater::calcTransientCraterRadiusGravity(Impactor impactor){
        return pow(3 * calcTransientVolumeGravity(impactor) / M_PI, 1.0/3.0);
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
        for (size_t i = 0; i < (z_model_shell_radius.size()); i++) {
                ejectaVolume.insert(ejectaVolume.begin(), Kg * M_PI * (pow(z_model_shell_radius[i+1], 3) - pow(z_model_shell_radius[i], 3)) / 0.8);
        }

        // Ejecta velocity, launch distance and area (Richardson 2007, Eq. 29, Richardson 2009, Eq. 17 and 23-25):
        double ejectaVelocityGravity;
        double finalEjectaVelocity;
        double horizontalVelocity;
        double verticalVelocity;
        for (size_t i = 0; i < (z_model_shell_radius.size() - 1); i++) {
                ejectaVelocityGravity = Cvpg * sqrt(g * transientRadiusGravity) * pow(z_model_shell_radius[i] / transientRadiusGravity, -1/mu);
                finalEjectaVelocity = sqrt(pow(ejectaVelocityGravity, 2) - pow(Cvpg, 2) * g * z_model_shell_radius[i] - pow(Cvps, 2) * Ybar / targetDensity);
                horizontalVelocity = finalEjectaVelocity * cos(M_PI/180 * (55 - (20 * z_model_shell_radius[i]/transientRadius)));
                verticalVelocity = finalEjectaVelocity * sin(M_PI/180 * (55 - (20 * z_model_shell_radius[i]/transientRadius)));
                ejectaDistance.insert(ejectaDistance.begin(), z_model_shell_radius[i] + 2 * horizontalVelocity * verticalVelocity / g);
        }

        // Ejecta area and thickness (Richardson 2009 Eq. 25-27):
        double ejectaArea;
        for (size_t i = (z_model_shell_radius.size() - 1); i > 0; i--) {
                ejectaArea = M_PI * (pow(ejectaDistance[i], 2) - pow(ejectaDistance[i - 1], 2));
                ejectaThickness.insert(ejectaThickness.begin(), ejectaVolume[i]/ejectaArea);
        }
}
