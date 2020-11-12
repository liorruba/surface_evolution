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
        calculateProducedMelt(impactor);

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
        calculateProducedMelt(impactor);

        // Number of secondaries
        numberOfSecondaries = pow(0.05 * finalRadius, slope_secondaries) * pow(resolution, -slope_secondaries);
}

// Third constructor: predetermined crater radius
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

///////////////////
// Melt Production:
///////////////////
void Crater::calculateProducedMelt(Impactor impactor) {
        double mixture_density = iceDensity * ice_fraction +
                                 regolithDensity * (1-ice_fraction);

        // Average sound speeds and slope of Hugoniot
        double c_mix = (c_ice * iceDensity + c_regolith * regolithDensity) /
                       (iceDensity + regolithDensity);
        double s_mix = (s_ice * iceDensity + -s_regolith * regolithDensity) /
                       (iceDensity + regolithDensity);

        // The effective impact velocity
        double Ustar = -c_ice * iceDensity
                       + sqrt( pow(c_ice, 2) *
                               pow(iceDensity, 2) + 2 * c_mix *
                               iceDensity * mixture_density * s_ice * impactor.velocity +
                               iceDensity * mixture_density * s_ice * s_mix *
                               pow(impactor.velocity, 2) ) / (iceDensity * s_ice);

        // The log10 of the mass of the melt + vapor / impactor mass
        double log10_M_melt_vapor_impactor_mass = -0.53 + 0.0017 * temperature + 0.7 *
                                                  log10(sin(impactAngle * M_PI / 180)) -
                                                  0.46 * porosity + 3 * (0.554 + 0.07 * porosity) / 2 *
                                                  log10(pow(Ustar, 2) / meltEnergy);

        // The log10 of the mass of vapor / impactor mass
        double log10_M_vapor_impactor_mass = -1.71 + 0.0011 * temperature + 0.6 *
                                             log10(sin(impactAngle * M_PI / 180))
                                             + 0.27 * porosity + 3 * (0.65 - 0.1 * porosity)/2 *
                                             log10(pow(Ustar, 2) / meltEnergy);


        meltVolume = 4 * M_PI / 3 * pow(impactor.radius, 3) * (
                pow(10, log10_M_melt_vapor_impactor_mass) -
                pow(10, log10_M_vapor_impactor_mass)
                );

        // Use the volume to calculate the height of a paraboloid or spher. cap:
        if (craterProfileType == 1) {
                std::cout << "ERROR: radius profile == 1 is not yet implemented" << std::endl;
                exit(0);
        }

        else if (craterProfileType == 2) {

                double sphereRadius = (pow(6000,2) +
                                       pow(1200,2)) / 2 / 1200;

                meltVolume = 12 * pow(1000, 3);
                std::complex<double> imag_unit(0, 1);
                std::complex<double> buff_radical = pow(3. * meltVolume / 2. / M_PI - pow(sphereRadius, 3.), 2.) - pow(sphereRadius, 6.);
                std::complex<double> buff_cube_root = pow(sqrt(buff_radical) - 3. * meltVolume / 2. / M_PI + pow(sphereRadius, 3.), 1./3.);

                meltHeight_c = sphereRadius - buff_cube_root / 2. -
                               pow(sphereRadius, 2.) / 2. / buff_cube_root -
                               sqrt(3.) * (buff_cube_root - pow(sphereRadius, 2.) /
                                           buff_cube_root) * imag_unit / 2.;

                meltHeight = meltHeight_c.real();
        }

        else {
                meltHeight = 0;
        }
}
