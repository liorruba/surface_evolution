// Class for crater object
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"
#include "../include/crater.hpp"

long Crater::zModelWarnings = 0;

// First constructor: randomize impact location
Crater::Crater(Impactor impactor) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness(), ejectaShells() {
        xLocation = randU(-regionWidth/2, regionWidth/2);
        yLocation = randU(-regionWidth/2, regionWidth/2);
        initializeFromImpactor(impactor, depthToDiameter);
}

// Second constructor: predetermined impact location
Crater::Crater(Impactor impactor, double _xLocation, double _yLocation) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness(), ejectaShells() {
        xLocation = _xLocation;
        yLocation = _yLocation;
        initializeFromImpactor(impactor, depthToDiameter);
}

// Third constructor: predetermined impact location and ejected mass (ghost craters)
Crater::Crater(Impactor impactor, double _xLocation, double _yLocation, Layer _ejectedMass) : ejectedMass(_ejectedMass), ejectaDistance(), ejectaThickness(), ejectaShells() {
        xLocation = _xLocation;
        yLocation = _yLocation;
        isGhost = true;
        initializeFromImpactor(impactor, depthToDiameter);
}

// Fourth constructor: predetermined location and depth/diameter ratio (secondary craters formed by
// ejecta fragments: the radius comes from the scaling with the fragment's mass and landing speed).
Crater::Crater(Impactor impactor, double _xLocation, double _yLocation, double depthRatio) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness(), ejectaShells() {
        xLocation = _xLocation;
        yLocation = _yLocation;
        initializeFromImpactor(impactor, depthRatio);
}

// Fifth constructor: predetermined crater radius and depth/diameter ratio. This type of crater has no ejecta.
Crater::Crater(double _xLocation, double _yLocation, double _finalRadius, double depthRatio) : ejectedMass(Layer(0,0,0,0)), ejectaDistance(), ejectaThickness(), ejectaShells() {
        xLocation = _xLocation;
        yLocation = _yLocation;
        finalRadius = _finalRadius;
        transientRadius = finalRadius / 1.18;
        transientRadiusGravity = transientRadius;
        finalDepth = depthRatio * 2 * finalRadius;
        finalDepth_init = finalDepth;
        rimHeight = rimToDiameter * 2 * finalRadius;
        floorElevation = 0;
}

// Shared part of the impactor-based constructors:
void Crater::initializeFromImpactor(const Impactor &impactor, double depthRatio) {
        projectileDiameter = 2 * impactor.radius;
        projectileVelocity = impactor.velocity;
        projectileDensity = impactor.density;
        transientRadius = calcTransientCraterRadius(impactor);
        transientRadiusGravity = calcTransientCraterRadiusGravity(impactor);
        finalRadius = calcFinalCraterRadius();
        finalDepth = depthRatio * 2 * finalRadius;
        finalDepth_init = finalDepth;
        rimHeight = rimToDiameter * 2 * finalRadius;
        floorElevation = 0;
        calcEjectaThickness(impactor);
}

///////////////////
// Crater physical parameters
///////////////////
// Transient crater volume:
double Crater::calcTransientVolume(const Impactor &impactor) {
        double buff1 = (g * impactor.radius / pow(impactor.velocity,2.0)) * pow(targetDensity/impactor.density, -1.0/3.0);
        double buff2 = pow(Ybar/targetDensity/pow(impactor.velocity, 2.0), (2.0 + mu)/2.0);

        return k1 * (impactor.mass/targetDensity) * pow(buff1 + k2 * buff2, -3 * mu / (2 + mu));
}

// Transient crater volume (gravity regime):
double Crater::calcTransientVolumeGravity(const Impactor &impactor) {
        double buff1 = pow(g * impactor.radius / pow(impactor.velocity, 2), (-3 * mu / (2 + mu)));
        double buff2 = pow(targetDensity / impactor.density, mu / (2 + mu));

        return k1 * (impactor.mass/targetDensity) * buff1 * buff2;
}

// Transient crater radius:
double Crater::calcTransientCraterRadius(const Impactor &impactor) {
        return pow(3 * calcTransientVolume(impactor) / M_PI, 1.0/3.0);
}

// Transient crater radius (gravity regime):
double Crater::calcTransientCraterRadiusGravity(const Impactor &impactor) {
        return pow(3 * calcTransientVolumeGravity(impactor) / M_PI, 1.0/3.0);
}

// Final crater radius:
double Crater::calcFinalCraterRadius() const {
        return 1.18 * transientRadius;
}

// Final crater diameter of any impactor (without constructing a crater):
double Crater::finalCraterDiameter(const Impactor &impactor) {
        return 2 * 1.18 * calcTransientCraterRadius(impactor);
}

// Inverse scaling by bisection in log(radius): the final diameter grows monotonically with the impactor radius.
double Crater::impactorRadiusForCraterDiameter(double craterDiameter, double velocity, double density) {
        double low = log(1e-6), high = log(1e6);   // m
        for (int iteration = 0; iteration < 60; iteration++) {
                const double middle = 0.5 * (low + high);
                if (finalCraterDiameter(Impactor(exp(middle), velocity, density)) < craterDiameter)
                        low = middle;
                else
                        high = middle;
        }
        return exp(0.5 * (low + high));
}

//////////////////
// Ejecta profile:
//////////////////
// Builds the (distance, thickness) table of the ejecta blanket from the Z-model: the transient
// cavity is divided into concentric launch annuli; the material of each annulus is ejected at the
// velocity of its inner radius and lands in an annulus on the surface (Richardson 2007, 2009).
void Crater::calcEjectaThickness(const Impactor &impactor){
        ejectaDistance.clear();
        ejectaThickness.clear();

        const int nShells = std::max(3, (int) numberOfZModelShells);
        const double innerLaunchRadius = 0.1;   // m
        if (transientRadius <= innerLaunchRadius) {
                return;   // Too small for the Z-model; no ejecta table.
        }
        std::vector<double> shellRadius = linspace(innerLaunchRadius, transientRadius, nShells);

        const double Ctg = 0.85; // Richardson 2009, Eq. 20
        const double Cvpg = sqrt(2) / Ctg * (mu/(mu + 1));
        const double transitionStr = targetDensity * pow(impactor.velocity, 2) * pow((g * impactor.radius/pow(impactor.velocity, 2)) * pow(impactor.density/targetDensity,1.0/3.0), 2/(2+mu)); // Richardson 2007, Eq. 18
        const double Cvps = Cvpg * sqrt(targetDensity * g * transientRadiusGravity / (Ybar + transitionStr)) * pow(transientRadiusGravity/transientRadius, 1/mu);
        const double Kg = pow(Cvpg, 2);

        // Ejecta velocity and landing distance of each launch radius (Richardson 2007, Eq. 29;
        // Richardson 2009, Eq. 17 and 23-25). Inner streamtubes are faster and land farther out.
        std::vector<double> landingDistance(nShells), ejectionSpeed(nShells);
        for (int i = 0; i < nShells; i++) {
                const double r = shellRadius[i];
                const double ejectaVelocityGravity = Cvpg * sqrt(g * transientRadiusGravity) * pow(r / transientRadiusGravity, -1/mu);
                const double velocitySquared = pow(ejectaVelocityGravity, 2) - pow(Cvpg, 2) * g * r - pow(Cvps, 2) * Ybar / targetDensity;
                const double finalEjectaVelocity = sqrt(std::max(velocitySquared, 0.0));
                ejectionSpeed[i] = finalEjectaVelocity;
                const double launchAngle = M_PI/180 * (55 - (20 * r/transientRadius));   // 55 deg at the center to 35 deg at the rim
                const double horizontalVelocity = finalEjectaVelocity * cos(launchAngle);
                const double verticalVelocity = finalEjectaVelocity * sin(launchAngle);
                landingDistance[i] = r + 2 * horizontalVelocity * verticalVelocity / g;
        }

        // Ejecta volume of each launch annulus (Richardson 2009, Eq. 22) spread over its landing
        // annulus (Richardson 2009, Eq. 25-27). The table is assembled from the outermost launch
        // annulus (landing nearest the rim) inward, so distances ascend.
        std::vector< std::pair<double, double> > table;   // (distance, thickness)
        int skipped = 0;
        ejectaShells.clear();
        ejectaShells.reserve(nShells - 1);
        for (int i = nShells - 2; i >= 0; i--) {
                const double ejectaVolume = Kg * M_PI * (pow(shellRadius[i+1], 3) - pow(shellRadius[i], 3)) / 0.8;
                const double outer = landingDistance[i];
                const double inner = landingDistance[i+1];
                const double ejectaArea = M_PI * (pow(outer, 2) - pow(inner, 2));
                ejectaShells.push_back(EjectaShell{shellRadius[i], ejectionSpeed[i], std::min(inner, outer), std::max(inner, outer), ejectaVolume});
                if (ejectaArea <= 0) {
                        skipped++;
                        continue;
                }
                table.emplace_back(0.5 * (inner + outer), ejectaVolume / ejectaArea);
        }
        if (skipped > 0) {
                // Logged for the first few craters only: with millions of craters this line would dominate the log.
                if (++zModelWarnings <= 10) {
                        addLogEntry("WARNING: " + std::to_string(skipped) + " Z-model annuli had a non-increasing landing distance and were skipped.", false);
                }
                if (zModelWarnings == 10) {
                        addLogEntry("WARNING: further Z-model annuli warnings are not logged; the total is reported at the end of the run.", false);
                }
        }
        std::sort(table.begin(), table.end());

        for (const std::pair<double, double> &entry : table) {
                ejectaDistance.push_back(entry.first);
                ejectaThickness.push_back(entry.second);
        }
}
