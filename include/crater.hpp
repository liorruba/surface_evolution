#pragma once
#include <vector>
#include "layer.hpp"
#include "impactor.hpp"

// One launch annulus of the Z-model: its material leaves at `speed` and lands between
// `landingInner` and `landingOuter` from the crater center (Richardson 2009). The secondary-crater
// model draws ejecta fragments from these shells.
struct EjectaShell {
        double launchRadius;   // inner radius of the launch annulus, m
        double speed;          // ejection (= landing) speed, m/s
        double landingInner;   // nearest landing distance from the crater center, m
        double landingOuter;   // farthest landing distance, m
        double volume;         // ejected volume of the annulus, m^3
};

// Class for crater object
class Crater {
public:
        static long zModelWarnings;   // craters whose Z-model table had skipped annuli (logged for the first few only)

double xLocation;
double yLocation;
double transientRadius;
double transientRadiusGravity;
double finalRadius;
double finalDepth;        // rim-to-floor depth, as opposed to transient depth
double finalDepth_init;   // depth at formation
double rimHeight;
double floorElevation;    // surface elevation at the crater center after formation
double projectileDiameter = 0;   // the impactor (0 for craters given by radius): drives the seismic shaking
double projectileVelocity = 0;
double projectileDensity = 0;
bool isGhost = false;            // periodic copy of a primary: formed, but does not shake or spawn secondaries again
Layer ejectedMass;                    // composition of the excavated material (thickness = summed excavated thickness)
std::vector<double> ejectaDistance;   // distances from the crater center, ascending, m
std::vector<double> ejectaThickness;  // ejecta blanket thickness at those distances, m
std::vector<EjectaShell> ejectaShells; // Z-model launch annuli (speed, landing ring, volume), by ascending launch radius

Crater(Impactor impactor);                                                         // random location
Crater(Impactor impactor, double xLocation, double yLocation);                     // given location
Crater(Impactor impactor, double xLocation, double yLocation, Layer ejectedMass);  // ghost crater: inherits the ejecta composition
Crater(Impactor impactor, double xLocation, double yLocation, double depthRatio);  // given depth/diameter (secondary craters from ejecta fragments)
Crater(double xLocation, double yLocation, double finalRadius, double depthRatio); // given radius and depth/diameter, no ejecta

// Crater scaling (Holsapple 1993 pi-scaling, gravity and strength regimes) for any impactor:
static double calcTransientVolume(const Impactor &impactor);
static double calcTransientVolumeGravity(const Impactor &impactor);
static double calcTransientCraterRadius(const Impactor &impactor);
static double calcTransientCraterRadiusGravity(const Impactor &impactor);
static double finalCraterDiameter(const Impactor &impactor);
// Inverse scaling: impactor radius that makes a crater of the given final diameter at the given speed and density.
static double impactorRadiusForCraterDiameter(double craterDiameter, double velocity, double density);

private:
void initializeFromImpactor(const Impactor &impactor, double depthRatio);
double calcFinalCraterRadius() const;
void calcEjectaThickness(const Impactor &impactor);
};
