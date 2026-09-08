#pragma once
#include <vector>
#include "layer.hpp"
#include "impactor.hpp"

// Class for crater object
class Crater {
public:
double xLocation;
double yLocation;
double transientRadius;
double transientRadiusGravity;
double finalRadius;
double finalDepth;        // rim-to-floor depth, as opposed to transient depth
double finalDepth_init;   // depth at formation
double rimHeight;
double floorElevation;    // surface elevation at the crater center after formation
int numberOfSecondaries;
Layer ejectedMass;                    // composition of the excavated material (thickness = summed excavated thickness)
std::vector<double> ejectaDistance;   // distances from the crater center, ascending, m
std::vector<double> ejectaThickness;  // ejecta blanket thickness at those distances, m

Crater(Impactor impactor);                                                         // random location
Crater(Impactor impactor, double xLocation, double yLocation);                     // given location
Crater(Impactor impactor, double xLocation, double yLocation, Layer ejectedMass);  // ghost crater: inherits the ejecta composition
Crater(double xLocation, double yLocation, double finalRadius, double depthRatio); // given radius and depth/diameter, no ejecta (secondaries)

private:
void initializeFromImpactor(const Impactor &impactor);
double calcTransientVolume(const Impactor &impactor) const;
double calcTransientVolumeGravity(const Impactor &impactor) const;
double calcTransientCraterRadius(const Impactor &impactor) const;
double calcTransientCraterRadiusGravity(const Impactor &impactor) const;
double calcFinalCraterRadius() const;
void calcEjectaThickness(const Impactor &impactor);
};
