#include <complex>

// Class for crater object
class Crater {
public:
double xLocation;
double yLocation;
double transientRadius;
double transientRadiusGravity;
double finalRadius;
double finalDepth;   // final depth, as opposed to transient depth
double finalDepth_init;   // initial final depth
double rimHeight;
double floorElevation;
int numberOfSecondaries;
Layer ejectedMass;
std::vector<double> ejectaDistance;
std::vector<double> ejectaThickness;

Crater(Impactor impactor);
Crater(Impactor impactor, double xLocation, double yLocation);
Crater(Impactor impactor, double xLocation, double yLocation, Layer ejectedMass);
Crater(double finalRadius, double xLocation, double yLocation);

private:
double calcTransientVolume(Impactor impactor);
double calcTransientVolumeGravity(Impactor impactor);
double calcTransientCraterRadius(Impactor impactor);
double calcTransientCraterRadiusGravity(Impactor impactor);
double calcFinalCraterRadius();
void calcEjectaThickness(Impactor impactor);
};
