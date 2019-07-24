// Structs
// A quick map to store variables:
typedef struct{
  std::string name;
  double value;
} var;

// Globals
// Simulation parameters:
extern double regionWidth; // km
extern double resolution; // 1/km
extern double endTime; // Ma
extern double printTimeStep; // Time step for printing data in Ma.
extern bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.

// Crater formation variables:
extern double depthToDiameter; // Dimensionless ratio
extern double rimToDiameter; // Dimensionless ratio, rim height to diameter
extern double outerRimExponent; // The exponent of the rim height decrease power law
extern double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
extern int craterProfile; // Chosen crater profile
extern int ejectaSpread; // The spread of the ejecta in crater radii
extern double slope_secondaries; // The spread of the ejecta in crater radii

// Impactor distribution variables:
extern double slope_b; // Slope of the impactor CDF
extern double moonEarthFluxRatio; // Moon-earth impactor flux ratio due to cross-section
extern double minimumImpactorDiameter; // The smallest impactor in the distribution
extern double fluxConstant_c; // The Flux of impactors > 1 m, Ma^-1 m^-2
extern double impactorDensity; // The impactor density in kg m^-3
extern double meanImpactVelocity; // The impact velocity, kg s^-1

// Surface physical properties:
extern double g;
extern double k1;
extern double Ybar;
extern double mu;
extern double targetDensity;

// Read config file:
std::vector<var> readConfig();

// Get variable from list:
// double setVariable(vector<var> varList, char * varName);

////////////////
// START MAIN //
////////////////
int main();
