// Globals
// Simulation parameters:
extern double regionWidth; // km
extern double resolution; // 1/km
extern double endTime; // Ma
extern double printTimeStep; // Time step for printing data in Ma.
extern double initialThickness; // Initial thickness of subsurface layer.
extern double latitude; // Latitude (for shadow calculation)
extern bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.
extern bool isEmplaceSecondaries; // Should emplace secondaries? Computationally extensive.
extern bool runTests; // Run tests?
extern int randomSeed; // Random number generator seed
extern bool isPrintSubsurface; // Print subsurface? May be disk space demanding 
extern double depthToIntegrate; // Depth to integrate when printing integrated subsurface
extern double downsamplingResolution; // Downsample the results to produce smaller files

// Crater formation variables:
extern double depthToDiameter; // Dimensionless ratio
extern double rimToDiameter; // Dimensionless ratio, rim height to diameter
extern double rimDropoffExponent; // The exponent of the rim height decrease power law
extern double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
extern int craterProfileType; // Chosen crater profile
extern int ejectaSpread; // The spread of the ejecta in crater radii
extern double ejectaVolatileRetention; // The fraction of volatiles that remain in the caterr ejecta
extern double ejectaSootRetention; // The fraction of soot that remain in the caterr ejecta
extern double slope_secondaries; // The spread of the ejecta in crater radii
extern double iceDensity; // The density of ice
extern double regolithDensity; // The density of ice
extern double sootDensity; // The density of "soot"
extern double c_ice; // The speed of sound in ice
extern double c_regolith; // The speed of sound in regolith
extern double s_ice; // slope of the shock / particle velocity relation (ice)
extern double s_regolith; // slope of the shock / particle velocity relation (regolith)
extern double ice_fraction; // The fraction of ice in the regolith-ice mixture
extern double porosity; // The impact target porosity
extern double iceEmplacementInterval; // The time between two episodic ice emplacements
extern double iceEmplacementThickness; // The thickness of ice deposited in each episodic event
extern double sublimationInterval; // The time between two erosion "events"
extern double sublimationThickness; // The thickness of sublimated layer


// Impactor distribution variables:
extern double slope_b; // Slope of the impactor CDF
extern double moonEarthFluxRatio; // Moon-earth impactor flux ratio due to cross-section
extern double minimumImpactorDiameter; // The smallest impactor in the distribution
extern double fluxConstant_c; // The Flux of impactors > 1 m, Ma^-1 m^-2
extern double impactorDensity; // The impactor density in kg m^-3
extern double meanImpactVelocity; // The impact velocity, kg s^-1
extern double impactAngle; // The impact angle, degrees
extern double angleOfRepose; // The regolith angle of repose, degrees

// Surface physical properties:
extern double g;
extern double k1;
extern double Ybar;
extern double mu;
extern double targetDensity;
extern double seismicEfficiency;
extern double Q_factor;
extern double prim_seis_freq;
extern double seis_mean_free;
extern double seis_wave_vel;

// Seismic diffusivity parameters
extern double Cs;
extern double Ki_a;
extern double Ki_b;
extern double Ki_c;
extern double Ki_d;

////////////////
// START MAIN //
////////////////
int main();
