#pragma once
// Globals, defined in src/regolit_main.cpp and read from config/config.cfg.

// Simulation parameters:
extern double regionWidth; // Width of the (square) simulated region, m
extern double resolution; // Pixel size, m
extern double endTime; // Ma
extern double printTimeStep; // Time step for printing data in Ma.
extern double initialThickness; // Initial thickness of the basement regolith layer, m.
extern double latitude; // Latitude (for shadow calculation; in development)
extern bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.
extern bool isEmplaceSecondaries; // Should emplace secondaries? Computationally extensive.
extern bool runTests; // Run tests?
extern int randomSeed; // Random number generator seed
extern int isPrintSubsurface; // Print the layer stacks: 0 never, 1 every print step, 2 final step only
extern double depthToIntegrate; // Depth to integrate when printing integrated subsurface
extern double downsamplingResolution; // Downsample the results to produce smaller files

// Crater formation variables:
extern double depthToDiameter; // Dimensionless ratio, crater depth to diameter
extern double rimToDiameter; // Dimensionless ratio, rim height to diameter
extern double rimDropoffExponent; // The exponent of the rim height decrease power law
extern double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
extern int craterProfileType; // Chosen crater profile
extern int ejectaSpread; // The spread of the ejecta in crater radii
extern double ejectaVolatileRetention; // The fraction of volatiles that remain in the crater ejecta
extern double ejectaSootRetention; // The fraction of soot that remain in the crater ejecta
extern double minimumLayerThickness; // Deposits thinner than this (m) are mixed into the surface layer
extern double slope_secondaries; // Slope of the secondary crater size distribution
extern double iceDensity; // The density of ice
extern double regolithDensity; // The density of regolith
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
extern double earthFluxRatioCoefficient; // Moon-earth impactor flux ratio due to cross-section
extern double minimumImpactorDiameter; // The smallest impactor in the distribution
extern double fluxConstant_c; // The flux of impactors > 1 m, Ma^-1 m^-2
extern double impactorDensity; // The impactor density in kg m^-3
extern double meanImpactVelocity; // The impact velocity, m s^-1
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
