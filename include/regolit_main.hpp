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
extern double secondaryLargestFraction; // Largest secondary radius as a fraction of the primary radius
extern double secondaryDepthToDiameter; // Depth to diameter ratio of secondary craters
extern double secondaryMinimumVelocity; // Landing speed below which ejecta only builds the blanket, m/s
extern double secondaryVelocityExponent; // Largest fragment shrinks with ejection speed as v^-exponent
extern double maximumSecondariesPerPrimary; // Budget: only the largest N secondaries of a primary are formed
extern bool isEmplaceDistantSecondaries; // Sample primaries outside the domain and form the fragments they send in
extern double secondaryMaximumRange; // Distance beyond the domain edge out to which distant primaries are sampled, m
extern int histogramBinsPerDecade; // Log bins per decade of the size and depth histograms
extern double testCraterDiameter; // Test mode: form one crater of this final diameter (m) instead of the random population (0 = off)
extern double testCraterX; // Test crater center, m from the domain center
extern double testCraterY;
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
extern double k2; // Strength-regime scaling constant (Holsapple 1993 K2; 1 reproduces the original model)
extern double Ybar;
extern double mu;
extern double targetDensity;
extern double seismicEfficiency; // Fraction of the impact energy radiated as seismic energy (Richardson 2005, 2009)
extern bool isVolatiles;         // Ice and soot processes (deposition, sublimation, loss from ejecta); off: passive tracers only
extern bool isSeismicShaking;    // Seismic shaking and local slope collapse after every impact
extern double Q_factor;          // Seismic quality factor
extern double prim_seis_freq;    // Primary seismic frequency, Hz
extern double seis_mean_free;    // Mean free path for seismic scattering, m
extern double seis_wave_vel;     // Seismic P-wave velocity, m/s
extern double seismicAccelerationThreshold; // Regolith moves where the peak acceleration exceeds this times g
extern double Cs, Ki_a, Ki_b, Ki_c, Ki_d;    // Diffusion dose K_i = Cs v^a D^b g^c / l^d (Richardson 2009, Eq. 32)
