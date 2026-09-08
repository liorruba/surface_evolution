#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <string>
#include <cstring>
#include <ctime>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../include/regolit_main.hpp"
#include "../include/log.hpp"
#include "../include/utility.hpp"
#include "../include/histogram.hpp"
#include "../include/impactor.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/grid.hpp"

//////////////////////////////
// DECLARE INPUT PARAMETERS //
//////////////////////////////
// Simulation parameters:
double regionWidth; // m
double resolution; // m / pixel
double endTime; // Ma
double printTimeStep; // Time step for printing data in Ma.
double initialThickness; // Initial thickness of subsurface layer.
double latitude; // Latitude (for shadow calculation; in development)
bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.
bool isEmplaceSecondaries; // Should emplace secondaries? Computationally extensive.
bool runTests; // Should run tests? 
int randomSeed; // Random number generator seed
int isPrintSubsurface; // Print the layer stacks: 0 never, 1 every print step, 2 final step only
double depthToIntegrate; // Depth to integrate when printing integrated subsurface
double downsamplingResolution; // Downsample the results to produce smaller files

// Crater formation variables:
double depthToDiameter; // Dimensionless ratio, crater depth to diameter
double rimToDiameter; // Dimensionless ratio, rim height to diameter
double rimDropoffExponent; // The exponent of the rim height decrease power law
double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
int craterProfileType; // Chosen crater profile
int ejectaSpread; // The spread of the ejecta in crater radii
double ejectaVolatileRetention; // The fraction of volatiles that remain in the crater ejecta
double ejectaSootRetention; // The fraction of soot that remain in the crater ejecta
double minimumLayerThickness; // Deposits thinner than this (m) are mixed into the surface layer
double slope_secondaries; // Slope of the secondary crater size distribution
double secondaryLargestFraction; // Largest secondary radius as a fraction of the primary radius
double secondaryDepthToDiameter; // Depth to diameter ratio of secondary craters
double iceDensity; // The density of ice
double regolithDensity; // The density of regolith
double sootDensity; // The density of "soot"
double c_ice; // The speed of sound in ice
double c_regolith; // The speed of sound in regolith
double s_ice; // slope of the shock / particle velocity relation (ice)
double s_regolith; // slope of the shock / particle velocity relation (regolith)
double ice_fraction; // The fraction of ice in the regolith-ice mixture
double porosity; // The impact target porosity
double sublimationInterval; // The time between two erosion "events"
double sublimationThickness; // The thickness of sublimated layer
double iceEmplacementInterval; // The time between two episodic ice emplacements
double iceEmplacementThickness; // The thickness of ice deposited in each episodic event

// Impactor distribution variables:
double slope_b; // Slope of the impactor CDF
double earthFluxRatioCoefficient; // Moon-earth impactor flux ratio due to cross-section
double minimumImpactorDiameter; // The smallest impactor in the distribution
double fluxConstant_c; // The Flux of impactors > 1 m, Ma^-1 m^-2
double impactorDensity; // The impactor density in kg m^-3
double meanImpactVelocity; // The impact velocity, m s^-1
double impactAngle;
double angleOfRepose; // The regolith angle of repose (deg)

// Surface physical properties:
double g;
double k1;
double k2; // Strength-regime scaling constant (Holsapple 1993 K2)
double Ybar;
double mu;
double targetDensity;
double seismicEfficiency;
double Q_factor;
double prim_seis_freq;
double seis_mean_free;
double seis_wave_vel;

// Seismic diffusivity parameters
double Cs;
double Ki_a;
double Ki_b;
double Ki_c;
double Ki_d;

////////////////
// START MAIN //
////////////////
int main() {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "******************************************************************************" << std::endl;
        std::cout << "******************************************************************************" << std::endl;
        std::cout << "**                                                                          **" << std::endl;
        std::cout << "** REGOLIT: REworking and Gardening Of Lunar Impacted Terrains. Version 1.0 **" << std::endl;
        std::cout << "**                                                                          **" << std::endl;
        std::cout << "******************************************************************************" << std::endl;
        std::cout << "******************************************************************************" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;

        // Prepare directories; create if needed:
        std::cout << "Creating output directory." << std::endl;
        if (std::filesystem::exists("./output")) {
                std::cout << "Clearing existing output." << std::endl;
                std::filesystem::remove_all("./output");
        }
        std::filesystem::create_directories("./output");

        if (std::filesystem::exists("./log/log.txt")) {
                std::cout << "Clearing existing logs." << std::endl;
                std::filesystem::remove("./log/log.txt");
        }
        std::filesystem::create_directories("./log");

        // Create log:
        std::cout << "Creating a new log file..." << std::endl;
        createLogFile("log/log.txt");
        // Read config:
        std::vector<var> varList = readConfig();

        // Set variables from parameters:
        // Simulation variables:
        regionWidth = setVariable(varList, "regionWidth");
        resolution = setVariable(varList, "resolution");
        endTime = setVariable(varList, "endTime");
        printTimeStep = setVariable(varList, "printTimeStep");
        initialThickness = setVariable(varList, "initialThickness");
        latitude = setVariable(varList, "latitude");
        isEmplaceEjecta = setVariable(varList, "isEmplaceEjecta");
        isEmplaceSecondaries = setVariable(varList, "isEmplaceSecondaries");
        runTests = setVariable(varList, "runTests");
        randomSeed = (int) setVariable(varList, "randomSeed");
        downsamplingResolution = (double) setVariable(varList, "downsamplingResolution");
        isPrintSubsurface = (int) setVariable(varList, "isPrintSubsurface");
        depthToIntegrate = (double) setVariable(varList, "depthToIntegrate");

        // Crater formation variables:
        depthToDiameter = setVariable(varList, "depthToDiameter");
        rimToDiameter = setVariable(varList, "rimToDiameter");
        rimDropoffExponent = setVariable(varList, "rimDropoffExponent");
        numberOfZModelShells = setVariable(varList, "numberOfZModelShells");
        craterProfileType = (int) setVariable(varList, "craterProfileType");
        ejectaSpread = (int) setVariable(varList, "ejectaSpread");
        ejectaVolatileRetention = (double) setVariable(varList, "ejectaVolatileRetention");
        ejectaSootRetention = (double) setVariable(varList, "ejectaSootRetention");
        minimumLayerThickness = setVariableOptional(varList, "minimumLayerThickness", 0.0);
        slope_secondaries = setVariable(varList, "slope_secondaries");
        secondaryLargestFraction = setVariableOptional(varList, "secondaryLargestFraction", 0.05);
        secondaryDepthToDiameter = setVariableOptional(varList, "secondaryDepthToDiameter", setVariable(varList, "depthToDiameter"));
        iceDensity = setVariable(varList, "iceDensity");
        regolithDensity = setVariable(varList, "regolithDensity");
        sootDensity = setVariable(varList, "sootDensity");
        c_ice = setVariable(varList, "c_ice");
        c_regolith = setVariable(varList, "c_regolith");
        s_ice = setVariable(varList, "s_ice");
        s_regolith = setVariable(varList, "s_regolith");
        ice_fraction = setVariable(varList, "ice_fraction");
        porosity = setVariable(varList, "porosity");
        iceEmplacementInterval = (double) setVariable(varList, "iceEmplacementInterval");
        iceEmplacementThickness = (double) setVariable(varList, "iceEmplacementThickness");
        sublimationInterval = (double) setVariable(varList, "sublimationInterval");
        sublimationThickness = (double) setVariable(varList, "sublimationThickness");

        // Impactor distribution variables:
        slope_b = setVariable(varList, "slope_b");
        earthFluxRatioCoefficient = setVariable(varList, "earthFluxRatioCoefficient");
        minimumImpactorDiameter = setVariable(varList, "minimumImpactorDiameter");
        fluxConstant_c = setVariable(varList, "fluxConstant_c");
        impactorDensity = setVariable(varList, "impactorDensity");
        meanImpactVelocity = setVariable(varList, "meanImpactVelocity");
        impactAngle = setVariable(varList, "impactAngle");

        // Surface physical properties:
        g = setVariable(varList, "g");
        k1 = setVariable(varList, "k1");
        k2 = setVariableOptional(varList, "k2", 1.0);
        Ybar = setVariable(varList, "Ybar");
        mu = setVariable(varList, "mu");
        targetDensity = setVariable(varList, "targetDensity");
        angleOfRepose = setVariable(varList, "angleOfRepose");

        // Sanity checks on the parameters that would otherwise fail silently:
        if (craterProfileType != 1 && craterProfileType != 2) {
                addLogEntry("ERROR: craterProfileType must be 1 (parabolic) or 2 (bowl-shaped).", true);
                return EXIT_FAILURE;
        }
        if (isPrintSubsurface < 0 || isPrintSubsurface > 2) {
                addLogEntry("ERROR: isPrintSubsurface must be 0 (never), 1 (every print step) or 2 (final step only).", true);
                return EXIT_FAILURE;
        }
        if (endTime <= 0 || printTimeStep <= 0 || regionWidth <= 0 || resolution <= 0) {
                addLogEntry("ERROR: endTime, printTimeStep, regionWidth and resolution must be positive.", true);
                return EXIT_FAILURE;
        }
        if (downsamplingResolution < resolution) {
                addLogEntry("WARNING: downsamplingResolution is finer than resolution; output is written at the grid resolution.", true);
                downsamplingResolution = resolution;
        }

        // Initialize the random number generator seed:
        srand48(randomSeed);

#ifdef _OPENMP
        addLogEntry("Slope relaxation runs on " + std::to_string(omp_get_max_threads()) + " threads (set OMP_NUM_THREADS to change).", true);
#endif

        int printIndex = 1; // The index appended to the output file names (elevation_01.out, ...).

        ////////////////////////////
        // Simulation parameters  //
        ////////////////////////////
        // Read initial layers and pixel index files:
        std::vector< std::vector<double> > initLayersList = readLayers();
        std::vector<int8_t> pixelIndexMatrix = readPixelIndex();

        // Generate grid:
        Grid grid(initLayersList, pixelIndexMatrix);

        // Number of impactors larger than minimumImpactorDiameter per Ma: N / (A t) = c D^-b.
        addLogEntry("Calculating number of craters to be created...", true);
        const double impactorsPerMa = fluxConstant_c * pow(minimumImpactorDiameter, -slope_b) * grid.area * earthFluxRatioCoefficient;
        const long totalNumberOfImpactors = (long) ceil(impactorsPerMa * endTime);
        const long numberOfCratersInTimestep = std::max(1L, (long) ceil(impactorsPerMa * printTimeStep));
        const long numberOfCratersInSublimationPeriod = sublimationInterval > 0 ? std::max(1L, (long) ceil(impactorsPerMa * sublimationInterval)) : 0;
        const long numberOfCratersInDepositionEvent = iceEmplacementInterval > 0 ? std::max(1L, (long) ceil(impactorsPerMa * iceEmplacementInterval)) : 0;

        addLogEntry("Number of craters in simulation: " + std::to_string(totalNumberOfImpactors) + ".", true);

        // Show a warning if the crater records will take a lot of memory:
        if (sizeof(CraterRecord) * totalNumberOfImpactors > 1e9)
                addLogEntry("WARNING: memory taken by the crater records exceeds 1 GB.", true);

        // Craters and impactors histograms:
        Histogram cratersHistogram(minimumImpactorDiameter * 10, regionWidth, 20); // Crater histogram from 10*minimumImpactorDiameter to regionWidth meters
        Histogram impactorsHistogram(minimumImpactorDiameter, 1e4, 20); // Impactor histogram from minimumImpactorDiameter m to 10 km
        Histogram cratersDepthHistogram(minimumImpactorDiameter, 1e4, 20); // Depth histogram of the visible craters at the end of the run

        //////////////////////
        // Start simulation //
        //////////////////////
        addLogEntry("Running simulation...", true);

        for (long i = 0; i < totalNumberOfImpactors; ++i) {
                // Randomize a new impactor:
                Impactor impactor;

                // Add diameter to crater histogram:
                impactorsHistogram.add(2 * impactor.radius);
                // Create a crater instance:
                Crater crater(impactor);
                // Record diameter in crater histogram:
                cratersHistogram.add(2 * crater.finalRadius);

                // Form a crater on the grid:
                if (crater.finalRadius > 0) {
                        grid.formCrater(crater);

                        // "Ghost" craters: the domain is periodic. A crater whose footprint (cavity, rim
                        // and ejecta blanket) crosses an edge of the grid is repeated on the opposite side,
                        // and in the opposite corner if it crosses two edges. The ghosts inherit the
                        // composition of the material excavated by the primary.
                        const double footprint = ejectaSpread * crater.finalRadius;
                        const double halfWidth = regionWidth / 2;
                        double xShift = 0, yShift = 0;
                        if (fabs(crater.xLocation) > halfWidth - footprint)
                                xShift = crater.xLocation > 0 ? -regionWidth : regionWidth;
                        if (fabs(crater.yLocation) > halfWidth - footprint)
                                yShift = crater.yLocation > 0 ? -regionWidth : regionWidth;

                        if (xShift != 0) {
                                Crater ghost(impactor, crater.xLocation + xShift, crater.yLocation, crater.ejectedMass);
                                grid.formCrater(ghost);
                        }
                        if (yShift != 0) {
                                Crater ghost(impactor, crater.xLocation, crater.yLocation + yShift, crater.ejectedMass);
                                grid.formCrater(ghost);
                        }
                        if (xShift != 0 && yShift != 0) {
                                Crater ghost(impactor, crater.xLocation + xShift, crater.yLocation + yShift, crater.ejectedMass);
                                grid.formCrater(ghost);
                        }

                        // Secondary craters: radii follow N(>r) = (r_max / r)^slope between one pixel and
                        // secondaryLargestFraction of the primary radius, within 4 primary radii of the primary.
                        if (isEmplaceSecondaries && crater.numberOfSecondaries > 0) {
                                const double largestSecondaryRadius = secondaryLargestFraction * crater.finalRadius;
                                const double truncation = 1 - pow(resolution / largestSecondaryRadius, slope_secondaries);
                                addLogEntry("Primary diameter: " + std::to_string(2*crater.finalRadius) + ". Number of secondaries: " + std::to_string(crater.numberOfSecondaries) + ".", false);
                                for (long j = 0; j < crater.numberOfSecondaries; j++) {
                                        double secondaryxLocation = randU(crater.xLocation - 4 * crater.finalRadius, crater.xLocation + 4 * crater.finalRadius);
                                        double secondaryyLocation = randU(crater.yLocation - 4 * crater.finalRadius, crater.yLocation + 4 * crater.finalRadius);
                                        double secondaryRadius = resolution * pow(1 - randU(0, 1) * truncation, -1 / slope_secondaries);
                                        Crater secondaryCrater(secondaryxLocation, secondaryyLocation, secondaryRadius, secondaryDepthToDiameter);
                                        grid.formCrater(secondaryCrater);
                                }
                        }
                }

                // Sublimate material every sublimation period:
                if (numberOfCratersInSublimationPeriod > 0 && i % numberOfCratersInSublimationPeriod == 0){
                        grid.sublimateIce();
                }

                // Deposit ice periodically
                if (numberOfCratersInDepositionEvent > 0 && i % numberOfCratersInDepositionEvent == 0){
                        if (iceEmplacementThickness > 0){ 
                                grid.depositLayer(Layer(iceEmplacementThickness, 0, 1, 0));
                        }
                }

                progressBar(i, totalNumberOfImpactors);

                // Print progress to a file:
                if (i % numberOfCratersInTimestep == 0) {
                        // Let slopes above the angle of repose fail:
                        grid.thresholdSlopes(angleOfRepose);
                        // Print the surface, the integrated subsurface and (optionally) the full subsurface:
                        grid.printSurface(printIndex, false);
                        grid.printIntegratedSubsurface(depthToIntegrate, printIndex);

                        if (isPrintSubsurface == 1) {
                                grid.printSubsurface(printIndex);
                        }

                        // Print to log:
                        addLogEntry("Progress: " + std::to_string((double) i / (double) totalNumberOfImpactors * 100) + "%.", false);
                        printIndex++;
                }

        }
        // Print craters histogram to file:
        addLogEntry("Finished running. Saving histograms and data.", true);
        addLogEntry("Printing crater histogram file.", true);
        cratersHistogram.print("./output/craters_histogram.txt");
        impactorsHistogram.print("./output/impactor_histogram.txt");
        grid.sublimateIce(); // Sublimate ice one last time
        grid.thresholdSlopes(angleOfRepose);
        grid.printSurface(printIndex, true);
        grid.printIntegratedSubsurface(depthToIntegrate, printIndex);
        
        if (isPrintSubsurface >= 1) {
                grid.printSubsurface(printIndex);
        }

        grid.printExistingCraters();
        grid.printExistingCratersToHistogram(20);
        addLogEntry("Number of visible craters at the end of the simulation: " + std::to_string(grid.numberOfVisibleCraters()) + ".", true);

        // Depth histogram of the craters that are still visible, with their current (degraded) depths:
        for (const CraterRecord &record : grid.craterRecords()) {
                if (record.isVisible) {
                        cratersDepthHistogram.add(record.finalDepth);
                }
        }
        cratersDepthHistogram.print("./output/depth_histogram.txt");

        addLogEntry("Simulation has ended.", true);
        return EXIT_SUCCESS;
}
