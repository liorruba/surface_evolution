#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS
#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING 1

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
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <inttypes.h>
#include <dirent.h>
#include <vector>
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
double regionWidth; // km
double resolution; // 1/km
double endTime; // Ma
double printTimeStep; // Time step for printing data in Ma.
double initialThickness; // Initial thickness of subsurface layer.
double latitude; // Latitude (for shadow calculation; in development)
bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.
bool isEmplaceSecondaries; // Should emplace ejecta? Computationally extensive.
bool runTests; // Should run tests? 
int randomSeed; // Random number generator seed
bool isPrintSubsurface; // Print the full subsurface?
double depthToIntegrate; // Depth to integrate when printing integrated subsurface
double downsamplingResolution; // Downsample the results to produce smaller files

// Crater formation variables:
double depthToDiameter; // Dimensionless ratio, crater depth to diameter
double rimToDiameter; // Dimensionless ratio, rim height to diameter
double rimDropoffExponent; // The exponent of the rim height decrease power law
double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
int craterProfileType; // Chosen crater profile
int ejectaSpread; // The spread of the ejecta in crater radii
double ejectaVolatileRetention; // The fraction of volatiles that remain in the caterr ejecta
double ejectaSootRetention; // The fraction of soot that remain in the caterr ejecta
double slope_secondaries; // The spread of the ejecta in crater radii
double iceDensity; // The density of ice
double regolithDensity; // The density of ice
double sootDensity; // The density of "soot"
double c_ice; // The speed of sound in ice
double c_regolith; // The speed of sound in regolith
double s_ice; // slope of the shock / particle velocity relation (ice)
double s_regolith; // slope of the shock / particle velocity relation (regolith)
double ice_fraction; // The fraction of ice in the regolith-ice mixture
double porosity; // The impact target porosity
double sublimationInterval; // The time between two erosion "events"
double sublimationThickness; // The thickness of sublimated layer
double iceEmplacementInterval; // The time between two erosion "events"
double iceEmplacementThickness; // The thickness of sublimated layer

// Impactor distribution variables:
double slope_b; // Slope of the impactor CDF
double earthFluxRatioCoefficient; // Moon-earth impactor flux ratio due to cross-section
double minimumImpactorDiameter; // The smallest impactor in the distribution
double fluxConstant_c; // The Flux of impactors > 1 m, Ma^-1 m^-2
double impactorDensity; // The impactor density in kg m^-3
double meanImpactVelocity; // The impact velocity, kg s^-1
double impactAngle;
double angleOfRepose; // The regolith angle of repose (deg)

// Surface physical properties:
double g;
double k1;
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
        DIR * outputdir = opendir("./output");
        if (outputdir) {
                closedir(outputdir);
                std::cout << "Clearing existing output." << std::endl;
                std::filesystem::remove_all("./output/");
                system("mkdir ./output");
        }
        else {
                system("mkdir ./output");
        }

        DIR * logdir = opendir("./log");
        if (logdir) {
                closedir(logdir);
                std::cout << "Clearing existing logs." << std::endl;
                std::filesystem::remove_all("./log/log.txt");
                system("mkdir ./log");
        }
        else {
                system("mkdir ./log");
        }

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
        isPrintSubsurface = (bool) setVariable(varList, "isPrintSubsurface");
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
        slope_secondaries = setVariable(varList, "slope_secondaries");
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
        Ybar = setVariable(varList, "Ybar");
        mu = setVariable(varList, "mu");
        targetDensity = setVariable(varList, "targetDensity");
        angleOfRepose = setVariable(varList, "angleOfRepose");

        // Initialize the random number generator seed:
        srand48(randomSeed);

        ////////////////////////////////////////////
        ////////////////////////////////////////////
        // Declare crater related varialbes (DO NO TOUCH!):
        double xGhost; // Ghost craters coordinates
        double yGhost; // Ghost craters coordinates
        int printIndex = 1; // The index to add to the end of the output matrix file (zmat_1.txt).
        ////////////////////////////////////////////
        ////////////////////////////////////////////

        ////////////////////////////
        // Simulation parameters  //
        ////////////////////////////
        // Read initial layers and pixel index files:
        std::vector< std::vector<double> > initLayersList = readLayers();
        std::vector<int8_t> pixelIndexMatrix = readPixelIndex();

        // Generate grid:
        Grid grid = Grid(initLayersList, pixelIndexMatrix);

        // Total number of craters to be created:
        addLogEntry("Calculating number of caters to be created...", true);
        long totalNumberOfImpactors = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * endTime * grid.area * earthFluxRatioCoefficient); // total number of impactors to generate larger than minimumDiameter: N/At = cD^-b.
        long numberOfCratersInTimestep = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * printTimeStep * grid.area * earthFluxRatioCoefficient); // number of impactors to generate larger than minimumDiameter: N/At = cD^-b in some time interval.
        long numberOfCratersInSublimationPeriod = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * sublimationInterval * grid.area * earthFluxRatioCoefficient); // number of impactors in 10 Ma, the time period after which material gets sublimated
        long numberOfCratersInDepositionEvent = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * iceEmplacementInterval * grid.area * earthFluxRatioCoefficient); // number of impactors in 10 Ma, the time period after which material gets sublimated

        char logEntry[50];
        sprintf(logEntry, "Number of craters in simulation: %ld.", totalNumberOfImpactors);
        addLogEntry(logEntry, true);

        // Show warning if totalNumberOfImpactors is too high. 
        // If the memory taken by the craters exceeds 1 GB, throw a memory warning:
        if (sizeof(Crater) * totalNumberOfImpactors > 1e9)
                addLogEntry("WARNING: memory taken by craters exceeds 1 GB. Press ENTER to continue.", true);

        // Craters and impactors histograms:
        Histogram cratersHistogram(minimumImpactorDiameter * 10, regionWidth, 20); // Crater histogram from 10*minimumImpactorDiameter to regionWidth meters
        Histogram impactorsHistogram(minimumImpactorDiameter, 1e4, 20); // Impactor histogram from minimumImpactorDiameter m to 10 km
        Histogram cratersDepthHistogram(minimumImpactorDiameter, 1e4, 20); // Impactor histogram from minimumImpactorDiameter m to 10 km

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

                        if (isEmplaceEjecta && !crater.ejectedMass.isEmpty()) {
                                grid.emplaceEjecta(crater);
                        }

                        
                        // "Ghost" craters:
                        // If the crater exceeds the grid, wrap around it by creating a ghost crater.
                        // If a corner:
                        if ( (fabs(crater.xLocation) > regionWidth/2 - crater.finalRadius) && (fabs(crater.yLocation) > regionWidth/2 - crater.finalRadius) ) {
                                // Calculate ghost crater location as sgn(x) * (region_width - x);
                                xGhost = (-crater.xLocation/fabs(crater.xLocation)) * (regionWidth - crater.xLocation * (crater.xLocation/fabs(crater.xLocation)));
                                yGhost = (-crater.yLocation/fabs(crater.yLocation)) * (regionWidth - crater.yLocation * (crater.yLocation/fabs(crater.yLocation)));
                                Crater ghost1 = Crater(impactor, xGhost, yGhost, crater.ejectedMass);
                                Crater ghost2 = Crater(impactor, xGhost, crater.yLocation, crater.ejectedMass);
                                Crater ghost3 = Crater(impactor, crater.xLocation, yGhost, crater.ejectedMass);


                                grid.formCrater(ghost1);
                                if (isEmplaceEjecta && !ghost1.ejectedMass.isEmpty()) {
                                        grid.emplaceEjecta(ghost1);
                                }

                                grid.formCrater(ghost2);
                                if (isEmplaceEjecta && !ghost2.ejectedMass.isEmpty()) {
                                        grid.emplaceEjecta(ghost2);
                                }

                                grid.formCrater(ghost3);
                                if (isEmplaceEjecta && !ghost3.ejectedMass.isEmpty()) {
                                        grid.emplaceEjecta(ghost3);
                                }
                        }

                        // If a side:
                        if ( (fabs(crater.xLocation) > regionWidth/2 - crater.finalRadius/2) && (fabs(crater.yLocation) < regionWidth/2 - crater.finalRadius/2)) {
                                // Calculate ghost crater location as sgn(x) * (region_width - x);
                                xGhost = (-crater.xLocation/fabs(crater.xLocation)) * (regionWidth - crater.xLocation * (crater.xLocation/fabs(crater.xLocation)));
                                Crater ghost2 = Crater(impactor, xGhost, crater.yLocation, crater.ejectedMass);

                                // Create one ghost crater:
                                grid.formCrater(ghost2);
                                if (isEmplaceEjecta && !ghost2.ejectedMass.isEmpty()) {
                                        grid.emplaceEjecta(ghost2);
                                }
                        }

                        // If another side:
                        if ( (fabs(crater.xLocation) < regionWidth/2 - crater.finalRadius/2) && (fabs(crater.yLocation) > regionWidth/2 - crater.finalRadius/2)) {
                                // Calculate ghost crater location as sgn(x) * (region_width - x);
                                yGhost = (-crater.yLocation/fabs(crater.yLocation)) * (regionWidth - crater.yLocation * (crater.yLocation/fabs(crater.yLocation)));
                                Crater ghost3 = Crater(impactor, crater.xLocation, yGhost, crater.ejectedMass);

                                // Create one ghost crater:
                                grid.formCrater(ghost3);
                                if (isEmplaceEjecta && !ghost3.ejectedMass.isEmpty()) {
                                        grid.emplaceEjecta(ghost3);
                                }
                        }

                        // Secondary craters:
                        if (isEmplaceSecondaries) {
                                // Form a secondary crater within 2 crater diameters from the primary:
                                if (crater.numberOfSecondaries > 0) {
                                        char logEntry[400];
                                        sprintf(logEntry, "Primary diameter: %f. Number of secondaries: %d.", 2*crater.finalRadius, crater.numberOfSecondaries);
                                        addLogEntry(logEntry, true);
                                }
                                for (long j = 0; j < crater.numberOfSecondaries; j++) {
                                        double secondaryxLocation = randU(crater.xLocation - 4 * crater.finalRadius, crater.xLocation + 4 * crater.finalRadius);
                                        double secondaryyLocation = randU(crater.yLocation - 4 * crater.finalRadius, crater.yLocation + 4 * crater.finalRadius);
                                        double secondaryRadius = resolution * pow(randU(0,1), -1/slope_secondaries); // Set impactor radius from the cumulative distribution, meters
                                        Crater secondaryCrater(secondaryxLocation, secondaryyLocation, secondaryRadius);
                                        grid.formCrater(secondaryCrater);
                                }
                        }
                }

                // Sublimate material every sublimation period:
                if (i%numberOfCratersInSublimationPeriod == 0){
                        grid.sublimateIce();
                }

                // Deposit ice periodically
                if (i%numberOfCratersInDepositionEvent == 0){
                        if (iceEmplacementThickness > 0){ 
                                grid.depositLayer(Layer(iceEmplacementThickness, 0, 1, 0));
                        }
                }

                // Update crater depths, after topography has changed:
                // grid.updateExistingCratersDepth(crater);
                // cratersDepthHistogram.add(crater.finalDepth);

                progressBar(i, totalNumberOfImpactors);

                // Print progress to a file:
                if (i%numberOfCratersInTimestep == 0) {
                        // Threshold slopes
                        grid.thresholdSlopes(angleOfRepose);
                        // Pring z matrix:
                        grid.printSurface(printIndex, false);
                        grid.printIntegratedSubsurface(depthToIntegrate, printIndex);

                        if (isPrintSubsurface) {
                                grid.printSubsurface(printIndex);
                        }

                        // Print to log:
                        char logEntry[100];
                        sprintf(logEntry, "Progress: %0.2f%%.",(double) i/(double) totalNumberOfImpactors * 100);
                        addLogEntry(logEntry, false);
                        printIndex++;
                }

        }
        // Print craters histogram to file:
        addLogEntry("Finished running. Saving histograms and data.", true);
        addLogEntry("Printing crater histogram file.", true);
        cratersHistogram.print("./output/craters_histogram.txt");
        impactorsHistogram.print("./output/impactor_histogram.txt");
        cratersDepthHistogram.print("./output/depth_histogram.txt");
        grid.sublimateIce(); // Sublimate ice one last time
        grid.printSurface(printIndex, true);
        grid.printIntegratedSubsurface(depthToIntegrate, printIndex);
        
        if (isPrintSubsurface) {
                grid.printSubsurface(printIndex);
        }

        addLogEntry("Simulation has ended.", true);
}
