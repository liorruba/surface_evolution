#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS
#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING 1

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <cstring>
#include <ctime>
#include <sstream>
#include <fstream>
#include <string>
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
#include "../include/utility.hpp"
#include "../include/tests.hpp"

//////////////////////////////
// DECLARE INPUT PARAMETERS //
//////////////////////////////
// Simulation parameters:
double regionWidth; // km
double resolution; // 1/km
double endTime; // Ma
double printTimeStep; // Time step for printing data in Ma.
double initialThickness; // Initial thickness of subsurface layer.
double latitude; // Latitude (for shadow calculation)
bool isEmplaceEjecta; // Should emplace ejecta? Computationally extensive.
bool isEmplaceSecondaries; // Should emplace ejecta? Computationally extensive.
bool runTests; // Should emplace ejecta? Computationally extensive.
bool isProduceMelt; // Should produce melt?

// Crater formation variables:
double depthToDiameter; // Dimensionless ratio, crater depth to diameter
double rimToDiameter; // Dimensionless ratio, rim height to diameter
double rimDropoffExponent; // The exponent of the rim height decrease power law
double numberOfZModelShells; // Number of shells in z model (for ejecta calc.)
int craterProfileType; // Chosen crater profile
int ejectaSpread; // The spread of the ejecta in crater radii
double slope_secondaries; // The spread of the ejecta in crater radii
double iceDensity; // The density of ice
double regolithDensity; // The density of ice
double c_ice; // The speed of sound in ice
double c_regolith; // The speed of sound in regolith
double s_ice; // slope of the shock / particle velocity relation (ice)
double s_regolith; // slope of the shock / particle velocity relation (regolith)
double meltEnergy; // The melt energy coefficient (Kraus et al., 2011, 3.1.2)
double ice_fraction; // The fraction of ice in the regolith-ice mixture
double porosity; // The impact target porosity
double temperature; // The subsurface constant temperature

// Impactor distribution variables:
double slope_b; // Slope of the impactor CDF
double earthFluxRatioCoefficient; // Moon-earth impactor flux ratio due to cross-section
double minimumImpactorDiameter; // The smallest impactor in the distribution
double fluxConstant_c; // The Flux of impactors > 1 m, Ma^-1 m^-2
double impactorDensity; // The impactor density in kg m^-3
double meanImpactVelocity; // The impact velocity, kg s^-1
double impactAngle;

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

/////////////
//Functions//
/////////////
// Read config file:
std::vector<var> readConfig(){
        std::ifstream configFile;
        configFile.open("./config/config.cfg");
        std::vector<var> varList; // A varlist vector to store variables.
        var tempVar; // Temporary var struct used as buffer.
        addLogEntry("Initializing config file");

        if (!configFile) {
                addLogEntry("Cannot read config file.");
                exit(EXIT_FAILURE);
        }

        // Count how many variables in config file:
        std::string line;
        while (std::getline(configFile, line)) {
                // Check if line is commented out:
                if (line[0] == '/' && line[1] == '/')
                        continue;

                // Check if line is empty:
                if (line.empty())
                        continue;

                std::istringstream iss(line);
                std::string name; double value;
                if (!(iss >> name >> value)) {
                        addLogEntry("Cannot read variable from config file.");
                        exit(EXIT_FAILURE);
                }

                tempVar.name = name;
                tempVar.value = value;
                varList.push_back(tempVar);
        }

        return varList;
}

// Get variable from list:
double setVariable(std::vector<var> varList, std::string varName){
        size_t i;

        for (i = 0; i < varList.size(); i++) {
                if (varList[i].name == varName) {
                        char logEntry[100];
                        sprintf(logEntry, "Getting variable %s with value %f.", varList[i].name.c_str(), varList[i].value);
                        addLogEntry(logEntry);
                        return varList[i].value;
                }
        }
// If variable was not found:
        return -1;
}

////////////////
// START MAIN //
////////////////
int main() {
        // Prepare directories:
        DIR * outputdir = opendir("./output");
        if (outputdir) {
                closedir(outputdir);
        }
        else {
                system("mkdir ./output");
        }

        DIR * logdir = opendir("./log");
        if (logdir) {
                closedir(logdir);
        }
        else {
                system("mkdir ./log");
        }

        // Create log:
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
        isProduceMelt = setVariable(varList, "isProduceMelt");

        // Crater formation variables:
        depthToDiameter = setVariable(varList, "depthToDiameter");
        rimToDiameter = setVariable(varList, "rimToDiameter");
        rimDropoffExponent = setVariable(varList, "rimDropoffExponent");
        numberOfZModelShells = setVariable(varList, "numberOfZModelShells");
        craterProfileType = (int) setVariable(varList, "craterProfileType");
        ejectaSpread = (int) setVariable(varList, "ejectaSpread");
        slope_secondaries = setVariable(varList, "slope_secondaries");
        iceDensity = setVariable(varList, "iceDensity");
        regolithDensity = setVariable(varList, "regolithDensity");
        c_ice = setVariable(varList, "c_ice");
        c_regolith = setVariable(varList, "c_regolith");
        s_ice = setVariable(varList, "s_ice");
        s_regolith = setVariable(varList, "s_regolith");
        meltEnergy = setVariable(varList, "meltEnergy");
        ice_fraction = setVariable(varList, "ice_fraction");
        porosity = setVariable(varList, "porosity");
        temperature = setVariable(varList, "temperature");

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

        // Initialize the random number generator seed:
        srand48((long) get_time());

        ///////////////
        // Run tests //
        ///////////////
        if (runTests) {
                // If tests did not pass:
                if (tests())
                        std::cout << "All tests passed successfully." << std::endl;
                else
                        std::cout << "Not all tests passed successfully." << std::endl;

                return 0;
        }

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
        // Generate grid:
        Grid grid = Grid();

        // Total number of craters to be created:
        long totalNumberOfImpactors = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * endTime * grid.area * earthFluxRatioCoefficient); // total number of impactors to generate larger than minimumDiameter: N/At = cD^-b.
        long numberOfCratersInTimestep = ceil(fluxConstant_c * pow(minimumImpactorDiameter,-slope_b) * printTimeStep * grid.area * earthFluxRatioCoefficient); // number of impactors to generate larger than minimumDiameter: N/At = cD^-b in some time interval.

        // Throw a warning if totalNumberOfImpactors is too high. If the memory
        // taken by the craters exceeds 1 GB, throw a memory warning:
        if (sizeof(Crater) * totalNumberOfImpactors > 1e9)
                std::cout << "WARNING: memory taken by craters exceeds 1 GB" << std::endl;

        // Craters and impactors histograms:
        Histogram cratersHistogram(minimumImpactorDiameter * 10, regionWidth, 20); // Crater histogram from 10*minimumImpactorDiameter to regionWidth meters
        Histogram impactorsHistogram(minimumImpactorDiameter, 1e4, 20); // Impactor histogram from minimumImpactorDiameter m to 10 km
        Histogram cratersDepthHistogram(minimumImpactorDiameter, 1e4, 20); // Impactor histogram from minimumImpactorDiameter m to 10 km

        //////////////////////
        // Start simulation //
        //////////////////////
        addLogEntry("Starting simulation:");
        char logEntry[50];
        sprintf(logEntry, "Number of craters in simulation: %ld.", totalNumberOfImpactors);
        addLogEntry(logEntry);

        for (long i = 0; i < totalNumberOfImpactors; i++) {
                Impactor impactor1(50);
                Impactor impactor2(50);
                Impactor impactor3(300);
                Crater crater1(impactor1, 0, -500);
                // Crater crater2(impactor2, 0, 500);
                // Crater crater3(impactor3, 0, 0);

                grid.formCrater(crater1);
                std::cout << crater1.meltHeight << std::endl;
                grid.emplaceEjecta(crater1);
                grid.updateExistingCratersDepth(crater1);
                //
                // grid.formCrater(crater2);
                // grid.emplaceEjecta(crater2);
                // grid.updateExistingCratersDepth(crater2);

                // grid.formCrater(crater3);
                // grid.emplaceEjecta(crater3);
                // grid.updateExistingCratersDepth(crater3);
                // std::cout << crater1.finalRadius << ", " << crater2.finalRadius << ", " << crater3.finalRadius << "\n";
                // std::cout << crater1.finalDepth << ", " << crater2.finalDepth << ", " << crater3.finalDepth << "\n";

                grid.printExistingCratersToHistogram(20);
                grid.printExistingCraters();

                break;
                // Randomize a new impactor:
                std::cout << i << std::endl;
                Impactor impactor;

                // Add diameter to crater histogram:
                impactorsHistogram.add(2 * impactor.radius);
                // Create a crater instance:
                Crater crater(impactor);
                // Record diameter in crater histogram:
                cratersHistogram.add(2 * crater.finalRadius);

                // Form a crater on the grid:
                grid.formCrater(crater);

                if (crater.finalRadius > 2 * resolution && isEmplaceEjecta) {
                        grid.emplaceEjecta(crater);
                }

                // "Ghost" craters:
                // If the crater exceeds the grid, wrap around it by creating a ghost crater.
                // If a corner:
                if ( (fabs(crater.xLocation) > regionWidth/2 - crater.finalRadius) && (fabs(crater.yLocation) > regionWidth/2 - crater.finalRadius) ) {
                        // Calculate ghost crater location as sgn(x) * (region_width - x);
                        xGhost = (-crater.xLocation/fabs(crater.xLocation)) * (regionWidth - crater.xLocation * (crater.xLocation/fabs(crater.xLocation)));
                        yGhost = (-crater.yLocation/fabs(crater.yLocation)) * (regionWidth - crater.yLocation * (crater.yLocation/fabs(crater.yLocation)));

                        Crater ghost1 = Crater(impactor, xGhost, yGhost);
                        Crater ghost2 = Crater(impactor, xGhost, crater.yLocation);
                        Crater ghost3 = Crater(impactor, crater.xLocation, yGhost);
                        grid.formCrater(ghost1);
                        grid.emplaceEjecta(ghost1);
                        grid.formCrater(ghost2);
                        grid.emplaceEjecta(ghost2);
                        grid.formCrater(ghost3);
                        grid.emplaceEjecta(ghost3);
                }

                // If a side:
                if ( (fabs(crater.xLocation) > regionWidth/2 - crater.finalRadius/2) && (fabs(crater.yLocation) < regionWidth/2 - crater.finalRadius/2)) {
                        // Calculate ghost crater location as sgn(x) * (region_width - x);
                        xGhost = (-crater.xLocation/fabs(crater.xLocation)) * (regionWidth - crater.xLocation * (crater.xLocation/fabs(crater.xLocation)));
                        Crater ghost2 = Crater(impactor, xGhost, crater.yLocation);

                        // Create one ghost crater:
                        grid.formCrater(ghost2);
                        grid.emplaceEjecta(ghost2);
                }

                // If another side:
                if ( (fabs(crater.xLocation) < regionWidth/2 - crater.finalRadius/2) && (fabs(crater.yLocation) > regionWidth/2 - crater.finalRadius/2)) {
                        // Calculate ghost crater location as sgn(x) * (region_width - x);
                        yGhost = (-crater.yLocation/fabs(crater.yLocation)) * (regionWidth - crater.yLocation * (crater.yLocation/fabs(crater.yLocation)));
                        Crater ghost3 = Crater(impactor, crater.xLocation, yGhost);

                        // Create one ghost crater:
                        grid.formCrater(ghost3);
                        grid.emplaceEjecta(ghost3);
                }

                // Secondary craters:
                // if (isEmplaceSecondaries) {
                //         // Form a secondary crater within 2 crater diameters from the primary:
                //         if (crater.numberOfSecondaries > 0) {
                //                 char logEntry[50];
                //                 sprintf(logEntry, "Primary diameter: %f. Number of secondaries: %d.", 2*crater.finalRadius, crater.numberOfSecondaries);
                //                 addLogEntry(logEntry);
                //         }
                //         for (long j = 0; j < crater.numberOfSecondaries; j++) {
                //                 double secondaryxLocation = randU(crater.xLocation - 4 * crater.finalRadius, crater.xLocation + 4 * crater.finalRadius);
                //                 double secondaryyLocation = randU(crater.yLocation - 4 * crater.finalRadius, crater.yLocation + 4 * crater.finalRadius);
                //                 double secondaryRadius = resolution * pow(randU(0,1), -1/slope_secondaries); // Set impactor radius from the cumulative distribution, meters
                //                 Crater secondaryCrater(secondaryxLocation, secondaryyLocation, secondaryRadius);
                //                 grid.formCrater(secondaryCrater);
                //         }
                // }

                // // Update crater depths, after topography has changed:
                // grid.updateExistingCratersDepth(crater);
                // cratersDepthHistogram.add(crater.finalDepth);

                // Print progress to a file:
                if (i%numberOfCratersInTimestep == 0) {
                        // Pring z matrix:
                        grid.printSurface(printIndex);

                        // Print to log:
                        char logEntry[100];
                        sprintf(logEntry, "Progress: %0.2f%%.",(double) i/(double) totalNumberOfImpactors * 100);
                        addLogEntry(logEntry);
                        printIndex++;
                }

        }
        // Print craters histogram to file:
        addLogEntry("Printing crater histogram file.");
        cratersHistogram.print("./output/craters_histogram.txt");
        impactorsHistogram.print("./output/impactor_histogram.txt");
        cratersDepthHistogram.print("./output/depth_histogram.txt");
        grid.printSurface(0);
        grid.printGrid(0);

        addLogEntry("Simulation has ended.");
}
