#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>
#include<time.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<sys/types.h>
#include<inttypes.h>
#include "../include/log.h"
#include "../include/grid.h"
#include "../include/histogram.h"
#include "../include/utility.h"

// Structs
// A quick map to store variables:
typedef struct{
  char name[50];
  double value;
} var;

typedef struct{
  var * vars;
  int numberOfVars;
} varlist;

/////////////
//Constants//
/////////////
// Simulation parameters:
double regionWidth;  // km
double resolution; // km
double endTime; // Ma
double depthToDiameter;
double minimumDiameter; // Crater minimum diameter, km.
double printTimeStep; // Time step for printing data in Ma.

// Crater production constants:
double b = 3.209;
double c = 1.467e-6;

/////////////
//Functions//
/////////////
// Read config file:
varlist readConfig(){
  FILE * configFile = fopen("config/config.cfg","r");
  varlist varList;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  int numberOfVars = 0;
  int i = 0; // Counter to count the line currently iterating.

  addLogEntry("Initializing config file");

  if (configFile == NULL){
    addLogEntry("Cannot read config file.");
    exit(EXIT_FAILURE);
  }

  // Count how many variables in config file:
  while ((read = getline(&line, &len, configFile)) != -1) {
    // Check if line is commented out:
    if (line[0] == '/' && line[1] == '/')
      continue;

    // Check if line is empty:
    if (isEmpty(line))
      continue;

    numberOfVars++;
  }

  varList.numberOfVars = numberOfVars;
  char logEntry[100]; sprintf(logEntry, "Found %d variables in configuration file.", numberOfVars);
  addLogEntry(logEntry);

  // Read variables into hash:
  rewind(configFile);
  varList.vars = (var *)malloc(sizeof(var) * numberOfVars);
  while ((read = getline(&line, &len, configFile)) != -1) {
    // Check if line is commented out:
    if (line[0] == '/' && line[1] == '/')
      continue;

    // Check if line is empty:
    if (isEmpty(line))
      continue;

    sscanf(line, "%s %lf\n", varList.vars[i].name, &varList.vars[i].value);
    i++;
  }

  return varList;
}

// Get variable from list:
double setVariable(varlist varList, char * varName){
  int i;

  for (i = 0; i < varList.numberOfVars; i++) {
    if (!strcmp(varList.vars[i].name, varName)){
      char logEntry[100];
      sprintf(logEntry, "Getting variable %s with value %f.", varList.vars[i].name, varList.vars[i].value);
      addLogEntry(logEntry);
      break;
    }
  }

  return varList.vars[i].value;
}

void createCraterInZmat(int gridSize, double *craterDiameter, double xPosition, double yPosition, double **xmat, double **ymat, double **zmat){
	// double tic = get_time();
	double distFromCraterCenter = 0; // For later use.
	double newDepth = 0; // The new crater depth.

	double craterRadius = *craterDiameter/2; // Calculate the crater radius

	// Find approximate inital i and j:
	int iInit = floor( ((xPosition + regionWidth/2) - 2 * craterRadius) / resolution ); if(iInit < 0) iInit = 0;
	int iFinal = ceil( ((xPosition + regionWidth/2) + 2 * craterRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
	int jInit = floor( ((yPosition + regionWidth/2) - 2 * craterRadius) / resolution ); if(jInit < 0) jInit = 0;
	int jFinal = ceil( ((yPosition + regionWidth/2) + 2 * craterRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;
	// int iInit = 0, iFinal = gridSize, jInit = 0, jFinal = gridSize;

	for (int i = iInit; i < iFinal; i++){
		for (int j = jInit; j < jFinal; j++){
			distFromCraterCenter = sqrt(pow(xmat[j][i] - xPosition, 2) + pow(ymat[j][i] - yPosition, 2));
      // printf("Pos: %f, %f %f %f\n",xPosition, yPosition, xmat[i][j], ymat[i][j]);
			// if inside the crater:
			if ( distFromCraterCenter <= craterRadius ){
				// printf("Pos: %f, %f; idx: %d %d %d %d; diam: %f\n",xPosition, yPosition, iInit, iFinal, jInit, jFinal,craterRadius*2);
				newDepth = 2 * depthToDiameter / craterRadius * pow(distFromCraterCenter,2) - 2 * depthToDiameter * craterRadius;

				if (newDepth < zmat[i][j]){
						zmat[i][j] = newDepth;
				}
			}
		}
	}
}

////////////////
// START MAIN //
////////////////

int main() {
  // Prepare directories:
  system("mkdir ./output");
  system("mkdir ./log");

  // Create log:
  createLogFile();

  // Read config:
  varlist varList = readConfig();

  // Set variables:
  regionWidth = setVariable(varList, "regionWidth");
  resolution = setVariable(varList, "resolution");
  endTime = setVariable(varList, "endTime");
  depthToDiameter = setVariable(varList, "depthToDiameter");
  minimumDiameter = setVariable(varList, "minimumDiameter");
  printTimeStep = setVariable(varList, "printTimeStep");

  // Crater production constants:
  b = setVariable(varList, "b");
  c = setVariable(varList, "c");

	// Initialize the random number generator seed:
	srand48((long) get_time());
  // Declare crater related varialbes:
  // (DO NO TOUCH)
  double diameter;
  double xPosition;
  double yPosition;
  double quantile;
  int zMatrixIndex = 1; // The index to add to the end of the output matrix file (zmat_1.txt).

	// Generate grid:
	// Create x and y equally spaced vectors:
	int gridSize = regionWidth/resolution;
	int lengthx = gridSize;
	int lengthy = gridSize;
	double * x = linspace(-regionWidth/2, regionWidth/2, lengthx);
	double * y = linspace(-regionWidth/2, regionWidth/2, lengthy);
  // Create meshgrid matrices:
	double ** xmat = meshgridX(x,y,lengthx,lengthy); // Create xmat, the x coordinates matrix.
	double ** ymat = meshgridY(x,y,lengthx,lengthy); // Create ymat, the y coordinates matrix.
	double ** zmat = zeros(lengthx,lengthy); // Initialize zmat, the height matrix.
	free(x);
	free(y);

	// Simulated area:
	double area = pow(regionWidth, 2);

  // Total number of craters to be created:
	long totalNumberOfCraters = ceil(c * pow(minimumDiameter,-b) * endTime * area); // total number of craters to generate larger than minimumDiameter: N/At = cD^-b.
  long numberOfCratersInTimestep = ceil(c * pow(minimumDiameter,-b) * printTimeStep * area); // number of craters to generate larger than minimumDiameter: N/At = cD^-b in 1e4 years.
  // Craters histogram:
  histogram cratersHistogram = createHistogram(minimumDiameter, regionWidth, 10);

  // Start simulation:
  addLogEntry("Starting simulation:");
	for (long i = 0; i < totalNumberOfCraters; i++){

    xPosition = randU(-regionWidth/2, regionWidth/2); // Randomize the x position of the crater.
    yPosition = randU(-regionWidth/2, regionWidth/2); // Randomize the y position of the crater.
    quantile = randU(0,1);
    diameter = minimumDiameter * pow(quantile, -1/b);

    cratersHistogram = addToHistogram(&diameter, cratersHistogram);
		createCraterInZmat(gridSize, &diameter, xPosition, yPosition, xmat, ymat, zmat);

    // Print progress to a file:
    if (i%numberOfCratersInTimestep == 0){
      // Pring z matrix:
      char fileName[50];
      sprintf(fileName, "./output/zmat_%d.txt",zMatrixIndex); zMatrixIndex++;
      printMatrixToFile(zmat, gridSize, fileName);

      // Print to log:
      char logEntry[100];
      sprintf(logEntry, "Progress:%f.",(double) i/(double) totalNumberOfCraters);
      addLogEntry(logEntry);
    }
	}

	// Print x and y matrices to file:
  printMatrixToFile(xmat, gridSize, "./output/xmat.txt");
  printMatrixToFile(ymat, gridSize, "./output/ymat.txt");
  printMatrixToFile(zmat, gridSize, "./output/zmat.txt");

  // Print craters histogram to file:
  addLogEntry("Creating histogram file.");
  FILE *cratersHistogramFile;
  cratersHistogramFile = fopen("./output/craters_histogram.txt","w+");
  if (cratersHistogramFile == NULL){
    addLogEntry("Cannot create histogram file.");
  }
  else{
    addLogEntry("Histogram file successfully created.");
  }

  for (int i = 0; i < cratersHistogram.length; i++){
    fprintf(cratersHistogramFile, "%f\t",cratersHistogram.bins[i]);
    fprintf(cratersHistogramFile, "%d\n",cratersHistogram.counts[i]);
  }
  fclose(cratersHistogramFile);

	// Free memory
  addLogEntry("Freeing memory.");
	free(xmat);
	free(ymat);
	free(zmat);
  addLogEntry("Simulation has ended.");
}
