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
#include <inttypes.h>

// Structs
// A histogram is a 2-D array in which the first column is bins and the second column is counts.
typedef struct{
  double * bins;
  int * counts;
  int length;
} histogram;

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
uint64_t get_time()
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e9 + t.tv_nsec;
}

// Log file:
void createLogFile(){
  FILE * logFile = fopen("./log/log.txt","w+"); // Create the log file. Remove if exists.
  if (logFile == NULL){
    printf("ERROR: Cannot create log file.");
    exit(EXIT_FAILURE);
  }
  else{
    fclose(logFile);
  }
}

// Function to check if string is empty:
int isEmpty(const char *s) {
  while (*s != '\0') {
    if (!isspace((unsigned char)*s))
      return 0;
    s++;
  }
  return 1;
}

// Add a new log entry
void addLogEntry(char *str){
  // Logging variables:
  time_t now = time(0);
  struct tm * timenow = localtime(&now);
  char timeInString[50];

  FILE * logFile = fopen("./log/log.txt","a+"); // Append to log file.
  strftime(timeInString, sizeof(timeInString), "%Y-%m-%d %H:%M:%S", timenow); // Format time string.

  fprintf(logFile, "%s\t%s\n", timeInString, str); // Write to log file.
  fclose(logFile);
}

// Read config file:
varlist readConfig(){
  FILE * configFile = fopen("./config.cfg","r");
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

// Generate a uniformly distributed random number:
double randU(double low, double high)
{
	return (high - low) * drand48() + low;
}

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
double * linspace(double low, double high, int numberOfElements) {
	double *linVec = (double *)malloc(sizeof(double) * (numberOfElements));
	double res = (high - low) / (double) (numberOfElements - 1);

	for (int i = 0; i < numberOfElements; i++) {
		linVec[i] = low + i * res;
	}

	return linVec;
}

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
double * logspace(double low, double high, int numberOfElements, double base) {
  double * logVec = (double *)malloc(sizeof(double) * (numberOfElements));
  double * exponents = linspace(low, high, numberOfElements);

  for (int i = 0; i < numberOfElements; i++){
    logVec[i] = pow(base, exponents[i]);
  }

  return logVec;
}

// This function creates an array xmat to serve as the x coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridX(double *x, double *y, int lengthx, int lengthy){
	// Allocate memory to the array:
	double **xmat = (double **)malloc(sizeof(double) * (lengthy));
	for (int i = 0; i < lengthy; i++)
		xmat[i] = (double *)malloc(sizeof(double) * (lengthx));

	// Populate the array:
	for (int i = 0; i < lengthx; i++) {
		for (int j = 0; j < lengthy; j++) {
			xmat[i][j] = x[j];
		}
	}

	return xmat;
}

// This function creates an array ymat to serve as the y coordinate of a 2D grid.
// The first element of the input array must be the array length.
double ** meshgridY(double *x, double *y, int lengthx, int lengthy){
	// Allocate memory to the array:
	double **ymat = (double **)malloc(sizeof(double) * (lengthx));
	for (int i = 0; i < lengthx; i++)
		ymat[i] = (double *)malloc(sizeof(double) * (lengthy));

	// Populate the array:
	for (int i = 0; i < lengthx; i++) {
		 for (int j = 0; j < lengthy; j++) {
			ymat[j][i] = y[j];
		}
	}

	return ymat;
}

// This function creates a matrix populated by zeros, whose size is nxm.
double ** zeros(int n, int m){
	double **zmat = (double **)malloc(sizeof(double) * (n));
	for (int i = 0; i < n; i++)
		zmat[i] = (double *)malloc(sizeof(double) * (m));

	// Populate with zeros:
	for (int j = 0; j < n; j++) {
		 for (int i = 0; i < m; i++) {
			zmat[i][j] = 0;
		}
	}

	return zmat;
}

histogram createHistogram(double minBin, double maxBin, int numOfBins){
  histogram hist;

  // Calulcate the exponent for the logvec function via change of base:
  double minBinExponent = log(minBin)/log(sqrt(2));
  double maxBinExponent = log(maxBin)/log(sqrt(2));

  // Create the bins and vals arrays:
  double * bins = logspace(minBinExponent, maxBinExponent, numOfBins, sqrt(2));

  // Set the histogram length:
  hist.length = numOfBins;
  // Allocate memory to the data array in histogram:
  hist.bins = (double *)malloc(sizeof(double) * (numOfBins));
  hist.counts = (int *)malloc(sizeof(int) * (numOfBins));

  for (int i = 0; i < numOfBins; i++){
    hist.bins[i] = bins[i];
    hist.counts[i] = 0;
  }

  return hist;
}

// This function recieves a value and an existing histogram and returns a histogram.
histogram addToHistogram(double *val, histogram hist){
  if ((*val < hist.bins[0]) || (*val > hist.bins[hist.length - 1])){
    // TODO: INCREASE HISTOGRAM SIZE.
  }
  else {
    // Iterate over histogram:
    for (int i = 1; i < hist.length; i++){
      // If in bin:
      if (*val < hist.bins[i]){
        hist.counts[i-1]++;
        break;
      }
    }
  }
  return hist;
}

// Generate list of crater diameters, write to file:
double * getDiameterList(int totalNumberOfCraters, double minimumDiameter){
	double * diameterList = (double *)malloc(sizeof(double) * (totalNumberOfCraters));
	double quantile;

	// Generate a file and open for writing:
	FILE * fileCraterList;
	fileCraterList = fopen("./diameterList.txt","w+");

	// Calculate the probability to find a crater with n>N. Choose n from a uniform distribution and calculate the quantile:
	for (int i = 0; i < totalNumberOfCraters; i++){
		quantile = randU(0,1);
		diameterList[i] = minimumDiameter * pow(quantile, -1/b);
	}

  //Write to a file:
  for (int i = 0; i < totalNumberOfCraters; i++){
    fprintf(fileCraterList, "%f\n", diameterList[i]);
  }

	fclose(fileCraterList);
	return diameterList;
}

void createCraterInZmat(int gridSize, double *craterDiameter, double **xmat, double **ymat, double **zmat){
	// double tic = get_time();
	double distFromCraterCenter = 0; // For later use.
	double newDepth = 0; // The new crater depth.

	double xPosition = randU(-regionWidth/2, regionWidth/2); // Randomize the x position of the crater.
	double yPosition = randU(-regionWidth/2, regionWidth/2); // Randomize the y position of the crater.

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

void printMatrixToFile(double ** matrix, int gridSize, char * fileName){
  FILE * file = fopen(fileName,"w+");

	for (int i = 0; i < gridSize; i++){
		for (int j = 0; j < gridSize; j++){
			fprintf(file, "%f ",matrix[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);
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
  double quantile;
  double p99; // 99th percentile.
  double median; // median depth
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
    quantile = randU(0,1);
    diameter = minimumDiameter * pow(quantile, -1/b);
    cratersHistogram = addToHistogram(&diameter, cratersHistogram);
		createCraterInZmat(gridSize, &diameter, xmat, ymat, zmat);

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
