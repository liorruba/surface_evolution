// This file contains various functions to help manage the code.
#define _XOPEN_SOURCE 700
#define __STDC_FORMAT_MACROS
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cctype>
#include<ctime>
#include<cinttypes>
#include<vector>
#include<fstream>
#include<sstream>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"


// Get current time
uint64_t get_time()
{
      struct timespec t;
      clock_gettime(CLOCK_MONOTONIC, &t);
      return t.tv_sec * 1e9 + t.tv_nsec;
}

// Generate a uniformly distributed random number:
double randU(double low, double high)
{
	return (high - low) * drand48() + low;
}

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
std::vector<double> linspace(double low, double high, int numberOfElements) {
	std::vector<double> linVec(numberOfElements);
	double res = (high - low) / (double) (numberOfElements - 1);

	for (int i = 0; i < numberOfElements; i++) {
		linVec[i] = low + i * res;
	}

	return linVec;
}

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
std::vector<double> logspace(double low, double high, int numberOfElements, double base) {
      std::vector<double> logVec(numberOfElements);
      std::vector<double> exponents = linspace(low, high, numberOfElements);

      for (int i = 0; i < numberOfElements; i++){
            logVec[i] = pow(base, exponents[i]);
    }

    return logVec;
}

// A simple progress bar 
void progressBar(long progress, long total) {
        std::string bar;

        // Set the bar resolution
        long res = 50;

        // Prevent division by zero:
        if (total < res)
                res = total - 1;

        for (long i = 0; i < res; i++) {

                if (i < (progress / (total / res)))
                        bar += "=";
                else
                        bar += " ";
        }
        
        std::cout << "\r" "[" << bar << "] " << 100 * progress / total << "%";
        std::fflush(stdout);
        if (progress == total - 1)
                std::cout << std::endl;
}

// Linear interpolation. Evaluates y(reqX) by interpolation the vector y(x).
double linearInterp(std::vector<double> x, std::vector<double> y, double reqX){
      int i = 0;
      double xPrev = x[0], xNext = x[1], yPrev = y[0], yNext = y[1];

  // Check if in range:
      if ((reqX < x.front()) || (reqX > x.back())){
            return 0;
    }

    while (!((reqX >= xPrev) & (reqX <= xNext))){
            xPrev = x[i]; xNext = x[i+1];
            yPrev = y[i]; yNext = y[i+1];
            i++;
    }

  // Linearly interpolate between points:
    double t = (reqX - xPrev)/(xNext - xPrev);

    return (1 - t) * yPrev + t * yNext;
}


// Calculate the norm of a vector
double vecNorm(std::vector<double> vec){
      return sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));
}

// Angle between xy plane and vector
double xyPlaneVecAngle(std::vector<double> vec){
      return M_PI/2 - asin(fabs(vec[2])/vecNorm(vec));
}

////////////
// Functions for reading configuation files
////////////
// Read config file:
std::vector<var> readConfig(){
        std::ifstream configFile;
        configFile.open("./config/config.cfg");
        std::vector<var> varList; // A varlist vector to store variables.
        var tempVar; // Temporary var struct used as buffer.
        addLogEntry("Initializing config file", true);

        if (!configFile) {
                addLogEntry("Cannot read config file.", true);
                exit(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(configFile, line)) {
                // Check if line is commented out:
                if (line[0] == '/' && line[1] == '/')
                        continue;

                // Check if line is empty:
                if (line.empty())
                        continue;

                // Read strings as stream
                std::istringstream iss(line);
                // Define buffer vars to store var name and value:
                std::string name; 
                double value;
                if (!(iss >> name >> value)) {
                        addLogEntry("Cannot read variable from config file.", true);
                        exit(EXIT_FAILURE);
                }

                tempVar.name = name;
                tempVar.value = value;
                varList.push_back(tempVar);
        }

        return varList;
}


// This function reads the layers file and stores the user input 
// layers into a vector, whose elements are the pixel index, layer thickness, 
// regolith fraction, ice fraction and soot fraction
std::vector< std::vector<double> > readLayers(){
        std::ifstream layersFile;
        std::vector< std::vector<double> > layersList;

        layersFile.open("./config/layers.cfg");
        addLogEntry("Initializing layers file", true);

        if (!layersFile) {
                addLogEntry("Cannot read layers file.", true);
                exit(EXIT_FAILURE);
        }

        std::string line;
        while (std::getline(layersFile, line)) {
                // Check if line is commented out:
                if (line[0] == '/' && line[1] == '/')
                        continue;

                // Check if line is empty:
                if (line.empty())
                        continue;

                // Read strings as stream
                std::istringstream iss(line);

                std::vector<double> tempLayer; // Temporary layer used as buffer.
                int pixelIndex;
                double thickness;
                double regolithFrac;
                double iceFrac;
                double sootFrac;

                // First, read layer index:
                if (!(iss >> pixelIndex >> thickness >> regolithFrac >> iceFrac >> sootFrac)) {
                        addLogEntry("Cannot read values from layers file.", true);
                        exit(EXIT_FAILURE);
                }
                tempLayer.push_back(pixelIndex);
                tempLayer.push_back(thickness);
                tempLayer.push_back(regolithFrac);
                tempLayer.push_back(iceFrac);
                tempLayer.push_back(sootFrac);

                layersList.push_back(tempLayer);
        }

        return layersList;
}

//
// The pixel index matrix is used to represent spatial patterns 
// in the initial conditions of the input surface (e.g. a cold spot)
//
// This function reads input values to the pixel index matrix,
// which is represented as a 1-d vector of 8 bit ints, with size 
// (1, gridSize^2).
//
std::vector<int8_t> readPixelIndex(){
        int gridSize = regionWidth / resolution;
        std::ifstream pxIdxFile("./config/pixelIndex.cfg", std::ios::binary);
        std::vector<int8_t> pixelsIndexMatrix;

        pxIdxFile.open("./config/pxIdxFile.cfg");

        // If the pixel index file doesn't exist, create a matrix of 1's
        if (!pxIdxFile) {
                addLogEntry("Cannot read pixel indexes file. Creating a uniform subsurface grid...", true);

                for (int i = 0; i < pow(gridSize, 2); ++i) {
                        pixelsIndexMatrix.push_back(1);

                        progressBar(i, gridSize);
                }
        }
        // If it exists, read the file:
        else {
                int8_t value;

                // Read file into pixelsIndexMatrix (represented as a 1-d vector):  
                int sz = 0;
                while (pxIdxFile.read((char*) &value, sizeof(int8_t))) {
                        pixelsIndexMatrix.push_back(value);
                        sz++;
                }
                std::cout << sz << std::endl;

                // Check the file has the same size as set by the user:
                if (gridSize != (int) sqrt(sz)) {
                        addLogEntry("Size of the input pixel matrix must equal the grid size, given by regionWidth/resolution. Terminating.", false);
                        throw std::runtime_error("Size of the input pixel matrix must equal the grid size, given by regionWidth/resolution.");
                }
        }

        pxIdxFile.close();
        return pixelsIndexMatrix;
}

// Get variable from list:
double setVariable(std::vector<var> varList, std::string varName){
        size_t i;

        for (i = 0; i < varList.size(); i++) {
                if (varList[i].name == varName) {
                        char logEntry[100];
                        sprintf(logEntry, "Getting variable %s with value %f.", varList[i].name.c_str(), varList[i].value);
                        addLogEntry(logEntry, false);
                        return varList[i].value;
                }
        }
// If variable was not found:
        return -1;
}

// Convert 2-d to linear index
int getLinearIndex(int i, int j, int numCols) {
      return i * numCols + j;
}
