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
#include "../include/utility.hpp"
#include "../include/log.hpp"

#include<alloca.h>
#include<Eigen/Dense>
#include <Eigen/Geometry>


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
  if (reqX < x.front() || reqX > x.back()){
    return 0;
  }

  while (!(reqX >= xPrev & reqX <= xNext)){
    xPrev = x[i]; xNext = x[i+1];
    yPrev = y[i]; yNext = y[i+1];
    i++;
  }

  // Linearly interpolate between points:
  double t = (reqX - xPrev)/(xNext - xPrev);

  return (1 - t) * yPrev + t * yNext;
}


// Rotate vector around x axis
Eigen::Vector3d rotate_vector_y(Eigen::Vector3d v, double alpha) {
  // Create rotation matrix
  Eigen::AngleAxisd rot_y(alpha, Eigen::Vector3d::UnitY());

  // Rotate vector
  Eigen::Vector3d rotated_v = rot_y * v;

  return rotated_v;
}

// Rotate vector around z axis
Eigen::Vector3d rotate_vector_z(Eigen::Vector3d v, double alpha) {
  // Create rotation matrix
  Eigen::AngleAxisd rot_z(alpha, Eigen::Vector3d::UnitZ());

  // Rotate vector
  Eigen::Vector3d rotated_v = rot_z * v;

  return rotated_v;
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

// Read layers file:
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
                double thickness;
                double regolithFrac;
                double iceFrac;
                double sootFrac;
                if (!(iss >> thickness >> regolithFrac >> iceFrac >> sootFrac)) {
                        addLogEntry("Cannot read values from layers file.", true);
                        exit(EXIT_FAILURE);
                }
                tempLayer.push_back(thickness);
                tempLayer.push_back(regolithFrac);
                tempLayer.push_back(iceFrac);
                tempLayer.push_back(sootFrac);

                layersList.push_back(tempLayer);
        }

        return layersList;
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
