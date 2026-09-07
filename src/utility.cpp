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
#include <tuple>
#include <algorithm>
#include <stdexcept>
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
std::vector<double> linspace(double low, double high, int numberOfElements) {
	if (numberOfElements <= 0) {
		return std::vector<double>();
	}
	if (numberOfElements == 1) {
		return std::vector<double>(1, low);
	}

	std::vector<double> linVec(numberOfElements);
	double res = (high - low) / (double) (numberOfElements - 1);

	for (int i = 0; i < numberOfElements; i++) {
		linVec[i] = low + i * res;
	}

	return linVec;
}

// Creates a logarithmic spaced vector of length numberOfElements over a range of exponents defined by [low,high].
std::vector<double> logspace(double low, double high, int numberOfElements, double base) {
	std::vector<double> logVec(numberOfElements);
	std::vector<double> exponents = linspace(low, high, numberOfElements);

	for (int i = 0; i < numberOfElements; i++){
		logVec[i] = pow(base, exponents[i]);
	}

	return logVec;
}

// A simple progress bar:
void progressBar(long progress, long total) {
	if (total <= 1) {
		return;
	}

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

// Linear interpolation. Evaluates y(reqX) from the table y(x), x ascending. Returns 0 outside the table.
double linearInterp(const std::vector<double> &x, const std::vector<double> &y, double reqX){
	if (x.size() < 2 || y.size() != x.size()) {
		return 0;
	}
	if ((reqX < x.front()) || (reqX > x.back())) {
		return 0;
	}

	// First table point strictly to the right of reqX:
	size_t k = std::upper_bound(x.begin(), x.end(), reqX) - x.begin();
	if (k == 0) {
		k = 1;
	}
	if (k >= x.size()) {
		k = x.size() - 1;
	}

	const double dx = x[k] - x[k - 1];
	if (dx <= 0) {
		return y[k];
	}
	const double t = (reqX - x[k - 1]) / dx;

	return (1 - t) * y[k - 1] + t * y[k];
}

// Number of cells per bin for the requested output resolution (1 = no downsampling).
static int binSizeFor(double bin_resolution) {
	if (bin_resolution < resolution) {
		throw std::runtime_error("The downsampling resolution cannot be finer than the grid resolution.");
	}
	return std::max(1, (int) std::lround(bin_resolution / resolution));
}

// Simple 1-d grid with average pooling
std::vector<double> bin_1d_vector(const std::vector<double> &input_vector, double bin_resolution) {
	const int bin_size = binSizeFor(bin_resolution);
	const int rows = input_vector.size();
	const int new_rows = rows / bin_size;
	std::vector<double> output_vector(new_rows, 0);

	for (int i = 0; i < new_rows; ++i) {
		for (int ii = i * bin_size; ii < (i + 1) * bin_size; ++ii) {
			output_vector[i] += input_vector[ii];
		}
		output_vector[i] /= bin_size;
	}

	return output_vector;
}

// Simple 2-d grid with average pooling
std::vector< std::vector<double> > bin_2d_vector(const std::vector< std::vector<double> > &input_vector, double bin_resolution) {
	const int bin_size = binSizeFor(bin_resolution);
	const int rows = input_vector.size();
	const int cols = rows > 0 ? (int) input_vector[0].size() : 0;
	const int new_rows = rows / bin_size;
	const int new_cols = cols / bin_size;
	std::vector<std::vector<double>> output_vector(new_rows, std::vector<double>(new_cols, 0));

	for (int i = 0; i < new_rows; ++i) {
		for (int j = 0; j < new_cols; ++j) {
			for (int ii = i * bin_size; ii < (i + 1) * bin_size; ++ii) {
				for (int jj = j * bin_size; jj < (j + 1) * bin_size; ++jj) {
					output_vector[i][j] += input_vector[ii][jj];
				}
			}
			output_vector[i][j] /= (double) bin_size * bin_size;
		}
	}

	return output_vector;
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
		if (line.compare(0, 2, "//") == 0)
			continue;

		// Check if line is empty:
		if (line.find_first_not_of(" \t\r") == std::string::npos)
			continue;

		// Read strings as stream
		std::istringstream iss(line);
		// Define buffer vars to store var name and value:
		std::string name; 
		double value;
		if (!(iss >> name >> value)) {
			addLogEntry("Cannot read variable from config file: " + line, true);
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
		if (line.compare(0, 2, "//") == 0)
			continue;

		// Check if line is empty:
		if (line.find_first_not_of(" \t\r") == std::string::npos)
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
			addLogEntry("Cannot read values from layers file: " + line, true);
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
// This function computes the cumulative integral of the input vector using the trapezoid method
//
std::vector<double> cumtrapz(const std::vector<double>& x, const std::vector<double>& y) {
	// Allocate a vector to hold the cumulative integral values
	int sz = y.size();
	std::vector<double> result(sz);

	// Calculate the step size
	double step = x[1] - x[0];

	// Calculate the first element of the result vector
	result[0] = 0.0;

	// Calculate the remaining elements of the result vector
	for (int i = 1; i < sz; i++) {
		result[i] = result[i - 1] + (y[i - 1] + y[i]) / 2.0 * step;
	}

	return result;
}

//
// The pixel index matrix is used to represent spatial patterns 
// in the initial conditions of the input surface (e.g. a cold spot)
//
// This function reads input values to the pixel index matrix,
// which is represented as a 1-d vector of 8 bit ints, with size 
// (1, gridSize^2). If the file does not exist an empty vector is
// returned and the grid is initialized uniformly.
//
std::vector<int8_t> readPixelIndex(){
	int gridSize = (int) std::lround(regionWidth / resolution);
	std::ifstream pxIdxFile("./config/pixelIndex.cfg", std::ios::binary);
	std::vector<int8_t> pixelsIndexMatrix;

	// If the pixel index file doesn't exist, the grid is uniform:
	if (!pxIdxFile) {
		addLogEntry("Cannot read pixel indexes file. Creating a uniform subsurface grid...", true);
		return pixelsIndexMatrix;
	}

	// If it exists, read the file:
	int8_t value;
	addLogEntry("Reading pixel index matrix...", true);
	// Read file into pixelsIndexMatrix (represented as a 1-d vector):  
	while (pxIdxFile.read((char*) &value, sizeof(int8_t))) {
		pixelsIndexMatrix.push_back(value);
	}

	// Check the file has the same size as set by the user:
	if ((long) pixelsIndexMatrix.size() != (long) gridSize * gridSize) {
		addLogEntry("Size of the input pixel matrix must equal the grid size, given by (regionWidth/resolution)^2. Terminating.", true);
		throw std::runtime_error("Size of the input pixel matrix must equal the grid size, given by (regionWidth/resolution)^2.");
	}

	pxIdxFile.close();
	return pixelsIndexMatrix;
}

// Get variable from list:
double setVariable(std::vector<var> varList, std::string varName){
	size_t i;

	for (i = 0; i < varList.size(); i++) {
		if (varList[i].name == varName) {
			addLogEntry("Getting variable " + varList[i].name + " with value " + std::to_string(varList[i].value) + ".", false);
			return varList[i].value;
		}
	}

	// If variable was not found, stop: a silent default would corrupt the run.
	addLogEntry("ERROR: variable '" + varName + "' is missing from config/config.cfg.", true);
	exit(EXIT_FAILURE);
}

// Get an optional variable from list:
double setVariableOptional(std::vector<var> varList, std::string varName, double defaultValue){
	for (size_t i = 0; i < varList.size(); i++) {
		if (varList[i].name == varName) {
			addLogEntry("Getting variable " + varList[i].name + " with value " + std::to_string(varList[i].value) + ".", false);
			return varList[i].value;
		}
	}
	addLogEntry("Variable '" + varName + "' is not in config/config.cfg; using the default " + std::to_string(defaultValue) + ".", true);
	return defaultValue;
}

// Convert 2-d to linear index
int getLinearIndex(int i, int j, int numCols) {
	return i * numCols + j;
}

// Zero-padded index for output file names. The width is taken from the largest index that can
// occur (the number of print steps plus the final print), so the padding never underflows.
std::string formatOutputIndex(int index) {
	const int maxIndex = (int) std::ceil(endTime / printTimeStep) + 2;
	const size_t width = std::to_string(maxIndex).length();
	std::string index_str = std::to_string(index);
	if (index_str.length() < width) {
		index_str.insert(0, width - index_str.length(), '0');
	}
	return index_str;
}
