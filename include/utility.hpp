#pragma once
#include <cstdint>
#include <string>
#include <vector>

// A quick map to store variables:
typedef struct {
        std::string name;
        double value;
} var;

// Get current time
uint64_t get_time();

// Generate a uniformly distributed random number:
double randU(double low, double high);

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
std::vector<double> linspace(double low, double high, int numberOfElements);

// Creates a logarithmic spaced vector of length numberOfElements over a range of exponents defined by [low,high].
std::vector<double> logspace(double low, double high, int numberOfElements, double base);

// Linear interpolation of the table y(x) at reqX. x must be ascending; returns 0 outside the table.
double linearInterp(const std::vector<double> &x, const std::vector<double> &y, double reqX);

// Average pooling for printing files with reduced resolution. Trailing cells that do not fill a
// whole bin are dropped, so both functions shrink a dimension of n cells to floor(n / binSize).
std::vector<double> bin_1d_vector(const std::vector<double> &input_vector, double bin_resolution);
std::vector< std::vector<double> > bin_2d_vector(const std::vector<std::vector<double>> &input_vector, double bin_resolution);

double vecNorm(std::vector<double> vec);

// Angle between vector and xy plane:
double xyPlaneVecAngle(std::vector<double> vec);

// Cumulative integration
std::vector<double> cumtrapz(const std::vector<double>& x, const std::vector<double>& y);

// Read config and layers files:
std::vector<var> readConfig();
std::vector< std::vector<double> > readLayers();
// Reads config/pixelIndex.cfg; returns an empty vector when the file does not exist (uniform grid).
std::vector<int8_t> readPixelIndex();
// Returns the value of a config variable; exits with an error if it is missing.
double setVariable(std::vector<var> varList, std::string varName);
// Returns the value of a config variable, or defaultValue (with a log entry) if it is missing.
double setVariableOptional(std::vector<var> varList, std::string varName, double defaultValue);

// A simple progress bar
void progressBar(long progress, long total);

// Convert 2-d to linear index
int getLinearIndex(int i, int j, int numCols);

// Zero-padded index used in output file names, wide enough for the last time step.
std::string formatOutputIndex(int index);
