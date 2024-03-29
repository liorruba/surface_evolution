// A quick map to store variables:
typedef struct {
        std::string name;
        double value;
} var;

// This file contains various functions to help manage the code.
// Function to check if string is empty:
int isEmpty(const char *s);

// Get current time
uint64_t get_time();

// Generate a uniformly distributed random number:
double randU(double low, double high);

// Creates a linearly spaced vector of length numberOfElements over a range defined by [low,high].
// The first element is the length of the array.
std::vector<double> linspace(double low, double high, int numberOfElements);

// Creates a logarithmic spaced vector of length numberOfElements over a range of expoenents defined by [low,high].
// The first element is the length of the array.
std::vector<double> logspace(double low, double high, int numberOfElements, double base);

// Print a matrix "object" to file
void printMatrixToFile(double ** matrix, int gridSize, char * fileName);

// Linear interpolation:
double linearInterp(std::vector<double> x, std::vector<double> y, double reqX);

// 2-d binning for printing a file with reduced resolution
std::vector<double> bin_1d_vector(const std::vector<double> &input_vector, double bin_resolution);
std::vector< std::vector<double> > bin_2d_vector(const std::vector<std::vector<double>> &input_vector, double bin_resolution);

// Surface normal
std::vector<double> surfNormal(std::vector<double> ptA, std::vector<double> ptB, std::vector<double> ptC);

double vecNorm(std::vector<double> vec);

// Angle between vector and xy plane:
double xyPlaneVecAngle(std::vector<double> vec);

// Cumulative integration
std::vector<double> cumtrapz(const std::vector<double>& x, const std::vector<double>& y);

// Round up to the nearest grid value
int roundUp(int numToRound, int multiple);

// Read config and layers files:
std::vector<var> readConfig();
std::vector< std::vector<double> > readLayers();
std::vector<int8_t> readPixelIndex();
double setVariable(std::vector<var> varList, std::string varName);

// A simple progress bar
void progressBar(long progress, long total);

// Convert 2-d to linear index
int getLinearIndex(int i, int j, int numCols);