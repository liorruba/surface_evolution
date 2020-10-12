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

// Surface normal
std::vector<double> surfNormal(std::vector<double> ptA, std::vector<double> ptB, std::vector<double> ptC);

double vecNorm(std::vector<double> vec);

// Angle between vector and xy plane:
double xyPlaneVecAngle(std::vector<double> vec);

// Round up to the nearest grid value
int roundUp(int numToRound, int multiple);
