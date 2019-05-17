// Get current time
uint64_t get_time();

// Function to check if string is empty:
int isEmpty(const char *s);

// Add a new log entry
void addLogEntry(char *str);

// Read config file:
varlist readConfig();

// Get variable from list:
double setVariable(varlist varList, char * varName);

// Generate a uniformly distributed random number:
double randU(double low, double high);


// This function creates a matrix populated by zeros, whose size is nxm.
double ** zeros(int n, int m);

histogram createHistogram(double minBin, double maxBin, int numOfBins);

// This function recieves a value and an existing histogram and returns a histogram.
histogram addToHistogram(double *val, histogram hist);

// Generate list of crater diameters, write to file:
double * getDiameterList(int totalNumberOfCraters, double minimumDiameter);

void createCraterInZmat(int gridSize, double *craterDiameter, double xPosition, double yPosition, double **xmat, double **ymat, double **zmat);

void printMatrixToFile(double ** matrix, int gridSize, char * fileName);

////////////////
// START MAIN //
////////////////

int main();
