// Function to check if string is empty:
int isEmpty(const char *s);

// Read config file:
varlist readConfig();

// Get variable from list:
double setVariable(varlist varList, char * varName);

// Generate a uniformly distributed random number:
double randU(double low, double high);

void createCraterInZmat(int gridSize, double *craterDiameter, double xPosition, double yPosition, double **xmat, double **ymat, double **zmat);

////////////////
// START MAIN //
////////////////

int main();
