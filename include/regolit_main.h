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

// Globals
// Simulation parameters:
extern double regionWidth;  // km
extern double resolution; // km
extern double endTime; // Ma
extern double depthToDiameter;
extern double minimumDiameter; // Crater minimum diameter, km.
extern double printTimeStep; // Time step for printing data in Ma.
extern int craterProfile;

// Crater production constants:
extern double b;
extern double c;

// Read config file:
varlist readConfig();

// Get variable from list:
double setVariable(varlist varList, char * varName);

////////////////
// START MAIN //
////////////////
int main();
