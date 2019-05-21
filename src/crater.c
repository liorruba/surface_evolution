#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../include/regolit_main.h"

// Crater profiles:
double craterParabolicProfile(double craterRadius, double depthToDiameter, double distanceFromCraterCenter){
  return 2 * depthToDiameter / craterRadius * pow(distanceFromCraterCenter,2) - 2 * depthToDiameter * craterRadius;
}

double craterBowlShapedProfile(double craterRadius, double depthToDiameter, double distanceFromCraterCenter){
  double craterDepth = 2 * craterRadius * depthToDiameter;
  double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;

  return sphereRadius - craterDepth - sqrt(pow(sphereRadius,2) - pow(distanceFromCraterCenter,2));
}

// Form a new crater:
void createCrater(int gridSize, double *craterDiameter, double xPosition, double yPosition, double **xmat, double **ymat, double **zmat, int craterProfile){
	double distanceFromCraterCenter = 0; // For later use.
	double newDepth = 0; // The new crater depth.

	double craterRadius = *craterDiameter/2; // Calculate the crater radius

	// Find approximate inital i and j:
	int iInit = floor( ((xPosition + regionWidth/2) - 2 * craterRadius) / resolution ); if(iInit < 0) iInit = 0;
	int iFinal = ceil( ((xPosition + regionWidth/2) + 2 * craterRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
	int jInit = floor( ((yPosition + regionWidth/2) - 2 * craterRadius) / resolution ); if(jInit < 0) jInit = 0;
	int jFinal = ceil( ((yPosition + regionWidth/2) + 2 * craterRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;

	for (int i = iInit; i < iFinal; i++){
		for (int j = jInit; j < jFinal; j++){
			distanceFromCraterCenter = sqrt(pow(xmat[j][i] - xPosition, 2) + pow(ymat[j][i] - yPosition, 2));
			// if inside the crater:
			if (distanceFromCraterCenter <= craterRadius){
        if (craterProfile == 1){
          newDepth = craterParabolicProfile(craterRadius, depthToDiameter, distanceFromCraterCenter);
        }
        if (craterProfile == 2){
          newDepth = craterBowlShapedProfile(craterRadius, depthToDiameter, distanceFromCraterCenter);
        }
        if (craterProfile == 3){
          // TODO: Create Fasset crater
          newDepth = craterBowlShapedProfile(craterRadius, depthToDiameter, distanceFromCraterCenter);
        }

				if (newDepth < zmat[j][i]){
          zmat[j][i] = newDepth;
				}
			}
		}
	}
}
