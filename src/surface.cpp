// This class defines the surface properties:
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include <fstream>
#include<string>
#include "../include/impactor.hpp"
#include "../include/crater.hpp"
#include "../include/utility.hpp"
#include "../include/surface.hpp"
#include "../include/regolit_main.hpp"

// A constructor for the surface class. Creates a surface with (size * size) elements.
Surface::Surface(double regionWidth, double resolution){
  gridSize = regionWidth / resolution;
  x = linspace(-regionWidth/2, regionWidth/2, gridSize);
  y = linspace(-regionWidth/2, regionWidth/2, gridSize);
  elevation = initMatrix(); // Initialize the elevation matrix.
  // std::vector< std::vector<double> > ice(gridSize,std::vector<double>(gridSize)); // Initialize the ice matrix.
}

// Initialize a matrix (2d vector):
std::vector< std::vector<double> > Surface::initMatrix(){
  std::vector< std::vector<double> > matrix(gridSize,std::vector<double>(gridSize));
  return matrix;
}

// Form a new crater:
void Surface::formCrater(Crater crater){
	double distanceFromCraterCenter = 0; // For later use.
	double newDepth = 0; // The new crater depth.

	// Find approximate inital i and j:
	int iInit = floor( ((crater.xPosition + regionWidth/2) - 2 * crater.finalRadius) / resolution ); if(iInit < 0) iInit = 0;
	int iFinal = ceil( ((crater.xPosition + regionWidth/2) + 2 * crater.finalRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
	int jInit = floor( ((crater.yPosition + regionWidth/2) - 2 * crater.finalRadius) / resolution ); if(jInit < 0) jInit = 0;
	int jFinal = ceil( ((crater.yPosition + regionWidth/2) + 2 * crater.finalRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;
  int iCenter = crater.xPosition + (gridSize)/2;
  int jCenter = crater.yPosition + (gridSize)/2;
  if (iCenter <= 0) iCenter = 0;
  if (iCenter >= gridSize) iCenter = gridSize - 1;
  if (jCenter <= 0) jCenter = 0;
  if (jCenter >= gridSize) jCenter = gridSize - 1;
  double refElevation = elevation[jCenter][iCenter];

	for (int i = iInit; i < iFinal; i++){
		for (int j = jInit; j < jFinal; j++){
			distanceFromCraterCenter = sqrt(pow(x[i] - crater.xPosition, 2) + pow(y[j] - crater.yPosition, 2));
			// if inside the crater:
			if (distanceFromCraterCenter <= crater.finalRadius){
        if (craterProfile == 1){
          newDepth = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter, refElevation, 0);
        }
        if (craterProfile == 2){
          newDepth = craterBowlShapedDepthProfile(crater.finalRadius, distanceFromCraterCenter, refElevation, 0);
        }
        elevation[j][i] = newDepth;
			}
		}
	}
}

// Emplace crater rim and ejecta:
void Surface::emplaceEjecta(Crater crater){
  double distanceFromCraterCenter = 0; // For later use.

  // Find approximate inital i and j:
  int iInit = floor( ((crater.xPosition + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(iInit < 0) iInit = 0;
  int iFinal = ceil( ((crater.xPosition + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
  int jInit = floor( ((crater.yPosition + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(jInit < 0) jInit = 0;
  int jFinal = ceil( ((crater.yPosition + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;

  for (int i = iInit; i < iFinal; i++){
    for (int j = jInit; j < jFinal; j++){
      distanceFromCraterCenter = sqrt(pow(x[i] - crater.xPosition, 2) + pow(y[j] - crater.yPosition, 2));
      // if inside the ejecta blanket:
      if (distanceFromCraterCenter <= ejectaSpread * crater.finalRadius && distanceFromCraterCenter > crater.finalRadius ){
        elevation[j][i] += rimHeight(crater.finalRadius, distanceFromCraterCenter);
        elevation[j][i] += linearInterp(crater.ejectaDistance, crater.ejectaThickness, distanceFromCraterCenter);
      }
    }
  }
}

// Print surface to file:
void Surface::print(int index){
  std::ofstream xFile;
  std::ofstream yFile;
  std::ofstream elevationFile;

  std::string elevationFileName = "./output/elevation_" + std::to_string(index) + ".txt";
  elevationFile.open(elevationFileName, std::ios_base::out);

  if (index != 0){
    // For the last print, print x,y vectors:
    std::string xFileName = "./output/x.txt";
    std::string yFileName = "./output/y.txt";
    xFile.open(xFileName, std::ios_base::out);
    yFile.open(yFileName, std::ios_base::out);
    for (int i = 0; i < gridSize; i++){
      xFile << x[i] << ",";
      yFile << y[i] << ",";
    }
  }

	for (int i = 0; i < gridSize; i++){
		for (int j = 0; j < gridSize; j++){
      elevationFile << elevation[i][j] << ",";
		}
    elevationFile << "\n";
	}
}

///////////////////
// Crater profiles:
///////////////////
// Parabolic:
double Surface::craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter, double refElevation, double normAngle){
  normAngle = 0;
  return depthToDiameter * 2 * craterRadius * (pow(distanceFromCraterCenter/craterRadius,2) - 1) + rimToDiameter * 2 * craterRadius + refElevation;
}

// Bowl shaped:
double Surface::craterBowlShapedDepthProfile(double craterRadius, double distanceFromCraterCenter, double refElevation, double normAngle){
  double craterDepth = 2 * craterRadius * depthToDiameter;
  double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;
  normAngle = 0;
  return sphereRadius - craterDepth - sqrt(pow(sphereRadius,2) - pow(distanceFromCraterCenter,2)) + 2 * craterRadius * rimToDiameter + refElevation;
}

// Rim height dropoff above the surface:
double Surface::rimHeight(double craterRadius, double distanceFromCraterCenter){
  return (rimToDiameter * craterRadius) * pow(distanceFromCraterCenter, -outerRimExponent);
}
