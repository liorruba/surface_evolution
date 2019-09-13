// This class defines the surface properties:
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<fstream>
#include<string>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/utility.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/grid.hpp"

// A constructor for the surface class. Creates a surface with (size * size) elements.
Grid::Grid(double regionWidth, double resolution){
  gridSize = regionWidth / resolution;
  x = linspace(-regionWidth/2, regionWidth/2, gridSize);
  y = linspace(-regionWidth/2, regionWidth/2, gridSize);
  area = pow(gridSize, 2);
  subsurfColumns = initializeSubsurface();
}

// Initialize a matrix (2d vector):
std::vector< std::vector<SubsurfColumn> > Grid::initializeSubsurface(){
  std::vector< std::vector<SubsurfColumn> > buffMat(gridSize,std::vector<SubsurfColumn>(gridSize));
  return buffMat;
}

// Form a new crater:
void Grid::formCrater(Crater &crater){
	double distanceFromCraterCenter = 0; // For later use.
	double materialToRemove = 0; // The new crater depth.

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

	for (int i = iInit; i < iFinal; i++){
		for (int j = jInit; j < jFinal; j++){
			distanceFromCraterCenter = sqrt(pow(x[i] - crater.xPosition, 2) + pow(y[j] - crater.yPosition, 2));
			// if inside the crater:
			if ((distanceFromCraterCenter <= crater.finalRadius)){
        if (craterProfile == 1){
          materialToRemove = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0);
        }
        if (craterProfile == 2){
          materialToRemove = craterBowlShapedDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0);
        }
        if (fabs(materialToRemove) > fabs(initialThickness - subsurfColumns[j][i].get_surfaceElevation())) {
          double removalDepthRelativeToSurface = materialToRemove - fabs(initialThickness - subsurfColumns[j][i].get_surfaceElevation());
          crater.ejectedMass.consolidate(subsurfColumns[j][i].integrateColumnComposition(removalDepthRelativeToSurface));
          subsurfColumns[j][i].removeMaterial(removalDepthRelativeToSurface);
        }
			}
		}
	}
}

// Emplace crater rim and ejecta:
void Grid::emplaceEjecta(Crater &crater){
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
      if (distanceFromCraterCenter <= ejectaSpread * crater.finalRadius && distanceFromCraterCenter > crater.finalRadius) {
        // Changes subsurface composition based on ejected mass.
        subsurfColumns[j][i].addLayer(
          Layer(linearInterp(crater.ejectaDistance, crater.ejectaThickness, distanceFromCraterCenter),
          crater.ejectedMass.regolithFraction,
          crater.ejectedMass.iceFraction,
          crater.ejectedMass.sootFraction)
        );
      }
    }
  }
}

///////////////////
// Crater profiles:
///////////////////
// Parabolic:
double Grid::craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle){
  normAngle = 0;
  if (isEmplaceEjecta)
    return depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2)) - rimToDiameter * 2 * craterRadius;
  else
    return depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2));
}

// Bowl shaped:
double Grid::craterBowlShapedDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle){
  double craterDepth = 2 * craterRadius * depthToDiameter;
  double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;
  normAngle = 0;

  if (isEmplaceEjecta)
    return -sphereRadius + craterDepth + sqrt(pow(sphereRadius,2) - pow(distanceFromCraterCenter,2)) - rimToDiameter * 2 * craterRadius;
  else
    return -sphereRadius + craterDepth + sqrt(pow(sphereRadius,2) - pow(distanceFromCraterCenter,2));

}

// Rim height dropoff above the surface:
double Grid::rimHeight(double craterRadius, double distanceFromCraterCenter){
  return (rimToDiameter * 2 * craterRadius) * pow(distanceFromCraterCenter, -outerRimExponent);
}

/////////
// Print:
/////////
// Print surface to file:
void Grid::printSurface(int index){
  std::ofstream xFile;
  std::ofstream yFile;
  std::ofstream elevationFile;

  std::string elevationFileName = "./output/elevation_" + std::to_string(index) + ".txt";
  elevationFile.open(elevationFileName, std::ios_base::out);

  if (index == 0){
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
      elevationFile << subsurfColumns[i][j].get_surfaceElevation() << ",";
		}
    elevationFile << "\n";
	}
}

// Print subsurface to file:
void Grid::printGrid(int index){

  SubsurfColumn col =  SubsurfColumn();
  std::fstream outputFile;
  std::string output_filename = "./output/output_" + std::to_string(index) + ".out";

  outputFile.open(output_filename, std::ios_base::binary|std::ios_base::app);

  for (size_t i = 0; i < subsurfColumns.size(); i++) {
    for (size_t j = 0; j < subsurfColumns.size(); j++) {
      SubsurfColumn col = subsurfColumns[i][j];
      // Prepare dummy layer whose first element is the number of layers in
      // column and second element is the surface elevation.
      Layer dummyLayer = Layer(col.subsurfLayers.size(), col.get_surfaceElevation());

      outputFile.write((char*)&dummyLayer, sizeof(Layer));
      outputFile.write((char*)&col.subsurfLayers[0],col.subsurfLayers.size() * sizeof(Layer));
    }
  }

  outputFile.close();

  // col_file.write((char*)&l, sizeof(Layer));
  // col_file.write((char*)&layersToPrint[0],layersToPrint.size() * sizeof(Layer));
  // col_file.write((char*)&l, sizeof(Layer));
  // col_file.write((char*)&layersToPrint[0],layersToPrint.size() * sizeof(Layer));
  //
  // // col_file.write((char*)col.subsurfLayers.size(),sizeof(unsigned long));
  // // col_file.write((char*)&col2.subsurfLayers[0],col2.subsurfLayers.size() * sizeof(Layer));
  // col_file.close();
  // Set data cubes:
  // std::vector<std::vector<std::vector<double> > > thickness_cube;
  // std::vector<std::vector<std::vector<double> > > regolith_cube;
  // std::vector<std::vector<std::vector<double> > > ice_cube;
  // std::vector<std::vector<std::vector<double> > > soot_cube;
  //
  // std::fstream thickness_file;
  // std::fstream regolith_file;
  // std::fstream ice_file;
  // std::fstream soot_file;
  //
  // std::string thickness_filename = "./output/thickness_" + std::to_string(index) + ".cube";
  // std::string regolith_filename = "./output/regolith_frac_" + std::to_string(index) + ".cube";
  // std::string ice_filename = "./output/ice_frac_" + std::to_string(index) + ".cube";
  // std::string soot_filename = "./output/soot_frac_" + std::to_string(index) + ".cube";
  //
  // for (int i = 0; i < )
  //
  // std::cout << sizeof(subsurfColumns) << std::endl;
  // outputFile.open(filename, std::ios_base::out|std::ios_base::binary);
  // outputFile.write((char*)&subsurfColumns,sizeof(subsurfColumns));
  // outputFile.close();
}
