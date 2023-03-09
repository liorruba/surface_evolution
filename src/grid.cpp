// This class defines the surface properties:
#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING 1

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/utility.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/grid.hpp"
#include "../include/histogram.hpp"
#include "../include/log.hpp"

// Boost packages
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/indexable.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// A constructor for the grid class. Creates a grid with (size * size) elements.
Grid::Grid(std::vector< std::vector<double> > _initLayersList, std::vector<int8_t> _pxIdxMat){
	gridSize = regionWidth / resolution;
	x = linspace(-regionWidth/2, regionWidth/2, gridSize);
	y = linspace(-regionWidth/2, regionWidth/2, gridSize);
	area = pow(regionWidth, 2);
	pixelIndexMatrix = _pxIdxMat;
	initLayersList = _initLayersList;
	subsurfColumns = initializeSubsurface();
}

// Initialize a matrix of columns (2d vector):
std::vector< std::vector<SubsurfColumn> > Grid::initializeSubsurface(){
	addLogEntry("Generating grid (this may take some time, depending on the grid size)", true);
	SubsurfColumn buffCol = SubsurfColumn(); // for storing the subsurface grid
	std::vector<SubsurfColumn> buffColVecByIdx; // for storing the subsurface grid BY PIXEL INDEX (a vector of columns)
	std::vector< std::vector<SubsurfColumn> > buffMat;
	std::vector< std::vector<double> >::iterator row;

	// Create a buffer subsurface column, and populate it with the
	// requested layers:
	addLogEntry("Creating initial subsurface layering for a single column.", true);

	int prev = (*initLayersList.begin())[0];
	for (row = initLayersList.begin(); row != initLayersList.end(); ++row) {
		if ((*row)[0] != prev){
			prev = (*row)[0];
			buffColVecByIdx.push_back(buffCol);
			buffCol = SubsurfColumn();       
		}
		Layer buffLayer = Layer((*row)[1], (*row)[2], (*row)[3], (*row)[4]);
		buffCol.addLayer(buffLayer);
	}
	buffColVecByIdx.push_back(buffCol);

	addLogEntry("Finished creating column.", true);
	addLogEntry("Populating grid by dulicating columns...", true);

	// Populate the grid with the buffer subsurface column:
	int idx; // linear index for i,j
	for (int i = 0; i < gridSize; ++i){
		std::vector<SubsurfColumn> buffVec;
		for (int j = 0; j < gridSize; ++j){
			
			idx = getLinearIndex(i, j, gridSize);

			buffVec.push_back(buffColVecByIdx[pixelIndexMatrix[idx]]);
		}

		buffMat.push_back(buffVec);

		progressBar(i, gridSize);
	}

	addLogEntry("Finished creating grid.", true);

	return buffMat;
}

// Form a new crater:
void Grid::formCrater(Crater &crater){
	double distanceFromCraterCenter = 0; // For later use.
	double craterProfile = 0; // The maximum depth of the new cavity.
	double elevationAtCenter = getSurfaceElevationAtPoint(crater.xLocation, crater.yLocation);

	// Add crater to vector:
	cratersDict[cratersDict.size()] = &crater;

	// Add crater location and index to r-tree:
	craters_rtree.insert({cratersDict.size() - 1, {crater.xLocation, crater.yLocation}});

	// Find approximate inital i and j (up to 2R to include crater rim dropoff):
	int iInit = floor( ((crater.xLocation + regionWidth/2) - ejectaSpread/2 * crater.finalRadius) / resolution ); 
	int iFinal = ceil( ((crater.xLocation + regionWidth/2) + ejectaSpread/2 * crater.finalRadius) / resolution );
	int jInit = floor( ((crater.yLocation + regionWidth/2) - ejectaSpread/2 * crater.finalRadius) / resolution );
	int jFinal = ceil( ((crater.yLocation + regionWidth/2) + ejectaSpread/2 * crater.finalRadius) / resolution ); 

	//////////
	/// The next piece of code iterates over the subgrid radially outwards, to first compute
	/// the material inside the crater for it to be emplaced in the ejecta
	/// NOTE: Some of this code was written by OpenGPT! :-)
	//////////
	// Trim the grid edges
	if(iInit < 0) iInit = 0;
	if(iFinal > gridSize) iFinal = gridSize;
	if(jInit < 0) jInit = 0;
	if(jFinal > gridSize) jFinal = gridSize;

	// Create a vector of tuples containing the subgrid indexes (i, j) and their radial distances from the center
	std::vector<std::tuple<int, int, double>> indexes;
	
	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {

		  // Calculate the radial distance from the center of the subgrid
			double radialDistance = sqrt(pow(x[i] - crater.xLocation, 2) + pow(y[j] - crater.yLocation, 2));

		  // Add the tuple to the vector
			indexes.emplace_back(i, j, radialDistance);
		}
	}
	
	// Sort the vector of tuples by the radial distance
	std::sort(indexes.begin(), indexes.end(), [](const std::tuple<int, int, double>& a, 
		const std::tuple<int, int, double>& b) { return std::get<2>(a) < std::get<2>(b); });


	// Loop over the sorted vector of tuples
	for (const std::tuple<int, int, double>& index : indexes) {

		// Extract the index variables from the tuple                
		int i = std::get<0>(index);
		int j = std::get<1>(index);
		
		distanceFromCraterCenter = sqrt(pow(x[i] - crater.xLocation, 2) + pow(y[j] - crater.yLocation, 2));

		// if inside the crater:
		if (distanceFromCraterCenter <= crater.finalRadius) {
			if (craterProfileType == 1) {
				if (isEmplaceEjecta) {
					craterProfile = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter)
					- crater.rimHeight - linearInterp(crater.ejectaDistance, crater.ejectaThickness, crater.finalRadius);
				}
				else {
					craterProfile = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter)
					- crater.rimHeight;
				}

			}
			if (craterProfileType == 2) {
				if (isEmplaceEjecta) {
					craterProfile = craterSphericalDepthProfile(crater.finalRadius, distanceFromCraterCenter)
					- crater.rimHeight - linearInterp(crater.ejectaDistance, crater.ejectaThickness, crater.finalRadius);
				}
				else {
					craterProfile = craterSphericalDepthProfile(crater.finalRadius, distanceFromCraterCenter)
					- crater.rimHeight;
				}
			}

			craterProfile = elevationAtCenter - craterProfile;
			// craterProfile = getSurfaceElevationAtPoint(x[i], y[j]) - craterProfile;
			

			double leftToRemove = subsurfColumns.at(j).at(i).getSurfaceElevation() - craterProfile;

			// If the thickness of the material left to remove is positive, remove material
			if (leftToRemove > 0) {
				// std::cout << x[i] << " " << y[j] << " " << getSurfaceElevationAtPoint(x[i], y[j]) << " " << leftToRemove << std::endl;
				crater.ejectedMass.consolidate(subsurfColumns.at(j).at(i).integrateColumnComposition(leftToRemove));
				subsurfColumns.at(j).at(i).removeMaterial(leftToRemove);
				
				if ((crater.ejectedMass.iceFraction > 0.001) || (crater.ejectedMass.sootFraction > 0.001)) {
					// Adjust the composition of the ejected material, by removing ice lost through heating:
					crater.ejectedMass.changeComposition(crater.ejectedMass.regolithFraction, 
						crater.ejectedMass.iceFraction * ejectaVolatileRetention, 
						crater.ejectedMass.sootFraction * ejectaSootRetention);
				}
			}       
			else if (!crater.ejectedMass.isEmpty()) {
				// If the thickness of the left to remove is negative, add material				
				subsurfColumns.at(j).at(i).addLayer(Layer(fabs(leftToRemove),
					crater.ejectedMass.regolithFraction,
					crater.ejectedMass.iceFraction,
					crater.ejectedMass.sootFraction));
			}

			
		}
		
		// If outside the crater, add rim dropoff
		if ((distanceFromCraterCenter > crater.finalRadius) && (!crater.ejectedMass.isEmpty())) {
			double rimDropoff = crater.rimHeight * pow((distanceFromCraterCenter / crater.finalRadius), -rimDropoffExponent);
			subsurfColumns.at(j).at(i).addLayer(Layer(rimDropoff,
				crater.ejectedMass.regolithFraction,
				crater.ejectedMass.iceFraction,
				crater.ejectedMass.sootFraction));
		}


	}

	// Set crater formation elevation:
	crater.floorElevation = getSurfaceElevationAtPoint(crater.xLocation,
		crater.yLocation);
}


// Emplace ejecta:
void Grid::emplaceEjecta(Crater &crater){
	double distanceFromCraterCenter = 0; // For later use.
	// Find approximate inital i and j:
	int iInit = floor( ((crater.xLocation + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(iInit < 0) iInit = 0;
	int iFinal = ceil( ((crater.xLocation + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
	int jInit = floor( ((crater.yLocation + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(jInit < 0) jInit = 0;
	int jFinal = ceil( ((crater.yLocation + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;

	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {
			distanceFromCraterCenter = sqrt(pow(x.at(i) - crater.xLocation, 2) + pow(y[j] - crater.yLocation, 2));

			// if inside the ejecta blanket:
			if (distanceFromCraterCenter <= ejectaSpread * crater.finalRadius && distanceFromCraterCenter > crater.finalRadius) {
				// Changes subsurface composition based on ejected mass.
				subsurfColumns.at(j).at(i).addLayer(
					Layer(linearInterp(crater.ejectaDistance, crater.ejectaThickness, distanceFromCraterCenter),
						crater.ejectedMass.regolithFraction,
						crater.ejectedMass.iceFraction,
						crater.ejectedMass.sootFraction
						)
					);
			}
		}
	}
}


// Get the surface elevation at a specific point on the surface:
double Grid::getSurfaceElevationAtPoint(double pt_x, double pt_y){
	// Get indices of point (x,y):
	int i = std::lower_bound(x.begin(), x.end(), pt_x) - x.begin();
	int j = std::lower_bound(y.begin(), y.end(), pt_y) - y.begin();

	// Ghost craters are a special case, since floor(location) equals
	// gridSize, causing a segmentation fault:
	if (i == gridSize) i -= 1;
	if (j == gridSize) j -= 1;

	return subsurfColumns.at(j).at(i).getSurfaceElevation();
}

// Update the stored depth of an existing crater
void Grid::updateExistingCratersDepth(Crater &crater) {
	double dist;
	// query point
	point query_pt(crater.xLocation, crater.yLocation);
	// buffer for NN craters location
	point neigborCraterLocation;

	// iterate over nearest Values. Iterate over all nearest craters, but break
	// when the distance from the new impact point exceeds the crater radius.
	for ( rtree_t::const_query_iterator
		it = craters_rtree.qbegin(bgi::nearest(query_pt, cratersDict.size()));
		it != craters_rtree.qend();
		++it ) {

		// Calculate distance
		neigborCraterLocation = it->location;
	dist = bg::distance(query_pt, neigborCraterLocation);

		// Break if the distance from crater to new crater center exceeds
		// the crater radius.
		// Only update up to 2 times the ejecta spread distance, for efficiency.
		// This is ok due to the fast dropoff in the ejecta.
	if (dist > 2 * crater.finalRadius)
		break;

		// The first NN will always point at the newly formed crater, so
		// ignore it:
	if (dist > 0) {
		double prevFloorElevation = cratersDict[it->index]
		->floorElevation;

		double newFloorElevation = getSurfaceElevationAtPoint(
			bg::get<0>(neigborCraterLocation),
			bg::get<1>(neigborCraterLocation)
			);

			// Update crater depth
		cratersDict[it->index]->finalDepth -=
		(newFloorElevation - prevFloorElevation);

			// Update crater center elevation
		cratersDict[it->index]->floorElevation = newFloorElevation;

			// If the depth of the crater changed by more than 50% (shallower or deeper)
			// remove the crater from the tree
		if ((fabs(cratersDict[it->index]->finalDepth - cratersDict[it->index]->finalDepth_init)
			/ (cratersDict[it->index]->finalDepth_init)) > 0.5) {

				// Remove crater from the map
			cratersDict.erase(it->index);

				// Remove crater from the tree
		craters_rtree.remove(*it);

	}
}
}
}

///////////////////
// Crater profiles:
///////////////////
// Parabolic:
double Grid::craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter){
	return depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2));

}

// Bowl shaped (hemispherical):
double Grid::craterSphericalDepthProfile(double craterRadius, double distanceFromCraterCenter){
	double craterDepth = 2 * craterRadius * depthToDiameter;
	double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;

	return -sphereRadius + craterDepth + sqrt(pow(sphereRadius,2) - pow(distanceFromCraterCenter,2));
}

///////////////////
// Simple sublimation/accumulation:
///////////////////
void Grid::sublimateIce() {
	return;
	// TODO: ADD NORBERT'S MODEL
}

void Grid::depositLayer(Layer layer) {
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j) {
			// If the top layer has some ice, sublimate:
			if (subsurfColumns[i][j].subsurfLayers.back().iceFraction > 0){
				subsurfColumns[i][j].addLayer(layer);
			}
		}
	}
}


//
// Compute the gradient of the surface topopgraphy
//
std::vector < std::vector<double> > Grid::computeGradient(std::vector < std::vector<double> > Z) {
	std::vector< std::vector<double> > slopeMap;
	double dx, dy;

	size_t subgridSz = Z.size();

	for (size_t i = 0; i < subgridSz; ++i) {
		std::vector<double> buff;

		for (size_t j = 0; j < subgridSz; ++j) {
			
			if (i == subgridSz - 1) {
				dx = (Z[i][j] - Z[i - 1][j]) / resolution;
			} else {
				dx = (Z[i + 1][j] - Z[i][j]) / resolution;
			}

			if (j == subgridSz - 1) {
				dy = (Z[i][j] - Z[i][j - 1]) / resolution;
			} else {
				dy = (Z[i][j + 1] - Z[i][j]) / resolution;
			}

			buff.push_back(sqrt(dx * dx + dy * dy));
		}
		slopeMap.push_back(buff);
	}

	return slopeMap;
}

///////////////////
// Linear slope diffusion:
///////////////////
// Helper function for the linear slope diffusion model
bool Grid::compareMaxSlope(std::vector< std::vector<double> > slopeMap, double maxSlope) {
	size_t sz = slopeMap.size();

	for (size_t i = 0; i < sz; ++i) {
		for (size_t j = 0; j < sz; ++j) {
			if (slopeMap[i][j] > maxSlope) {
				return true;
			}
		}
	}

	return false;
}

std::vector< std::vector<double> > Grid::linearSlopeDiffusion(std::vector< std::vector<double> > Z, double maxSlope, double K) {
	// If K == -1, choose the diffusivity based on the courant condition:
	if (K == -1) {
		K = 0.5 * resolution;
	}

	double Kbuff = K / (resolution * resolution);
	size_t sz_Z = Z.size();

	while (true) {
		// Compute the gradient magnitude
		std::vector< std::vector<double> > slopeMap = computeGradient(Z);
		bool converged = true;

		// Apply diffusion with constant boundary conditions
		if (compareMaxSlope(slopeMap, maxSlope)) {
			for (size_t i = 1; i < sz_Z - 1; ++i) {
				for (size_t j = 1; j < sz_Z - 1; ++j) {
					double buffii = (Z[i+1][j] - 2 * Z[i][j] + Z[i-1][j]);
					double buffjj = (Z[i][j+1] - 2 * Z[i][j] + Z[i][j-1]);
					Z[i][j] = Z[i][j] + Kbuff * (buffii + buffjj);
					converged = false;
				}
			}
		}
		
		if (converged) {
			break;
		}
	}

	return Z;
}

///////////////////
// Surface modification below the angle of repose (in deg):
///////////////////
void Grid::thresholdSlopes(double angleOfRepose) {
	// Compute the slope of repose
	double slopeOfRepose = tan(M_PI * angleOfRepose / 180.0);

	// Create a grid subset based on the initial and final indexes:
	std::vector< std::vector<double> > subGrid;

	for (size_t i = 0; i < subGrid.size(); ++i) {
		std::vector<double> buff;

		for (size_t j = 0; j < subGrid.size(); ++j) {
			buff.push_back(subsurfColumns[i][j].getSurfaceElevation());
		}

		subGrid.push_back(buff);
	}

	// Diffuse the slopes until they reach the angle of repose:
	subGrid = linearSlopeDiffusion(subGrid, slopeOfRepose, -1);

	// Reproduce the elevation
	for (size_t i = 0; i < subGrid.size(); i++) {
		for (size_t j = 0; j < subGrid.size(); j++) {
			Layer buffLayer = subsurfColumns[i][j].subsurfLayers.back();
			double topoDiff = subsurfColumns[i][j].getSurfaceElevation() - subGrid[i][j];

			if (topoDiff >= 0) {
				subsurfColumns[i][j].removeMaterial(topoDiff);
			}

			else {
				subsurfColumns[i][j].addLayer(Layer(fabs(topoDiff), buffLayer.regolithFraction, buffLayer.iceFraction, buffLayer.sootFraction));
			}
		}
	}
}


/////////////////
// Ice Stability:
/////////////////
// The temperature depndent surface loss rate based on WMB 1962
// Returns loss rate in kg Ma^-1
// ********************
// ***IN DEVELOPMENT***
// ********************
// }

/////////
// Print:
/////////
// Print surface to file:

void Grid::printSurface(int index, bool isfinal) {
	char logEntry[100];
	sprintf(logEntry, "Writing surface elevation and composition to file for time step: %d.", index);
	addLogEntry(logEntry, false);
	std::ofstream elevationFile;
	std::ofstream regolithFractionFile;
	std::ofstream iceFractionFile;
	std::ofstream sootFractionFile;

	std::string index_str = std::string(
		std::to_string(int(endTime / printTimeStep)).length() -
		std::to_string(index).length(), '0'
		) + std::to_string(index);

	std::string elevationFileName = "./output/elevation_" + index_str + ".out";
	std::string regolithFractionFileName = "./output/regolithFraction_" + index_str + ".out";
	std::string iceFractionFileName = "./output/iceFraction_" + index_str + ".out";
	std::string sootFractionFileName = "./output/sootFraction_" + index_str + ".out";
	
	elevationFile.open(elevationFileName, std::ios_base::binary|std::ios_base::app);
	regolithFractionFile.open(regolithFractionFileName, std::ios_base::binary|std::ios_base::app);
	iceFractionFile.open(iceFractionFileName, std::ios_base::binary|std::ios_base::app);
	sootFractionFile.open(sootFractionFileName, std::ios_base::binary|std::ios_base::app);

	// Put data in matrices:
	std::vector< std::vector<double> > elevationMatrix;
	std::vector< std::vector<double> > regFractionMatrix;
	std::vector< std::vector<double> > iceFractionMatrix;
	std::vector< std::vector<double> > sootFractionMatrix;

	for (int i = 0; i < gridSize; ++i) {
		std::vector<double> elevationBuffVector;
		std::vector<double> regFractionBuffVector;
		std::vector<double> iceFractionBuffVector;
		std::vector<double> sootFractionBuffVector;

		for (int j = 0; j < gridSize; ++j) {
			Layer buffLayer = subsurfColumns.at(j).at(i).subsurfLayers.back();

			elevationBuffVector.push_back(subsurfColumns.at(j).at(i).getSurfaceElevation());
			regFractionBuffVector.push_back(buffLayer.regolithFraction);
			iceFractionBuffVector.push_back(buffLayer.iceFraction);
			sootFractionBuffVector.push_back(buffLayer.sootFraction);
		}

		elevationMatrix.push_back(elevationBuffVector);
		regFractionMatrix.push_back(regFractionBuffVector);
		iceFractionMatrix.push_back(iceFractionBuffVector);
		sootFractionMatrix.push_back(sootFractionBuffVector);
	}

	// If the input bin size is > resolution, re-bin the data to save a low-resolution version of the file:
	int griddedBinSize = gridSize;

	if (downsamplingResolution > resolution) {
		addLogEntry("The resampling resolution is greater than the resolution. Downsampling surface...", false);
		elevationMatrix = bin_2d_vector(elevationMatrix, downsamplingResolution);
		regFractionMatrix = bin_2d_vector(regFractionMatrix, downsamplingResolution);
		iceFractionMatrix = bin_2d_vector(iceFractionMatrix, downsamplingResolution);
		sootFractionMatrix = bin_2d_vector(sootFractionMatrix, downsamplingResolution);

		griddedBinSize = (int) gridSize / (downsamplingResolution / resolution);
	}

	for (int i = 0; i < griddedBinSize; ++i) {
		for (int j = 0; j < griddedBinSize; ++j) {
			elevationFile.write(reinterpret_cast<const char*>(&elevationMatrix[i][j]), sizeof(double));
			regolithFractionFile.write(reinterpret_cast<const char*>(&regFractionMatrix[i][j]), sizeof(double));
			iceFractionFile.write(reinterpret_cast<const char*>(&iceFractionMatrix[i][j]), sizeof(double));
			sootFractionFile.write(reinterpret_cast<const char*>(&sootFractionMatrix[i][j]), sizeof(double));
		}
	}

	// If it is the final time step, also write the horizontal coordinates:
	if (isfinal) {
		addLogEntry("Writing x,y coordinates to file...", false);
		std::ofstream xFile;
		std::ofstream yFile;

		std::string xFileName = "./output/x.out";
		std::string yFileName = "./output/y.out";

		xFile.open(xFileName, std::ios_base::binary|std::ios_base::app);
		yFile.open(yFileName, std::ios_base::binary|std::ios_base::app);

		std::vector<double> x_binned = bin_1d_vector(x, downsamplingResolution);
		std::vector<double> y_binned = bin_1d_vector(y, downsamplingResolution);

		xFile.write(reinterpret_cast<const char*>(x_binned.data()), x_binned.size() * sizeof(double));
		yFile.write(reinterpret_cast<const char*>(y_binned.data()), y_binned.size() * sizeof(double));

		xFile.close();
		yFile.close();
	}

	elevationFile.close();
	regolithFractionFile.close();
	iceFractionFile.close();
	sootFractionFile.close();
}

// Integrate the subsurface composition down to some depth and print it
void Grid::printIntegratedSubsurface(double depth, int index){
	std::ofstream regolithFractionFile;
	std::ofstream iceFractionFile;
	std::ofstream sootFractionFile;
	
	std::string index_str = std::string(
		std::to_string(int(endTime / printTimeStep)).length() -
		std::to_string(index).length(), '0'
		) + std::to_string(index);

	std::string regolithFractionFileName = "./output/depthRegolithFraction_" + index_str + ".out";
	std::string iceFractionFileName = "./output/depthIceFraction_" + index_str + ".out";
	std::string sootFractionFileName = "./output/depthSootFraction_" + index_str + ".out";
	
	regolithFractionFile.open(regolithFractionFileName, std::ios_base::binary|std::ios_base::app);
	iceFractionFile.open(iceFractionFileName, std::ios_base::binary|std::ios_base::app);
	sootFractionFile.open(sootFractionFileName, std::ios_base::binary|std::ios_base::app);

	// Put data in matrices:
	std::vector< std::vector<double> > regFractionMatrix;
	std::vector< std::vector<double> > iceFractionMatrix;
	std::vector< std::vector<double> > sootFractionMatrix;

	for (int i = 0; i < gridSize; ++i) {
		std::vector<double> regFractionBuffVector;
		std::vector<double> iceFractionBuffVector;
		std::vector<double> sootFractionBuffVector;

		for (int j = 0; j < gridSize; ++j) {
			Layer buffLayer = subsurfColumns.at(j).at(i).integrateColumnComposition(depth);

			regFractionBuffVector.push_back(buffLayer.regolithFraction);
			iceFractionBuffVector.push_back(buffLayer.iceFraction);
			sootFractionBuffVector.push_back(buffLayer.sootFraction);
		}

		regFractionMatrix.push_back(regFractionBuffVector);
		iceFractionMatrix.push_back(iceFractionBuffVector);
		sootFractionMatrix.push_back(sootFractionBuffVector);
	}

	// If the input bin size is > resolution, re-bin the data to save a low-resolution version of the file:
	int griddedBinSize = gridSize;

	if (downsamplingResolution > resolution) {
		addLogEntry("The resampling resolution is greater than the resolution. Downsampling integrated subsurface....", false);
		regFractionMatrix = bin_2d_vector(regFractionMatrix, downsamplingResolution);
		iceFractionMatrix = bin_2d_vector(iceFractionMatrix, downsamplingResolution);
		sootFractionMatrix = bin_2d_vector(sootFractionMatrix, downsamplingResolution);

		griddedBinSize = (int) gridSize / (downsamplingResolution / resolution);
	}

	for (int i = 0; i < griddedBinSize; ++i) {
		for (int j = 0; j < griddedBinSize; ++j) {
			regolithFractionFile.write(reinterpret_cast<const char*>(&regFractionMatrix[i][j]), sizeof(double));
			iceFractionFile.write(reinterpret_cast<const char*>(&iceFractionMatrix[i][j]), sizeof(double));
			sootFractionFile.write(reinterpret_cast<const char*>(&sootFractionMatrix[i][j]), sizeof(double));
		}
	}

	regolithFractionFile.close();
	iceFractionFile.close();
	sootFractionFile.close();
}

// Print subsurface to file:
void Grid::printSubsurface(int index){
	SubsurfColumn col =  SubsurfColumn();
	std::fstream outputFile;

	std::string index_str = std::string(
		std::to_string(int(endTime / printTimeStep)).length() -
		std::to_string(index).length(), '0'
		) + std::to_string(index);

	std::string output_filename = "./output/subsurface_" + index_str + ".out";

	outputFile.open(output_filename, std::ios_base::binary|std::ios_base::app);

	for (size_t i = 0; i < subsurfColumns.size(); ++i) {
		for (size_t j = 0; j < subsurfColumns.size(); ++j) {
			SubsurfColumn col = subsurfColumns.at(j).at(i);
			// Prepare dummy layer whose first element is the number of layers in
			// column and second element is the surface elevation.
			Layer dummyLayer = Layer(col.subsurfLayers.size(), col.getSurfaceElevation());

			outputFile.write((char*)&dummyLayer, sizeof(Layer));
			outputFile.write((char*)&col.subsurfLayers[0],col.subsurfLayers.size() * sizeof(Layer));
		}
	}

	outputFile.close();
}

// Print existing craters to histogram:
void Grid::printExistingCratersToHistogram(double bins){
	addLogEntry("Printing craters histrogram.", false);
	Histogram hist(minimumImpactorDiameter * 10, regionWidth, bins);

	for (auto cm : cratersDict) {
		// Add the crater diameter to the histogram
		hist.add(2 * cm.second->finalRadius);
	}

	hist.print("./output/existing_craters_histogram.txt");
}

// Print existing craters to histogram:
void Grid::printExistingCraters(){
	addLogEntry("Printing craters to file.", false);

	std::ofstream craterFile;

	std::string elevationFileName = "./output/existing_craters.txt";
	craterFile.open(elevationFileName, std::ios_base::out);

	// Print titles:
	craterFile << "diameter, depth, initial_depth\n";

	for (auto cm : cratersDict) {
		// Add the crater diameter to the histogram
		craterFile << 2 * cm.second->finalRadius << ","
		<< cm.second->finalDepth << ","
		<< cm.second->finalDepth_init << "\n";
	}
	craterFile.close();
}
