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
Grid::Grid(std::vector< std::vector<double> > _initLayersList, std::vector< std::vector<int> > _pxIdxMat){
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
        for (int i = 0; i < gridSize; ++i){
                std::vector<SubsurfColumn> buffVec;
                for (int j = 0; j < gridSize; ++j){
                        
                        buffVec.push_back(buffColVecByIdx[pixelIndexMatrix[i][j] - 1]);
                }

                buffMat.push_back(buffVec);

                progressBar(i, gridSize);
        }

        addLogEntry("Finished creating grid.", true);
        // for (int i = 0; i < gridSize; ++i){
        //         std::cout << "hi" << std::endl;
        //         for (int j = 0; j < gridSize; ++j){
                        
        //                 buffMat[i][j].subsurfLayers.back().print();
        //         }
        // }

        return buffMat;
}

// Form a new crater:
void Grid::formCrater(Crater &crater){
        double distanceFromCraterCenter = 0; // For later use.
        double craterProfile = 0; // The maximum depth of the new cavity.
        // double elevationAtCenter = getSurfaceElevationAtPoint(crater.xLocation, crater.yLocation);

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

        // Compute the slope the craters is formed on
        // auto [craterSlope, craterAspect] = calculateSlope(crater);

        // Create a vector of tuples containing the subgrid indexes (i, j) and their radial distances from the center
        std::vector<std::tuple<int, int, double>> indexes;
        
        for (int i = iInit; i < iFinal; ++i) {
            for (int j = jInit; j < jFinal; ++j) {
                  
                  // Calculate the radial distance from the center of the subgrid
                  double radialDistance = sqrt(pow(x.at(i) - crater.xLocation, 2) + pow(y.at(j) - crater.yLocation, 2));

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
                
                distanceFromCraterCenter = sqrt(pow(x.at(i) - crater.xLocation, 2) + pow(y.at(j) - crater.yLocation, 2));

                // if inside the crater:
                if (distanceFromCraterCenter <= crater.finalRadius) {
                        if (craterProfileType == 1) {
                                if (isEmplaceEjecta) {
                                        craterProfile = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0)
                                                        - crater.rimHeight - linearInterp(crater.ejectaDistance, crater.ejectaThickness, crater.finalRadius);
                                }
                                else {
                                        craterProfile = craterParabolicDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0)
                                                        - crater.rimHeight;
                                }

                        }
                        if (craterProfileType == 2) {
                                if (isEmplaceEjecta) {
                                        craterProfile = craterSphericalDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0)
                                                        - crater.rimHeight - linearInterp(crater.ejectaDistance, crater.ejectaThickness, crater.finalRadius);
                                }
                                else {
                                        craterProfile = craterSphericalDepthProfile(crater.finalRadius, distanceFromCraterCenter, 0)
                                                        - crater.rimHeight;
                                }
                        }

                        craterProfile = getSurfaceElevationAtPoint(x.at(i), y.at(j)) - craterProfile;

                        double leftToRemove = subsurfColumns.at(j).at(i).getSurfaceElevation() - craterProfile;
                        
                        // If the thickness of the material left to remove is positive, remove material
                        if (leftToRemove > 0) {
                                crater.ejectedMass.consolidate(subsurfColumns.at(j).at(i).integrateColumnComposition(leftToRemove));
                                subsurfColumns.at(j).at(i).removeMaterial(leftToRemove);
                                
                                if (crater.ejectedMass.iceFraction > 0.01) {
                                        // std::cout << "before " << crater.ejectedMass.iceFraction << std::endl;
                                        // Adjust the composition of the ejected material, by removing ice lost through heating:
                                        crater.ejectedMass.changeComposition(crater.ejectedMass.regolithFraction, 
                                                crater.ejectedMass.iceFraction * ejectaVolatileRetention, 
                                                crater.ejectedMass.sootFraction);
                                        // std::cout << "after "<< crater.ejectedMass.iceFraction << std::endl;
                                }
                        }       
                        // If the thickness of the left to remove is negative, add material
                        else if (!crater.ejectedMass.isEmpty()) {
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

// Compute the slope and aspect the crater will rest on:
std::tuple<double, double> Grid::calculateSlope(const Crater crater){
        // Three points define a surface:
        double x1 = crater.xLocation;
        double y1 = crater.yLocation;
        double z1 = getSurfaceElevationAtPoint(x1, y1);

        double x2 = crater.xLocation + 2 * crater.finalRadius;
        double y2 = crater.yLocation - 2 * crater.finalRadius;
        double z2 = getSurfaceElevationAtPoint(x2, y2);

        double x3 = crater.xLocation - 2 * crater.finalRadius;
        double y3 = crater.yLocation - 2 * crater.finalRadius;
        double z3 = getSurfaceElevationAtPoint(x3, y3);

        Eigen::Vector3d p1(x1, y1, z1);
        Eigen::Vector3d p2(x2, y2, z2);
        Eigen::Vector3d p3(x3, y3, z3);

        Eigen::Vector3d v1 = p2 - p1;
        Eigen::Vector3d v2 = p3 - p1;

        Eigen::Vector3d surfNormal = v1.cross(v2);
        surfNormal.normalize();

        Eigen::Vector3d up(0, 0, 1); 
        Eigen::Vector3d east(1, 0, 0);

        double slope = acos(surfNormal.dot(up)); 
        double aspect = acos(surfNormal.dot(east));

        // if slope > 90, wrap angle:
        if (slope > (M_PI / 2)) {
                slope = M_PI - slope;
        }
        return {slope, aspect};
}

///////////////////
// Crater profiles:
///////////////////
// Parabolic:
double Grid::craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle){
        normAngle = 0; // To do
        return depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2));

}

// Bowl shaped:
double Grid::craterSphericalDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle){
        double craterDepth = 2 * craterRadius * depthToDiameter;
        double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;
        normAngle = 0; // To do

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

/////////////////
// Ice Stability:
/////////////////
// The temperature depndent surface loss rate based on WMB 1962
// Returns loss rate in kg Ma^-1

/////////
// ********************
// ***IN DEVELOPMENT***
// ********************
// Calculate Permanent Shadow: permanent shadow is accounted for by calculating the noontime shadow.
// Since there is no preferred directionallity, the solar azimuth is always assumed to be zero.
// This greatly simplifies calculations, as the Sun is always in the direction of the grid, so that
// no interpolation is needed. In Future versions this will be improved to account for more accurate
// permanent shadows.
/////////
// Calculate transient shadow:
bool Grid::calculatePermanentShadow(int faceti, int facetj, double solarZenithAngle) {
        // GOT TO ITERATE OVER ALL SURFACE, AS CHANGING FACETS CAN AFFECT *OTHER* FACET SHADOW

        double heightDiff;
        double horDistanceDiff;
        double angleBetweenFacets;
        double facetElevation = subsurfColumns[facetj][faceti].getSurfaceElevation();

        for (int i = faceti + 1; i < gridSize; ++i) {
                heightDiff = facetElevation - subsurfColumns.at(facetj).at(i).getSurfaceElevation();
                horDistanceDiff = fabs(x.at(facetj) - x.at(i));
                angleBetweenFacets = heightDiff / horDistanceDiff;

                if (angleBetweenFacets > tan( (90 - solarZenithAngle) / 180 * M_PI )) {
                        return true;
                }
        }
        return false;
}

/////////
// Print:
/////////
// Print surface to file:
void Grid::printSurface(int index){
        std::ofstream xFile;
        std::ofstream yFile;
        std::ofstream elevationFile;

        std::string index_str = std::string(
                std::to_string(int(endTime / printTimeStep)).length() -
                std::to_string(index).length(), '0'
                ) + std::to_string(index);

        std::string elevationFileName = "./output/elevation_" + index_str + ".txt";
        elevationFile.open(elevationFileName, std::ios_base::out);

        if (index == 0) {
                // For the last print, print x,y vectors:
                std::string xFileName = "./output/x.txt";
                std::string yFileName = "./output/y.txt";
                xFile.open(xFileName, std::ios_base::out);
                yFile.open(yFileName, std::ios_base::out);
                for (int i = 0; i < gridSize; ++i) {
                        xFile << x.at(i) << ",";
                        yFile << y.at(i) << ",";
                }
        }

        for (int i = 0; i < gridSize; ++i) {
                for (int j = 0; j < gridSize; ++j) {
                        elevationFile << subsurfColumns.at(j).at(i).getSurfaceElevation() << ",";
                }
                elevationFile << "\n";
        }
}

// Print shadow matrix to file:
void Grid::printShadow(int index){
        std::ofstream xFile;
        std::ofstream yFile;
        std::ofstream shadowFile;

        std::string elevationFileName = "./output/shadow_" + std::to_string(index) + ".txt";
        shadowFile.open(elevationFileName, std::ios_base::out);

        if (index == 0) {
                // For the last print, print x,y vectors:
                std::string xFileName = "./output/x.txt";
                std::string yFileName = "./output/y.txt";
                xFile.open(xFileName, std::ios_base::out);
                yFile.open(yFileName, std::ios_base::out);
                for (int i = 0; i < gridSize; i++) {
                        xFile << x.at(i) << ",";
                        yFile << y.at(i) << ",";
                }
        }

        for (int i = 0; i < gridSize; ++i) {
                for (int j = 0; j < gridSize; ++j) {
                        shadowFile << subsurfColumns.at(j).at(i).isPermShadow << ",";
                }
                shadowFile << "\n";
        }
}

// Print subsurface to file:
void Grid::printGrid(int index){

        SubsurfColumn col =  SubsurfColumn();
        std::fstream outputFile;
        std::string output_filename = "./output/output_" + std::to_string(index) + ".out";

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
