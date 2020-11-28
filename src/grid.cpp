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

// A constructor for the surface class. Creates a surface with (size * size) elements.
Grid::Grid(){
        gridSize = regionWidth / resolution;
        x = linspace(-regionWidth/2, regionWidth/2, gridSize);
        y = linspace(-regionWidth/2, regionWidth/2, gridSize);
        area = pow(regionWidth, 2);
        subsurfColumns = initializeSubsurface();
}

// Initialize a matrix of columns (2d vector):
std::vector< std::vector<SubsurfColumn> > Grid::initializeSubsurface(){
        std::vector< std::vector<SubsurfColumn> > buffMat(gridSize,std::vector<SubsurfColumn>(gridSize));
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
        int iInit = floor( ((crater.xLocation + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(iInit < 0) iInit = 0;
        int iFinal = ceil( ((crater.xLocation + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(iFinal > gridSize) iFinal = gridSize;
        int jInit = floor( ((crater.yLocation + regionWidth/2) - ejectaSpread * crater.finalRadius) / resolution ); if(jInit < 0) jInit = 0;
        int jFinal = ceil( ((crater.yLocation + regionWidth/2) + ejectaSpread * crater.finalRadius) / resolution ); if(jFinal > gridSize) jFinal = gridSize;

        for (int i = iInit; i < iFinal; i++) {
                for (int j = jInit; j < jFinal; j++) {
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

                                // Melt production
                                if (isProduceMelt) {
                                        // std::cout << craterProfile << std::endl;
                                        // std::cout << crater.meltHeight - crater.rimHeight << std::endl;
                                        if (craterProfile > crater.meltHeight - crater.rimHeight) {
                                                std::cout << crater.meltHeight - crater.rimHeight << std::endl;
                                                craterProfile = crater.meltHeight - crater.rimHeight;
                                        }
                                }

                                craterProfile = elevationAtCenter - craterProfile;
                                double leftToRemove = subsurfColumns.at(j).at(i).getSurfaceElevation() - craterProfile;

                                // If the thickness of the materialToRemove is positive, remove material
                                if (leftToRemove >= 0) {
                                        crater.ejectedMass.consolidate(subsurfColumns.at(j).at(i).integrateColumnComposition(leftToRemove));
                                        subsurfColumns.at(j).at(i).removeMaterial(leftToRemove);
                                }
                                // If the thickness of the materialToRemove is negative, add material
                                else {
                                        subsurfColumns.at(j).at(i).addLayer(Layer(fabs(leftToRemove),1,0,0));
                                }
                        }

                        // If outside the crater, add rim dropoff
                        if (distanceFromCraterCenter > crater.finalRadius) {
                                double rimDropoff = crater.rimHeight * pow((distanceFromCraterCenter / crater.finalRadius), -rimDropoffExponent);

                                // Add a layer of pure regolith (1,0,0) as the rim dropoff.
                                subsurfColumns.at(j).at(i).addLayer(Layer(rimDropoff,1,0,0));
                        }
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

        for (int i = iInit; i < iFinal; i++) {
                for (int j = jInit; j < jFinal; j++) {
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

double Grid::calculateSlope(const Crater crater){
        // Three points define a surface:
        double x1 = crater.xLocation;
        double y1 = crater.yLocation;
        double z1 = getSurfaceElevationAtPoint(x1, y1);

        double x2 = crater.xLocation + crater.finalRadius;
        double y2 = crater.yLocation - crater.finalRadius;
        double z2 = getSurfaceElevationAtPoint(x2, y2);

        double x3 = crater.xLocation - crater.finalRadius;
        double y3 = crater.yLocation - crater.finalRadius;
        double z3 = getSurfaceElevationAtPoint(x3, y3);

        point3 v1(x1 - x3, y1 - y3, z1 - z3);
        point3 v2(x1 - x2, y1 - y2, z1 - z2);

        point3 v3 = bg::cross_product(v1, v2);
        return xyPlaneVecAngle(std::vector<double> {v3.get<0>(), v3.get<1>(), v3.get<2>()});
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

        for (int i = faceti + 1; i < gridSize; i++) {
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

        std::string elevationFileName = "./output/elevation_" + std::to_string(index) + ".txt";
        elevationFile.open(elevationFileName, std::ios_base::out);

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

        for (int i = 0; i < gridSize; i++) {
                for (int j = 0; j < gridSize; j++) {
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

        for (int i = 0; i < gridSize; i++) {
                for (int j = 0; j < gridSize; j++) {
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

        for (size_t i = 0; i < subsurfColumns.size(); i++) {
                for (size_t j = 0; j < subsurfColumns.size(); j++) {
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
        addLogEntry("Printing existing craters histrogram.");
        Histogram hist(minimumImpactorDiameter * 10, regionWidth, bins);

        for (auto cm : cratersDict) {
                // Add the crater diameter to the histogram
                hist.add(2 * cm.second->finalRadius);
        }

        hist.print("./output/existing_craters_histogram.txt");
}

// Print existing craters to histogram:
void Grid::printExistingCraters(){
        addLogEntry("Printing existing craters to file.");

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
