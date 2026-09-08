#pragma once
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "layer.hpp"
#include "impactor.hpp"
#include "crater.hpp"
#include "subsurf_column.hpp"
#include "spatial_index.hpp"

// Lightweight record of a formed crater, kept to track crater degradation over time.
struct CraterRecord {
        double x;
        double y;
        double finalRadius;
        double finalDepth;       // current rim-to-floor depth
        double finalDepth_init;  // depth at formation
        double floorElevation;   // current surface elevation at the crater center
        bool isVisible;          // false once the depth changed by more than 50% either way
};

// Reference plane z = z0 + sx * (x - x0) + sy * (y - y0) fitted to the pre-impact surface.
struct ReferencePlane {
        double x0, y0, z0, sx, sy;
        double at(double x, double y) const { return z0 + sx * (x - x0) + sy * (y - y0); }
};

////////////////////////
// Grid class definition
////////////////////////
class Grid {

public:
double area;
std::vector<double> x;   // cell-center coordinates, m
std::vector<double> y;
std::vector< std::vector<SubsurfColumn> > subsurfColumns;   // indexed [j (y index)][i (x index)]
// Elevation changes from settling that are still too small to be applied to the columns (|change|
// below minimumLayerThickness), indexed j * gridSize + i. The true surface is column + pending;
// the pending part is applied once it accumulates past the threshold, and flushed before output.
std::vector<double> pendingElevation;

Grid(std::vector< std::vector<double> > _initLayersList, std::vector<int8_t> _pxIdxMat);
// Forms a crater: carves the cavity, adds the rim and (if enabled) the ejecta blanket, records the
// crater and updates the depths of older craters affected by it.
void formCrater(Crater &crater);
void thresholdSlopes(double angleOfRepose);
// Settle the surface in the (periodically wrapped) square of the given half-width around (xc, yc):
// an optional downslope-diffusion dose K(l) [m^2] as a function of the distance l from (xc, yc)
// (periodic minimum-image distances when periodicDistance, plain distances otherwise), then the
// collapse of slopes above slopeOfRepose (rise over run). The changes are applied to the columns.
void updateCraterDepthsWithin(double x, double y, double reach);
void refreshCraterDepths();   // re-read the floor elevation of every visible crater (after whole-domain shakes)
void flushPendingElevation();   // apply every pending settling change (before output)
static double settleTimers[6];  // seconds spent in settleRegion: copy, dose, diffusion, relaxation, apply, crater update (for the log)
static void applyElevationChange(SubsurfColumn &column, double change);
void settleRegion(double xc, double yc, double halfWidth, const std::function<double(double)> *dose, double slopeOfRepose, bool periodicDistance);
void printSurface(int index, bool isfinal);
void printSubsurface(int index);
void printIntegratedSubsurface(double depth, int index);
void printExistingCratersToHistogram(double bins);
void printExistingCraters();
void sublimateIce();
void depositLayer(Layer layer);
size_t numberOfVisibleCraters() const;
const std::vector<CraterRecord>& craterRecords() const { return craters; }

private:
int gridSize;
std::vector< std::vector<double> > initLayersList;
std::vector<int8_t> pixelIndexMatrix;
std::vector<CraterRecord> craters;
BucketGrid craterIndex;

std::vector< std::vector<SubsurfColumn> > initializeSubsurface();
void carveCavity(Crater &crater);
void emplaceRimDropoff(const Crater &crater);
void emplaceEjecta(const Crater &crater);
void registerCrater(const Crater &crater);
void updateExistingCratersDepth(const Crater &crater);
ReferencePlane fitReferencePlane(const Crater &crater);
void footprintIndexRange(double xc, double yc, double halfSize, int &iInit, int &iFinal, int &jInit, int &jFinal) const;
double cavityDepthProfile(double craterRadius, double craterDepth, double distanceFromCraterCenter) const;
double craterParabolicDepthProfile(double craterRadius, double craterDepth, double distanceFromCraterCenter) const;
double craterSphericalDepthProfile(double craterRadius, double craterDepth, double distanceFromCraterCenter) const;
double getSurfaceElevationAtPoint(double x, double y) const;
std::vector<double> surfaceElevationMap() const;
void relaxSlopes(std::vector<double> &z, double maxSlope) const;
static void relaxBuffer(std::vector<double> &z, long ni, long nj, bool periodicI, bool periodicJ, double maxSlope, double cellSize);
static void diffuseBuffer(std::vector<double> &z, const std::vector<double> &dose, long ni, long nj, bool periodicI, bool periodicJ, double cellSize);
void writeMatrix(const std::string &fileName, const std::vector< std::vector<double> > &matrix);
};
