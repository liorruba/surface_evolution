#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING 1

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// This class defines the surface properties:
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// For the r-tree struct
typedef bg::model::point<double, 2, bg::cs::cartesian> point;

// For later, when doing cross product
typedef bg::model::point<double, 3, bg::cs::cartesian> point3;

///////////////
// R-TREE SETUP
///////////////
// Struct for r-tree which stores index of insertion and location, to be later
// accessed craters_vector vector.
struct CraterRef {
        size_t index;
        point location;
};

template <>
struct bgi::indexable<CraterRef>
{
        typedef point result_type;
        point operator()(const CraterRef& c) const {
                return c.location;
        }
};

struct my_equal {
        using result_type = bool;
        bool operator() (CraterRef const& v1, CraterRef const& v2) const {
                return v1.index == v2.index;
        }
};

// typedef CraterRef value;
typedef bgi::rtree< CraterRef, bgi::quadratic<16>, bgi::indexable<CraterRef>,  my_equal> rtree_t;

////////////////////////
// Grid class definition
////////////////////////
class Grid {

public:   
double area;
std::vector<double> x;
std::vector<double> y;
std::vector< std::vector<SubsurfColumn> > subsurfColumns;

Grid(std::vector< std::vector<double> > _initLayersList, std::vector<int8_t> _pxIdxMat);
void formCrater(Crater &crater);
void emplaceEjecta(Crater &crater);
void thresholdSlopes(double angleOfRepose);
void printSurface(int index, bool isfinal);
void printSubsurface(int index);
void printIntegratedSubsurface(double depth, int index);
void updateExistingCratersDepth(Crater &crater);
void printExistingCratersToHistogram(double bins);
void printExistingCraters();
void sublimateIce();
void depositLayer(Layer layer);
// bool calculatePermanentShadow(int faceti, int facetj, double solarZenith);

private:
int gridSize;
long numberOfCratersOnGrid;
double craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter);
double craterSphericalDepthProfile(double craterRadius, double distanceFromCraterCenter);
double rimHeight(double craterRadius, double distanceFromCraterCenter);
double getSurfaceElevationAtPoint(double x, double y);
bool compareMaxSlope(std::vector< std::vector<double> > slopeMap, double maxSlope);
std::vector < std::vector<double> > computeGradient(std::vector < std::vector<double> > Z);
std::vector< std::vector<double> > linearSlopeDiffusion(std::vector< std::vector<double> > Z, double maxSlope, double K);
std::tuple<double, double> calculateSlope(const Crater crater);
std::vector< std::vector<double> > initLayersList;
std::vector<int8_t> pixelIndexMatrix;

std::vector< std::vector<SubsurfColumn> > initializeSubsurface();
std::map<size_t, Crater *> cratersDict;

rtree_t craters_rtree;
};
