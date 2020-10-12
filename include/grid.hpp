#define BOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING 1

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
// #include <boost/geometry/index/indexable.hpp>
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
    point operator()(const CraterRef& c) const { return c.location; }
};

typedef CraterRef value;
typedef bgi::rtree< value, bgi::quadratic<16> > rtree_t;

////////////////////////
// Grid class definition
////////////////////////
class Grid {

public:
  Grid();
  double area;
  std::vector<double> x;
  std::vector<double> y;
  std::vector< std::vector<SubsurfColumn> > subsurfColumns;
  void formCrater(Crater &crater);
  void emplaceEjecta(Crater &crater);
  void printSurface(int index);
  void printShadow(int index);
  void printGrid(int index);
  void updateExistingCratersDepth(Crater &crater);
  bool calculatePermanentShadow(int faceti, int facetj, double solarZenith);

private:
  int gridSize;
  double craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle);
  double craterSphericalDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle);
  double rimHeight(double craterRadius, double distanceFromCraterCenter);
  double getSurfaceElevationAtPoint(double x, double y);
  double calculateSlope(const Crater crater);
  std::vector< std::vector<SubsurfColumn> > initializeSubsurface();
  std::vector<Crater *> cratersVector;

  rtree_t rtree;
};
