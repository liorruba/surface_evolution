// This class defines the surface properties:
class Surface {
public:
  Surface(double regionWidth, double resolution);
  double area;
  std::vector<double> x;
  std::vector<double> y;
  std::vector< std::vector<SubsurfColumn> > subsurfColumns;
  void formCrater(Crater crater);
  void emplaceEjecta(Crater crater);
  void print(int index);


private:
  int gridSize;
  double craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle);
  double craterBowlShapedDepthProfile(double craterRadius, double distanceFromCraterCenter, double normAngle);
  double rimHeight(double craterRadius, double distanceFromCraterCenter);
  std::vector< std::vector<SubsurfColumn> > initializeSubsurface();
};
