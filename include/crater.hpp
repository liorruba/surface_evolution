// Class for crater object
class Crater {
public:
  double xLocation;
  double yLocation;
  double transientRadius;
  double transientRadiusGravity;
  double finalRadius;
  double finalDepth;
  double finalDepth_init;
  double rimHeight;
  double floorElevation;
  int numberOfSecondaries;
  Layer ejectedMass;
  std::vector<double> ejectaDistance;
  std::vector<double> ejectaThickness;

  Crater(Impactor impactor);
  Crater(Impactor impactor, double xLocation, double yLocation);
  Crater(double _finalRadius, double _xLocation, double _yLocation);

private:
  double calcTransientVolume(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientVolumeGravity(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientCraterRadius(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientCraterRadiusGravity(double impactorRadius, double impactorMass, double impactVelocity);
  double calcFinalCraterRadius();
  void calcEjectaThickness(Impactor impactor);
};
