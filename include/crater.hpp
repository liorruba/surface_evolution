// Class for crater object
class Crater {
public:
  double xPosition;
  double yPosition;
  double transientRadius;
  double transientRadiusGravity;
  double finalRadius;
  int numberOfSecondaries;
  std::vector<double> ejectaDistance;
  std::vector<double> ejectaThickness;

  Crater(Impactor impactor);
  Crater(Impactor impactor, double xPosition, double yPosition);
  Crater(double _finalRadius, double _xPosition, double _yPosition);

private:
  double calcTransientVolume(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientVolumeGravity(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientCraterRadius(double impactorRadius, double impactorMass, double impactVelocity);
  double calcTransientCraterRadiusGravity(double impactorRadius, double impactorMass, double impactVelocity);
  double calcFinalCraterRadius();
  void calcEjectaThickness(Impactor impactor);
};
