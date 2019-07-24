// Class for impactor objects
class Impactor {
public:
  double radius;
  double mass;
  double velocity;

  Impactor(); // Initialize an impactor from a cumulative distribution.
  Impactor(double radius); // Initialize an impactor with given parameters.

private:
  double calcMass(double radius);
};
