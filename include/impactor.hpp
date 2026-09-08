#pragma once

// Class for impactor objects
class Impactor {
public:
  double radius;
  double mass;
  double velocity;
  double density;   // kg/m^3; the global impactorDensity unless given (ejecta fragments are target material)

  Impactor(); // Initialize an impactor from a cumulative distribution.
  Impactor(double radius); // Initialize an impactor with given radius (mean impact velocity, impactor density).
  Impactor(double radius, double velocity, double density); // Fully specified (used for ejecta fragments).

private:
  double calcMass(double radius) const;
};
