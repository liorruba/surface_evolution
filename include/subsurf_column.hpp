#pragma once
#include <vector>
#include "layer.hpp"

// Class for subsurf column object. Layers are stored bottom-up: subsurfLayers.back() is the surface.
// The bottom layer is the basement: it is never removed and is treated as extending indefinitely
// downward, so a column always has a composition even if a crater digs below the initial stack.
class SubsurfColumn {
public:
  std::vector<Layer> subsurfLayers;
  bool isPermShadow;

  SubsurfColumn();
  double getSurfaceElevation() const;
  void addLayer(Layer layer);
  void removeMaterial(double depthToRemove);
  Layer integrateColumnComposition(double depthToIntegrate);
  void print(bool isNiceInterface = true);

private:
  double surfaceElevation;
};
