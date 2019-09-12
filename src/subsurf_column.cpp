// Class for column object
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/utility.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/subsurf_column.hpp"

SubsurfColumn::SubsurfColumn() {
  subsurfLayers.push_back(Layer(initialThickness, 1, 0, 0));

  surfaceElevation = initialThickness;
}

double SubsurfColumn::get_surfaceElevation() {
  return surfaceElevation;
}

////
// Add material to column:
void SubsurfColumn::addLayer(Layer newLayer) {
  // If new layer composition is the same as the topmost layer, consolidate:

  if (newLayer.compareComposition(subsurfLayers.back())) {
    subsurfLayers.back().consolidate(newLayer);
  }
  // if new layer is smaller than the threshold, also consolidate with topmost layer:
  else if (newLayer.thickness < resolution / initialThickness) {
    subsurfLayers.back().consolidate(newLayer);
  }
  // Else, add the layer on top:
  else {
    subsurfLayers.push_back(newLayer);
  }

  surfaceElevation += newLayer.thickness;
}

////
// Remove material from column:
void SubsurfColumn::removeMaterial(double depthToRemove) {
  // Change the surface elevation:
  surfaceElevation -= depthToRemove;

  if (surfaceElevation < 0) {
    removeAllLayers();
    surfaceElevation = 0;
    return;
  }

  // Remove layers to some depth:
  while ((depthToRemove >= subsurfLayers.back().thickness) && (depthToRemove > 0)) {
    depthToRemove -= subsurfLayers.back().thickness;
    subsurfLayers.pop_back();
  }

  subsurfLayers.back().shrink(depthToRemove);
}

////
// Integrate column composition, return single normalized layer:
Layer SubsurfColumn::integrateColumnComposition(double depthToIntergrate) {
  Layer buffLayer = Layer(0,0,0,0);
  // std::cout << std::endl;
  // subsurfLayers.begin() -> print();
  // subsurfLayers.end() -> print();
  // std::cout << std::endl;
  //
  std::vector<Layer>::iterator it = subsurfLayers.end();
  // (*it).print();
  // return buffLayer;

  it = subsurfLayers.end();
  it--;
  while(depthToIntergrate > (it -> thickness) && (depthToIntergrate > 0)){
    buffLayer.consolidate(*it);
    depthToIntergrate -= (it -> thickness);
    it--;
  }
  Layer remainingLayer = *it;
  remainingLayer.shrink(fabs((it -> thickness) - depthToIntergrate));
  buffLayer.consolidate(remainingLayer);
  return buffLayer;
  //
  // // Consolidate with the remaining material:
  // Layer remainingLayer = subsurfLayers.back();
  // remainingLayer.shrink(fabs(subsurfLayers.back().thickness - depthToIntergrate));
  // buffLayer.consolidate(remainingLayer);


}

// Print layers in column:
void SubsurfColumn::print() {

  for (long i = (subsurfLayers.size() - 1); i >= 0; i--)
  {
    std::cout << "***" << " Layer " << i << " ***" << std::endl;
    subsurfLayers[i].print();
  }
  std::cout << "Number of layers in column: " << subsurfLayers.size() << std::endl;
  std::cout << "Surface elevation: " << surfaceElevation << std::endl;
}

// Remove all the layers in the subsurface column:
void SubsurfColumn::removeAllLayers() {
  subsurfLayers.clear();
}
