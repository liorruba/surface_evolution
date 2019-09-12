// Class for layer object
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/utility.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/log.hpp"

Layer::Layer(double _thickness, double _regolithFraction, double _iceFraction, double _sootFraction){
  thickness = _thickness;
  regolithFraction = _regolithFraction;
  iceFraction = _iceFraction;
  sootFraction = _sootFraction;

  normalizeComposition(); // Normalize layer composition.
}

// Normlize later composition:
void Layer::normalizeComposition() {
  double sum = regolithFraction + iceFraction + sootFraction;
  if (sum == 0)
    return;
  else {
    regolithFraction /= sum;
    iceFraction /= sum;
    sootFraction /= sum;
  }
}

// Compare the composition of two layers:
bool Layer::compareComposition(Layer layer) {
  if (regolithFraction ==  layer.regolithFraction && iceFraction == layer.iceFraction && sootFraction == layer.sootFraction){
    return true;
  }
  return false;
}

// Consolidate two layers:
void Layer::consolidate(Layer otherLayer) {
  regolithFraction = (thickness * regolithFraction + otherLayer.thickness * otherLayer.regolithFraction);
  iceFraction = (thickness * iceFraction + otherLayer.thickness * otherLayer.iceFraction);
  sootFraction = (thickness * sootFraction + otherLayer.thickness * otherLayer.sootFraction);
  thickness += otherLayer.thickness;

  normalizeComposition();
}

// Shrink a layer:
void Layer::shrink(double thicknessToRemove) {
  thickness -= thicknessToRemove;
  if (thickness < 0) {
    thickness = 0;
  }
}

// Print layer
void Layer::print() {
  std::cout << "T: " << thickness << " R: " << regolithFraction << " I: " << iceFraction << " S: " << sootFraction << std::endl;
}
