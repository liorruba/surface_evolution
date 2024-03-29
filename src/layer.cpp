// Class for layer object
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
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

        normalizeComposition(); // Normalize layer composition to 1.
}

// Second constructor for a dummy layer used in printing the column. This layer is not normalized.
Layer::Layer(double dummy_numOfLayers, double dummy_elevation){
        thickness = dummy_numOfLayers;
        regolithFraction = dummy_elevation;
        iceFraction = -1;
        sootFraction = -1;
}

// Normalize layer composition:
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
        double eps = 1e-2;

        if (
                (fabs(regolithFraction - layer.regolithFraction) < eps) && 
                (fabs(iceFraction - layer.iceFraction) < eps) && 
                (fabs(sootFraction - layer.sootFraction) < eps)
                ) {
                return true;
        }
        return false;
}

// Change layer composition
void Layer::changeComposition(double _regolithFraction, double _iceFraction, double _sootFraction){
        regolithFraction = _regolithFraction;
        iceFraction = _iceFraction;
        sootFraction = _sootFraction;

        normalizeComposition();
}

// Consolidate two layers:
void Layer::consolidate(Layer otherLayer) {
        regolithFraction = (thickness * regolithFraction + otherLayer.thickness * otherLayer.regolithFraction);
        iceFraction = (thickness * iceFraction + otherLayer.thickness * otherLayer.iceFraction);
        sootFraction = (thickness * sootFraction + otherLayer.thickness * otherLayer.sootFraction);
        thickness += otherLayer.thickness;

        normalizeComposition();
}

// Check if layer is empty
bool Layer::isEmpty() {
        if (regolithFraction == 0 && iceFraction == 0 && sootFraction == 0) {
                return true;
        }
        return false;
}

// Shrink a layer:
void Layer::shrink(double thicknessToRemove) {
        thickness -= thicknessToRemove;
        if (thickness < 0) {
                thickness = 0;
        }
}

// Print layer
void Layer::print(bool isNiceInterface) {
        if  (isNiceInterface) {
                std::cout << "T: " << thickness << " R: " << regolithFraction << " I: " << iceFraction << " S: " << sootFraction << std::endl;
        }
        else {
                std::cout << thickness << ", " << regolithFraction << ", " << iceFraction << ", " << sootFraction << std::endl;
        }
}
