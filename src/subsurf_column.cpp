// Class for column object
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cassert>
#include <stdexcept>
#include "../include/regolit_main.hpp"
#include "../include/impactor.hpp"
#include "../include/utility.hpp"
#include "../include/layer.hpp"
#include "../include/crater.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/log.hpp"

// The empty constructor returns a sursurface column with just regolith
SubsurfColumn::SubsurfColumn() {
        subsurfLayers.push_back(Layer(initialThickness, 1, 0, 0));
        // Initially set the surface elevation to the initial thickness
        surfaceElevation = initialThickness;
}


double SubsurfColumn::getSurfaceElevation() {
        return surfaceElevation;
}

////
// Add material to column:
void SubsurfColumn::addLayer(Layer newLayer) {
        // Thickness cannot be a zero or a negative number
        if (newLayer.thickness <= 0){
                throw std::invalid_argument("Cannot add layer with thickness <= 0.");
        }
        // Layer cannot be empty
        if (newLayer.isEmpty()){
                throw std::invalid_argument("Cannot add layer with null composition.");
        }

        // If new layer composition is the same as the topmost layer, consolidate:
        if (newLayer.compareComposition(subsurfLayers.back())) {
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
        assert(depthToRemove >= 0);

        // Change the surface elevation:
        surfaceElevation -= depthToRemove;

        if (surfaceElevation < 0) {
                removeAllLayers();
                surfaceElevation = 0;
                return;
        }

        // Remove layers to some depth:
        while ((depthToRemove >= subsurfLayers.back().thickness)) {
                depthToRemove -= subsurfLayers.back().thickness;
                subsurfLayers.pop_back();
        }

        subsurfLayers.back().shrink(depthToRemove);
}

////
// Integrate column composition, return single normalized layer:
Layer SubsurfColumn::integrateColumnComposition(double depthToIntergrate) {
        Layer buffLayer = Layer(0,0,0,0);

        std::vector<Layer>::iterator it = subsurfLayers.end();

        it = subsurfLayers.end();
        it--;

        while(depthToIntergrate > (it->thickness) && (depthToIntergrate > 0)) {
                buffLayer.consolidate(*it);
                depthToIntergrate -= (it->thickness);
                it--;
        }
        Layer remainingLayer = *it;
        remainingLayer.shrink(fabs((it->thickness) - depthToIntergrate));
        buffLayer.consolidate(remainingLayer);

        if ((buffLayer.regolithFraction == 0) && (buffLayer.sootFraction == 0) && (buffLayer.iceFraction == 0)){
                throw std::invalid_argument("Integrated layer composition is null.");
        }

        return buffLayer;
}

// Print layers in column:
void SubsurfColumn::print(bool isNiceInterface) {

        if (isNiceInterface) {
                for (long i = (subsurfLayers.size() - 1); i >= 0; i--)
                {
                        std::cout << "***" << " Layer " << i << " ***" << std::endl;
                        subsurfLayers[i].print(isNiceInterface);
                }
                std::cout << "Number of layers in column: " << subsurfLayers.size() << std::endl;
                std::cout << "Surface elevation: " << surfaceElevation << std::endl;
        }
        else {
                for (long i = (subsurfLayers.size() - 1); i >= 0; i--)
                {
                        std::cout << i << ", " << std::endl;
                        subsurfLayers[i].print(isNiceInterface);
                        std::cout << std::endl;
                }
        }
}

// Remove all the layers in the subsurface column:
void SubsurfColumn::removeAllLayers() {
        subsurfLayers.clear();
}
