// Class for column object
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/layer.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/log.hpp"

// The empty constructor returns a subsurface column with just regolith (the basement layer).
SubsurfColumn::SubsurfColumn() {
        subsurfLayers.push_back(Layer(initialThickness, 1, 0, 0));
        // Initially set the surface elevation to the initial thickness
        surfaceElevation = initialThickness;
        isPermShadow = false;
}


double SubsurfColumn::getSurfaceElevation() const {
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

        // If the new layer has the same composition as the topmost layer, or is thinner than the
        // minimum layer thickness (a veneer that should not redefine the surface composition), mix it
        // into the topmost layer. Mass is conserved either way:
        if (newLayer.thickness < minimumLayerThickness || newLayer.compareComposition(subsurfLayers.back())) {
                subsurfLayers.back().consolidate(newLayer);
        }

        // Else, add the layer on top:
        else {
                subsurfLayers.push_back(newLayer);
        }

        surfaceElevation += newLayer.thickness;
}

////
// Remove material from the top of the column. The bottom (basement) layer is never removed, so
// the column always keeps a composition; the surface elevation is lowered by the full amount.
void SubsurfColumn::removeMaterial(double depthToRemove) {
        if (depthToRemove <= 0) {
                return;
        }

        // Change the surface elevation:
        surfaceElevation -= depthToRemove;

        // Peel off whole layers from the top, keeping the basement:
        while (subsurfLayers.size() > 1 && depthToRemove >= subsurfLayers.back().thickness) {
                depthToRemove -= subsurfLayers.back().thickness;
                subsurfLayers.pop_back();
        }

        subsurfLayers.back().shrink(depthToRemove);
}

////
// Integrate the column composition from the surface down to some depth and return it as a single
// normalized layer whose thickness is the integrated depth. The basement layer is treated as
// extending indefinitely downward.
Layer SubsurfColumn::integrateColumnComposition(double depthToIntegrate) {
        Layer buffLayer = Layer(0, 0, 0, 0);

        if (depthToIntegrate <= 0) {
                const Layer &top = subsurfLayers.back();
                return Layer(0, top.regolithFraction, top.iceFraction, top.sootFraction);
        }

        for (size_t k = subsurfLayers.size(); k-- > 0;) {
                const Layer &layer = subsurfLayers[k];
                double take = (k == 0) ? depthToIntegrate : std::min(depthToIntegrate, layer.thickness);

                if (take > 0) {
                        buffLayer.consolidate(Layer(take, layer.regolithFraction, layer.iceFraction, layer.sootFraction));
                        depthToIntegrate -= take;
                }

                if (depthToIntegrate <= 0) {
                        break;
                }
        }

        if (buffLayer.isEmpty()) {
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
