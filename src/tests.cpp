#include <iostream>
#include <cmath>
#include <vector>
#include "../include/regolit_main.hpp"
#include "../include/tests.hpp"
#include "../include/layer.hpp"
#include "../include/subsurf_column.hpp"
#include "../include/impactor.hpp"
#include "../include/crater.hpp"
#include "../include/grid.hpp"
#include "../include/utility.hpp"

bool tests() {
        bool all_test_passed = 1;

        //////////////////////
        // TEST LAYER CLASS //
        //////////////////////

        // Create a layer:
        Layer testLayer = Layer(1, 1, 0, 0);
        Layer testLayer_iden = testLayer;
        Layer testLayer_diff = Layer(2, 0, 1, 1);

        // Try to read layer parameters:
        if (testLayer.thickness == 1 &&
            testLayer.regolithFraction == 1 &&
            testLayer.iceFraction == 0 &&
            testLayer.sootFraction == 0)
        {
                std::cout << "Test passed: creating a layer." << std::endl;
        }
        else {
                std::cout << "Test failed: creating a layer." << std::endl;
                all_test_passed = 0;
        }

        // Compare two layers:
        if (testLayer.compareComposition(testLayer_iden) == true &&
            testLayer.compareComposition(testLayer_diff) == false) {
                std::cout << "Test passed: comparing two layers." << std::endl;
        }
        else {
                std::cout << "Test failed: comparing two layers." << std::endl;
                all_test_passed = 0;
        }

        // Consolidate two layers:
        testLayer.consolidate(testLayer_diff);
        double eps = 1e-3;

        if (fabs(testLayer.thickness - 3) < eps &&
            fabs(testLayer.regolithFraction - 1.0/3.0) < eps &&
            fabs(testLayer.iceFraction - 1.0/3.0) < eps &&
            fabs(testLayer.sootFraction - 1.0/3.0) < eps) {
                std::cout << "Test passed: consolidate layers." << std::endl;
        }
        else {
                std::cout << "Test failed: consolidate layers." << std::endl;
                all_test_passed = 0;
        }

        // Shrink a layer:
        testLayer.shrink(2);
        if (fabs(testLayer.thickness - 1) < eps) {
                std::cout << "Test passed: shrink a layer." << std::endl;
        }
        else {
                std::cout << "Test failed: shrink a layer." << std::endl;
                all_test_passed = 0;
        }

        ///////////////////////////////
        // TEST SUBSURF COLUMN CLASS //
        ///////////////////////////////
        SubsurfColumn column = SubsurfColumn();

        // Add layer tests:
        // Same layer, should consolidate:
        column.addLayer(testLayer);
        if (fabs(column.getSurfaceElevation() - (initialThickness + 1)) < eps && fabs(column.subsurfLayers.size() - 2) < eps) {
                std::cout << "Test passed: add a layer." << std::endl;
        }
        else {
                std::cout << "Test failed: add a layer." << std::endl;
                all_test_passed = 0;
        }

        // Add another identical layer, should consolidate:
        column.addLayer(testLayer);
        if (fabs(column.getSurfaceElevation() - (initialThickness + 2)) < eps && fabs(column.subsurfLayers.size() - 2) < eps) {
                std::cout << "Test passed: add an identical layer." << std::endl;
        }
        else {
                std::cout << "Test failed: add an identical layer." << std::endl;
                all_test_passed = 0;
        }

        // Add a different layer but smaller than the threshold:
        column.addLayer(Layer(1.0/101.0, 1, 0, 0));
        if (fabs(column.getSurfaceElevation() -  (initialThickness + 2.01)) < eps && fabs(column.subsurfLayers.size() - 2) < eps) {
                std::cout << "Test passed: add a layer thinner than threshold." << std::endl;
        }
        else {
                std::cout << "Test failed: add a layer thinner than threshold." << std::endl;
                all_test_passed = 0;
        }

        // Remove material:
        // First add some layers:
        column.addLayer(Layer(1, 1, 0, 0));
        column.addLayer(Layer(2, 1, 1, 0));
        column.addLayer(Layer(.3, 1, 0, 1));
        column.removeMaterial(2.5);
        if (fabs(column.getSurfaceElevation() - (initialThickness + 2.81)) < eps && fabs(column.subsurfLayers.size() - 3) < eps) {
                std::cout << "Test passed: remove layer." << std::endl;
        }
        else {
                std::cout << "Test failed: remove layer." << std::endl;
                all_test_passed = 0;
        }

        // Integrate composition:
        SubsurfColumn column2 = SubsurfColumn();
        column2.addLayer(Layer(1, 0, 0, 1));
        column2.addLayer(Layer(2, 0, 1, 0));
        column2.addLayer(Layer(1, 1, 0, 0));
        Layer l = column2.integrateColumnComposition(10);

        if (fabs(l.thickness - 10) < eps) {
                std::cout << "Test passed: integrate column." << std::endl;
        }
        else {
                std::cout << "Test failed: integrate column." << std::endl;
                all_test_passed = 0;
        }

        ///////////////////////////
        // TEST GRID + SHADOWING //
        ///////////////////////////
        std::vector<var> varList = readConfig();
        std::vector< std::vector<double> > initLayersList = readLayers();

        Grid testGrid = Grid(initLayersList);
        Crater ctr = Crater(0, 0, 100);
        testGrid.formCrater(ctr);
        testGrid.printSurface(0);
        testGrid.printShadow(0);

        //Return tests outcome:
        return all_test_passed;
}
