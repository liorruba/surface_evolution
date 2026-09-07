// This class defines the surface properties:
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <initializer_list>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"
#include "../include/histogram.hpp"
#include "../include/grid.hpp"

namespace {
// Elevation changes smaller than this (m) are not worth a layer.
const double kElevationTolerance = 1e-9;
}

// A constructor for the grid class. Creates a grid with (size * size) elements.
Grid::Grid(std::vector< std::vector<double> > _initLayersList, std::vector<int8_t> _pxIdxMat){
	gridSize = (int) std::lround(regionWidth / resolution);
	if (gridSize < 1) {
		addLogEntry("ERROR: regionWidth / resolution must be at least 1.", true);
		exit(EXIT_FAILURE);
	}
	if (std::fabs(gridSize * resolution - regionWidth) > 1e-9 * regionWidth) {
		addLogEntry("WARNING: regionWidth is not a multiple of resolution; the grid spans " + std::to_string(gridSize * resolution) + " m.", true);
	}

	// Cell centers. The domain is periodic (see the ghost craters in main), so the cells are spaced
	// exactly one resolution apart, from -regionWidth/2 + resolution/2 to regionWidth/2 - resolution/2.
	x.resize(gridSize);
	y.resize(gridSize);
	for (int i = 0; i < gridSize; ++i) {
		x[i] = -regionWidth / 2 + (i + 0.5) * resolution;
		y[i] = x[i];
	}

	area = pow(regionWidth, 2);
	pixelIndexMatrix = _pxIdxMat;
	initLayersList = _initLayersList;
	subsurfColumns = initializeSubsurface();

	// Spatial index of crater centers, about 8 x 8 cells per bucket:
	craterIndex = BucketGrid(regionWidth / 2, std::max(1, gridSize / 8));
}

// Initialize a matrix of columns (2d vector):
std::vector< std::vector<SubsurfColumn> > Grid::initializeSubsurface(){
	addLogEntry("Generating grid (this may take some time, depending on the grid size)", true);

	if (initLayersList.empty()) {
		addLogEntry("ERROR: config/layers.cfg does not define any layer.", true);
		exit(EXIT_FAILURE);
	}

	// Build one template column per block of layers.cfg. A block is a run of consecutive lines with
	// the same class index in the first column; the k-th block (k = 0, 1, ...) is used for the cells
	// whose entry in the pixel index matrix is k.
	addLogEntry("Creating initial subsurface layering for each column class.", true);
	std::vector<SubsurfColumn> templates;
	SubsurfColumn buffCol = SubsurfColumn();
	double prev = initLayersList.front()[0];
	for (const std::vector<double> &row : initLayersList) {
		if (row[0] != prev){
			prev = row[0];
			templates.push_back(buffCol);
			buffCol = SubsurfColumn();
		}
		buffCol.addLayer(Layer(row[1], row[2], row[3], row[4]));
	}
	templates.push_back(buffCol);

	// Without a pixel index file the grid is uniform. The second block is used when it exists (this
	// is how the model has always behaved: the first block is then only a reference), otherwise the
	// first and only block.
	const size_t numberOfCells = (size_t) gridSize * gridSize;
	if (pixelIndexMatrix.empty()) {
		const int defaultClass = templates.size() > 1 ? 1 : 0;
		addLogEntry("No pixel index file: using layer class " + std::to_string(defaultClass) + " for the whole grid.", true);
		pixelIndexMatrix.assign(numberOfCells, (int8_t) defaultClass);
	}
	else if (pixelIndexMatrix.size() != numberOfCells) {
		addLogEntry("ERROR: the pixel index matrix has " + std::to_string(pixelIndexMatrix.size()) + " entries but the grid has " + std::to_string(numberOfCells) + " cells.", true);
		exit(EXIT_FAILURE);
	}

	addLogEntry("Populating grid by duplicating columns...", true);
	std::vector< std::vector<SubsurfColumn> > buffMat;
	buffMat.reserve(gridSize);
	for (int j = 0; j < gridSize; ++j) {        // rows (y)
		std::vector<SubsurfColumn> buffVec;
		buffVec.reserve(gridSize);
		for (int i = 0; i < gridSize; ++i) {    // columns (x)
			const int cls = pixelIndexMatrix[getLinearIndex(j, i, gridSize)];
			if (cls < 0 || cls >= (int) templates.size()) {
				addLogEntry("ERROR: pixel class " + std::to_string(cls) + " is not defined in config/layers.cfg (" + std::to_string(templates.size()) + " classes defined).", true);
				exit(EXIT_FAILURE);
			}
			buffVec.push_back(templates[cls]);
		}
		buffMat.push_back(buffVec);
		progressBar(j, gridSize);
	}

	addLogEntry("Finished creating grid.", true);
	return buffMat;
}

// Index range [iInit, iFinal) x [jInit, jFinal) of the cells whose centers may lie within halfSize
// of (xc, yc), clipped to the grid.
void Grid::footprintIndexRange(double xc, double yc, double halfSize, int &iInit, int &iFinal, int &jInit, int &jFinal) const {
	iInit = std::max(0, (int) std::floor((xc - halfSize + regionWidth / 2) / resolution));
	iFinal = std::min(gridSize, (int) std::ceil((xc + halfSize + regionWidth / 2) / resolution) + 1);
	jInit = std::max(0, (int) std::floor((yc - halfSize + regionWidth / 2) / resolution));
	jFinal = std::min(gridSize, (int) std::ceil((yc + halfSize + regionWidth / 2) / resolution) + 1);
}

// Get the surface elevation of the cell containing point (x, y); points outside the grid are
// clamped to the nearest edge cell (needed for ghost craters).
double Grid::getSurfaceElevationAtPoint(double pt_x, double pt_y) const {
	int i = (int) std::floor((pt_x + regionWidth / 2) / resolution);
	int j = (int) std::floor((pt_y + regionWidth / 2) / resolution);
	i = std::min(std::max(i, 0), gridSize - 1);
	j = std::min(std::max(j, 0), gridSize - 1);

	return subsurfColumns[j][i].getSurfaceElevation();
}

// Least-squares plane through the pre-impact surface inside the rim. The crater is shaped
// relative to this plane, so a crater on a slope (or on the wall of an older crater) is tilted
// with the terrain instead of being cut as a horizontal bowl with cliffs at the rim. Falls back to
// the mean elevation when there are too few cells for a fit, and to the elevation of the nearest
// cell when the footprint is entirely off the grid.
ReferencePlane Grid::fitReferencePlane(const Crater &crater) {
	ReferencePlane plane = {crater.xLocation, crater.yLocation, 0.0, 0.0, 0.0};
	const double R = crater.finalRadius;

	int iInit, iFinal, jInit, jFinal;
	footprintIndexRange(crater.xLocation, crater.yLocation, R, iInit, iFinal, jInit, jFinal);

	// Sums for the fit, in coordinates relative to the crater center:
	double n = 0, Sx = 0, Sy = 0, Sz = 0, Sxx = 0, Sxy = 0, Syy = 0, Sxz = 0, Syz = 0;
	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {
			const double dx = x[i] - crater.xLocation;
			const double dy = y[j] - crater.yLocation;
			if (dx * dx + dy * dy > R * R) {
				continue;
			}
			const double z = subsurfColumns[j][i].getSurfaceElevation();
			n += 1;
			Sx += dx; Sy += dy; Sz += z;
			Sxx += dx * dx; Sxy += dx * dy; Syy += dy * dy;
			Sxz += dx * z; Syz += dy * z;
		}
	}

	if (n == 0) {
		plane.z0 = getSurfaceElevationAtPoint(crater.xLocation, crater.yLocation);
		return plane;
	}

	const double meanX = Sx / n, meanY = Sy / n, meanZ = Sz / n;
	plane.z0 = meanZ;
	if (n < 3) {
		return plane;
	}

	// Centered second moments; the fit is degenerate when the cells are collinear.
	const double Cxx = Sxx - n * meanX * meanX;
	const double Cyy = Syy - n * meanY * meanY;
	const double Cxy = Sxy - n * meanX * meanY;
	const double Cxz = Sxz - n * meanX * meanZ;
	const double Cyz = Syz - n * meanY * meanZ;
	const double det = Cxx * Cyy - Cxy * Cxy;
	if (det <= 1e-9 * (Cxx + Cyy) * (Cxx + Cyy)) {
		return plane;
	}

	plane.sx = (Cxz * Cyy - Cyz * Cxy) / det;
	plane.sy = (Cyz * Cxx - Cxz * Cxy) / det;
	// The fitted plane passes through the centroid; express it at the crater center:
	plane.z0 = meanZ - plane.sx * meanX - plane.sy * meanY;
	return plane;
}

// Carve the crater cavity. The target surface inside the rim is the crater profile hung from the
// reference plane, raised by the rim height and by the ejecta thickness at the rim so that the
// interior joins the rim and the ejecta blanket continuously. Two passes: the first computes the
// target of every cell and integrates the composition of everything that has to be removed; the
// second applies the changes, filling cells that lie below the target (e.g. an older cavity under
// the new floor) with the ejected material. The volatile/soot loss is applied once to the bulk
// ejecta, not per cell.
void Grid::carveCavity(Crater &crater){
	const double R = crater.finalRadius;
	const ReferencePlane plane = fitReferencePlane(crater);
	const double rimEjecta = isEmplaceEjecta ? linearInterp(crater.ejectaDistance, crater.ejectaThickness, R) : 0.0;

	int iInit, iFinal, jInit, jFinal;
	footprintIndexRange(crater.xLocation, crater.yLocation, R, iInit, iFinal, jInit, jFinal);

	// Pass 1: targets and excavated composition.
	struct CellChange {
		int i;
		int j;
		double delta;   // target minus current elevation; negative = excavate
	};
	std::vector<CellChange> changes;
	Layer excavated = Layer(0, 0, 0, 0);

	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {
			const double r = std::hypot(x[i] - crater.xLocation, y[j] - crater.yLocation);
			if (r > R) {
				continue;
			}
			const double target = plane.at(x[i], y[j]) - cavityDepthProfile(R, r) + crater.rimHeight + rimEjecta;
			const double delta = target - subsurfColumns[j][i].getSurfaceElevation();
			if (delta < -kElevationTolerance) {
				excavated.consolidate(subsurfColumns[j][i].integrateColumnComposition(-delta));
			}
			changes.push_back({i, j, delta});
		}
	}

	// Ice and soot lost from the ejecta by impact heating:
	if (!excavated.isEmpty()) {
		excavated.changeComposition(excavated.regolithFraction,
			excavated.iceFraction * ejectaVolatileRetention,
			excavated.sootFraction * ejectaSootRetention);
	}
	crater.ejectedMass.consolidate(excavated);

	// Pass 2: apply.
	for (const CellChange &change : changes) {
		SubsurfColumn &column = subsurfColumns[change.j][change.i];
		if (change.delta < -kElevationTolerance) {
			column.removeMaterial(-change.delta);
		}
		else if (change.delta > kElevationTolerance && !crater.ejectedMass.isEmpty()) {
			column.addLayer(Layer(change.delta,
				crater.ejectedMass.regolithFraction,
				crater.ejectedMass.iceFraction,
				crater.ejectedMass.sootFraction));
		}
	}
}

// Outside the rim the raised rim decays as a power law of the distance from the center.
void Grid::emplaceRimDropoff(const Crater &crater){
	const double R = crater.finalRadius;
	const double outerRadius = 0.5 * ejectaSpread * R;

	int iInit, iFinal, jInit, jFinal;
	footprintIndexRange(crater.xLocation, crater.yLocation, outerRadius, iInit, iFinal, jInit, jFinal);

	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {
			const double r = std::hypot(x[i] - crater.xLocation, y[j] - crater.yLocation);
			if (r <= R || r > outerRadius) {
				continue;
			}
			const double rimDropoff = crater.rimHeight * pow(r / R, -rimDropoffExponent);
			if (rimDropoff > kElevationTolerance) {
				subsurfColumns[j][i].addLayer(Layer(rimDropoff,
					crater.ejectedMass.regolithFraction,
					crater.ejectedMass.iceFraction,
					crater.ejectedMass.sootFraction));
			}
		}
	}
}

// Emplace the ejecta blanket (Z-model thickness table) between the rim and ejectaSpread radii.
void Grid::emplaceEjecta(const Crater &crater){
	if (crater.ejectaDistance.size() < 2) {
		return;
	}
	const double R = crater.finalRadius;
	const double outerRadius = ejectaSpread * R;

	int iInit, iFinal, jInit, jFinal;
	footprintIndexRange(crater.xLocation, crater.yLocation, outerRadius, iInit, iFinal, jInit, jFinal);

	for (int i = iInit; i < iFinal; ++i) {
		for (int j = jInit; j < jFinal; ++j) {
			const double r = std::hypot(x[i] - crater.xLocation, y[j] - crater.yLocation);
			if (r <= R || r > outerRadius) {
				continue;
			}
			const double thickness = linearInterp(crater.ejectaDistance, crater.ejectaThickness, r);
			if (thickness > kElevationTolerance) {
				subsurfColumns[j][i].addLayer(Layer(thickness,
					crater.ejectedMass.regolithFraction,
					crater.ejectedMass.iceFraction,
					crater.ejectedMass.sootFraction));
			}
		}
	}
}

// Record a crater and add it to the spatial index. Craters centered outside the grid are ghost
// images of a recorded crater and are not recorded again.
void Grid::registerCrater(const Crater &crater){
	if (std::fabs(crater.xLocation) > regionWidth / 2 || std::fabs(crater.yLocation) > regionWidth / 2) {
		return;
	}
	craters.push_back({crater.xLocation, crater.yLocation, crater.finalRadius,
		crater.finalDepth, crater.finalDepth_init, crater.floorElevation, true});
	craterIndex.insert(craters.size() - 1, crater.xLocation, crater.yLocation);
}

// Update the stored depth of the older craters affected by a new crater: their floor elevation is
// re-read from the topography and the change is applied to their depth. A crater whose depth has
// changed by more than 50% either way (filled, or obliterated by a deeper crater) is no longer
// visible and leaves the index.
void Grid::updateExistingCratersDepth(const Crater &crater) {
	// The rim dropoff and the ejecta are negligible beyond two radii.
	const double reach = 2 * crater.finalRadius;

	for (size_t id : craterIndex.candidatesWithin(crater.xLocation, crater.yLocation, reach)) {
		CraterRecord &record = craters[id];
		if (std::hypot(record.x - crater.xLocation, record.y - crater.yLocation) > reach) {
			continue;
		}

		const double newFloorElevation = getSurfaceElevationAtPoint(record.x, record.y);
		record.finalDepth -= (newFloorElevation - record.floorElevation);
		record.floorElevation = newFloorElevation;

		if (std::fabs(record.finalDepth - record.finalDepth_init) / record.finalDepth_init > 0.5) {
			record.isVisible = false;
			craterIndex.remove(id, record.x, record.y);
		}
	}
}

size_t Grid::numberOfVisibleCraters() const {
	size_t count = 0;
	for (const CraterRecord &record : craters) {
		if (record.isVisible) {
			count++;
		}
	}
	return count;
}

// Form a new crater:
void Grid::formCrater(Crater &crater){
	if (crater.finalRadius <= 0) {
		return;
	}

	carveCavity(crater);

	if (!crater.ejectedMass.isEmpty()) {
		emplaceRimDropoff(crater);
		if (isEmplaceEjecta) {
			emplaceEjecta(crater);
		}
	}

	// Set crater formation elevation:
	crater.floorElevation = getSurfaceElevationAtPoint(crater.xLocation, crater.yLocation);

	// Older craters first (the new one is not in the index yet), then record the new crater:
	updateExistingCratersDepth(crater);
	registerCrater(crater);
}

///////////////////
// Crater profiles: depth below the reference plane, from depthToDiameter * D at the center to 0 at the rim.
///////////////////
double Grid::cavityDepthProfile(double craterRadius, double distanceFromCraterCenter) const {
	if (craterProfileType == 1) {
		return craterParabolicDepthProfile(craterRadius, distanceFromCraterCenter);
	}
	return craterSphericalDepthProfile(craterRadius, distanceFromCraterCenter);
}

// Parabolic:
double Grid::craterParabolicDepthProfile(double craterRadius, double distanceFromCraterCenter) const {
	return depthToDiameter * 2 * craterRadius * (1 - pow(distanceFromCraterCenter/craterRadius,2));
}

// Bowl shaped (spherical cap):
double Grid::craterSphericalDepthProfile(double craterRadius, double distanceFromCraterCenter) const {
	double craterDepth = 2 * craterRadius * depthToDiameter;
	double sphereRadius = (pow(craterRadius,2) + pow(craterDepth,2)) / 2 / craterDepth;

	return -sphereRadius + craterDepth + sqrt(std::max(0.0, pow(sphereRadius,2) - pow(distanceFromCraterCenter,2)));
}

///////////////////
// Simple sublimation/accumulation:
///////////////////
void Grid::sublimateIce() {
	return;
	// TODO: ADD NORBERT'S MODEL
}

// Deposit a layer on every cell whose surface already contains ice:
void Grid::depositLayer(Layer layer) {
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j) {
			if (subsurfColumns[i][j].subsurfLayers.back().iceFraction > 0){
				subsurfColumns[i][j].addLayer(layer);
			}
		}
	}
}

///////////////////
// Slopes above the angle of repose
///////////////////
// Surface elevation as a flat array indexed [j * gridSize + i]:
std::vector<double> Grid::surfaceElevationMap() const {
	std::vector<double> z((size_t) gridSize * gridSize);
	for (int j = 0; j < gridSize; ++j) {
		for (int i = 0; i < gridSize; ++i) {
			z[(size_t) j * gridSize + i] = subsurfColumns[j][i].getSurfaceElevation();
		}
	}
	return z;
}

// Move material downslope wherever the surface is steeper than maxSlope (rise over run) until no
// cell-to-cell slope exceeds it. Only interfaces steeper than the threshold transport material, so
// terrain below the angle of repose is left untouched. The transport is written in flux form (it
// conserves volume) on a periodic domain, consistent with the ghost craters, and the step size is
// within the explicit stability limit of the five-point stencil.
void Grid::relaxSlopes(std::vector<double> &z, double maxSlope) const {
	const size_t n = gridSize;
	const double maxRise = maxSlope * resolution;
	const double relaxation = 0.25;
	const int maxIterations = 100000;
	std::vector<double> dz(z.size());

	for (int iteration = 0; iteration < maxIterations; ++iteration) {
		std::fill(dz.begin(), dz.end(), 0.0);
		bool anySteep = false;

		for (size_t j = 0; j < n; ++j) {
			for (size_t i = 0; i < n; ++i) {
				const size_t k = j * n + i;
				const size_t kRight = j * n + (i + 1) % n;
				const size_t kDown = ((j + 1) % n) * n + i;

				for (size_t kn : {kRight, kDown}) {
					const double diff = z[k] - z[kn];
					const double excess = std::fabs(diff) - maxRise;
					if (excess > kElevationTolerance) {
						anySteep = true;
						const double flux = relaxation * excess * (diff > 0 ? 1.0 : -1.0);   // from k to kn
						dz[k] -= flux;
						dz[kn] += flux;
					}
				}
			}
		}

		if (!anySteep) {
			return;
		}
		for (size_t k = 0; k < z.size(); ++k) {
			z[k] += dz[k];
		}
	}

	addLogEntry("WARNING: slope relaxation did not converge within " + std::to_string(maxIterations) + " iterations.", true);
}

// Surface modification below the angle of repose (in deg). The elevation changes are applied to
// the columns; material that arrives in a cell takes the composition of that cell's surface layer.
void Grid::thresholdSlopes(double angleOfRepose) {
	const double slopeOfRepose = tan(M_PI * angleOfRepose / 180.0);

	std::vector<double> z = surfaceElevationMap();
	relaxSlopes(z, slopeOfRepose);

	for (int j = 0; j < gridSize; ++j) {
		for (int i = 0; i < gridSize; ++i) {
			SubsurfColumn &column = subsurfColumns[j][i];
			const double topoDiff = z[(size_t) j * gridSize + i] - column.getSurfaceElevation();

			if (topoDiff < -kElevationTolerance) {
				column.removeMaterial(-topoDiff);
			}
			else if (topoDiff > kElevationTolerance) {
				const Layer top = column.subsurfLayers.back();
				column.addLayer(Layer(topoDiff, top.regolithFraction, top.iceFraction, top.sootFraction));
			}
		}
	}
}

/////////
// Print:
/////////
// Write a matrix of doubles row by row in binary:
void Grid::writeMatrix(const std::string &fileName, const std::vector< std::vector<double> > &matrix) {
	std::ofstream file(fileName, std::ios_base::binary);
	for (const std::vector<double> &row : matrix) {
		file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
	}
}

// Print surface elevation and surface composition to file:
void Grid::printSurface(int index, bool isfinal) {
	addLogEntry("Writing surface elevation and composition to file for time step: " + std::to_string(index) + ".", false);
	const std::string index_str = formatOutputIndex(index);

	// Put data in matrices (indexed [x][y], as in the output files):
	std::vector< std::vector<double> > elevationMatrix(gridSize, std::vector<double>(gridSize));
	std::vector< std::vector<double> > regFractionMatrix(gridSize, std::vector<double>(gridSize));
	std::vector< std::vector<double> > iceFractionMatrix(gridSize, std::vector<double>(gridSize));
	std::vector< std::vector<double> > sootFractionMatrix(gridSize, std::vector<double>(gridSize));

	for (int i = 0; i < gridSize; ++i) {
		for (int j = 0; j < gridSize; ++j) {
			const SubsurfColumn &column = subsurfColumns[j][i];
			const Layer &top = column.subsurfLayers.back();

			elevationMatrix[i][j] = column.getSurfaceElevation();
			regFractionMatrix[i][j] = top.regolithFraction;
			iceFractionMatrix[i][j] = top.iceFraction;
			sootFractionMatrix[i][j] = top.sootFraction;
		}
	}

	// If the requested output resolution is coarser than the grid, re-bin the data:
	if (downsamplingResolution > resolution) {
		addLogEntry("Downsampling the surface to " + std::to_string(downsamplingResolution) + " m/pixel...", false);
		elevationMatrix = bin_2d_vector(elevationMatrix, downsamplingResolution);
		regFractionMatrix = bin_2d_vector(regFractionMatrix, downsamplingResolution);
		iceFractionMatrix = bin_2d_vector(iceFractionMatrix, downsamplingResolution);
		sootFractionMatrix = bin_2d_vector(sootFractionMatrix, downsamplingResolution);
	}

	writeMatrix("./output/elevation_" + index_str + ".out", elevationMatrix);
	writeMatrix("./output/regolithFraction_" + index_str + ".out", regFractionMatrix);
	writeMatrix("./output/iceFraction_" + index_str + ".out", iceFractionMatrix);
	writeMatrix("./output/sootFraction_" + index_str + ".out", sootFractionMatrix);

	// If it is the final time step, also write the horizontal coordinates:
	if (isfinal) {
		addLogEntry("Writing x,y coordinates to file...", false);
		std::vector<double> x_out = x;
		std::vector<double> y_out = y;
		if (downsamplingResolution > resolution) {
			x_out = bin_1d_vector(x, downsamplingResolution);
			y_out = bin_1d_vector(y, downsamplingResolution);
		}

		std::ofstream xFile("./output/x.out", std::ios_base::binary);
		std::ofstream yFile("./output/y.out", std::ios_base::binary);
		xFile.write(reinterpret_cast<const char*>(x_out.data()), x_out.size() * sizeof(double));
		yFile.write(reinterpret_cast<const char*>(y_out.data()), y_out.size() * sizeof(double));
	}
}

// Integrate the subsurface composition down to some depth and print it
void Grid::printIntegratedSubsurface(double depth, int index){
	const std::string index_str = formatOutputIndex(index);

	std::vector< std::vector<double> > regFractionMatrix(gridSize, std::vector<double>(gridSize));
	std::vector< std::vector<double> > iceFractionMatrix(gridSize, std::vector<double>(gridSize));
	std::vector< std::vector<double> > sootFractionMatrix(gridSize, std::vector<double>(gridSize));

	for (int i = 0; i < gridSize; ++i) {
		for (int j = 0; j < gridSize; ++j) {
			const Layer buffLayer = subsurfColumns[j][i].integrateColumnComposition(depth);

			regFractionMatrix[i][j] = buffLayer.regolithFraction;
			iceFractionMatrix[i][j] = buffLayer.iceFraction;
			sootFractionMatrix[i][j] = buffLayer.sootFraction;
		}
	}

	if (downsamplingResolution > resolution) {
		addLogEntry("Downsampling the integrated subsurface to " + std::to_string(downsamplingResolution) + " m/pixel...", false);
		regFractionMatrix = bin_2d_vector(regFractionMatrix, downsamplingResolution);
		iceFractionMatrix = bin_2d_vector(iceFractionMatrix, downsamplingResolution);
		sootFractionMatrix = bin_2d_vector(sootFractionMatrix, downsamplingResolution);
	}

	writeMatrix("./output/depthRegolithFraction_" + index_str + ".out", regFractionMatrix);
	writeMatrix("./output/depthIceFraction_" + index_str + ".out", iceFractionMatrix);
	writeMatrix("./output/depthSootFraction_" + index_str + ".out", sootFractionMatrix);
}

// Print subsurface to file: for every column, a header "layer" holding the number of layers and
// the surface elevation, followed by the layers bottom-up as raw Layer structs.
void Grid::printSubsurface(int index){
	const std::string index_str = formatOutputIndex(index);
	std::ofstream outputFile("./output/subsurface_" + index_str + ".out", std::ios_base::binary);

	for (int i = 0; i < gridSize; ++i) {
		for (int j = 0; j < gridSize; ++j) {
			const SubsurfColumn &col = subsurfColumns[j][i];
			Layer dummyLayer = Layer(col.subsurfLayers.size(), col.getSurfaceElevation());

			outputFile.write(reinterpret_cast<const char*>(&dummyLayer), sizeof(Layer));
			outputFile.write(reinterpret_cast<const char*>(col.subsurfLayers.data()), col.subsurfLayers.size() * sizeof(Layer));
		}
	}
}

// Print the visible craters to a histogram:
void Grid::printExistingCratersToHistogram(double bins){
	addLogEntry("Printing craters histogram.", false);
	Histogram hist(minimumImpactorDiameter * 10, regionWidth, (int) bins);

	for (const CraterRecord &record : craters) {
		if (record.isVisible) {
			hist.add(2 * record.finalRadius);
		}
	}

	hist.print("./output/existing_craters_histogram.txt");
}

// Print the visible craters to file:
void Grid::printExistingCraters(){
	addLogEntry("Printing the list of visible craters to file.", false);

	std::ofstream craterFile("./output/existing_craters.txt");
	craterFile << "x, y, diameter, depth, initial_depth\n";

	for (const CraterRecord &record : craters) {
		if (!record.isVisible) {
			continue;
		}
		craterFile << record.x << "," << record.y << ","
		<< 2 * record.finalRadius << ","
		<< record.finalDepth << ","
		<< record.finalDepth_init << "\n";
	}
}
