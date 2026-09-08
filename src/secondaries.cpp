// Secondary craters as ejecta fragments. See include/secondaries.hpp for the model.
#include <cmath>
#include <algorithm>
#include <string>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"
#include "../include/impactor.hpp"
#include "../include/secondaries.hpp"
#include "../include/seismic.hpp"

SecondaryModel::SecondaryModel(double _regionWidth) : halfWidth(_regionWidth / 2), rings(), speedTable(), smallestFragmentTable() {
        buildFragmentTable();
}

// The smallest fragment that still makes a crater spanning two pixels depends only on the landing
// speed; it is tabulated once (log-log) instead of solving the inverse scaling per shell.
void SecondaryModel::buildFragmentTable() {
        const double smallestCraterDiameter = 2 * resolution;
        const double lowSpeed = std::max(1.0, secondaryMinimumVelocity), highSpeed = 1e5;
        const int n = 240;
        for (int i = 0; i < n; i++) {
                const double speed = lowSpeed * pow(highSpeed / lowSpeed, (double) i / (n - 1));
                speedTable.push_back(log(speed));
                smallestFragmentTable.push_back(log(Crater::impactorRadiusForCraterDiameter(smallestCraterDiameter, speed, targetDensity)));
        }
}

double SecondaryModel::smallestFragmentRadius(double speed) const {
        return exp(linearInterp(speedTable, smallestFragmentTable, log(std::max(speed, exp(speedTable.front())))));
}

// Speed of the ejecta landing where the secondary field begins, about three primary radii out
// (the annulus whose landing ring contains 3 R; the nearest annulus if none does).
double SecondaryModel::anchorSpeed(const Crater &primary) const {
        const double target = 3 * primary.finalRadius;
        double best = primary.ejectaShells.front().speed, bestDistance = INFINITY;
        for (const EjectaShell &shell : primary.ejectaShells) {
                if (shell.landingInner <= target && target <= shell.landingOuter)
                        return shell.speed;
                const double gap = std::min(fabs(shell.landingInner - target), fabs(shell.landingOuter - target));
                if (gap < bestDistance) {
                        bestDistance = gap;
                        best = shell.speed;
                }
        }
        return best;
}

// Largest fragment of the annulus ejected at `speed`: at the anchor speed it makes the largest
// secondary (secondaryLargestFraction of the primary diameter); faster fragments are smaller by
// (v / v_anchor)^-exponent (spallation scaling, Melosh 1984; Vickery 1986, 1987); slower ones are capped.
double SecondaryModel::largestFragmentRadius(const Crater &primary, double anchor, double speed) const {
        const double atAnchor = Crater::impactorRadiusForCraterDiameter(secondaryLargestFraction * 2 * primary.finalRadius, anchor, targetDensity);
        return speed <= anchor ? atAnchor : atAnchor * pow(speed / anchor, -secondaryVelocityExponent);
}

// Smallest impactor whose primary still sends a fragment that makes a two-pixel crater to a place
// where the ejecta lands at `landingSpeed` (bisection in log(radius); monotonic in the impactor size).
double SecondaryModel::smallestImpactorRadius(double landingSpeed) const {
        const double smallestCraterDiameter = 2 * resolution;
        double low = log(1e-3), high = log(1e5);
        for (int iteration = 0; iteration < 50; iteration++) {
                const double middle = 0.5 * (low + high);
                const Crater primary(Impactor(exp(middle)), 0, 0);
                double craterDiameter = 0;
                if (!primary.ejectaShells.empty()) {
                        const double fragment = largestFragmentRadius(primary, anchorSpeed(primary), landingSpeed);
                        craterDiameter = Crater::finalCraterDiameter(Impactor(fragment, landingSpeed, targetDensity));
                }
                if (craterDiameter < smallestCraterDiameter)
                        low = middle;
                else
                        high = middle;
        }
        return exp(high);
}

void SecondaryModel::sampleFragments(const Crater &primary, std::vector<Fragment> &out, double azimuthCenter, double azimuthHalfWidth) const {
        const double primaryDiameter = 2 * primary.finalRadius;
        const double largestSecondaryDiameter = secondaryLargestFraction * primaryDiameter;
        if (largestSecondaryDiameter < 2 * resolution || primary.ejectaShells.empty())
                return;
        const double vMin = secondaryMinimumVelocity;

        double fastVolume = 0;
        for (const EjectaShell &shell : primary.ejectaShells)
                if (shell.speed >= vMin)
                        fastVolume += shell.volume;
        if (fastVolume <= 0)
                return;
        const double anchor = anchorSpeed(primary);
        const double blanketEdge = ejectaSpread * primary.finalRadius;   // the grid emplaces the blanket out to here

        // Distances from the primary at which the domain square can be hit:
        const double px = primary.xLocation, py = primary.yLocation;
        const double nearest = std::hypot(std::max(fabs(px) - halfWidth, 0.0), std::max(fabs(py) - halfWidth, 0.0));
        const double farthest = std::hypot(fabs(px) + halfWidth, fabs(py) + halfWidth);

        struct Plan {
                const EjectaShell *shell;
                double smallest, largest, expected, overlapInner, overlapOuter;
        };
        std::vector<Plan> plans;
        double totalExpected = 0;   // all fragments of this primary above the two-pixel size, anywhere
        for (const EjectaShell &shell : primary.ejectaShells) {
                if (shell.speed < vMin)
                        continue;
                const double largest = largestFragmentRadius(primary, anchor, shell.speed);
                const double smallest = smallestFragmentRadius(shell.speed);
                if (largest <= smallest)
                        continue;
                const double count = shell.volume / fastVolume * pow(largest / smallest, slope_secondaries);
                totalExpected += count;
                const double overlapInner = std::max(shell.landingInner, nearest);
                const double overlapOuter = std::min(shell.landingOuter, farthest);
                if (overlapOuter <= overlapInner)
                        continue;
                const double ringArea = pow(shell.landingOuter, 2) - pow(shell.landingInner, 2);
                const double areaFraction = ringArea > 0 ? (pow(overlapOuter, 2) - pow(overlapInner, 2)) / ringArea : 1.0;
                plans.push_back(Plan{&shell, smallest, largest, count * areaFraction * (azimuthHalfWidth / M_PI), overlapInner, overlapOuter});
        }
        if (plans.empty())
                return;

        // Budget: only the largest maximumSecondariesPerPrimary fragments of a primary are formed. The
        // size floor is raised so that the expected population fits; smaller fragments are treated
        // as part of the ejecta blanket.
        double floorFactor = 1.0;
        if (totalExpected > maximumSecondariesPerPrimary)
                floorFactor = pow(totalExpected / maximumSecondariesPerPrimary, 1.0 / slope_secondaries);

        for (const Plan &plan : plans) {
                const double smallest = plan.smallest * floorFactor;
                if (smallest >= plan.largest)
                        continue;
                const long n = randPoisson(plan.expected * pow(floorFactor, -slope_secondaries));
                const double truncation = 1 - pow(smallest / plan.largest, slope_secondaries);
                double volumeBudget = plan.shell->volume;   // the fragments cannot carry more than their shell
                for (long k = 0; k < n; k++) {
                        const double radius = smallest * pow(1 - randU(0, 1) * truncation, -1 / slope_secondaries);
                        const double volume = 4.0 / 3.0 * M_PI * pow(radius, 3);
                        if (volume > volumeBudget)
                                break;
                        volumeBudget -= volume;
                        const double distance = sqrt(pow(plan.overlapInner, 2) + randU(0, 1) * (pow(plan.overlapOuter, 2) - pow(plan.overlapInner, 2)));
                        const double azimuth = azimuthCenter + randU(-azimuthHalfWidth, azimuthHalfWidth);
                        const double x = px + distance * cos(azimuth), y = py + distance * sin(azimuth);
                        if (fabs(x) > halfWidth || fabs(y) > halfWidth)
                                continue;
                        // Buried in the primary's blanket? (crater depth vs. blanket thickness at the landing distance)
                        const double blanket = distance < blanketEdge ? linearInterp(primary.ejectaDistance, primary.ejectaThickness, distance) : 0.0;
                        if (blanket > 0) {
                                const double craterDiameter = Crater::finalCraterDiameter(Impactor(radius, plan.shell->speed, targetDensity));
                                if (secondaryDepthToDiameter * craterDiameter <= blanket)
                                        continue;
                        }
                        out.push_back(Fragment{x, y, radius, plan.shell->speed});
                }
        }
}

long SecondaryModel::formSecondaries(const Crater &primary, Grid &grid) {
        std::vector<Fragment> fragments;
        sampleFragments(primary, fragments);
        for (const Fragment &fragment : fragments) {
                Crater secondary(Impactor(fragment.radius, fragment.speed, targetDensity), fragment.x, fragment.y, secondaryDepthToDiameter);
                grid.formCrater(secondary);
        }
        if (!fragments.empty())
                primariesWithSecondaries++;
        secondariesFormed += fragments.size();
        return fragments.size();
}

// Square annuli of doubling half-width around the domain, out to secondaryMaximumRange from its edge.
void SecondaryModel::initializeDistantPrimaries(long totalNumberOfImpactors, double endTime) {
        rings.clear();
        const double outermost = halfWidth + secondaryMaximumRange;
        double inner = halfWidth;
        while (inner < outermost) {
                const double outer = std::min(2 * inner, outermost);
                Ring ring;
                ring.innerHalfWidth = inner;
                ring.outerHalfWidth = outer;
                // A fragment must fly at least to the ring's nearest edge; the speed that takes it there
                // (45-degree ballistics, flat surface) bounds the fragment size and thus the primary size.
                const double nearestDistance = inner - halfWidth;
                const double landingSpeed = std::max(secondaryMinimumVelocity, sqrt(g * nearestDistance));
                ring.smallestImpactorRadius = smallestImpactorRadius(landingSpeed);
                const double smallestPrimary = Crater::finalCraterDiameter(Impactor(ring.smallestImpactorRadius));
                const double area = 4 * (outer * outer - inner * inner);
                ring.expectedCount = fluxConstant_c * pow(2 * ring.smallestImpactorRadius, -slope_b) * earthFluxRatioCoefficient * area * endTime;
                ring.meanGap = ring.expectedCount > 0 ? totalNumberOfImpactors / ring.expectedCount : INFINITY;
                ring.nextEvent = randExponential(ring.meanGap);
                rings.push_back(ring);
                addLogEntry("Distant primaries " + std::to_string((long) (inner - halfWidth)) + "-" + std::to_string((long) (outer - halfWidth)) +
                            " m beyond the edge: craters larger than " + std::to_string((long) smallestPrimary) + " m, about " +
                            std::to_string((long) ring.expectedCount) + " expected over the run.", true);
                inner = outer;
        }
}

long SecondaryModel::processDistantPrimaries(long impactorIndex, Grid &grid) {
        long formed = 0;
        std::vector<Fragment> fragments;
        for (Ring &ring : rings) {
                while (ring.nextEvent <= impactorIndex) {
                        ring.nextEvent += randExponential(ring.meanGap);
                        distantPrimariesSampled++;
                        // Uniform position in the square annulus (rejection from the outer square):
                        double x, y;
                        do {
                                x = randU(-ring.outerHalfWidth, ring.outerHalfWidth);
                                y = randU(-ring.outerHalfWidth, ring.outerHalfWidth);
                        } while (fabs(x) < ring.innerHalfWidth && fabs(y) < ring.innerHalfWidth);
                        const double impactorRadius = ring.smallestImpactorRadius * pow(randU(0, 1), -1 / slope_b);
                        Crater primary(Impactor(impactorRadius), x, y);
                        // Only the azimuths that face the domain (its circumscribed circle) are drawn:
                        const double distance = std::hypot(x, y);
                        const double circumscribed = halfWidth * M_SQRT2;
                        const double azimuthHalfWidth = distance <= circumscribed ? M_PI : asin(circumscribed / distance);
                        fragments.clear();
                        sampleFragments(primary, fragments, atan2(-y, -x), azimuthHalfWidth);
                        for (const Fragment &fragment : fragments) {
                                Crater secondary(Impactor(fragment.radius, fragment.speed, targetDensity), fragment.x, fragment.y, secondaryDepthToDiameter);
                                grid.formCrater(secondary);
                        }
                        formed += fragments.size();
                        SeismicShaking::shakeFromOutside(grid, x, y, primary.finalRadius, 2 * impactorRadius, meanImpactVelocity, impactorDensity);
                }
        }
        distantSecondariesFormed += formed;
        return formed;
}
