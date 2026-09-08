// Impact-induced seismic shaking. See include/seismic.hpp.
#include <cmath>
#include <algorithm>
#include <functional>
#include "../include/regolit_main.hpp"
#include "../include/utility.hpp"
#include "../include/log.hpp"
#include "../include/seismic.hpp"

long SeismicShaking::shakes = 0;
long SeismicShaking::wholeDomainShakes = 0;
long SeismicShaking::distantShakes = 0;

// Richardson 2009, Eqs. 35-36: E_s = (eta/12) pi rho_i v^2 d^3 exp(-2 pi f l^2 / (K_s pi^2 Q)), K_s = v_s l_s / 3.
double SeismicShaking::energyAtDistance(double d, double v, double rho, double l) {
        const double seismicDiffusivity = seis_wave_vel * seis_mean_free / 3.0;
        const double injected = seismicEfficiency / 12.0 * M_PI * rho * v * v * d * d * d;
        return injected * exp(-2 * M_PI * prim_seis_freq * l * l / (seismicDiffusivity * M_PI * M_PI * Q_factor));
}

// Richardson 2009, Eqs. 33-34: energy density in a 1 m thick hemispherical shell, a_s = 2 pi f sqrt(2 eps / rho_t).
double SeismicShaking::accelerationAtDistance(double d, double v, double rho, double l) {
        const double energyDensity = energyAtDistance(d, v, rho, l) / (2 * M_PI * l * l * 1.0);
        return 2 * M_PI * prim_seis_freq * sqrt(2 * energyDensity / targetDensity);
}

// The acceleration decreases monotonically with distance: bisection in log(l).
double SeismicShaking::range(double d, double v, double rho) {
        const double threshold = seismicAccelerationThreshold * g;
        double low = log(0.5 * resolution), high = log(1e8);
        if (accelerationAtDistance(d, v, rho, exp(low)) < threshold)
                return 0;
        for (int iteration = 0; iteration < 60; iteration++) {
                const double middle = 0.5 * (low + high);
                if (accelerationAtDistance(d, v, rho, exp(middle)) >= threshold)
                        low = middle;
                else
                        high = middle;
        }
        return exp(low);
}

// Richardson 2009, Eq. 32.
double SeismicShaking::dose(double d, double v, double l, double craterRadius) {
        return Cs * pow(v, Ki_a) * pow(d, Ki_b) * pow(g, Ki_c) / pow(std::max(l, craterRadius), Ki_d);
}

void SeismicShaking::settleAfterImpact(Grid &grid, const Crater &crater) {
        const double slopeOfRepose = tan(M_PI * angleOfRepose / 180.0);
        double seismicRange = 0;
        if (isSeismicShaking && crater.projectileDiameter > 0)
                seismicRange = range(crater.projectileDiameter, crater.projectileVelocity, crater.projectileDensity);
        const double collapseRadius = 1.5 * crater.finalRadius + resolution;   // Richardson 2009: slopes within 1.5 R always checked
        const double halfWidth = std::max(seismicRange, collapseRadius);
        if (seismicRange >= resolution) {
                shakes++;
                if (seismicRange >= regionWidth / 2)
                        wholeDomainShakes++;
                const double d = crater.projectileDiameter, v = crater.projectileVelocity, R = crater.finalRadius, limit = seismicRange;
                std::function<double(double)> k = [d, v, R, limit](double l) { return l <= limit ? dose(d, v, l, R) : 0.0; };
                grid.settleRegion(crater.xLocation, crater.yLocation, halfWidth, &k, slopeOfRepose, true);
        } else {
                grid.settleRegion(crater.xLocation, crater.yLocation, halfWidth, nullptr, slopeOfRepose, true);
        }
}

void SeismicShaking::shakeFromOutside(Grid &grid, double x, double y, double craterRadius, double d, double v, double rho) {
        if (!isSeismicShaking)
                return;
        const double seismicRange = range(d, v, rho);
        const double halfDomain = regionWidth / 2;
        const double nearest = std::hypot(std::max(fabs(x) - halfDomain, 0.0), std::max(fabs(y) - halfDomain, 0.0));
        if (seismicRange <= nearest)
                return;
        distantShakes++;
        const double slopeOfRepose = tan(M_PI * angleOfRepose / 180.0);
        const double limit = seismicRange;
        std::function<double(double)> k = [d, v, craterRadius, limit](double l) { return l <= limit ? dose(d, v, l, craterRadius) : 0.0; };
        // The whole domain, with plain (not periodic) distances from the external impact point:
        grid.settleRegion(x, y, halfDomain, &k, slopeOfRepose, false);
}
