#pragma once
// Impact-induced seismic shaking (Richardson et al. 2005; Richardson 2009, Section 2.6).
//
// Every impact injects a fraction seismicEfficiency of its kinetic energy as seismic energy that
// spreads in a thin hemispherical shell and is attenuated by scattering (Richardson 2009, Eqs. 33-36).
// Where the peak seismic acceleration exceeds seismicAccelerationThreshold times g, the regolith
// receives one dose of downslope diffusion K_i = Cs v^a D^b g^c / l^d (Eq. 32, D the impactor
// diameter as fitted in Richardson et al. 2005, Fig. 14B) and slopes above the angle of repose
// collapse. The affected area is settled locally right after each impact; there is no global
// relaxation at output steps.
#include "grid.hpp"
#include "crater.hpp"

class SeismicShaking {
public:
        // Seismic energy (J) reaching distance l from an impactor of diameter d (m), speed v (m/s) and density rho.
        static double energyAtDistance(double d, double v, double rho, double l);
        // Peak seismic acceleration (m/s^2) at distance l.
        static double accelerationAtDistance(double d, double v, double rho, double l);
        // Distance (m) out to which the acceleration exceeds the threshold; 0 if it never does.
        static double range(double d, double v, double rho);
        // Downslope diffusion dose K_i (m^2) at distance l from the impact (l floored at the crater radius).
        static double dose(double d, double v, double l, double craterRadius);
        static double dosePrefactor(double d, double v);   // Cs v^a D^b g^c, evaluated once per impact

        // Settle the surface after a crater formed inside the domain: seismic diffusion within its
        // range (periodic distances) and slope collapse, at least within 1.5 crater radii.
        static void settleAfterImpact(Grid &grid, const Crater &crater);
        // Shaking by a primary that formed outside the domain (distant primaries of the secondary model).
        static void shakeFromOutside(Grid &grid, double x, double y, double craterRadius, double d, double v, double rho);

        static long shakes;               // impacts whose seismic range exceeded a cell
        static long wholeDomainShakes;    // impacts that shook the entire domain
        static long distantShakes;        // distant primaries that shook the domain
};
