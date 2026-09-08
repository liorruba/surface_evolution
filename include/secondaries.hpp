#pragma once
// Secondary craters as ejecta fragments (see README, "Secondary craters").
//
// Every primary launches its ejecta in the Z-model annuli stored in Crater::ejectaShells. Material
// that lands slower than secondaryMinimumVelocity only builds the ejecta blanket (already emplaced
// by the grid). Faster material is treated as a population of fragments with a cumulative
// size-frequency distribution N(>L) = (L_max / L)^slope_secondaries per annulus (weighted by the
// annulus volume). The largest fragment is anchored where the secondary field begins, at the annulus
// landing about three primary radii out: there it makes a crater of secondaryLargestFraction times
// the primary diameter (Allen 1979). Faster annuli have smaller largest fragments,
// L_max(v) = L_anchor (v / v_anchor)^-secondaryVelocityExponent (spallation scaling); slower ones are
// capped at L_anchor. Each fragment forms a crater at its landing point with the same pi-scaling as
// the primaries (fragment mass at the target density, landing speed), but only where that crater
// is deeper than the primary's own ejecta blanket at that distance: nearer the rim the fragments are
// buried in the blanket, which is why secondary fields begin beyond the continuous ejecta.
//
// Fragments that leave the domain are dropped. With isEmplaceDistantSecondaries the primaries that
// form outside the domain (out to secondaryMaximumRange) are sampled too, and the fragments they
// send into the domain are formed, so the domain receives the background of distant secondaries.
#include <vector>
#include "crater.hpp"
#include "grid.hpp"

struct Fragment {
        double x;        // landing point, m
        double y;
        double radius;   // fragment radius, m
        double speed;    // landing speed, m/s
};

class SecondaryModel {
public:
        explicit SecondaryModel(double regionWidth);

        // Fragments of `primary` that land inside the domain. Only fragments landing within
        // [azimuthCenter - azimuthHalfWidth, azimuthCenter + azimuthHalfWidth] are drawn (the
        // expected counts are scaled accordingly), which keeps distant primaries cheap.
        void sampleFragments(const Crater &primary, std::vector<Fragment> &out, double azimuthCenter = 0, double azimuthHalfWidth = M_PI) const;

        // Form the secondaries of a primary (already formed on the grid). Returns the number formed.
        long formSecondaries(const Crater &primary, Grid &grid);

        // Distant primaries: rings of increasing distance around the domain, each with the smallest
        // primary that can still deliver a resolvable secondary and its formation rate.
        void initializeDistantPrimaries(long totalNumberOfImpactors, double endTime);
        // Called once per domain impactor; forms the secondaries of the distant primaries due by then.
        long processDistantPrimaries(long impactorIndex, Grid &grid);

        long secondariesFormed = 0;          // from primaries inside the domain
        long distantPrimariesSampled = 0;
        long distantSecondariesFormed = 0;
        long primariesWithSecondaries = 0;

private:
        struct Ring {
                double innerHalfWidth;        // square annulus [inner, outer] half-widths around the domain center, m
                double outerHalfWidth;
                double smallestImpactorRadius; // impactors below this cannot send a resolvable fragment into the domain
                double meanGap;               // mean spacing of ring primaries, in domain-impactor units
                double nextEvent;             // impactor index at which the next ring primary forms
                double expectedCount;
        };
        double halfWidth;
        std::vector<Ring> rings;
        std::vector<double> speedTable, smallestFragmentTable;   // smallest crater-forming fragment vs landing speed

        void buildFragmentTable();
        double smallestFragmentRadius(double speed) const;
        double anchorSpeed(const Crater &primary) const;
        double largestFragmentRadius(const Crater &primary, double anchor, double speed) const;
        double smallestImpactorRadius(double landingSpeed) const;
};
