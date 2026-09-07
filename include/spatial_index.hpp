#pragma once
// A uniform bucket grid ("spatial hash") for point locations in a square domain centered on the
// origin. It replaces the Boost.Geometry R-tree: crater centers are roughly uniformly distributed
// over a fixed domain and the only query is "which craters lie within a distance of this point",
// which a bucket grid answers by visiting the few buckets that overlap the query square.
// Insert and remove are O(1); points outside the domain (ghost craters) fall into the edge buckets.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

class BucketGrid {
public:
        BucketGrid() : halfWidth_(0), bucketSize_(1), nBuckets_(1), buckets_(1) {}

        BucketGrid(double halfWidth, int nBucketsPerSide)
                : halfWidth_(halfWidth),
                  bucketSize_(2.0 * halfWidth / std::max(1, nBucketsPerSide)),
                  nBuckets_(std::max(1, nBucketsPerSide)),
                  buckets_(static_cast<size_t>(nBuckets_) * nBuckets_) {}

        void insert(size_t id, double x, double y) {
                buckets_[bucketOf(x, y)].push_back(id);
        }

        void remove(size_t id, double x, double y) {
                std::vector<size_t> &bucket = buckets_[bucketOf(x, y)];
                bucket.erase(std::remove(bucket.begin(), bucket.end(), id), bucket.end());
        }

        // Ids of all points stored in buckets overlapping the square of half-size `radius` around
        // (x, y). This is a superset of the points within `radius`; callers check exact distances.
        std::vector<size_t> candidatesWithin(double x, double y, double radius) const {
                std::vector<size_t> out;
                const int i0 = axisIndex(x - radius), i1 = axisIndex(x + radius);
                const int j0 = axisIndex(y - radius), j1 = axisIndex(y + radius);
                for (int i = i0; i <= i1; ++i) {
                        for (int j = j0; j <= j1; ++j) {
                                const std::vector<size_t> &bucket = buckets_[static_cast<size_t>(i) * nBuckets_ + j];
                                out.insert(out.end(), bucket.begin(), bucket.end());
                        }
                }
                return out;
        }

private:
        double halfWidth_;
        double bucketSize_;
        int nBuckets_;
        std::vector< std::vector<size_t> > buckets_;

        int axisIndex(double v) const {
                int k = static_cast<int>(std::floor((v + halfWidth_) / bucketSize_));
                return std::min(std::max(k, 0), nBuckets_ - 1);
        }

        size_t bucketOf(double x, double y) const {
                return static_cast<size_t>(axisIndex(x)) * nBuckets_ + axisIndex(y);
        }
};
