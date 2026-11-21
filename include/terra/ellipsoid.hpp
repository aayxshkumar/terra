// Reference ellipsoids, geodetic <-> ECEF conversion, local frames, and the
// Vincenty inverse geodesic. Angles are radians, lengths metres.
#pragma once

#include "terra/vec.hpp"

namespace terra {

struct Geodetic {
  double lat = 0, lon = 0, h = 0;
};

struct Enu {
  dvec3 east, north, up;
};

class Ellipsoid {
 public:
  Ellipsoid(double semi_major, double flattening);

  static Ellipsoid wgs84() { return {6378137.0, 1.0 / 298.257223563}; }
  static Ellipsoid sphere(double radius) { return {radius, 0.0}; }

  double a() const { return a_; }
  double b() const { return b_; }
  double f() const { return f_; }
  double e2() const { return e2_; }

  dvec3 to_ecef(const Geodetic& g) const;
  // Iterative (Bowring-style); converges to ~1e-15 rad for any point not
  // within a few km of the centre.
  Geodetic to_geodetic(const dvec3& p) const;

  // Unit normal to the ellipsoid at a geodetic latitude / longitude.
  static dvec3 normal(double lat, double lon);
  static Enu enu(double lat, double lon);

  // Radii for scaled-space horizon tests: x and y by a, z by b.
  dvec3 radii() const { return {a_, a_, b_}; }

 private:
  double a_, f_, b_, e2_;
};

struct GeodesicResult {
  double distance = 0;  // metres along the ellipsoid
  double az1 = 0;       // forward azimuth at point 1, radians clockwise from north
  double az2 = 0;       // forward azimuth at point 2
  int iterations = 0;
  bool converged = false;  // false near antipodal points, where Vincenty fails
};

GeodesicResult vincenty_inverse(const Ellipsoid& e, double lat1, double lon1, double lat2, double lon2);

}  // namespace terra
