#include "terra/ellipsoid.hpp"

#include <cmath>

namespace terra {

Ellipsoid::Ellipsoid(double semi_major, double flattening)
    : a_(semi_major), f_(flattening), b_(semi_major * (1.0 - flattening)), e2_(flattening * (2.0 - flattening)) {}

dvec3 Ellipsoid::to_ecef(const Geodetic& g) const {
  const double s = std::sin(g.lat), c = std::cos(g.lat);
  const double n = a_ / std::sqrt(1.0 - e2_ * s * s);
  return {(n + g.h) * c * std::cos(g.lon), (n + g.h) * c * std::sin(g.lon), (n * (1.0 - e2_) + g.h) * s};
}

Geodetic Ellipsoid::to_geodetic(const dvec3& q) const {
  Geodetic g;
  const double p = std::hypot(q.x, q.y);
  g.lon = std::atan2(q.y, q.x);
  // lat = atan2(z + e2 N sin(lat), p) is well conditioned at the poles too.
  double lat = std::atan2(q.z, p * (1.0 - e2_));
  for (int k = 0; k < 30; ++k) {
    const double s = std::sin(lat);
    const double n = a_ / std::sqrt(1.0 - e2_ * s * s);
    const double next = std::atan2(q.z + e2_ * n * s, p);
    const bool done = std::abs(next - lat) < 1e-15;
    lat = next;
    if (done) break;
  }
  const double s = std::sin(lat), c = std::cos(lat);
  const double n = a_ / std::sqrt(1.0 - e2_ * s * s);
  g.lat = lat;
  g.h = p * c + q.z * s - a_ * a_ / n;  // stable form, no division by cos(lat)
  return g;
}

dvec3 Ellipsoid::normal(double lat, double lon) {
  const double c = std::cos(lat);
  return {c * std::cos(lon), c * std::sin(lon), std::sin(lat)};
}

Enu Ellipsoid::enu(double lat, double lon) {
  const double sl = std::sin(lat), cl = std::cos(lat), so = std::sin(lon), co = std::cos(lon);
  return {{-so, co, 0.0}, {-sl * co, -sl * so, cl}, {cl * co, cl * so, sl}};
}

// T. Vincenty (1975), "Direct and inverse solutions of geodesics on the
// ellipsoid with application of nested equations", Survey Review 23(176).
GeodesicResult vincenty_inverse(const Ellipsoid& e, double lat1, double lon1, double lat2, double lon2) {
  const double a = e.a(), b = e.b(), f = e.f();
  const double L = lon2 - lon1;
  const double u1 = std::atan((1 - f) * std::tan(lat1)), u2 = std::atan((1 - f) * std::tan(lat2));
  const double su1 = std::sin(u1), cu1 = std::cos(u1), su2 = std::sin(u2), cu2 = std::cos(u2);

  GeodesicResult r;
  double lambda = L, sin_sigma = 0, cos_sigma = 1, sigma = 0, cos2_alpha = 1, cos_2sm = 0;
  for (r.iterations = 1; r.iterations <= 200; ++r.iterations) {
    const double sl = std::sin(lambda), cl = std::cos(lambda);
    sin_sigma = std::hypot(cu2 * sl, cu1 * su2 - su1 * cu2 * cl);
    if (sin_sigma == 0) {  // coincident points
      r.converged = true;
      return r;
    }
    cos_sigma = su1 * su2 + cu1 * cu2 * cl;
    sigma = std::atan2(sin_sigma, cos_sigma);
    const double sin_alpha = cu1 * cu2 * sl / sin_sigma;
    cos2_alpha = 1 - sin_alpha * sin_alpha;
    cos_2sm = cos2_alpha != 0 ? cos_sigma - 2 * su1 * su2 / cos2_alpha : 0;  // equatorial line
    const double C = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha));
    const double prev = lambda;
    lambda = L + (1 - C) * f * sin_alpha *
                     (sigma + C * sin_sigma * (cos_2sm + C * cos_sigma * (-1 + 2 * cos_2sm * cos_2sm)));
    if (std::abs(lambda - prev) < 1e-12) {
      r.converged = true;
      break;
    }
  }
  const double u_sq = cos2_alpha * (a * a - b * b) / (b * b);
  const double A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));
  const double B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)));
  const double d_sigma =
      B * sin_sigma *
      (cos_2sm + B / 4 *
                     (cos_sigma * (-1 + 2 * cos_2sm * cos_2sm) -
                      B / 6 * cos_2sm * (-3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos_2sm * cos_2sm)));
  r.distance = b * A * (sigma - d_sigma);
  const double sl = std::sin(lambda), cl = std::cos(lambda);
  r.az1 = std::atan2(cu2 * sl, cu1 * su2 - su1 * cu2 * cl);
  r.az2 = std::atan2(cu1 * sl, -su1 * cu2 + cu1 * su2 * cl);
  const double two_pi = 2 * kPi;
  r.az1 = std::fmod(r.az1 + two_pi, two_pi);
  r.az2 = std::fmod(r.az2 + two_pi, two_pi);
  return r;
}

}  // namespace terra
