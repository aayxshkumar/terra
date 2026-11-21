// Geodesy against closed forms that don't go through the code under test.
#include <gtest/gtest.h>

#include <cmath>

#include "terra/ellipsoid.hpp"

using namespace terra;

TEST(Ellipsoid, Wgs84DerivedConstants) {
  const Ellipsoid e = Ellipsoid::wgs84();
  EXPECT_NEAR(e.b(), 6356752.314245, 1e-6);  // published semi-minor axis
  EXPECT_NEAR(e.e2(), 0.00669437999014, 1e-14);
}

TEST(Ellipsoid, AxisPointsAreExact) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const dvec3 eq = e.to_ecef({0, 0, 0});
  EXPECT_DOUBLE_EQ(eq.x, e.a());
  EXPECT_DOUBLE_EQ(eq.y, 0);
  const dvec3 np = e.to_ecef({kPi / 2, 1.0, 100});
  EXPECT_NEAR(np.x, 0, 1e-9);
  EXPECT_NEAR(np.z, e.b() + 100, 1e-6);
  const dvec3 east = e.to_ecef({0, kPi / 2, -50});
  EXPECT_NEAR(east.y, e.a() - 50, 1e-6);
}

TEST(Ellipsoid, InverseAtThePoleAndCentreLine) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const Geodetic g = e.to_geodetic({0, 0, e.b() + 1234.5});
  EXPECT_NEAR(g.lat, kPi / 2, 1e-12);
  EXPECT_NEAR(g.h, 1234.5, 1e-6);
  const Geodetic s = e.to_geodetic({0, 0, -(e.b() + 10)});
  EXPECT_NEAR(s.lat, -kPi / 2, 1e-12);
  EXPECT_NEAR(s.h, 10, 1e-6);
}

TEST(Ellipsoid, EnuIsOrthonormalAndRightHanded) {
  for (double lat : {-1.2, 0.0, 0.7, 1.5})
    for (double lon : {-3.0, 0.0, 2.0}) {
      const Enu f = Ellipsoid::enu(lat, lon);
      EXPECT_NEAR(length(f.east), 1, 1e-15);
      EXPECT_NEAR(dot(f.east, f.north), 0, 1e-15);
      EXPECT_NEAR(dot(f.north, f.up), 0, 1e-15);
      EXPECT_NEAR(length(cross(f.east, f.north) - f.up), 0, 1e-15);
    }
}

TEST(Ellipsoid, SphereGeodeticEqualsGeocentric) {
  const Ellipsoid s = Ellipsoid::sphere(1000);
  const dvec3 p = s.to_ecef({0.3, -1.1, 25});
  EXPECT_NEAR(length(p), 1025, 1e-9);
  EXPECT_NEAR(std::asin(p.z / length(p)), 0.3, 1e-15);
}

TEST(Vincenty, EquatorIsAnArcOfTheMajorCircle) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const auto r = vincenty_inverse(e, 0, 0, 0, kPi / 2);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(r.distance, e.a() * kPi / 2, 1e-5);  // lambda converges to 1e-12 rad ~ 6 um
  EXPECT_NEAR(r.az1, kPi / 2, 1e-12);
}

TEST(Vincenty, MeridianQuarterIsTheMeridianArc) {
  // Equator to pole along a meridian: published WGS84 quarter meridian.
  const auto r = vincenty_inverse(Ellipsoid::wgs84(), 0, 0, kPi / 2, 0);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(r.distance, 10001965.729, 1e-3);
}

TEST(Vincenty, SymmetricAndZeroForCoincidentPoints) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const auto ab = vincenty_inverse(e, 0.5, 0.1, -0.2, 1.3), ba = vincenty_inverse(e, -0.2, 1.3, 0.5, 0.1);
  EXPECT_NEAR(ab.distance, ba.distance, 1e-6);
  EXPECT_NEAR(vincenty_inverse(e, 0.4, 0.4, 0.4, 0.4).distance, 0, 1e-9);
}

TEST(Vincenty, ReportsFailureNearAntipodes) {
  const auto r = vincenty_inverse(Ellipsoid::wgs84(), 0, 0, 0.5e-3, kPi - 1e-3);
  EXPECT_FALSE(r.converged);
}
