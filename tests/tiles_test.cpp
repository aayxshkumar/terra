#include <gtest/gtest.h>

#include <cmath>

#include "terra/select.hpp"
#include "terra/tiles.hpp"

using namespace terra;

TEST(Tiles, LevelZeroIsTwoHemispheres) {
  EXPECT_EQ(columns(0), 2u);
  EXPECT_EQ(rows(0), 1u);
  const Bounds w = bounds({0, 0, 0}), e = bounds({0, 1, 0});
  EXPECT_DOUBLE_EQ(w.west, -kPi);
  EXPECT_DOUBLE_EQ(w.east, 0);
  EXPECT_DOUBLE_EQ(e.east, kPi);
  EXPECT_DOUBLE_EQ(w.north, kPi / 2);
  EXPECT_DOUBLE_EQ(w.south, -kPi / 2);
  EXPECT_NEAR(unit_area(w) + unit_area(e), 4 * kPi, 1e-14);
}

TEST(Tiles, KeyForIsRowMajorFromTheNorthWest) {
  EXPECT_EQ(key_for(kPi / 2 - 1e-9, -kPi + 1e-9, 3), (TileKey{3, 0, 0}));
  EXPECT_EQ(key_for(-kPi / 2, kPi, 3), (TileKey{3, 15, 7}));  // edges clamp inward
  EXPECT_EQ(key_for(0.1, 0.1, 1), (TileKey{1, 2, 0}));
}

TEST(Tiles, ParentOfChildRoundTrips) {
  const TileKey k{5, 37, 11};
  for (int q = 0; q < 4; ++q) EXPECT_EQ(parent(child(k, q)), k);
  EXPECT_EQ(child(k, 3), (TileKey{6, 75, 23}));
}

TEST(Tiles, SharedEdgesAreBitIdenticalAcrossLevels) {
  const TileKey k{7, 100, 40};
  EXPECT_EQ(bounds(k).east, bounds({7, 101, 40}).west);
  EXPECT_EQ(bounds(k).south, bounds({7, 100, 41}).north);
  EXPECT_EQ(bounds(child(k, 1)).east, bounds(k).east);
  EXPECT_EQ(bounds(child(child(k, 3), 3)).south, bounds(k).south);
}

TEST(Horizon, NadirVisibleAntipodeHidden) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const dvec3 cam = e.to_ecef({0.3, 0.2, 1000000});
  EXPECT_FALSE(behind_horizon(e, cam, e.to_ecef({0.3, 0.2, 0})));
  EXPECT_TRUE(behind_horizon(e, cam, e.to_ecef({-0.3, 0.2 + kPi, 0})));
}

TEST(Horizon, LowCameraKeepsTheTileBelowIt) {
  // From 10 km the visible cap is ~3 degrees wide, much smaller than a
  // level-2 tile: the tile underneath must still not be culled.
  const Ellipsoid e = Ellipsoid::wgs84();
  const dvec3 cam = e.to_ecef({0.5, 0.5, 10000});
  EXPECT_FALSE(horizon_culled(e, cam, key_for(0.5, 0.5, 2), 9000));
  EXPECT_TRUE(horizon_culled(e, cam, key_for(-0.5, 0.5 + kPi, 4), 9000));
}

TEST(Horizon, CameraBelowSurfaceCullsNothing) {
  const Ellipsoid e = Ellipsoid::wgs84();
  EXPECT_FALSE(behind_horizon(e, {0, 0, 0}, e.to_ecef({1, 1, 0})));
}

TEST(Frustum, SphereBehindCameraIsCulled) {
  Camera c;
  c.pos = {0, 0, 0};
  EXPECT_TRUE(frustum_culled(c, {{-100, 0, 0}, 10}));
  EXPECT_FALSE(frustum_culled(c, {{100, 0, 0}, 10}));
  EXPECT_FALSE(frustum_culled(c, {{-5, 0, 0}, 10}));  // camera inside the sphere
  EXPECT_TRUE(frustum_culled(c, {{10, 0, 500}, 10}));  // far above the top plane
}

TEST(Select, DistantCameraSeesAboutHalfTheGlobe) {
  const Ellipsoid e = Ellipsoid::wgs84();
  const Camera c = Camera::look_at({3e7, 0, 0}, {0, 0, 0}, {0, 0, 1}, 30, 400, 400);
  const Selection s = select_tiles(e, c, {});
  double area = 0;
  for (const TileKey& k : s.visible) area += unit_area(bounds(k));
  EXPECT_GT(area, 2 * kPi * 0.95);
  EXPECT_LT(area, 4 * kPi * 0.9);  // tiles straddling the limb are kept whole
  EXPECT_GT(s.horizon_culled, 0u);
}

TEST(Select, RefinesMoreWhenCloser) {
  const Ellipsoid e = Ellipsoid::wgs84();
  auto deepest = [&](double alt) {
    const Camera c = Camera::from_geodetic(e, {0.5, 0.5, alt}, 0, -90, 60, 320, 240);
    int d = 0;
    for (const TileKey& k : select_tiles(e, c, {}).visible) d = std::max(d, k.level);
    return d;
  };
  EXPECT_GT(deepest(10000), deepest(1000000));
}
