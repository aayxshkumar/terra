// Geographic (plate carree) tile pyramid, the same layout as osgEarth's
// "global-geodetic" profile: level 0 is two 180 x 180 degree tiles, and each
// level splits every tile into four. y = 0 is the northernmost row.
#pragma once

#include <compare>
#include <cstdint>

namespace terra {

struct TileKey {
  int level = 0;
  std::uint32_t x = 0, y = 0;
  auto operator<=>(const TileKey&) const = default;
};

struct Bounds {
  double west, south, east, north;  // radians
};

std::uint32_t columns(int level);  // 2 << level
std::uint32_t rows(int level);     // 1 << level
double tile_span(int level);       // radians, same in lon and lat

Bounds bounds(const TileKey& k);
TileKey child(const TileKey& k, int quadrant);  // 0 NW, 1 NE, 2 SW, 3 SE
TileKey parent(const TileKey& k);
TileKey key_for(double lat, double lon, int level);

// Area of the tile on the unit sphere (steradians); the whole pyramid level
// sums to 4 pi.
double unit_area(const Bounds& b);

}  // namespace terra
