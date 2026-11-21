#include "terra/tiles.hpp"

#include <algorithm>
#include <cmath>

#include "terra/vec.hpp"

namespace terra {

std::uint32_t columns(int level) { return 2u << level; }
std::uint32_t rows(int level) { return 1u << level; }
double tile_span(int level) { return kPi / static_cast<double>(1u << level); }

// Edges are computed from integer indices with one expression, so the east
// edge of a tile and the west edge of its neighbour are bit-identical, at
// every level (scaling by 2 commutes with rounding).
Bounds bounds(const TileKey& k) {
  const double s = tile_span(k.level);
  return {-kPi + k.x * s, kPi / 2 - (k.y + 1.0) * s, -kPi + (k.x + 1.0) * s, kPi / 2 - k.y * s};
}

TileKey child(const TileKey& k, int quadrant) {
  return {k.level + 1, 2 * k.x + static_cast<std::uint32_t>(quadrant & 1),
          2 * k.y + static_cast<std::uint32_t>(quadrant >> 1)};
}

TileKey parent(const TileKey& k) { return {k.level - 1, k.x / 2, k.y / 2}; }

TileKey key_for(double lat, double lon, int level) {
  const double s = tile_span(level);
  const auto clampi = [](double v, std::uint32_t n) {
    return static_cast<std::uint32_t>(std::clamp(std::floor(v), 0.0, static_cast<double>(n - 1)));
  };
  return {level, clampi((lon + kPi) / s, columns(level)), clampi((kPi / 2 - lat) / s, rows(level))};
}

double unit_area(const Bounds& b) { return (b.east - b.west) * (std::sin(b.north) - std::sin(b.south)); }

}  // namespace terra
