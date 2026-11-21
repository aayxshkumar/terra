#include "terra/terrain.hpp"

#include <algorithm>
#include <cmath>

#include "terra/ellipsoid.hpp"

namespace terra {

namespace {

std::uint32_t hash(std::int32_t x, std::int32_t y, std::int32_t z, std::uint32_t seed) {
  std::uint32_t h = seed * 0x9E3779B9u;
  h ^= static_cast<std::uint32_t>(x) * 0x85EBCA6Bu;
  h = (h << 13 | h >> 19) * 5u + 0xE6546B64u;
  h ^= static_cast<std::uint32_t>(y) * 0xC2B2AE35u;
  h = (h << 13 | h >> 19) * 5u + 0xE6546B64u;
  h ^= static_cast<std::uint32_t>(z) * 0x27D4EB2Fu;
  h ^= h >> 16;
  h *= 0x85EBCA6Bu;
  h ^= h >> 13;
  h *= 0xC2B2AE35u;
  h ^= h >> 16;
  return h;
}

double lattice(std::int32_t x, std::int32_t y, std::int32_t z, std::uint32_t seed) {
  return hash(x, y, z, seed) * (2.0 / 4294967295.0) - 1.0;
}

double fade(double t) { return t * t * t * (t * (t * 6 - 15) + 10); }

double lerp(double a, double b, double t) { return a + (b - a) * t; }

}  // namespace

double value_noise(const dvec3& p, std::uint32_t seed) {
  const double fx = std::floor(p.x), fy = std::floor(p.y), fz = std::floor(p.z);
  const auto ix = static_cast<std::int32_t>(fx), iy = static_cast<std::int32_t>(fy), iz = static_cast<std::int32_t>(fz);
  const double tx = fade(p.x - fx), ty = fade(p.y - fy), tz = fade(p.z - fz);
  double c[2][2];
  for (int dz = 0; dz < 2; ++dz)
    for (int dy = 0; dy < 2; ++dy)
      c[dz][dy] = lerp(lattice(ix, iy + dy, iz + dz, seed), lattice(ix + 1, iy + dy, iz + dz, seed), tx);
  return lerp(lerp(c[0][0], c[0][1], ty), lerp(c[1][0], c[1][1], ty), tz);
}

double planet_height(const PlanetParams& pp, double lat, double lon) {
  const dvec3 n = Ellipsoid::normal(lat, lon);
  // Continents: a few smooth octaves.
  double base = 0, amp = 1, freq = pp.continent_freq, norm = 0;
  for (int o = 0; o < 4; ++o) {
    base += amp * value_noise(n * freq + dvec3{17.1, -3.7, 9.2} * o, pp.seed + o);
    norm += amp;
    amp *= 0.5;
    freq *= 2.0;
  }
  base /= norm;
  // Mountains: ridged octaves normalised to [0, 1], raised only inside
  // low-frequency "range" patches and away from the coast.
  double ridge = 0, rnorm = 0;
  amp = 1;
  for (int o = 4; o < pp.octaves; ++o) {
    const double v = 1.0 - std::abs(value_noise(n * freq, pp.seed + o));
    ridge += amp * v * v;
    rnorm += amp;
    amp *= 0.5;
    freq *= 2.0;
  }
  ridge /= rnorm;
  const double land = base - pp.sea_level;
  double h = 0;
  if (land <= 0) {
    h = pp.ocean_depth * land;
  } else {
    const double range = std::clamp(2.5 * value_noise(n * 2.3 + dvec3{5.5, 1.3, -7.9}, pp.seed + 97), 0.0, 1.0);
    const double inland = std::clamp(land / 0.12, 0.0, 1.0);
    h = pp.lowland * land + pp.mountain * range * inland * std::pow(ridge, 6.0) + 120.0 * (ridge - 0.5) * inland;
  }
  return std::clamp(h, -kMaxTerrainHeight, kMaxTerrainHeight);
}

HeightFn planet(const PlanetParams& params) {
  return [params](double lat, double lon) { return planet_height(params, lat, lon); };
}

HeightFn flat() {
  return [](double, double) { return 0.0; };
}

}  // namespace terra
