// Procedural planet. Heights come from 3D noise evaluated on the unit
// sphere, so there are no seams at the antimeridian and no pinching at the
// poles. Integer hashing only: the field is identical on every platform, up
// to the libm sin/cos used to get the unit vector.
#pragma once

#include <cstdint>
#include <functional>

#include "terra/vec.hpp"

namespace terra {

// Height above the ellipsoid in metres (negative = sea floor) at a geodetic
// latitude / longitude in radians.
using HeightFn = std::function<double(double lat, double lon)>;

constexpr double kMaxTerrainHeight = 9000.0;  // every HeightFn here stays within +-this

struct PlanetParams {
  std::uint32_t seed = 7;
  int octaves = 14;           // finest wavelength ~ 6378 km / (1.6 * 2^13) ~ 0.5 km
  double continent_freq = 1.6;
  double sea_level = 0.02;    // noise value that maps to height 0
  double lowland = 2200.0;      // metres per unit of noise above sea level
  double mountain = 11000.0;    // peak scale inside ranges (ridge^6, so most is far lower)
  double ocean_depth = 12000.0; // metres per unit of noise below sea level
};

double value_noise(const dvec3& p, std::uint32_t seed);  // in [-1, 1]
double planet_height(const PlanetParams& params, double lat, double lon);
HeightFn planet(const PlanetParams& params = {});
HeightFn flat();

}  // namespace terra
