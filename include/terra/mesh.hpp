// Tile meshes: a (grid+1)^2 lattice of geodetic samples, stored as float
// offsets from a double-precision tile centre (relative-to-centre, so float
// jitter stays proportional to the tile size, not to the Earth's radius),
// plus a skirt hanging down from every edge to hide T-junction cracks
// between tiles of different levels.
#pragma once

#include <cstdint>
#include <vector>

#include "terra/ellipsoid.hpp"
#include "terra/terrain.hpp"
#include "terra/tiles.hpp"
#include "terra/vec.hpp"

namespace terra {

struct TileMesh {
  TileKey key;
  int grid = 0;
  dvec3 center;               // ECEF
  std::vector<fvec3> pos;     // offsets from center; lattice first, then skirt
  std::vector<fvec3> normal;  // unit, ECEF directions
  std::vector<fvec3> albedo;  // linear RGB in [0, 1]
  std::vector<std::uint32_t> index;
  std::size_t lattice_vertices = 0;
  double skirt = 0;  // metres
};

double skirt_depth(const Ellipsoid& e, const TileKey& k, int grid);

// Surface colour from height (negative = water), latitude and slope
// (cosine of the angle between the terrain and ellipsoid normals).
fvec3 surface_albedo(double height, double lat, double slope_cos);

TileMesh build_mesh(const Ellipsoid& e, const HeightFn& height, const TileKey& k, int grid);

// Geodetic sample positions of lattice vertex (i, j), j = 0 on the north edge.
double lattice_lat(const Bounds& b, int j, int grid);
double lattice_lon(const Bounds& b, int i, int grid);

}  // namespace terra
