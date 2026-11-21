// Level-of-detail tile selection: refine by screen-space error, cull by the
// view frustum and by the ellipsoid horizon, then balance the quadtree so
// that neighbouring tiles differ by at most one level.
#pragma once

#include <cstddef>
#include <map>
#include <vector>

#include "terra/camera.hpp"
#include "terra/ellipsoid.hpp"
#include "terra/tiles.hpp"

namespace terra {

struct SelectOptions {
  double max_sse_px = 2.0;  // refine while a tile's geometric error projects to more than this
  int max_level = 18;
  int grid = 16;            // mesh cells per tile edge
  bool frustum_cull = true;
  bool horizon_cull = true;
  bool balance = true;
  double max_height = 9000.0;  // terrain bound used by culling and bounding spheres
};

struct Selection {
  // Every leaf of the quadtree; together they tile the whole globe exactly
  // once. `true` means the leaf is culled (not drawn).
  std::map<TileKey, bool> leaves;
  std::vector<TileKey> visible;
  std::size_t tested = 0, frustum_culled = 0, horizon_culled = 0, balance_splits = 0;
};

struct Sphere {
  dvec3 center;
  double radius = 0;
};

Sphere tile_sphere(const Ellipsoid& e, const TileKey& k, double max_height);
double geometric_error(const Ellipsoid& e, const TileKey& k, int grid);
double screen_space_error(const Ellipsoid& e, const Camera& cam, const TileKey& k, int grid, double max_height);

// Ellipsoid horizon test in scaled space (the occluder is the ellipsoid
// itself; terrain does not occlude). True if p is hidden from the camera.
bool behind_horizon(const Ellipsoid& e, const dvec3& camera, const dvec3& p);
bool sphere_behind_horizon(const Ellipsoid& e, const dvec3& camera, const Sphere& s);
bool horizon_culled(const Ellipsoid& e, const dvec3& camera, const TileKey& k, double max_height);
bool frustum_culled(const Camera& cam, const Sphere& s);

Selection select_tiles(const Ellipsoid& e, const Camera& cam, const SelectOptions& opt);

// The leaf containing a point, found by walking levels 0..max_level.
const TileKey* find_leaf(const Selection& s, double lat, double lon, int max_level);

// Largest level difference between edge-adjacent leaves (sampled along
// every edge); 1 for a balanced tree.
int max_neighbour_level_gap(const Selection& s, int max_level);

}  // namespace terra
