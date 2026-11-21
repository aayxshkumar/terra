#include "terra/select.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <limits>

namespace terra {

namespace {

// Calls f(lat, lon) for points just outside each edge of k, at fractions ts
// along the edge. Edges on a pole have no neighbour and are skipped.
template <class F>
void outside_points(const TileKey& k, std::initializer_list<double> ts, F&& f) {
  const Bounds b = bounds(k);
  const double eps = tile_span(k.level) * 1e-3;
  auto wrap = [](double lon) { return lon > kPi ? lon - 2 * kPi : (lon < -kPi ? lon + 2 * kPi : lon); };
  for (double t : ts) {
    const double lat = b.south + (b.north - b.south) * t;
    const double lon = b.west + (b.east - b.west) * t;
    f(lat, wrap(b.east + eps));
    f(lat, wrap(b.west - eps));
    if (b.north < kPi / 2) f(b.north + eps, lon);
    if (b.south > -kPi / 2) f(b.south - eps, lon);
  }
}

}  // namespace

Sphere tile_sphere(const Ellipsoid& e, const TileKey& k, double max_height) {
  const Bounds b = bounds(k);
  constexpr int n = 8;
  dvec3 pts[2 * (n + 1) * (n + 1)];
  int count = 0;
  dvec3 sum;
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i)
      for (double h : {0.0, max_height}) {
        const double lat = b.north - (b.north - b.south) * j / n;
        const double lon = b.west + (b.east - b.west) * i / n;
        pts[count] = e.to_ecef({lat, lon, h});
        sum += pts[count++];
      }
  Sphere s;
  s.center = sum / static_cast<double>(count);
  for (int k2 = 0; k2 < count; ++k2) s.radius = std::max(s.radius, length(pts[k2] - s.center));
  // Samples are 1/8 of a tile apart; pad for the surface bulging between them.
  s.radius *= 1.02;
  return s;
}

double geometric_error(const Ellipsoid& e, const TileKey& k, int grid) { return tile_span(k.level) * e.a() / grid; }

double screen_space_error(const Ellipsoid& e, const Camera& cam, const TileKey& k, int grid, double max_height) {
  const Sphere s = tile_sphere(e, k, max_height);
  const double d = length(s.center - cam.pos) - s.radius;
  if (d <= 0) return std::numeric_limits<double>::infinity();
  return geometric_error(e, k, grid) * cam.focal_px() / d;
}

bool behind_horizon(const Ellipsoid& e, const dvec3& camera, const dvec3& p) {
  const dvec3 r = e.radii();
  const dvec3 cv{camera.x / r.x, camera.y / r.y, camera.z / r.z};
  const double vh_mag_sq = dot(cv, cv) - 1.0;
  if (vh_mag_sq <= 0) return false;  // camera inside the ellipsoid: no horizon to use
  const dvec3 ps{p.x / r.x, p.y / r.y, p.z / r.z};
  const dvec3 vt = ps - cv;
  const double vt_dot_vc = -dot(vt, cv);
  return vt_dot_vc > vh_mag_sq && vt_dot_vc * vt_dot_vc / dot(vt, vt) > vh_mag_sq;
}

// Conservative: the whole bounding sphere must be hidden. In scaled space
// the ellipsoid is the unit sphere; everything inside the horizon cone and
// beyond the plane of the horizon circle is behind it. A world sphere of
// radius r fits in a scaled sphere of radius r / b.
bool sphere_behind_horizon(const Ellipsoid& e, const dvec3& camera, const Sphere& s) {
  const dvec3 r = e.radii();
  const dvec3 cv{camera.x / r.x, camera.y / r.y, camera.z / r.z};
  const double d = length(cv);
  if (d <= 1.0) return false;
  const dvec3 c_hat = cv / d;
  const dvec3 sc{s.center.x / r.x, s.center.y / r.y, s.center.z / r.z};
  const double rs = s.radius / e.b();
  if (dot(sc, c_hat) + rs >= 1.0 / d) return false;  // reaches the camera side of the horizon plane
  const dvec3 to = sc - cv;
  const double dist = length(to);
  if (dist <= rs) return false;
  const double angle = std::acos(std::clamp(dot(to / dist, -c_hat), -1.0, 1.0));
  return angle + std::asin(rs / dist) < std::asin(1.0 / d);
}

bool horizon_culled(const Ellipsoid& e, const dvec3& camera, const TileKey& k, double max_height) {
  return sphere_behind_horizon(e, camera, tile_sphere(e, k, max_height));
}

bool frustum_culled(const Camera& cam, const Sphere& s) {
  const dvec3 c = cam.to_camera(s.center);
  const double r = s.radius;
  if (c.z + r < cam.near_m) return true;
  const double vy = 0.5 * cam.fov_y;
  const double vx = std::atan(std::tan(vy) * cam.width / cam.height);
  const double cx = std::cos(vx), sx = std::sin(vx), cy = std::cos(vy), sy = std::sin(vy);
  return c.x * cx - c.z * sx > r || -c.x * cx - c.z * sx > r || c.y * cy - c.z * sy > r || -c.y * cy - c.z * sy > r;
}

const TileKey* find_leaf(const Selection& s, double lat, double lon, int max_level) {
  for (int l = 0; l <= max_level; ++l) {
    auto it = s.leaves.find(key_for(lat, lon, l));
    if (it != s.leaves.end()) return &it->first;
  }
  return nullptr;
}

Selection select_tiles(const Ellipsoid& e, const Camera& cam, const SelectOptions& opt) {
  Selection sel;
  // 0 = drawn, 1 = outside the frustum, 2 = behind the horizon
  auto cull = [&](const TileKey& k) {
    if (k.level < 2) return 0;  // hemisphere-sized tiles: bounds too loose to be worth it
    ++sel.tested;
    if (opt.frustum_cull && frustum_culled(cam, tile_sphere(e, k, opt.max_height))) return 1;
    if (opt.horizon_cull && horizon_culled(e, cam.pos, k, opt.max_height)) return 2;
    return 0;
  };

  std::vector<TileKey> stack = {{0, 0, 0}, {0, 1, 0}};
  while (!stack.empty()) {
    const TileKey k = stack.back();
    stack.pop_back();
    const int c = cull(k);
    if (c != 0) {
      sel.leaves[k] = true;
      (c == 1 ? sel.frustum_culled : sel.horizon_culled)++;
      continue;
    }
    if (k.level < opt.max_level && screen_space_error(e, cam, k, opt.grid, opt.max_height) > opt.max_sse_px) {
      for (int q = 0; q < 4; ++q) stack.push_back(child(k, q));
    } else {
      sel.leaves[k] = false;
    }
  }

  // Restricted quadtree: split any leaf that is two or more levels coarser
  // than an edge neighbour, until none is.
  for (bool changed = opt.balance; changed;) {
    changed = false;
    std::vector<TileKey> keys;
    keys.reserve(sel.leaves.size());
    for (const auto& [k, culled] : sel.leaves) keys.push_back(k);
    for (const TileKey& k : keys) {
      if (!sel.leaves.count(k)) continue;
      outside_points(k, {0.5}, [&](double lat, double lon) {
        const TileKey* n = find_leaf(sel, lat, lon, opt.max_level);
        if (!n || n->level >= k.level - 1) return;
        const TileKey coarse = *n;
        sel.leaves.erase(coarse);
        for (int q = 0; q < 4; ++q) {
          const TileKey ch = child(coarse, q);
          sel.leaves[ch] = cull(ch) != 0;
        }
        ++sel.balance_splits;
        changed = true;
      });
    }
  }

  for (const auto& [k, culled] : sel.leaves)
    if (!culled) sel.visible.push_back(k);
  return sel;
}

int max_neighbour_level_gap(const Selection& s, int max_level) {
  int worst = 0;
  for (const auto& [k, culled] : s.leaves) {
    outside_points(k, {0.25, 0.5, 0.75}, [&](double lat, double lon) {
      if (const TileKey* n = find_leaf(s, lat, lon, max_level)) worst = std::max(worst, std::abs(n->level - k.level));
    });
  }
  return worst;
}

}  // namespace terra
