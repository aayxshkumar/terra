// terra: render the planet to PNG, and check the engine against geodesy
// references, closed-form geometry and brute force.
//
//   terra render [--lat D] [--lon D] [--alt M] [--heading D] [--pitch D] [--fov D]
//                [--size WxH] [--sse PX] [--grid N] [--flat] [--out file.png] [--stats]
//   terra --check
//   terra --bench
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "terra/camera.hpp"
#include "terra/ellipsoid.hpp"
#include "terra/mesh.hpp"
#include "terra/render.hpp"
#include "terra/select.hpp"
#include "terra/terrain.hpp"
#include "terra/tiles.hpp"

using namespace terra;

namespace {

double dms(double d, double m, double s) { return d + m / 60.0 + s / 3600.0; }

struct Row {
  std::string name;
  enum { Pass, Fail, Gap } result;
  double value, threshold;
  std::string measures;
};

Row judge(const std::string& name, double value, double threshold, const std::string& measures,
          bool known_gap = false) {
  const bool ok = value <= threshold;
  return {name, ok ? Row::Pass : (known_gap ? Row::Gap : Row::Fail), value, threshold, measures};
}

// Random camera somewhere around the planet: log-uniform altitude from
// 100 m to 20 000 km, any heading, pitch between straight down and level.
Camera random_camera(const Ellipsoid& e, std::mt19937_64& rng, int w = 640, int h = 360) {
  std::uniform_real_distribution<double> u(0, 1);
  const double lat = std::asin(2 * u(rng) - 1), lon = (2 * u(rng) - 1) * kPi;
  const double alt = 100.0 * std::pow(10.0, u(rng) * 5.3);
  return Camera::from_geodetic(e, {lat, lon, alt + kMaxTerrainHeight * (alt < 20000 ? 1 : 0)}, 360 * u(rng),
                               -90 + 90 * u(rng), 60, w, h);
}

std::vector<Row> run_checks() {
  std::vector<Row> rows;
  const Ellipsoid wgs = Ellipsoid::wgs84();
  std::mt19937_64 rng(20261121);
  std::uniform_real_distribution<double> u(0, 1);

  // 1. Geodetic -> ECEF, worked example from IOGP Guidance Note 7-2 (EPSG 9602).
  {
    const dvec3 p = wgs.to_ecef({deg(dms(53, 48, 33.82)), deg(dms(2, 7, 46.38)), 73.0});
    const dvec3 want{3771793.968, 140253.342, 5124304.349};
    rows.push_back(judge("ecef-iogp-example", length(p - want), 0.001,
                         "|ECEF - IOGP GN 7-2 worked example| (m)"));
  }
  // 2. ECEF -> geodetic -> ECEF round trip, poles and high orbits included.
  {
    double worst = 0;
    for (int k = 0; k < 200000; ++k) {
      double lat = std::asin(2 * u(rng) - 1);
      if (k % 1000 == 0) lat = (k % 2000 == 0 ? 1 : -1) * kPi / 2;
      const Geodetic g{lat, (2 * u(rng) - 1) * kPi, -10000.0 + u(rng) * u(rng) * 4.0e7};
      const dvec3 p = wgs.to_ecef(g);
      worst = std::max(worst, length(wgs.to_ecef(wgs.to_geodetic(p)) - p));
    }
    rows.push_back(judge("ecef-roundtrip", worst, 1e-6, "max position error over 200k points, h to 40 000 km (m)"));
  }
  // 3-4. Vincenty inverse on his own published example (Flinders Peak -> Buninyong).
  {
    const auto r = vincenty_inverse(wgs, deg(-dms(37, 57, 3.72030)), deg(dms(144, 25, 29.52440)),
                                    deg(-dms(37, 39, 10.15610)), deg(dms(143, 55, 35.38390)));
    rows.push_back(judge("vincenty-distance", std::abs(r.distance - 54972.271), 0.001,
                         "|s - 54 972.271 m| Flinders Peak to Buninyong (m)"));
    const double az1 = deg(dms(306, 52, 5.37)), back = deg(dms(127, 10, 25.07));
    const double err = std::max(std::abs(r.az1 - az1), std::abs(std::fmod(r.az2 + kPi, 2 * kPi) - back));
    rows.push_back(judge("vincenty-azimuths", to_deg(err) * 3600, 0.01, "max azimuth error vs published (arcsec)"));
  }
  // 5. Tile pyramid: children partition parents; key_for lands inside bounds.
  {
    int bad = 0;
    for (int k = 0; k < 20000; ++k) {
      const int level = static_cast<int>(rng() % 20);
      const double lat = std::asin(2 * u(rng) - 1), lon = (2 * u(rng) - 1) * kPi;
      const TileKey key = key_for(lat, lon, level);
      const Bounds b = bounds(key);
      if (lat < b.south || lat > b.north || lon < b.west || lon > b.east) ++bad;
      double sum = 0;
      for (int q = 0; q < 4; ++q) sum += unit_area(bounds(child(key, q)));
      if (std::abs(sum - unit_area(b)) > 1e-12 * unit_area(b) + 1e-300) ++bad;
    }
    rows.push_back(judge("tile-pyramid", bad, 0, "points outside their tile + children not summing to parent"));
  }
  // 6-7. LOD selection covers the globe exactly once, and is balanced.
  std::vector<Camera> cams;
  for (int k = 0; k < 40; ++k) cams.push_back(random_camera(wgs, rng));
  {
    double worst_area = 0;
    int worst_gap = 0, overlaps = 0;
    SelectOptions opt;
    for (const Camera& c : cams) {
      const Selection s = select_tiles(wgs, c, opt);
      double area = 0;
      for (const auto& [k, culled] : s.leaves) {
        area += unit_area(bounds(k));
        for (TileKey a = k; a.level > 0;) {
          a = parent(a);
          if (s.leaves.count(a)) ++overlaps;
        }
      }
      worst_area = std::max(worst_area, std::abs(area - 4 * kPi) / (4 * kPi));
      worst_gap = std::max(worst_gap, max_neighbour_level_gap(s, opt.max_level));
    }
    rows.push_back(judge("lod-coverage", worst_area + overlaps, 1e-12,
                         "|leaf area - 4 pi| / 4 pi (+ overlapping leaves), 40 random cameras"));
    rows.push_back(judge("lod-balanced", worst_gap, 1, "max level difference between edge neighbours"));
  }
  // 8. Horizon culling never removes a tile that has a visible point.
  {
    const HeightFn h = planet();
    std::size_t culled = 0, wrong = 0;
    for (std::size_t ci = 0; ci < 16; ++ci) {
      const Camera& c = cams[ci];
      const Selection s = select_tiles(wgs, c, SelectOptions{});
      for (const auto& [k, is_culled] : s.leaves) {
        if (!is_culled || k.level < 2 || !horizon_culled(wgs, c.pos, k, kMaxTerrainHeight)) continue;
        ++culled;
        const Bounds b = bounds(k);
        bool visible = false;
        for (int j = 0; j <= 32 && !visible; ++j)
          for (int i = 0; i <= 32 && !visible; ++i) {
            const double lat = lattice_lat(b, j, 32), lon = lattice_lon(b, i, 32);
            visible = !behind_horizon(wgs, c.pos, wgs.to_ecef({lat, lon, std::max(0.0, h(lat, lon))}));
          }
        wrong += visible ? 1 : 0;
      }
    }
    rows.push_back(judge("horizon-cull-safe", static_cast<double>(wrong), 0,
                         "horizon-culled tiles with a visible point (33x33 brute force), of " +
                             std::to_string(culled)));
  }
  // Low cameras looking straight down must see ground in every pixel: catches
  // over-eager culling that the sampled brute force above can miss.
  {
    RenderOptions ro;
    ro.sky = false;
    ro.shading = false;
    double worst = 1.0;
    for (double alt : {300.0, 10000.0, 150000.0})
      for (int k = 0; k < 4; ++k) {
        const double lat = std::asin(2 * u(rng) - 1), lon = (2 * u(rng) - 1) * kPi;
        const Camera c = Camera::from_geodetic(wgs, {lat, lon, alt + kMaxTerrainHeight}, 360 * u(rng), -90, 50, 160, 90);
        FrameStats st;
        render(wgs, planet(), c, ro, &st);
        worst = std::min(worst, st.pixels_covered / (160.0 * 90.0));
      }
    rows.push_back(judge("nadir-coverage", 1.0 - worst, 0, "1 - covered fraction, cameras at 300 m-150 km looking down"));
  }
  // 9-10. Cracks: same-level edges coincide; skirts are deeper than T-junction gaps.
  {
    const HeightFn h = planet();
    double same_level = 0, tj = 0;
    for (int t = 0; t < 300; ++t) {
      const int level = 3 + static_cast<int>(rng() % 12);
      const TileKey a{level, static_cast<std::uint32_t>(rng() % (columns(level) - 1)),
                      static_cast<std::uint32_t>(rng() % terra::rows(level))};
      const TileKey b{level, a.x + 1, a.y};
      const int g = 16;
      const TileMesh ma = build_mesh(wgs, h, a, g), mb = build_mesh(wgs, h, b, g);
      for (int j = 0; j <= g; ++j) {
        const dvec3 pa = ma.center + to_double(ma.pos[static_cast<std::size_t>(j * (g + 1) + g)]);
        const dvec3 pb = mb.center + to_double(mb.pos[static_cast<std::size_t>(j * (g + 1))]);
        same_level = std::max(same_level, length(pa - pb) / (tile_span(level) * wgs.a()));
      }
      // Coarse tile a against the two level+1 children of b along the shared edge.
      const Bounds ba = bounds(a);
      const double skirt = std::min(skirt_depth(wgs, a, g), skirt_depth(wgs, child(b, 0), g));
      auto at = [&](double lat) { return wgs.to_ecef({lat, ba.east, std::max(0.0, h(lat, ba.east))}); };
      for (int j = 0; j < g; ++j) {
        const double lat0 = lattice_lat(ba, j, g), lat1 = lattice_lat(ba, j + 1, g);
        const dvec3 p0 = at(lat0), p1 = at(lat1), mid = at(0.5 * (lat0 + lat1));
        const dvec3 d = p1 - p0;
        const double s = std::clamp(dot(mid - p0, d) / dot(d, d), 0.0, 1.0);
        tj = std::max(tj, length(mid - (p0 + d * s)) / skirt);
      }
    }
    rows.push_back(judge("edge-seams", same_level, 1e-6, "same-level shared edge gap / tile size (float offsets)"));
    rows.push_back(judge("skirt-covers-tjunction", tj, 1.0, "T-junction gap / skirt depth, coarse-fine pairs"));
  }
  // 11-12. Rendering a sphere: silhouette area and centre depth in closed form.
  {
    const double R = 6371000.0, d = 3 * R;
    const Ellipsoid sph = Ellipsoid::sphere(R);
    const Camera c = Camera::look_at({d, 0, 0}, {0, 0, 0}, {0, 0, 1}, 40, 960, 640);
    RenderOptions ro;
    ro.sky = false;
    ro.shading = false;
    ro.select.max_height = 0;
    FrameStats st;
    const Image img = render(sph, flat(), c, ro, &st);
    // Closed form, independent of Camera: half the image height spans tan(fov/2).
    const double r_px = 320.0 * std::tan(std::asin(R / d)) / std::tan(deg(20));
    rows.push_back(judge("sphere-silhouette", std::abs(st.pixels_covered / (kPi * r_px * r_px) - 1.0), 0.005,
                         "covered pixels vs pi r^2 of the projected sphere (relative)"));
    const float centre = img.inv_depth[static_cast<std::size_t>(320 * 960 + 480)];
    rows.push_back(judge("sphere-centre-depth", std::abs(1.0 / centre - (d - R)) / (d - R), 1e-4,
                         "depth at the image centre vs d - R (relative)"));
  }
  // 13. Known gap: plate carree tiles degenerate towards the poles.
  {
    const HeightFn h = flat();
    const int g = 16;
    const TileMesh m = build_mesh(wgs, h, {8, 0, 0}, g);
    double worst = 0;
    for (std::size_t t = 0; t < static_cast<std::size_t>(g * g * 6); t += 3) {
      const dvec3 p[3] = {to_double(m.pos[m.index[t]]), to_double(m.pos[m.index[t + 1]]),
                          to_double(m.pos[m.index[t + 2]])};
      const double e[3] = {length(p[1] - p[0]), length(p[2] - p[1]), length(p[0] - p[2])};
      const double lo = std::min({e[0], e[1], e[2]}), hi = std::max({e[0], e[1], e[2]});
      if (lo > 1e-6) worst = std::max(worst, hi / lo);
    }
    rows.push_back(judge("pole-slivers", worst, 20, "worst edge ratio in a polar level-8 tile (non-degenerate)", true));
  }
  // 14. Determinism: one scene, rendered twice.
  {
    const Camera c = Camera::from_geodetic(wgs, {deg(30), deg(40), 30000}, 20, -30, 60, 320, 180);
    const std::uint64_t h1 = image_hash(render(wgs, planet(), c, RenderOptions{}));
    const std::uint64_t h2 = image_hash(render(wgs, planet(), c, RenderOptions{}));
    rows.push_back(judge("determinism", h1 == h2 ? 0 : 1, 0, "same scene rendered twice: hashes differ (0 = identical)"));
  }
  return rows;
}

int check() {
  const auto rows = run_checks();
  std::printf("%-22s %-10s %12s %10s  %s\n", "check", "result", "value", "threshold", "measures");
  int failed = 0, gaps = 0;
  for (const Row& r : rows) {
    const char* res = r.result == Row::Pass ? "PASS" : r.result == Row::Fail ? "FAIL" : "KNOWN GAP";
    std::printf("%-22s %-10s %12.6g %10.3g  %s\n", r.name.c_str(), res, r.value, r.threshold, r.measures.c_str());
    failed += r.result == Row::Fail;
    gaps += r.result == Row::Gap;
  }
  std::printf("\n%d failed, %d known gap%s\n", failed, gaps, gaps == 1 ? "" : "s");
  return failed == 0 ? 0 : 1;
}

struct Args {
  double lat = 30.2, lon = -51.2, alt = 60000, heading = 35, pitch = -28, fov = 60, sse = 2.0;
  int w = 960, h = 540, grid = 16;
  bool flat_terrain = false, stats = false;
  std::string out = "terra.png";
};

int render_cmd(const Args& a) {
  const Ellipsoid wgs = Ellipsoid::wgs84();
  const Camera cam = Camera::from_geodetic(wgs, {deg(a.lat), deg(a.lon), a.alt}, a.heading, a.pitch, a.fov, a.w, a.h);
  RenderOptions ro;
  ro.select.max_sse_px = a.sse;
  ro.select.grid = a.grid;
  FrameStats st;
  const auto t0 = std::chrono::steady_clock::now();
  const Image img = render(wgs, a.flat_terrain ? flat() : planet(), cam, ro, &st);
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  if (!write_png(a.out, img)) {
    std::fprintf(stderr, "cannot write %s\n", a.out.c_str());
    return 1;
  }
  if (a.stats)
    std::printf("%s: %zu tiles, %zu triangles, %zu culled by frustum, %zu by horizon, %zu balance splits\n"
                "  select %.1f ms, meshes %.1f ms, raster %.1f ms, total %.1f ms\n",
                a.out.c_str(), st.tiles, st.triangles, st.frustum_culled, st.horizon_culled, st.balance_splits,
                st.select_ms, st.mesh_ms, st.raster_ms, ms);
  return 0;
}

int bench() {
  Args a;
  a.w = 1280;
  a.h = 720;
  a.stats = true;
  a.out = "bench.png";
  return render_cmd(a);
}

int usage() {
  std::fprintf(stderr,
               "usage: terra render [--lat D] [--lon D] [--alt M] [--heading D] [--pitch D] [--fov D]\n"
               "                    [--size WxH] [--sse PX] [--grid N] [--flat] [--out file.png] [--stats]\n"
               "       terra --check | --bench\n");
  return 2;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return usage();
  const std::string cmd = argv[1];
  if (cmd == "--check") return check();
  if (cmd == "--bench") return bench();
  if (cmd != "render") return usage();
  Args a;
  for (int k = 2; k < argc; ++k) {
    const std::string f = argv[k];
    auto num = [&](double& v) {
      if (k + 1 >= argc) return false;
      v = std::atof(argv[++k]);
      return true;
    };
    bool ok = true;
    if (f == "--lat") ok = num(a.lat);
    else if (f == "--lon") ok = num(a.lon);
    else if (f == "--alt") ok = num(a.alt);
    else if (f == "--heading") ok = num(a.heading);
    else if (f == "--pitch") ok = num(a.pitch);
    else if (f == "--fov") ok = num(a.fov);
    else if (f == "--sse") ok = num(a.sse);
    else if (f == "--grid" && k + 1 < argc) a.grid = std::atoi(argv[++k]);
    else if (f == "--size" && k + 1 < argc) ok = std::sscanf(argv[++k], "%dx%d", &a.w, &a.h) == 2;
    else if (f == "--out" && k + 1 < argc) a.out = argv[++k];
    else if (f == "--flat") a.flat_terrain = true;
    else if (f == "--stats") a.stats = true;
    else ok = false;
    if (!ok) return usage();
  }
  return render_cmd(a);
}
