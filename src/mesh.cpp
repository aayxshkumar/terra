#include "terra/mesh.hpp"

#include <algorithm>
#include <cmath>

namespace terra {

double lattice_lat(const Bounds& b, int j, int grid) {
  return std::clamp(b.north - (b.north - b.south) * j / grid, -kPi / 2, kPi / 2);
}

double lattice_lon(const Bounds& b, int i, int grid) { return b.west + (b.east - b.west) * i / grid; }

double skirt_depth(const Ellipsoid& e, const TileKey& k, int grid) {
  return 0.5 * tile_span(k.level) * e.a() / grid + 50.0;
}

fvec3 surface_albedo(double h, double lat, double slope_cos) {
  auto mix = [](fvec3 a, fvec3 b, double t) {
    const float u = static_cast<float>(std::clamp(t, 0.0, 1.0));
    return a * (1 - u) + b * u;
  };
  if (h <= 0) return mix({0.05f, 0.22f, 0.40f}, {0.01f, 0.06f, 0.20f}, -h / 4000.0);
  const double alat = std::abs(lat);
  const double snowline = 5200.0 - 4800.0 * std::pow(std::sin(alat), 2.0);
  const fvec3 sand{0.74f, 0.68f, 0.50f}, grass{0.18f, 0.36f, 0.13f}, dry{0.50f, 0.42f, 0.27f},
      rock{0.40f, 0.37f, 0.33f}, snow{0.93f, 0.94f, 0.96f};
  fvec3 c = h < 40 ? sand : mix(grass, dry, (h - 40) / 2500.0);
  c = mix(c, dry, 1.0 - std::abs(alat - deg(25)) / deg(12));  // subtropical dry band
  c = mix(c, rock, (0.93 - slope_cos) / 0.12);
  c = mix(c, snow, (h - snowline) / 300.0);
  return c;
}

TileMesh build_mesh(const Ellipsoid& e, const HeightFn& height, const TileKey& k, int grid) {
  const Bounds b = bounds(k);
  const int n = grid, ring = n + 3;  // one extra sample around the lattice for normals
  std::vector<dvec3> p(static_cast<std::size_t>(ring) * ring);
  std::vector<double> raw(static_cast<std::size_t>(n + 1) * (n + 1));
  for (int j = -1; j <= n + 1; ++j)
    for (int i = -1; i <= n + 1; ++i) {
      const double lat = lattice_lat(b, j, n), lon = lattice_lon(b, i, n);
      const double h = height(lat, lon);
      if (i >= 0 && j >= 0 && i <= n && j <= n) raw[static_cast<std::size_t>(j * (n + 1) + i)] = h;
      p[static_cast<std::size_t>((j + 1) * ring + (i + 1))] = e.to_ecef({lat, lon, std::max(h, 0.0)});
    }
  auto P = [&](int i, int j) { return p[static_cast<std::size_t>((j + 1) * ring + (i + 1))]; };

  TileMesh m;
  m.key = k;
  m.grid = n;
  m.center = e.to_ecef({0.5 * (b.north + b.south), 0.5 * (b.east + b.west), 0.0});
  m.skirt = skirt_depth(e, k, n);
  for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i) {
      const double lat = lattice_lat(b, j, n), lon = lattice_lon(b, i, n);
      const dvec3 up = Ellipsoid::normal(lat, lon);
      dvec3 nrm = cross(P(i + 1, j) - P(i - 1, j), P(i, j - 1) - P(i, j + 1));
      const double len = length(nrm);
      nrm = len > 1e-9 ? nrm / len : up;  // degenerate at the poles
      if (dot(nrm, up) < 0) nrm = -nrm;
      m.pos.push_back(to_float(P(i, j) - m.center));
      m.normal.push_back(to_float(nrm));
      m.albedo.push_back(surface_albedo(raw[static_cast<std::size_t>(j * (n + 1) + i)], lat, dot(nrm, up)));
    }
  m.lattice_vertices = m.pos.size();
  auto vid = [n](int i, int j) { return static_cast<std::uint32_t>(j * (n + 1) + i); };
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      m.index.insert(m.index.end(), {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)});
      m.index.insert(m.index.end(), {vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)});
    }

  // Skirt: walk the perimeter once, dropping a copy of each vertex.
  std::vector<std::uint32_t> perimeter;
  for (int i = 0; i < n; ++i) perimeter.push_back(vid(i, 0));
  for (int j = 0; j < n; ++j) perimeter.push_back(vid(n, j));
  for (int i = n; i > 0; --i) perimeter.push_back(vid(i, n));
  for (int j = n; j > 0; --j) perimeter.push_back(vid(0, j));
  const auto first_skirt = static_cast<std::uint32_t>(m.pos.size());
  for (std::uint32_t v : perimeter) {
    const int i = static_cast<int>(v) % (n + 1), j = static_cast<int>(v) / (n + 1);
    const dvec3 up = Ellipsoid::normal(lattice_lat(b, j, n), lattice_lon(b, i, n));
    m.pos.push_back(to_float(P(i, j) - up * m.skirt - m.center));
    m.normal.push_back(m.normal[v]);
    m.albedo.push_back(m.albedo[v]);
  }
  const auto count = static_cast<std::uint32_t>(perimeter.size());
  for (std::uint32_t s = 0; s < count; ++s) {
    const std::uint32_t a = perimeter[s], c = perimeter[(s + 1) % count];
    const std::uint32_t as = first_skirt + s, cs = first_skirt + (s + 1) % count;
    m.index.insert(m.index.end(), {a, c, cs, a, cs, as});
  }
  return m;
}

}  // namespace terra
