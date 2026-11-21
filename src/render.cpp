#include "terra/render.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>

#include "terra/mesh.hpp"

namespace terra {

Image::Image(int w, int h)
    : width(w), height(h), rgb(static_cast<std::size_t>(w) * h * 3, 0), inv_depth(static_cast<std::size_t>(w) * h, 0.f) {}

namespace {

double ms_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

std::uint8_t to_byte(double v) {
  // Linear to sRGB-ish gamma 2.2.
  return static_cast<std::uint8_t>(std::lround(255.0 * std::pow(std::clamp(v, 0.0, 1.0), 1.0 / 2.2)));
}

// Top-left fill rule for an edge a -> b of a triangle wound so its signed
// area is positive in y-down screen space.
bool top_left(const ScreenVertex& a, const ScreenVertex& b) {
  const double dx = b.x - a.x, dy = b.y - a.y;
  return dy < 0 || (dy == 0 && dx > 0);
}

double edge(const ScreenVertex& a, const ScreenVertex& b, double px, double py) {
  return (b.x - a.x) * (py - a.y) - (b.y - a.y) * (px - a.x);
}

struct CamVertex {
  dvec3 p, c;
};

}  // namespace

std::size_t raster_triangle(Image& img, ScreenVertex v0, ScreenVertex v1, ScreenVertex v2) {
  double area = edge(v0, v1, v2.x, v2.y);
  if (area == 0 || !std::isfinite(area)) return 0;
  if (area < 0) {
    std::swap(v1, v2);
    area = -area;
  }
  const int x0 = std::max(0, static_cast<int>(std::floor(std::min({v0.x, v1.x, v2.x}))));
  const int x1 = std::min(img.width - 1, static_cast<int>(std::ceil(std::max({v0.x, v1.x, v2.x}))));
  const int y0 = std::max(0, static_cast<int>(std::floor(std::min({v0.y, v1.y, v2.y}))));
  const int y1 = std::min(img.height - 1, static_cast<int>(std::ceil(std::max({v0.y, v1.y, v2.y}))));
  const bool tl0 = top_left(v1, v2), tl1 = top_left(v2, v0), tl2 = top_left(v0, v1);
  std::size_t written = 0;
  for (int y = y0; y <= y1; ++y) {
    const double py = y + 0.5;
    for (int x = x0; x <= x1; ++x) {
      const double px = x + 0.5;
      const double w0 = edge(v1, v2, px, py), w1 = edge(v2, v0, px, py), w2 = edge(v0, v1, px, py);
      if (w0 < 0 || w1 < 0 || w2 < 0) continue;
      if ((w0 == 0 && !tl0) || (w1 == 0 && !tl1) || (w2 == 0 && !tl2)) continue;
      const double b0 = w0 / area, b1 = w1 / area, b2 = w2 / area;
      const double iz = b0 * v0.inv_z + b1 * v1.inv_z + b2 * v2.inv_z;
      const std::size_t idx = static_cast<std::size_t>(y) * img.width + x;
      if (iz <= img.inv_depth[idx]) continue;
      img.inv_depth[idx] = static_cast<float>(iz);
      // Perspective-correct colour: interpolate c/z, divide by 1/z.
      const dvec3 c = (v0.color * (b0 * v0.inv_z) + v1.color * (b1 * v1.inv_z) + v2.color * (b2 * v2.inv_z)) / iz;
      img.rgb[idx * 3 + 0] = to_byte(c.x);
      img.rgb[idx * 3 + 1] = to_byte(c.y);
      img.rgb[idx * 3 + 2] = to_byte(c.z);
      ++written;
    }
  }
  return written;
}

std::size_t draw_camera_triangle(Image& img, const Camera& cam, const dvec3 p[3], const dvec3 c[3]) {
  // Sutherland-Hodgman against z >= near: a triangle becomes 0, 3 or 4 vertices.
  CamVertex in[3] = {{p[0], c[0]}, {p[1], c[1]}, {p[2], c[2]}};
  CamVertex out[4];
  int n = 0;
  const double near = cam.near_m;
  for (int k = 0; k < 3; ++k) {
    const CamVertex& a = in[k];
    const CamVertex& b = in[(k + 1) % 3];
    const bool ain = a.p.z >= near, bin = b.p.z >= near;
    if (ain) out[n++] = a;
    if (ain != bin) {
      const double t = (near - a.p.z) / (b.p.z - a.p.z);
      out[n++] = {a.p + (b.p - a.p) * t, a.c + (b.c - a.c) * t};
    }
  }
  if (n < 3) return 0;
  const double f = cam.focal_px(), cx = 0.5 * img.width, cy = 0.5 * img.height;
  ScreenVertex s[4];
  for (int k = 0; k < n; ++k) {
    const dvec3& q = out[k].p;
    s[k] = {cx + f * q.x / q.z, cy - f * q.y / q.z, 1.0 / q.z, out[k].c};
  }
  std::size_t written = raster_triangle(img, s[0], s[1], s[2]);
  if (n == 4) written += raster_triangle(img, s[0], s[2], s[3]);
  return written;
}

Image render(const Ellipsoid& e, const HeightFn& height, const Camera& cam, const RenderOptions& opt,
             FrameStats* stats) {
  FrameStats st;
  Image img(cam.width, cam.height);
  const Geodetic here = e.to_geodetic(cam.pos);
  const dvec3 local_up = Ellipsoid::normal(here.lat, here.lon);
  const double alt = std::max(1.0, here.h);
  // 0 near the ground, 1 in space: fades the sky to black and the haze out.
  const double space = std::clamp(std::log10(alt / 20000.0) / 0.8, 0.0, 1.0);
  const dvec3 horizon{0.72, 0.82, 0.93}, zenith{0.18, 0.38, 0.78}, black{0.005, 0.005, 0.02};

  if (opt.sky) {
    const double f = cam.focal_px();
    for (int y = 0; y < img.height; ++y)
      for (int x = 0; x < img.width; ++x) {
        const dvec3 dir = normalize(cam.right * ((x + 0.5 - 0.5 * img.width) / f) -
                                    cam.up * ((y + 0.5 - 0.5 * img.height) / f) + cam.fwd);
        const double s = std::clamp(dot(dir, local_up), 0.0, 1.0);
        const dvec3 sky = horizon + (zenith - horizon) * std::sqrt(s);
        const dvec3 c = sky + (black - sky) * space;
        const std::size_t idx = static_cast<std::size_t>(y) * img.width + x;
        img.rgb[idx * 3 + 0] = to_byte(c.x);
        img.rgb[idx * 3 + 1] = to_byte(c.y);
        img.rgb[idx * 3 + 2] = to_byte(c.z);
      }
  }

  auto t0 = std::chrono::steady_clock::now();
  const Selection sel = select_tiles(e, cam, opt.select);
  st.select_ms = ms_since(t0);
  st.tiles = sel.visible.size();
  st.frustum_culled = sel.frustum_culled;
  st.horizon_culled = sel.horizon_culled;
  st.balance_splits = sel.balance_splits;

  const double haze_km = 120.0 + 2000.0 * space;
  std::vector<dvec3> cpos, ccol;
  for (const TileKey& k : sel.visible) {
    t0 = std::chrono::steady_clock::now();
    const TileMesh m = build_mesh(e, height, k, opt.select.grid);
    st.mesh_ms += ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    cpos.resize(m.pos.size());
    ccol.resize(m.pos.size());
    for (std::size_t v = 0; v < m.pos.size(); ++v) {
      const dvec3 world = m.center + to_double(m.pos[v]);
      cpos[v] = cam.to_camera(world);
      if (!opt.shading) {
        ccol[v] = {1, 1, 1};
        continue;
      }
      const double light = 0.16 + 0.84 * std::max(0.0, dot(to_double(m.normal[v]), opt.sun));
      const dvec3 lit = to_double(m.albedo[v]) * light;
      const double haze = 1.0 - std::exp(-length(world - cam.pos) / (haze_km * 1000.0));
      ccol[v] = lit + (horizon * 0.9 - lit) * (haze * (1.0 - 0.85 * space));
    }
    for (std::size_t t = 0; t + 2 < m.index.size(); t += 3) {
      const dvec3 p[3] = {cpos[m.index[t]], cpos[m.index[t + 1]], cpos[m.index[t + 2]]};
      const dvec3 c[3] = {ccol[m.index[t]], ccol[m.index[t + 1]], ccol[m.index[t + 2]]};
      draw_camera_triangle(img, cam, p, c);
      ++st.triangles;
    }
    st.raster_ms += ms_since(t0);
  }
  for (float d : img.inv_depth) st.pixels_covered += d > 0 ? 1 : 0;
  if (stats) *stats = st;
  return img;
}

std::uint64_t image_hash(const Image& img) {
  std::uint64_t h = 1469598103934665603ull;
  for (std::uint8_t b : img.rgb) h = (h ^ b) * 1099511628211ull;
  return h;
}

std::uint32_t crc32(const std::uint8_t* data, std::size_t n, std::uint32_t crc) {
  crc = ~crc;
  for (std::size_t i = 0; i < n; ++i) {
    crc ^= data[i];
    for (int k = 0; k < 8; ++k) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
  }
  return ~crc;
}

// PNG with an uncompressed (stored) zlib stream: no zlib dependency.
// tools/pngshrink.py recompresses for the docs.
bool write_png(const std::string& path, const Image& img) {
  std::vector<std::uint8_t> raw;
  raw.reserve(static_cast<std::size_t>(img.height) * (img.width * 3 + 1));
  for (int y = 0; y < img.height; ++y) {
    raw.push_back(0);  // filter: none
    const auto row = img.rgb.begin() + static_cast<std::ptrdiff_t>(y) * img.width * 3;
    raw.insert(raw.end(), row, row + img.width * 3);
  }
  std::vector<std::uint8_t> z = {0x78, 0x01};
  for (std::size_t off = 0; off < raw.size() || off == 0;) {
    const std::size_t len = std::min<std::size_t>(65535, raw.size() - off);
    const bool last = off + len == raw.size();
    z.push_back(last ? 1 : 0);
    z.push_back(static_cast<std::uint8_t>(len));
    z.push_back(static_cast<std::uint8_t>(len >> 8));
    z.push_back(static_cast<std::uint8_t>(~len));
    z.push_back(static_cast<std::uint8_t>(~len >> 8));
    z.insert(z.end(), raw.begin() + static_cast<std::ptrdiff_t>(off), raw.begin() + static_cast<std::ptrdiff_t>(off + len));
    off += len;
    if (last) break;
  }
  std::uint32_t s1 = 1, s2 = 0;
  for (std::uint8_t b : raw) {
    s1 = (s1 + b) % 65521;
    s2 = (s2 + s1) % 65521;
  }
  const std::uint32_t adler = s2 << 16 | s1;
  for (int k = 3; k >= 0; --k) z.push_back(static_cast<std::uint8_t>(adler >> (8 * k)));

  std::ofstream out(path, std::ios::binary);
  if (!out) return false;
  auto be32 = [&](std::uint32_t v) {
    const char b[4] = {static_cast<char>(v >> 24), static_cast<char>(v >> 16), static_cast<char>(v >> 8),
                       static_cast<char>(v)};
    out.write(b, 4);
  };
  auto chunk = [&](const char* type, const std::vector<std::uint8_t>& data) {
    be32(static_cast<std::uint32_t>(data.size()));
    std::vector<std::uint8_t> td(type, type + 4);
    td.insert(td.end(), data.begin(), data.end());
    out.write(reinterpret_cast<const char*>(td.data()), static_cast<std::streamsize>(td.size()));
    be32(crc32(td.data(), td.size()));
  };
  out.write("\x89PNG\r\n\x1a\n", 8);
  const auto w = static_cast<std::uint32_t>(img.width), h = static_cast<std::uint32_t>(img.height);
  chunk("IHDR", {static_cast<std::uint8_t>(w >> 24), static_cast<std::uint8_t>(w >> 16), static_cast<std::uint8_t>(w >> 8),
                 static_cast<std::uint8_t>(w), static_cast<std::uint8_t>(h >> 24), static_cast<std::uint8_t>(h >> 16),
                 static_cast<std::uint8_t>(h >> 8), static_cast<std::uint8_t>(h), 8, 2, 0, 0, 0});
  chunk("IDAT", z);
  chunk("IEND", {});
  return static_cast<bool>(out);
}

}  // namespace terra
