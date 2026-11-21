// Software renderer: select tiles, build their meshes, shade per vertex, clip
// against the near plane in camera space, and rasterize with a 1/z depth
// buffer (float 1/z keeps relative precision from 1 m to the far horizon,
// so there is no far plane to tune).
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "terra/camera.hpp"
#include "terra/ellipsoid.hpp"
#include "terra/select.hpp"
#include "terra/terrain.hpp"

namespace terra {

struct Image {
  int width = 0, height = 0;
  std::vector<std::uint8_t> rgb;
  std::vector<float> inv_depth;  // 0 = nothing drawn
  Image() = default;
  Image(int w, int h);
};

struct RenderOptions {
  SelectOptions select;
  dvec3 sun = normalize(dvec3{1.0, 0.35, 0.45});  // ECEF direction towards the sun
  bool sky = true;
  bool shading = true;  // false: flat white, for coverage measurements
};

struct FrameStats {
  std::size_t tiles = 0, triangles = 0, pixels_covered = 0;
  std::size_t frustum_culled = 0, horizon_culled = 0, balance_splits = 0;
  double select_ms = 0, mesh_ms = 0, raster_ms = 0;
};

Image render(const Ellipsoid& e, const HeightFn& height, const Camera& cam, const RenderOptions& opt,
             FrameStats* stats = nullptr);

// Lower-level pieces, exposed for tests.
struct ScreenVertex {
  double x, y;   // pixels, y down
  double inv_z;  // 1 / camera depth
  dvec3 color;
};
// Fills pixels whose centres are inside the triangle, top-left rule on
// shared edges, depth-tested on inv_z. Returns pixels written.
std::size_t raster_triangle(Image& img, ScreenVertex v0, ScreenVertex v1, ScreenVertex v2);
// Camera-space triangle (x right, y up, z forward): clipped at z = near, then projected.
std::size_t draw_camera_triangle(Image& img, const Camera& cam, const dvec3 p[3], const dvec3 c[3]);

std::uint64_t image_hash(const Image& img);
std::uint32_t crc32(const std::uint8_t* data, std::size_t n, std::uint32_t crc = 0);
bool write_png(const std::string& path, const Image& img);

}  // namespace terra
