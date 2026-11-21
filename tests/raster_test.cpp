#include <gtest/gtest.h>

#include <cstring>
#include <fstream>
#include <vector>

#include "terra/render.hpp"

using namespace terra;

namespace {
ScreenVertex sv(double x, double y) { return {x, y, 1.0, {1, 1, 1}}; }
}  // namespace

TEST(Raster, MeshOnPixelCentresCoversEachPixelExactly1x) {
  // A fan of triangles whose edges all pass through pixel centres (the
  // hardest case for a fill rule). Each pixel must be written exactly once.
  const int n = 24;
  std::vector<int> hits(static_cast<std::size_t>(n) * n, 0);
  const ScreenVertex c = sv(12.5, 12.5);
  const ScreenVertex ring[] = {sv(2.5, 2.5), sv(12.5, 2.5), sv(22.5, 2.5), sv(22.5, 12.5),
                               sv(22.5, 22.5), sv(12.5, 22.5), sv(2.5, 22.5), sv(2.5, 12.5)};
  for (int k = 0; k < 8; ++k) {
    Image img(n, n);
    raster_triangle(img, c, ring[k], ring[(k + 1) % 8]);
    for (std::size_t p = 0; p < hits.size(); ++p) hits[p] += img.inv_depth[p] > 0;
  }
  int twice = 0, covered = 0;
  for (int h : hits) {
    twice += h > 1;
    covered += h > 0;
  }
  EXPECT_EQ(twice, 0);
  EXPECT_EQ(covered, 20 * 20);  // the 20x20 square of centres from 2.5 to 22.5, half-open
}

TEST(Raster, WindingDoesNotMatter) {
  Image a(16, 16), b(16, 16);
  EXPECT_EQ(raster_triangle(a, sv(1, 1), sv(15, 2), sv(4, 14)), raster_triangle(b, sv(1, 1), sv(4, 14), sv(15, 2)));
}

TEST(Raster, NearerWins) {
  Image img(8, 8);
  ScreenVertex far[3] = {sv(0, 0), sv(8, 0), sv(0, 8)};
  ScreenVertex near[3] = {sv(0, 0), sv(8, 0), sv(0, 8)};
  for (auto& v : near) v.inv_z = 2.0, v.color = {0, 0, 0};
  raster_triangle(img, near[0], near[1], near[2]);
  EXPECT_EQ(raster_triangle(img, far[0], far[1], far[2]), 0u);
  EXPECT_EQ(img.rgb[0], 0);
}

TEST(Raster, NearPlaneClipping) {
  Camera cam;
  cam.width = cam.height = 64;
  Image img(64, 64);
  const dvec3 white[3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  const dvec3 behind[3] = {{-1, -1, -5}, {1, -1, -5}, {0, 1, -5}};
  EXPECT_EQ(draw_camera_triangle(img, cam, behind, white), 0u);
  // A floor 1 m below the camera, running from behind it to 1 km ahead:
  // once clipped it must fill the bottom rows, and nothing above the horizon.
  const dvec3 floor[3] = {{-1000, -1, -10}, {1000, -1, -10}, {0, -1, 1000}};
  EXPECT_GT(draw_camera_triangle(img, cam, floor, white), 0u);
  EXPECT_GT(img.inv_depth[63 * 64 + 32], 0.f);
  EXPECT_GT(img.inv_depth[63 * 64 + 0], 0.f);
  for (int x = 0; x < 64; ++x) EXPECT_EQ(img.inv_depth[10 * 64 + x], 0.f);
}

TEST(Camera, FocalLengthFromVerticalFov) {
  Camera c;
  c.width = 640;
  c.height = 480;
  c.fov_y = deg(90);
  EXPECT_NEAR(c.focal_px(), 240.0, 1e-9);
}

TEST(Png, Crc32KnownAnswers) {
  const char* s = "123456789";
  EXPECT_EQ(crc32(reinterpret_cast<const std::uint8_t*>(s), 9), 0xCBF43926u);
  EXPECT_EQ(crc32(reinterpret_cast<const std::uint8_t*>("IEND"), 4), 0xAE426082u);
}

TEST(Png, WritesSignatureAndHeader) {
  Image img(3, 2);
  const std::string path = ::testing::TempDir() + "terra_test.png";
  ASSERT_TRUE(write_png(path, img));
  std::ifstream f(path, std::ios::binary);
  char buf[24];
  f.read(buf, 24);
  EXPECT_EQ(std::memcmp(buf, "\x89PNG\r\n\x1a\n", 8), 0);
  EXPECT_EQ(std::memcmp(buf + 12, "IHDR", 4), 0);
  EXPECT_EQ(static_cast<unsigned char>(buf[19]), 3);  // width
  EXPECT_EQ(static_cast<unsigned char>(buf[23]), 2);  // height
}
