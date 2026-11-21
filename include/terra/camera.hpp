// Pinhole camera in ECEF. Camera space: x = right, y = up, z = forward
// (depth). Square pixels.
#pragma once

#include "terra/ellipsoid.hpp"
#include "terra/vec.hpp"

namespace terra {

struct Camera {
  dvec3 pos, fwd{1, 0, 0}, right{0, 1, 0}, up{0, 0, 1};
  double fov_y = deg(60);
  int width = 960, height = 540;
  double near_m = 1.0;

  double focal_px() const;
  dvec3 to_camera(const dvec3& world) const;  // (right, up, forward) components of world - pos

  // heading: degrees clockwise from north; pitch: degrees above the local
  // horizon (-90 looks straight down).
  static Camera from_geodetic(const Ellipsoid& e, const Geodetic& where, double heading_deg, double pitch_deg,
                              double fov_y_deg, int width, int height);
  static Camera look_at(const dvec3& pos, const dvec3& target, const dvec3& up_hint, double fov_y_deg, int width,
                        int height);
};

}  // namespace terra
