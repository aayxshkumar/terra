#include "terra/camera.hpp"

#include <cmath>

namespace terra {

double Camera::focal_px() const { return 0.5 * height / std::tan(0.5 * fov_y); }

dvec3 Camera::to_camera(const dvec3& world) const {
  const dvec3 v = world - pos;
  return {dot(v, right), dot(v, up), dot(v, fwd)};
}

Camera Camera::from_geodetic(const Ellipsoid& e, const Geodetic& where, double heading_deg, double pitch_deg,
                             double fov_y_deg, int width, int height) {
  const Enu f = Ellipsoid::enu(where.lat, where.lon);
  const double h = deg(heading_deg), p = deg(pitch_deg);
  Camera c;
  c.pos = e.to_ecef(where);
  const dvec3 horizontal = f.east * std::sin(h) + f.north * std::cos(h);
  c.fwd = normalize(horizontal * std::cos(p) + f.up * std::sin(p));
  // Right stays horizontal whatever the pitch, so looking straight down is fine.
  c.right = normalize(f.east * std::cos(h) - f.north * std::sin(h));
  c.up = cross(c.right, c.fwd);
  c.fov_y = deg(fov_y_deg);
  c.width = width;
  c.height = height;
  return c;
}

Camera Camera::look_at(const dvec3& pos, const dvec3& target, const dvec3& up_hint, double fov_y_deg, int width,
                       int height) {
  Camera c;
  c.pos = pos;
  c.fwd = normalize(target - pos);
  c.right = normalize(cross(c.fwd, up_hint));
  c.up = cross(c.right, c.fwd);
  c.fov_y = deg(fov_y_deg);
  c.width = width;
  c.height = height;
  return c;
}

}  // namespace terra
