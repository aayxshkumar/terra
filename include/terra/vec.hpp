// Minimal 3-vector. World coordinates are always double (ECEF metres); only
// tile-local mesh offsets are stored as float.
#pragma once

#include <cmath>

namespace terra {

template <class T>
struct Vec3 {
  T x = 0, y = 0, z = 0;

  constexpr Vec3() = default;
  constexpr Vec3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator-() const { return {-x, -y, -z}; }
  constexpr Vec3 operator*(T s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(T s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) { return *this = *this + o; }
  Vec3& operator-=(const Vec3& o) { return *this = *this - o; }
  constexpr bool operator==(const Vec3&) const = default;
};

template <class T>
constexpr Vec3<T> operator*(T s, const Vec3<T>& v) { return v * s; }

template <class T>
constexpr T dot(const Vec3<T>& a, const Vec3<T>& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

template <class T>
constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template <class T>
T length(const Vec3<T>& v) { return std::sqrt(dot(v, v)); }

template <class T>
Vec3<T> normalize(const Vec3<T>& v) {
  const T l = length(v);
  return l > 0 ? v / l : v;
}

using dvec3 = Vec3<double>;
using fvec3 = Vec3<float>;

inline fvec3 to_float(const dvec3& v) {
  return {static_cast<float>(v.x), static_cast<float>(v.y), static_cast<float>(v.z)};
}
inline dvec3 to_double(const fvec3& v) { return {v.x, v.y, v.z}; }

constexpr double kPi = 3.14159265358979323846;
constexpr double deg(double d) { return d * kPi / 180.0; }
constexpr double to_deg(double r) { return r * 180.0 / kPi; }

}  // namespace terra
