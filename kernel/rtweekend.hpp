#pragma once

// Common Headers
#include "ray.hpp"
#include "vec3.hpp"
#include "hittable.hpp"
#include "random.hpp"

//constexpr value_t infinity_f() {
//  value_t v { 0 };
//  v.range() = -1LL;
//  v[value_t::width-1] = 0;
//  return v;
//}

//constexpr value_t pi_f() {
//  value_t v { 3.1415926535897932385 };
//  return v;
//}

// Utility Functions
inline value_t degrees_to_radians(value_t degrees) {
  return degrees * pi / 180;
}

inline value_t clamp(value_t x, value_t min, value_t max) {
  if (x < min) return min;
  if (x > max) return max;
  return x;
}

template<typename T>
constexpr const T& clamp(const T& v, const T& low, const T& high) {
  return v < low ? low : v > high ? high : v;
}
