#pragma once

#include "ap_fixed.h"
#include "hls_math.h"

#include "random.hpp"

using value_t = float;

// Constants
const value_t infinity = 10000.0;
const value_t pi = 3.1415926535897932385;

struct vec3 {
  //vec3() : e{0,0,0} {}
  vec3() {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=e complete
  }
  vec3(value_t e0, value_t e1, value_t e2) : e {e0, e1, e2} {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=e complete
  }

  value_t x() const { return e[0]; }
  value_t y() const { return e[1]; }
  value_t z() const { return e[2]; }

  vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
  value_t operator[](int i) const { return e[i]; }
  value_t& operator[](int i) { return e[i]; }

  vec3& operator+=(const vec3 &v) {
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];
    return *this;
  }

  vec3& operator*=(const value_t t) {
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;
    return *this;
  }

  vec3& operator/=(const value_t t) {
    //return *this *= 1/t;
    e[0] /= t;
    e[1] /= t;
    e[2] /= t;
    return *this;
  }

  value_t length() const {
    return hls::sqrt(length_squared());
  }

  value_t length_squared() const {
    return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
  }

//  template<typename Random>
//  inline static vec3 random(Random& r) {
//    return vec3(random_norm::gen(r), random_norm::gen(r), random_norm::gen(r));
//  };
//
//  template<typename Random>
//  inline static vec3 random(Random& r, value_t min, value_t max) {
//    return vec3(random_norm::gen(r, min, max), random_norm::gen(r, min, max), random_norm::gen(r, min, max));
//  };

  value_t e[3];
};

// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color

#ifndef __SYNTHESIS__
inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
  return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}
#endif

inline vec3 operator+(const vec3 &u, const vec3 &v) {
  return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
  return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
  return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(value_t t, const vec3 &v) {
  return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, value_t t) {
  return t * v;
}

inline vec3 operator/(vec3 v, value_t t) {
  //return (1/t) * v;
  v.e[0] /= t;
  v.e[1] /= t;
  v.e[2] /= t;
  return v;
}

inline value_t dot(const vec3 &u, const vec3 &v) {
  return u.e[0] * v.e[0]
    + u.e[1] * v.e[1]
    + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
  return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
              u.e[2] * v.e[0] - u.e[0] * v.e[2],
              u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
  auto len = v.length();
  assert(len>0);
  return len == 0 ? v : v / len;
}

inline vec3 reflect(const vec3& v, const vec3& n) {
  return v - 2 * dot(v, n) * n;
}

inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
  auto cos_theta = dot(-uv, n);
  vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
  vec3 r_out_parallel = -hls::sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
  return r_out_perp + r_out_parallel;
}

inline float schlick(float cosine, float r0) {
  return r0 + (1-r0) * hls::pow((1.0f - cosine), 5.0f);
}

struct random_uniform
{
  xorshift32 r;
  uint32_t raw;
  float raw_float;

  random_uniform(uint32_t seed): r(seed) {}

  float next() {
#pragma HLS INLINE
    raw = r.next();
    raw_float = float(raw) * xorshift32::max_float_inv;
    return raw_float;
  }

  float next(float max) {
#pragma HLS INLINE
    return next() * max;
  }

  float next(float min, float max) {
#pragma HLS INLINE
    return next(max - min) + min;
  }
};

struct random_in_unit_sphere
{
  random_uniform ru;
  random_uniform rv;
  random_uniform rw;

  random_in_unit_sphere():
    ru(12345),
    rv(12346),
    rw(12347)
  {}

  vec3 next() {
#pragma HLS INLINE
    vec3 vec;
    value_t u = ru.next();
    value_t v = rv.next();
    value_t w = hls::pow(rw.next(), value_t(1.0/3.0));
    vec[0] = w * (1 - 2 * u);
    value_t z = hls::sqrt(1 - vec[0] * vec[0]);
    vec[1] = w * z * hls::cos(2 * pi * v);
    vec[2] = w * z * hls::sin(2 * pi * v);
    return vec;
  }
};

struct random_unit_vector
{
  xorshift32 ra, rz;

  random_unit_vector():
    ra(22222),
    rz(22223)
  {}

  vec3 next() {
#pragma HLS INLINE
    auto a = random_norm<value_t>::gen(ra, 0, 2*pi);
    auto z = random_norm<value_t>::gen(rz, -1, 1);
    auto r = hls::sqrt(1 - z*z);
    return vec3(r * hls::cos(a), r * hls::sin(a), z);
  }
};
