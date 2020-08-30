#pragma once

#ifndef OPENCL
#include "vec3.hpp"
#else
struct vec3 {
  vec3() {}
  vec3(float e0, float e1, float e2) : e {e0, e1, e2} {}
  float e[3];
};
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color
#endif

struct material
{
  uint8_t type;

  color albedo;
  float fuzz;
  float ref_idx;

  // pre calc
  float ref_idx_inv;
  float r0;
  float r0_inv;

  void pre_calculation() {
    ref_idx_inv = 1.0f / ref_idx;
    r0 = (1 - ref_idx) / (1 + ref_idx);

    r0 = r0 * r0;

    r0_inv = (1 - ref_idx_inv) / (1 + ref_idx_inv);
    r0_inv = r0_inv * r0_inv;
  }
};

inline material make_lambertian(color albedo)
{
  material m;
  m.type = 0;
  m.albedo = albedo;
  return m;
}

inline material make_metal(color albedo, float fuzz)
{
  material m;
  m.type = 1;
  m.albedo = albedo;
  m.fuzz = fuzz;
  return m;
}

inline material make_dielectric(float ref_idx)
{
  material m;
  m.type = 2;
  m.albedo = color(1,1,1);
  m.ref_idx = ref_idx;
  m.pre_calculation();
  return m;
}

struct object
{
  vec3 center;
  float radius;

  material m;
};

inline object make_sphere(vec3 center, float radius, material m) {
  object o;
  o.center = center;
  o.radius = radius;
  o.m = m;
  return o;
}
