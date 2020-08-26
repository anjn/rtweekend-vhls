#pragma once

#include "ray.hpp"

const int MAX_OBJECTS = 2;

struct hit_record
{
  point3 p;
  vec3 normal;
  value_t t;
  bool front_face;

  inline void set_face_normal(const ray& r, const vec3& outward_normal) {
    front_face = dot(r.direction(), outward_normal) < 0;
    normal = front_face ? outward_normal :-outward_normal;
  }
};

struct hittable
{
  bool valid;

  point3 center;
  value_t radius;

  void set(const point3& c, value_t r) {
#pragma HLS INLINE
    valid = true;
    center = c;
    radius = r;
  }

  bool hit_sphere(const ray& r, value_t t_min, value_t t_max, hit_record& rec) const {
#pragma HLS INLINE
    vec3 oc = r.origin() - center;
    value_t a = r.direction().length_squared();
    value_t half_b = dot(oc, r.direction());
    value_t c = oc.length_squared() - radius*radius;
    value_t discriminant = half_b*half_b - a*c;

    if (discriminant > 0) {
      value_t root = hls::sqrt(discriminant);
      value_t temp = (-half_b - root)/a;
      if (temp < t_max && temp > t_min) {
        rec.t = temp;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        return true;
      }
      temp = (-half_b + root) / a;
      if (temp < t_max && temp > t_min) {
        rec.t = temp;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        return true;
      }
    }
    return false;
  }

  bool hit(const ray& r, value_t t_min, value_t t_max, hit_record& rec) const {
#pragma HLS INLINE
    return hit_sphere(r, t_min, t_max, rec);
  }
};

struct hittable_list
{
  hittable objects[MAX_OBJECTS];

  hittable_list() {
  }

  bool hit(const ray& r, value_t t_min, value_t t_max, hit_record& rec) const {
#pragma HLS INLINE
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects) {
#pragma HLS UNROLL
      if (object.valid && object.hit(r, t_min, closest_so_far, temp_rec)) {
        hit_anything = true;
        closest_so_far = temp_rec.t;
        rec = temp_rec;
      }
    }

    return hit_anything;
  }
};
