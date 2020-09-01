#pragma once

#include "ray.hpp"
#include "object.hpp"

const int MAX_OBJECTS = 8;

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

bool hit(const object& obj, const ray& r, value_t t_min, value_t t_max, hit_record& rec) {
#pragma HLS INLINE
  vec3 oc = r.origin() - obj.center;
  value_t a = r.direction().length_squared();
  value_t half_b = dot(oc, r.direction());
  value_t c = oc.length_squared() - obj.radius*obj.radius;
  value_t discriminant = half_b*half_b - a*c;

  if (discriminant > 0) {
    value_t root = hls::sqrt(discriminant);
    value_t temp = (-half_b - root)/a;
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.at(rec.t);
      vec3 outward_normal = (rec.p - obj.center) / obj.radius;
      rec.set_face_normal(r, outward_normal);
      return true;
    }
    temp = (-half_b + root) / a;
    if (temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.at(rec.t);
      vec3 outward_normal = (rec.p - obj.center) / obj.radius;
      rec.set_face_normal(r, outward_normal);
      return true;
    }
  }
  return false;
}
