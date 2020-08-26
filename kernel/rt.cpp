#include <iostream>

#include "ap_int.h"
#include "hls_stream.h"

#include "hlslib/xilinx/Stream.h"
#include "hlslib/xilinx/Simulation.h"

#include "rt.hpp"
#include "rtweekend.hpp"
#include "camera.hpp"

struct render_info
{
  uint16_t image_w;
  uint16_t image_h;
  uint16_t start_x;
  uint16_t start_y;
  uint16_t end_x;
  uint16_t end_y;
  uint16_t samples_per_pixel;

  uint16_t render_w;
  uint16_t render_h;
  uint32_t num_pixels;
  uint32_t num_rays;
};

struct ray_info
{
  uint32_t pixel_index;
  ray r;
  color attenuation;
  uint8_t hit_count;

  uint8_t object_index;
  hit_record rec;
  bool hit_anything;
  value_t closest_so_far;
};

struct object_info
{
  point3 center;
  value_t radius;
  color albedo;
  value_t fuzz;
};

void generate_ray(
  const render_info& p,
  hlslib::Stream<ray_info, 8>& ray_o
) {
  static xorshift32 rx(2463534242);
  static xorshift32 ry(1463534242);

  const float a = 256.0f / p.samples_per_pixel;

  camera cam;

samples:
  for (int s=0; s<p.samples_per_pixel; s++) {
scanlines:
    for (int y=0; y<p.render_h; y++) {
pixels:
      for (int x=0; x<p.render_w; x++) {
#pragma HLS PIPELINE II=1 rewind
        value_t rand_x = ((x + p.start_x) << 16) + (rx.next() & ((1<<16)-1));
        value_t rand_y = ((y + p.start_y) << 16) + (ry.next() & ((1<<16)-1));

        value_t u = rand_x / (p.image_w << 16);
        value_t v = rand_y / (p.image_h << 16);

        ray_info ri;
        ri.pixel_index = p.render_w * y + x;
        ri.r = cam.get_ray(u, v);
        ri.attenuation = color(a, a, a);
        ri.hit_count = 0;
        ri.object_index = 0;
        ri.hit_anything = false;
        ri.closest_so_far = infinity;

        ray_o.Push(ri);
      }
    }
  }
}

void merge_stream(
  hlslib::Stream<ray_info, 8>& ray_i,
  hlslib::Stream<ray_info, 1024>& ray_loop_i,
  hlslib::Stream<bool, 1>& done_i,
  hlslib::Stream<ray_info, 8>& ray_o
) {
#pragma HLS INLINE off // Prevent auto inlining
  while (true)
  {
#pragma HLS PIPELINE II=1
    ray_info ri;
    bool ri_valid = false;
    bool done;

    if (ray_loop_i.ReadNonBlocking(ri)) {
      ri_valid = true;
    } else if (ray_i.ReadNonBlocking(ri)) {
      ri_valid = true;
    } else if (done_i.ReadNonBlocking(done)) {
      break;
    }

    if (ri_valid) {
      ray_o.Push(ri);
    }
  }
}

void hit_sphere(
  const render_info& p,
  const hittable_list& world,
  hlslib::Stream<ray_info, 8>& ray_i,
  hlslib::Stream<bool, 1>& done_i,
  hlslib::Stream<ray_info, 8>& ray_o,
  hlslib::Stream<ray_info, 1024>& ray_loop_o
) {
  static random_in_unit_sphere rs;

  ray_info ri;
  bool done;

  while (true)
  {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=ray_i false
#pragma HLS DEPENDENCE variable=ray_o false
#pragma HLS DEPENDENCE variable=ray_loop_o false
#pragma HLS DEPENDENCE variable=ri false

    if (done_i.ReadNonBlocking(done)) {
      break;
    } else if (ray_i.ReadNonBlocking(ri)) {
      //bool loop = false;
      //
      //hit_record rec;
      //if (world.hit(ri.r, 0.001, infinity, rec)) {
      //  ri.r = ray(rec.p, rec.normal + rs.next());
      //  //ri.r = ray(rec.p, rec.normal);
      //  ri.attenuation *= 0.5;
      //  if (ri.hit_count < 8) loop = true;
      //  ri.hit_count++;
      //}

      bool loop = true;
      hit_record temp_rec;
      if (world.objects[ri.object_index].hit(ri.r, 0.001, ri.closest_so_far, temp_rec)) {
        ri.hit_anything = true;
        ri.closest_so_far = temp_rec.t;
        ri.rec = temp_rec;
      }

      if (ri.object_index == MAX_OBJECTS - 1) {
        loop = false;
        if (ri.hit_anything) {
          ri.r = ray(ri.rec.p, ri.rec.normal + rs.next());
          ri.attenuation *= 0.5;
          if (ri.hit_count < 8) loop = true;
          ri.hit_count++;
          // Reset
          ri.object_index = 0;
          ri.rec;
          ri.hit_anything = false;
          ri.closest_so_far = infinity;
        }
      } else {
        ri.object_index++;
      }

      if (loop) {
        ray_loop_o.Push(ri);
      } else {
        ray_o.Push(ri);
      }
    }
  }
}

void light(
  const render_info& p,
  hlslib::Stream<ray_info, 8>& ray_i,
  hlslib::Stream<ray_info, 8>& ray_o,
  hlslib::Stream<bool, 1> done_o[2]
) {
  for (int i=0; i<p.num_rays; i++) {
#pragma HLS PIPELINE II=1
    auto ri = ray_i.Pop();

    color rc { 0, 0, 0 };

    if (ri.hit_count < 8) {
      vec3 unit_direction = unit_vector(ri.r.direction());
      value_t t = value_t(0.5)*(unit_direction.y() + value_t(1.0));
      rc = (value_t(1.0)-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
      rc[0] *= ri.attenuation[0];
      rc[1] *= ri.attenuation[1];
      rc[2] *= ri.attenuation[2];
    }

    ri.attenuation = rc;

    ray_o.Push(ri);
  }

  done_o[0].Push(true);
  done_o[1].Push(true);
}

void buffer(
  const render_info& p,
  hlslib::Stream<ray_info, 8>& ray_i,
  hlslib::Stream<pixel_block, 256>& color_o
) {
  constexpr int num_urams = std::tuple_size<pixel_block>::value;
  constexpr int num_bytes = num_urams * 32 * 1024;
  constexpr int max_pixels = num_bytes / 4 / 3; // float
  static float buf_r[4096*2][num_urams] = {0};
  static float buf_g[4096*2][num_urams] = {0};
  static float buf_b[4096*2][num_urams] = {0};
#pragma HLS BIND_STORAGE variable=buf_r type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_g type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_b type=RAM_S2P impl=uram
#pragma HLS ARRAY_PARTITION variable=buf_r complete dim=2
#pragma HLS ARRAY_PARTITION variable=buf_g complete dim=2
#pragma HLS ARRAY_PARTITION variable=buf_b complete dim=2
#pragma HLS ARRAY_RESHAPE variable=buf_r cyclic factor=2 dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_g cyclic factor=2 dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_b cyclic factor=2 dim=1

accum:
  for (int i=0; i<p.num_rays; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=buf_r inter false
#pragma HLS DEPENDENCE variable=buf_g inter false
#pragma HLS DEPENDENCE variable=buf_b inter false
    auto ri = ray_i.Pop();

    int index = ri.pixel_index;
    int addr_a = index / (num_urams);
    int addr_b = index % (num_urams);
    buf_r[addr_a][addr_b] += ri.attenuation[0];
    buf_g[addr_a][addr_b] += ri.attenuation[1];
    buf_b[addr_a][addr_b] += ri.attenuation[2];
  }

output:
  for (int i=0; i<p.num_pixels/num_urams; i++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=buf_r inter false
#pragma HLS DEPENDENCE variable=buf_g inter false
#pragma HLS DEPENDENCE variable=buf_b inter false
    pixel_block c;
    for (int j=0; j<num_urams; j++) {
      c[j].r = clamp(int(buf_r[i][j]), 0, 255);
      c[j].g = clamp(int(buf_g[i][j]), 0, 255);
      c[j].b = clamp(int(buf_b[i][j]), 0, 255);
      c[j].a = 255;
      buf_r[i][j] = 0;
      buf_g[i][j] = 0;
      buf_b[i][j] = 0;
    }
    color_o.Push(c);
  }
}

void write_mem(
  const render_info& p,
  hlslib::Stream<pixel_block, 256>& color_i,
  pixel_block* image
) {
  constexpr auto block_w = std::tuple_size<pixel_block>::value;
  for (int y=p.start_y; y<p.end_y; y++) {
    for (int x=0; x<p.render_w/block_w; x++) {
#pragma HLS PIPELINE II=1 rewind
      image[y * p.image_w/block_w + p.start_x/block_w + x] = color_i.Pop();
    }
  }
}

void rt(
  const int image_w,
  const int image_h,
  const int start_x,
  const int start_y,
  const int end_x,
  const int end_y,
  const int samples_per_pixel,
  pixel_block* image
) {
#pragma HLS INTERFACE m_axi port=image offset=slave bundle=gmem
#pragma HLS INTERFACE s_axilite port=image_w
#pragma HLS INTERFACE s_axilite port=image_h
#pragma HLS INTERFACE s_axilite port=start_x
#pragma HLS INTERFACE s_axilite port=start_y
#pragma HLS INTERFACE s_axilite port=end_x
#pragma HLS INTERFACE s_axilite port=end_y
#pragma HLS INTERFACE s_axilite port=samples_per_pixel
#pragma HLS INTERFACE s_axilite port=image
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS INTERFACE ap_ctrl_hs port=return

#pragma HLS DATAFLOW

  render_info p;
  p.image_w = image_w;
  p.image_h = image_h;
  p.start_x = start_x;
  p.start_y = start_y;
  p.end_x = end_x;
  p.end_y = end_y;
  p.samples_per_pixel = samples_per_pixel;

  p.render_w = end_x - start_x;
  p.render_h = end_y - start_y;
  p.num_pixels = p.render_w * p.render_h;
  p.num_rays = p.num_pixels * samples_per_pixel;

  hittable_list world;
  world.objects[0].set(point3(0,0,-1), 0.5);
  world.objects[1].set(point3(0,-100.5,-1), 100);

  HLSLIB_DATAFLOW_INIT();

  hlslib::Stream<ray_info, 8>      ray_strm0;
  hlslib::Stream<ray_info, 1024>   ray_strm1;
  hlslib::Stream<ray_info, 8>      ray_strm2;
  hlslib::Stream<ray_info, 8>      ray_strm3;
  hlslib::Stream<ray_info, 8>      ray_strm4;
  hlslib::Stream<pixel_block, 256> color_strm0("pixel_strm");
  hlslib::Stream<bool, 1>          done_strm0[2];

  HLSLIB_DATAFLOW_FUNCTION(generate_ray, p, ray_strm0);
  HLSLIB_DATAFLOW_FUNCTION(merge_stream, ray_strm0, ray_strm1, done_strm0[0], ray_strm4);
  HLSLIB_DATAFLOW_FUNCTION(hit_sphere, p, world, ray_strm4, done_strm0[1], ray_strm2, ray_strm1);
  HLSLIB_DATAFLOW_FUNCTION(light, p, ray_strm2, ray_strm3, done_strm0);
  HLSLIB_DATAFLOW_FUNCTION(buffer, p, ray_strm3, color_strm0);
  HLSLIB_DATAFLOW_FUNCTION(write_mem, p, color_strm0, image);

  HLSLIB_DATAFLOW_FINALIZE();
}

