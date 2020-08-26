#include <iostream>
#include <tuple>

#include "ap_int.h"
#include "hls_stream.h"

#include "hlslib/xilinx/Stream.h"
#include "hlslib/xilinx/Simulation.h"
#include "hlslib/xilinx/Utility.h"

#include "rt.hpp"
#include "rtweekend.hpp"
#include "camera.hpp"

const int MAX_HIT = 16;
const uint16_t MAX_BUFFER_HEIGHT = 32;

struct render_info
{
  uint16_t image_w;
  uint16_t image_h;
  uint16_t start_x;
  uint16_t start_y;
  uint16_t end_x;
  uint16_t end_y;
  uint16_t samples_per_pixel;

  uint8_t cache_w_log2;
  uint8_t cache_h_log2;
  uint16_t render_w;
  uint16_t render_h;
  uint32_t num_pixels;
  uint32_t num_rays;
};

struct ray_info
{
  ray r;
  color attenuation;
  uint8_t hit_count;
};

struct hit_info
{
  uint8_t object_index;
  hit_record rec;
  bool hit_anything;
  value_t closest_so_far;
};

void init(hit_info& hi) {
#pragma HLS INLINE
  hi.object_index = 0;
  hi.hit_anything = false;
  hi.closest_so_far = infinity;
}

struct pixel_position
{
  uint16_t x;
  uint16_t y;
};

struct pixel_info
{
  color p;
  uint16_t sample_count;
};

struct object_info
{
  point3 center;
  value_t radius;
  color albedo;
  value_t fuzz;
};

using pos_stream            = hlslib::Stream<pixel_position, 8>;
using pix_stream            = hlslib::Stream<std::tuple<pixel_position, pixel_info>, 8>;
using pix_stream_fb         = hlslib::Stream<std::tuple<pixel_position, pixel_info>, 1024>;
using pix_ray_stream        = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info>, 8>;
using pix_ray_hit_stream    = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info, hit_info>, 8>;
using pix_ray_hit_stream_fb = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info, hit_info>, 512>;
using ray_hit_stream        = hlslib::Stream<std::tuple<ray_info, hit_info>, 8>;
using pix_blk_stream        = hlslib::Stream<pixel_block, 128>;
using line_stream           = hlslib::Stream<uint16_t, 8>;
using line_stream_fb        = hlslib::Stream<uint16_t, MAX_BUFFER_HEIGHT>;
using done_stream           = hlslib::Stream<bool>;

void generate_pixel(
  const render_info& p,
  line_stream_fb& y_i,
  pos_stream& pos_o
) {
scanlines:
  for (uint16_t y=0; y<p.render_h; y++) {
    uint16_t line = y;
    if (line >= MAX_BUFFER_HEIGHT) {
      line = y_i.Pop();
    }
pixels:
    for (uint16_t x=0; x<p.render_w; x++) {
#pragma HLS PIPELINE II=1 rewind
      pos_o.Push({x, line});
    }
  }
}

void merge_pixel(
  pos_stream& pos_i,
  pix_stream_fb& pix_loop_i,
  done_stream& done_i,
  pix_stream& pix_o
) {
#pragma HLS INLINE off // Prevent auto inlining
  while (true)
  {
#pragma HLS PIPELINE II=1
    pixel_position pos;
    pixel_info pix;
    bool input_valid = false;
    bool done;

    if (!pix_loop_i.IsEmpty())
    {
      input_valid = true;

      std::tie(pos, pix) = pix_loop_i.Pop();
    }
    else if (!pos_i.IsEmpty())
    {
      input_valid = true;

      pos = pos_i.Pop();
      // Init pixel_info
      pix.p[0] = 0.0f;
      pix.p[1] = 0.0f;
      pix.p[2] = 0.0f;
      pix.sample_count = 0;
    }
    else if (done_i.ReadNonBlocking(done))
    {
      break;
    }

    if (input_valid) {
      pix_o.Push({pos, pix});
    }
  }
}

void generate_ray(
  const render_info& p,
  pix_stream& pix_i,
  pix_ray_stream& ray_o
) {
  static xorshift32 rx(2463534242);
  static xorshift32 ry(1463534242);

  const float a = 256.0f / p.samples_per_pixel;

  camera cam;

  for (uint32_t i=0; i<p.num_rays; i++) {
//#pragma HLS PIPELINE
    //auto [pos, pix] = pix_i.Pop();
    pixel_position pos;
    pixel_info pix;
    std::tie(pos, pix) = pix_i.Pop();

    value_t rand_x = ((pos.x + p.start_x) << 16) + (rx.next() & ((1<<16)-1));
    value_t rand_y = ((pos.y + p.start_y) << 16) + (ry.next() & ((1<<16)-1));

    value_t u = rand_x / (p.image_w << 16);
    value_t v = rand_y / (p.image_h << 16);

    ray_info ri;
    ri.r = cam.get_ray(u, v);
    ri.attenuation = color(a, a, a);
    ri.hit_count = 0;

    ray_o.Push({pos, pix, ri});
  }
}

void merge_hit(
  pix_ray_stream& ray_i,
  pix_ray_hit_stream_fb& ray_loop_i,
  done_stream& done_i,
  ray_hit_stream& ray_o,
  pix_ray_hit_stream& pix_o
) {
#pragma HLS INLINE off // Prevent auto inlining
  while (true)
  {
#pragma HLS PIPELINE II=1
    pixel_position pos;
    pixel_info pix;
    ray_info ri;
    hit_info hi;
    bool input_valid = false;
    bool done;

    if (!ray_loop_i.IsEmpty())
    {
      input_valid = true;

      std::tie(pos, pix, ri, hi) = ray_loop_i.Pop();
    }
    else if (!ray_i.IsEmpty())
    {
      input_valid = true;

      std::tie(pos, pix, ri) = ray_i.Pop();
      // Init hit_info
      init(hi);
    }
    else if (done_i.ReadNonBlocking(done))
    {
      break;
    }

    if (input_valid) {
      ray_o.Push({ri, hi});
      pix_o.Push({pos, pix, ri, hi});
    }
  }
}

void hit_sphere(
  const render_info& p,
  const hittable_list& world,
  ray_hit_stream& ray_i,
  pix_ray_hit_stream& pix_i,
  done_stream& done_i,
  pix_ray_stream& ray_o,
  pix_ray_hit_stream_fb& ray_loop_o
) {
  static random_in_unit_sphere rs;

hit_sphere_main:
  while (true)
  {
    bool done;

    if (done_i.ReadNonBlocking(done)) {
      break;
    } else if (!ray_i.IsEmpty()) {

      bool hit;
      hit_record rec;
      {
        //auto [ri, hi] = ray_i.Pop();
        ray_info ri;
        hit_info hi;
        std::tie(ri, hi) = ray_i.Pop();
        hit = world.objects[hi.object_index].hit(ri.r, 0.001, hi.closest_so_far, rec);
      }

      pixel_position pos;
      pixel_info pix;
      ray_info ri;
      hit_info hi;
      std::tie(pos, pix, ri, hi) = pix_i.Pop();
      //auto [pos, pix, ri, hi] = pix_i.Pop();

      if (hit) {
        hi.hit_anything = true;
        hi.closest_so_far = rec.t;
        hi.rec = rec;
      }

      auto rs_v = rs.next();
      bool loop = true;

      if (hi.object_index == MAX_OBJECTS - 1) {
        loop = false;

        if (hi.hit_anything) {
          // Create new ray
          ri.r = ray(hi.rec.p, hi.rec.normal + rs_v);
          ri.attenuation *= 0.5;
          if (ri.hit_count < MAX_HIT) loop = true;
          ri.hit_count++;
          // Reset hit_info
          init(hi);
        }
      } else {
        hi.object_index++;
      }

      if (loop) {
        ray_loop_o.Push({pos, pix, ri, hi});
      } else {
        ray_o.Push({pos, pix, ri});
      }
    }
  }
}

void light(
  const render_info& p,
  pix_ray_stream& ray_i,
  pix_stream& pix_o,
  pix_stream_fb& pix_loop_o,
  done_stream done_o[3]
) {
  for (int i=0; i<p.num_rays; i++) {
#pragma HLS PIPELINE
    pixel_position pos;
    pixel_info pix;
    ray_info ri;
    std::tie(pos, pix, ri) = ray_i.Pop();
    //auto [pos, pix, ri] = ray_i.Pop();

    color rc { 0, 0, 0 };

    if (ri.hit_count < MAX_HIT) {
      vec3 unit_direction = unit_vector(ri.r.direction());
      value_t t = value_t(0.5)*(unit_direction.y() + value_t(1.0));
      rc = (value_t(1.0)-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
      rc[0] *= ri.attenuation[0];
      rc[1] *= ri.attenuation[1];
      rc[2] *= ri.attenuation[2];
    }

    pix.p += rc;
    pix.sample_count++;

    if (pix.sample_count == p.samples_per_pixel) {
      pix_o.Push({pos, pix});
    } else {
      pix_loop_o.Push({pos, pix});
    }
  }

  done_o[0].Push(true);
  done_o[1].Push(true);
  done_o[2].Push(true);
}

struct buffer_memory
{
  constexpr static int num_urams = std::tuple_size<pixel_block>::value;
  constexpr static int num_pixels = 4096 * 8 * num_urams;
  constexpr static int num_pixels_clog2 = hlslib::ConstLog2(num_pixels - 1); // 19

  // URAM: 64bit(=uint8 * 8) 4K words
  uint16_t buffer_w;
  uint16_t buffer_h;
  uint8_t buffer_w_clog2;

  buffer_memory(const render_info& p)
  {
#pragma HLS INLINE

    // width
    buffer_w_clog2 = 8;
calc_w_clog2:
    do {
      buffer_w_clog2++;
      buffer_w = 1 << buffer_w_clog2;
    } while (p.render_w > buffer_w);

    // height
    buffer_h = std::min(p.render_h, MAX_BUFFER_HEIGHT);
  }
};


template<typename T, int N>
struct fifo
{
  std::array<T, N> arr;
  uint16_t wp, rp;

  fifo() : wp(0), rp(0) {}

  bool empty() {
#pragma HLS INLINE
    return wp == rp;
  }

  inline uint16_t next_ptr(uint16_t p) {
    uint16_t np = p + 1;
    return np == N ? 0 : np;
  }

  bool full() {
#pragma HLS INLINE
    return next_ptr(wp) == rp;
  }

  void push(const T& v) {
#pragma HLS INLINE
    assert(!full());
    arr[wp] = v;
    wp = next_ptr(wp);
  }

  T pop() {
#pragma HLS INLINE
    assert(!empty());
    T v = arr[rp];
    rp = next_ptr(rp);
    return v;
  }
};

void buffer_write(
  const render_info& p,
  pix_stream& pix_i,
  line_stream& y_o,
  line_stream_fb& y_loop_o,
  pix_blk_stream& blk_o
) {
  uint8_t buf_r[4096*8][buffer_memory::num_urams];
  uint8_t buf_g[4096*8][buffer_memory::num_urams];
  uint8_t buf_b[4096*8][buffer_memory::num_urams];
#pragma HLS BIND_STORAGE variable=buf_r type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_g type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_b type=RAM_S2P impl=uram
#pragma HLS ARRAY_PARTITION variable=buf_r complete dim=2
#pragma HLS ARRAY_PARTITION variable=buf_g complete dim=2
#pragma HLS ARRAY_PARTITION variable=buf_b complete dim=2
#pragma HLS ARRAY_RESHAPE variable=buf_r cyclic factor=8 dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_g cyclic factor=8 dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_b cyclic factor=8 dim=1
#pragma HLS SHARED variable=buf_r
#pragma HLS SHARED variable=buf_g
#pragma HLS SHARED variable=buf_b

  uint16_t pixel_count[MAX_BUFFER_HEIGHT];

  buffer_memory buf(p);

init:
  for (uint16_t y=0; y<buf.buffer_h; y++) {
#pragma HLS PIPELINE
    pixel_count[y] = 0;
  }

  fifo<uint16_t, 256> fifo_read_y;

  bool read = false;
  uint16_t read_y = 0;
  uint16_t read_x = 0;
  uint16_t read_count = 0;

write:
  while (true) {
#pragma HLS PIPELINE
#pragma HLS DEPENDENCE variable=buf_r false
#pragma HLS DEPENDENCE variable=buf_g false
#pragma HLS DEPENDENCE variable=buf_b false

    if (read)
    {
      uint32_t buffer_y = read_y & (MAX_BUFFER_HEIGHT - 1);
      uint32_t index = buffer_y << buf.buffer_w_clog2;
      int addr_a = index / (buffer_memory::num_urams);

      pixel_block b;
read:
      for (int k=0; k<buffer_memory::num_urams; k++) {
#pragma HLS UNROLL
        b[k].r = buf_r[addr_a + read_x][k];
        b[k].g = buf_g[addr_a + read_x][k];
        b[k].b = buf_b[addr_a + read_x][k];
        b[k].a = 255;
      }
      blk_o.Push(b);

      if (read_x == p.render_w / buffer_memory::num_urams - 1) {
        read = false;
        uint16_t next_y = read_y + MAX_BUFFER_HEIGHT;
        if (next_y < p.render_h)
          y_loop_o.Push(next_y);
        if (read_count == p.render_h - 1) break;
        read_count++;
      }
      read_x++;
    }
    else if (!fifo_read_y.empty())
    {
      read = true;
      read_y = fifo_read_y.pop();
      read_x = 0;
    }

    if (!pix_i.IsEmpty())
    {
      pixel_position pos;
      pixel_info pix;
      std::tie(pos, pix) = pix_i.Pop();

      // Write to buffer
      uint32_t buffer_y = pos.y & (MAX_BUFFER_HEIGHT - 1);
      uint32_t index = (buffer_y << buf.buffer_w_clog2) | pos.x;
      int addr_a = index / (buffer_memory::num_urams);
      int addr_b = index % (buffer_memory::num_urams);

      assert(addr_a < 4096*8);

      buf_r[addr_a][addr_b] = clamp(int(pix.p[0]), 0, 255);
      buf_g[addr_a][addr_b] = clamp(int(pix.p[1]), 0, 255);
      buf_b[addr_a][addr_b] = clamp(int(pix.p[2]), 0, 255);

      // Count pixels
      auto count = pixel_count[buffer_y];

      if (count == p.render_w - 1) {
        y_o.Push(pos.y);
        fifo_read_y.push(pos.y);

        count = 0;
      } else {
        count++;
      }

      pixel_count[buffer_y] = count;
    }
  }
}

void write_mem(
  const render_info& p,
  line_stream& y_i,
  pix_blk_stream& blk_i,
  pixel_block* image
) {
  constexpr auto block_w = std::tuple_size<pixel_block>::value;
  for (uint16_t i=0; i<p.render_h; i++) {
    auto y = y_i.Pop();

#ifndef HLSLIB_SYNTHESIS
    std::cout << "Write line : " << y << std::endl;
#endif

    for (int x=0; x<p.render_w/block_w; x++) {
      image[y * p.image_w/block_w + p.start_x/block_w + x] = blk_i.Pop();
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

#define INST_STREAM(type, name) type name(#name)

  INST_STREAM(pos_stream,            stage0_pos);
  INST_STREAM(pix_stream,            stage1_pix);
  INST_STREAM(pix_ray_stream,        stage2_ray);
  INST_STREAM(ray_hit_stream,        stage3_ray);
  INST_STREAM(pix_ray_hit_stream,    stage3_pix);
  INST_STREAM(pix_ray_stream,        stage4_ray);
  INST_STREAM(pix_ray_hit_stream_fb, stage4_ray_loop);
  INST_STREAM(pix_stream,            stage5_pix);
  INST_STREAM(pix_stream_fb,         stage5_pix_loop);
  INST_STREAM(pix_blk_stream,        stage6_blk);
  INST_STREAM(line_stream,           stage6_line);
  INST_STREAM(line_stream_fb,        stage6_line_loop);

  done_stream done[3];

  HLSLIB_DATAFLOW_INIT();
  // Stage 0
  HLSLIB_DATAFLOW_FUNCTION(generate_pixel, p, stage6_line_loop, stage0_pos);
  // Stage 1
  HLSLIB_DATAFLOW_FUNCTION(merge_pixel, stage0_pos, stage5_pix_loop, done[0], stage1_pix);
  // Stage 2
  HLSLIB_DATAFLOW_FUNCTION(generate_ray, p, stage1_pix, stage2_ray);
  // Stage 3
  HLSLIB_DATAFLOW_FUNCTION(merge_hit, stage2_ray, stage4_ray_loop, done[1], stage3_ray, stage3_pix);
  // Stage 4
  HLSLIB_DATAFLOW_FUNCTION(hit_sphere, p, world, stage3_ray, stage3_pix, done[2], stage4_ray, stage4_ray_loop);
  // Stage 5
  HLSLIB_DATAFLOW_FUNCTION(light, p, stage4_ray, stage5_pix, stage5_pix_loop, done);
  // Stage 6
  HLSLIB_DATAFLOW_FUNCTION(buffer_write, p, stage5_pix, stage6_line, stage6_line_loop, stage6_blk);
  // Stage 7
  HLSLIB_DATAFLOW_FUNCTION(write_mem, p, stage6_line, stage6_blk, image);

  HLSLIB_DATAFLOW_FINALIZE();
}

