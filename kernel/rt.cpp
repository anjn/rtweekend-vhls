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
#include "shifter.hpp"

constexpr int MAX_HIT = 16;

constexpr int MEM_DATA_WIDTH = 32 * std::tuple_size<pixel_block>::value; // 512
constexpr int BUFFER_DATA_WIDTH = MEM_DATA_WIDTH / 4; // 128
constexpr int BUFFER_DATA_WORD = BUFFER_DATA_WIDTH / 8; // 16
constexpr int BUFFER_MIN_URAM_PER_COLOR = std::max(BUFFER_DATA_WIDTH / 64, 1); // 2

constexpr int BUFFER_WIDTH = 4096;
constexpr int BUFFER_WIDTH_CLOG2 = hlslib::ConstLog2(BUFFER_WIDTH-1);
constexpr int BUFFER_HEIGHT = std::max(16, 64/8 * BUFFER_MIN_URAM_PER_COLOR * 4096 / BUFFER_WIDTH); // 16
constexpr int BUFFER_SIZE = BUFFER_WIDTH * BUFFER_HEIGHT;

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

// Stream types
using pos_stream            = hlslib::Stream<pixel_position, 8>;
using pix_stream            = hlslib::Stream<std::tuple<pixel_position, pixel_info>, 8>;
using pix_stream_fb         = hlslib::Stream<std::tuple<pixel_position, pixel_info>, 1024>;
using pix_ray_stream        = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info>, 8>;
using pix_ray_hit_stream    = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info, hit_info>, 8>;
using pix_ray_hit_stream_fb = hlslib::Stream<std::tuple<pixel_position, pixel_info, ray_info, hit_info>, 512>;
using ray_hit_stream        = hlslib::Stream<std::tuple<ray_info, hit_info>, 8>;
using pix_blk_stream        = hlslib::Stream<pixel_block, 128>;
using line_stream           = hlslib::Stream<uint16_t, 8>;
using line_stream_fb        = hlslib::Stream<uint16_t, BUFFER_HEIGHT>;
using wr_req_stream         = hlslib::Stream<std::tuple<uint32_t, uint8_t, uint8_t, uint8_t>, 8>;
using rd_req_stream         = hlslib::Stream<uint32_t, 8>;
using done_stream           = hlslib::Stream<bool>;

void generate_pixel(
  const render_info& p,
  line_stream_fb& y_i,
  pos_stream& pos_o
) {
scanlines:
  for (uint16_t y=0; y<p.render_h; y++) {
    uint16_t line = y;
    if (line >= BUFFER_HEIGHT) {
      line = y_i.Pop();
    }

pixels:
    for (uint16_t x=0; x<p.render_w; x++) {
#pragma HLS PIPELINE
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
  done_stream& done_i,
  pix_ray_stream& ray_o
) {
  static xorshift32 rx(2463534242);
  static xorshift32 ry(1463534242);

  const float a = 256.0f / p.samples_per_pixel;

  camera cam;

  while (true) {
#pragma HLS PIPELINE
    if (!pix_i.IsEmpty())
    {
      auto [pos, pix] = pix_i.Pop();

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

    bool done;
    if (done_i.ReadNonBlocking(done)) {
      break;
    }
  }
}

void merge_hit(
  pix_ray_stream& ray_i,
  pix_ray_hit_stream_fb& ray_loop_i,
  done_stream& done_i,
  ray_hit_stream& ray_o,
  pix_ray_hit_stream_fb& pix_o
) {
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
  pix_ray_hit_stream_fb& pix_i,
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
        auto [ri, hi] = ray_i.Pop();
        hit = world.objects[hi.object_index].hit(ri.r, 0.001, hi.closest_so_far, rec);
      }

      auto [pos, pix, ri, hi] = pix_i.Pop();

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
  done_stream& done_i,
  pix_stream& pix_o,
  pix_stream_fb& pix_loop_o
) {
  while (true) {
#pragma HLS PIPELINE
    if (!ray_i.IsEmpty())
    {
      auto [pos, pix, ri] = ray_i.Pop();

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

    bool done;
    if (done_i.ReadNonBlocking(done)) {
      break;
    }
  }
}

void buffer_write(
  const render_info& p,
  pix_stream& pix_i,
  done_stream& done_i,
  wr_req_stream& wr_o,
  line_stream& y_o
) {

write:
  while (true) {
#pragma HLS PIPELINE
    if (!pix_i.IsEmpty())
    {
      auto [pos, pix] = pix_i.Pop();

      // Write to buffer
      uint32_t buffer_y = pos.y & (BUFFER_HEIGHT - 1);
      uint32_t index = (buffer_y << BUFFER_WIDTH_CLOG2) | pos.x;

      assert(index < BUFFER_SIZE);
      uint8_t w[3];
      w[0] = clamp(int(pix.p[0]), 0, 255);
      w[1] = clamp(int(pix.p[1]), 0, 255);
      w[2] = clamp(int(pix.p[2]), 0, 255);
      wr_o.Push({index, w[0], w[1], w[2]});

      y_o.Push(pos.y);
    }

    bool done;
    if (done_i.ReadNonBlocking(done)) {
      break;
    }
  }
}

void count_pixel(
  const render_info& p,
  line_stream& y_i,
  done_stream& done_i,
  line_stream& y_o
) {
  const int latency = 4;
  using y_shifter = shifter<uint16_t, latency+1>;

  y_shifter ys;

  uint16_t pixel_count[BUFFER_HEIGHT];

clear_counter:
  for (uint16_t y=0; y<BUFFER_HEIGHT; y++) {
#pragma HLS PIPELINE
    pixel_count[y] = 0;
  }

  while (true) {
#pragma HLS PIPELINE
#pragma HLS LATENCY max=latency
#pragma HLS DEPENDENCE variable=pixel_count inter false

    y_shifter::item y;

    if (!y_i.IsEmpty()) {
      y = ys.shift(y_i.Pop());
    } else {
      y = ys.shift();
    }

    if (y.valid) {
      int num = 1;
      for (int i=0; i<latency; i++) {
#pragma HLS UNROLL
        auto tmp_y = ys.arr[i];
        if (tmp_y.valid && tmp_y.data == y.data) {
          num++;
          ys.arr[i].valid = false;
        }
      }

      uint32_t buffer_y = y.data & (BUFFER_HEIGHT - 1);
      uint16_t count = pixel_count[buffer_y] + num;

      assert(count <= p.render_w);

      if (count == p.render_w) {
        y_o.Push(y.data);
        count = 0;
      }

      pixel_count[buffer_y] = count;
    }

    bool done;
    if (done_i.ReadNonBlocking(done)) {
      break;
    }
  }
}

void buffer(
  wr_req_stream& wr_i,
  rd_req_stream& rd_req_i,
  done_stream& done_i,
  pix_blk_stream& blk_o
) {
  uint8_t buf_r[BUFFER_SIZE];
  uint8_t buf_g[BUFFER_SIZE];
  uint8_t buf_b[BUFFER_SIZE];

#pragma HLS BIND_STORAGE variable=buf_r type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_g type=RAM_S2P impl=uram
#pragma HLS BIND_STORAGE variable=buf_b type=RAM_S2P impl=uram
#pragma HLS ARRAY_RESHAPE variable=buf_r cyclic factor=BUFFER_DATA_WORD dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_g cyclic factor=BUFFER_DATA_WORD dim=1
#pragma HLS ARRAY_RESHAPE variable=buf_b cyclic factor=BUFFER_DATA_WORD dim=1

  while (true) {
#pragma HLS PIPELINE
    if (!rd_req_i.IsEmpty())
    {
      auto addr = rd_req_i.Pop();

      pixel_block b;
      for (int k=0; k<BUFFER_DATA_WORD; k++) {
#pragma HLS UNROLL
        b[k].r = buf_r[(addr << 4) | k];
        b[k].g = buf_g[(addr << 4) | k];
        b[k].b = buf_b[(addr << 4) | k];
        b[k].a = 255;
      }
      blk_o.Push(b);
    }

    if (!wr_i.IsEmpty())
    {
      uint32_t addr;
      uint8_t p[3];
      std::tie(addr, p[0], p[1], p[2]) = wr_i.Pop();
      buf_r[addr] = p[0];
      buf_g[addr] = p[1];
      buf_b[addr] = p[2];
    }

    bool done;
    if (done_i.ReadNonBlocking(done)) {
      break;
    }
  }
}

void buffer_read(
  const render_info& p,
  line_stream& y_i,
  rd_req_stream& rd_req_o,
  line_stream& y_o,
  line_stream_fb& y_loop_o
) {
read:
  for (int i=0; i<p.render_h; i++) {
    auto y = y_i.Pop();
    auto buffer_y = y & (BUFFER_HEIGHT - 1);
    auto index = buffer_y << BUFFER_WIDTH_CLOG2;

    y_o.Push(y);
    
    uint16_t next_y = y + BUFFER_HEIGHT;
    if (next_y < p.render_h)
      y_loop_o.Push(next_y);

    for (int j=0; j<p.render_w/BUFFER_DATA_WORD; j++) {
#pragma HLS PIPELINE
      uint32_t addr =  index + BUFFER_DATA_WORD * j;
      rd_req_o.Push(addr >> 4);
    }
  }
}

void write_mem(
  const render_info& p,
  line_stream& y_i,
  pix_blk_stream& blk_i,
  done_stream done_o[5],
  pixel_block* image
) {
  constexpr auto block_w = std::tuple_size<pixel_block>::value;
  for (uint16_t i=0; i<p.render_h; i++) {
    auto y = y_i.Pop();

#ifndef HLSLIB_SYNTHESIS
    std::cout << "Write line : " << y << std::endl;
#endif

    for (int x=0; x<p.render_w/block_w; x++) {
#pragma HLS PIPELINE
      image[y * p.image_w/block_w + p.start_x/block_w + x] = blk_i.Pop();
    }
  }

  for (int i=0; i<8; i++) {
#pragma HLS UNROLL
    done_o[i].Push(true);
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
  INST_STREAM(pix_ray_hit_stream_fb, stage3_pix);
  INST_STREAM(pix_ray_stream,        stage4_ray);
  INST_STREAM(pix_ray_hit_stream_fb, stage4_ray_loop);
  INST_STREAM(pix_stream,            stage5_pix);
  INST_STREAM(pix_stream_fb,         stage5_pix_loop);
  INST_STREAM(wr_req_stream,         stage6_wr);
  INST_STREAM(line_stream,           stage6_line);
  INST_STREAM(line_stream,           stage7_line);
  INST_STREAM(pix_blk_stream,        stage8_blk);
  INST_STREAM(rd_req_stream,         stage9_rd);
  INST_STREAM(line_stream,           stage9_line);
  INST_STREAM(line_stream_fb,        stage9_line_loop);

  done_stream done[8];

  HLSLIB_DATAFLOW_INIT();
  // Stage 0
  HLSLIB_DATAFLOW_FUNCTION(generate_pixel, p, stage9_line_loop, stage0_pos);
  // Stage 1
  HLSLIB_DATAFLOW_FUNCTION(merge_pixel, stage0_pos, stage5_pix_loop, done[0], stage1_pix);
  // Stage 2
  HLSLIB_DATAFLOW_FUNCTION(generate_ray, p, stage1_pix, done[1], stage2_ray);
  // Stage 3
  HLSLIB_DATAFLOW_FUNCTION(merge_hit, stage2_ray, stage4_ray_loop, done[2], stage3_ray, stage3_pix);
  // Stage 4
  HLSLIB_DATAFLOW_FUNCTION(hit_sphere, p, world, stage3_ray, stage3_pix, done[3], stage4_ray, stage4_ray_loop);
  // Stage 5
  HLSLIB_DATAFLOW_FUNCTION(light, p, stage4_ray, done[4], stage5_pix, stage5_pix_loop);
  // Stage 6
  HLSLIB_DATAFLOW_FUNCTION(buffer_write, p, stage5_pix, done[5], stage6_wr, stage6_line);
  // Stage 7
  HLSLIB_DATAFLOW_FUNCTION(count_pixel, p, stage6_line, done[7], stage7_line);
  // Stage 8
  HLSLIB_DATAFLOW_FUNCTION(buffer, stage6_wr, stage9_rd, done[6], stage8_blk);
  // Stage 9
  HLSLIB_DATAFLOW_FUNCTION(buffer_read, p, stage7_line, stage9_rd, stage9_line, stage9_line_loop);
  // Stage 10
  HLSLIB_DATAFLOW_FUNCTION(write_mem, p, stage9_line, stage8_blk, done, image);

  HLSLIB_DATAFLOW_FINALIZE();
}

