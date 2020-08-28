#pragma once
#include <array>
#include <stdint.h>

template<typename T>
struct pixel
{
  T b;
  T g;
  T r;
  T a;
};

struct object
{
  float center_x;
  float center_y;
  float center_z;
  float radius;
};

constexpr int pixel_block_size = 4;
using pixel_type = float;
using pixel_block = std::array<pixel<pixel_type>, pixel_block_size>;

void rt(
  const int image_w,
  const int image_h,
  const int start_x,
  const int start_y,
  const int end_x,
  const int end_y,
  const int samples_per_pixel,
  const float output_ratio,
  const int num_objects,
  object* objects,
  pixel_block* image
);

