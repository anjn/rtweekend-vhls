#pragma once
#include <array>
#include <stdint.h>

#include "ap_int.h"

constexpr auto aspect_ratio = 16.0 / 9.0;
constexpr int SIM_WIDTH  = 384;
constexpr int SIM_HEIGHT = static_cast<int>(SIM_WIDTH / aspect_ratio);

struct pixel
{
  uint8_t b;
  uint8_t g;
  uint8_t r;
  uint8_t a;
};

constexpr int output_block_width = 16;

using pixel_block = std::array<pixel, output_block_width>;

void rt(
  const int image_w,
  const int image_h,
  const int start_x,
  const int start_y,
  const int end_x,
  const int end_y,
  const int samples_per_pixel,
  pixel_block* image
);

