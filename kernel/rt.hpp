#pragma once
#include <array>
#include <stdint.h>

#include "render_info.hpp"
#include "object.hpp"

template<typename T>
struct pixel
{
  T b;
  T g;
  T r;
  T a;
};

constexpr int pixel_block_size = 4;
using pixel_type = float;
using pixel_block = std::array<pixel<pixel_type>, pixel_block_size>;

void rt(
  uint32_t* render,
  object* objects,
  material* materials,
  pixel_block* image
);

