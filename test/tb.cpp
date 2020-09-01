#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "rt.hpp"

int main(int argc, char** argv)
{
  const int samples_per_pixel = 8;
  const int image_width = 480/2;
  const int image_height = 270/2;

  std::vector<pixel<pixel_type>> image(image_width * image_height);
  auto image_ptr = reinterpret_cast<pixel_block*>(image.data());

  for (auto& p: image) {
    p.r = 0;
    p.g = 0;
    p.b = 0;
  }

  std::vector<object> objects;
  std::vector<material> materials;
  objects.push_back(make_sphere({ 0,      0, -1}, 0.5)); materials.push_back(make_lambertian({0.8, 0.3, 0.3}));
  objects.push_back(make_sphere({ 0, -100.5, -1}, 100)); materials.push_back(make_lambertian({0.8, 0.8, 0.3}));
  objects.push_back(make_sphere({-1,      0, -1}, 0.4)); materials.push_back(make_metal({0.8, 0.8, 0.8}, 0.0));

  render_info p;
  p.image_w = image_width;
  p.image_h = image_height;
  p.start_x = 0;
  p.start_y = 0;
  p.end_x = image_width;
  p.end_y = image_height;
  p.samples_per_pixel = samples_per_pixel;
  p.output_ratio = 1.0f;
  p.num_objects = objects.size();

  p.pre_calculation();

  auto p_arr = to_array<uint32_t>(p);

  rt(p_arr.data(), objects.data(), materials.data(), image_ptr);

  // Output in ppm format
  std::ofstream ofs("../../../../out.ppm");
  ofs << "P3\n" << image_width << ' ' << image_height << "\n255\n";
  for (int y=image_height-1; y>=0; y--) {
    for (int x=0; x<image_width; x++) {
      auto p = image[image_width * y + x];
      p.r = std::clamp(std::sqrt(p.r / 256) * 256, 0.0f, 255.0f);
      p.g = std::clamp(std::sqrt(p.g / 256) * 256, 0.0f, 255.0f);
      p.b = std::clamp(std::sqrt(p.b / 256) * 256, 0.0f, 255.0f);
      ofs << int(p.r) << ' ' << int(p.g) << ' ' << int(p.b) << '\n';
    }
  }
}
