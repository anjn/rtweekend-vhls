#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "rt.hpp"

int main(int argc, char** argv)
{
  const int samples_per_pixel = 8*4;
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
  objects.push_back(make_sphere({ 0,      0, -1}, 0.5, make_lambertian({0.8, 0.3, 0.3})));
  objects.push_back(make_sphere({ 0, -100.5, -1}, 100, make_lambertian({0.8, 0.8, 0.3})));
  objects.push_back(make_sphere({-1,      0, -1}, 0.4, make_metal({0.8, 0.8, 0.8}, 0.0)));
  //objects.push_back(make_sphere({ 1,      0, -1}, 0.4, make_dielectric(1.5)));

  rt(image_width, image_height, 0, 0, image_width, image_height, samples_per_pixel, 1.0f, objects.size(), objects.data(), image_ptr);
  //rt(image_width, image_height, 0, 0, image_width, image_height, samples_per_pixel, 0.5f, image_ptr);

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
