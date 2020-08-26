#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "rt.hpp"

int main(int argc, char** argv)
{
  //const int samples_per_pixel = 8;
  //const int image_width = SIM_WIDTH;
  //const int image_height = SIM_HEIGHT;
  const int samples_per_pixel = 1;
  const int image_width = 480*4*2;
  const int image_height = 270*4*2;

  std::vector<pixel> image(image_width * image_height);
  auto image_ptr = reinterpret_cast<pixel_block*>(image.data());

  rt(image_width, image_height, 0, 0, image_width, image_height, samples_per_pixel, image_ptr);

  // Output in ppm format
  std::ofstream ofs("../../../../out.ppm");
  ofs << "P3\n" << image_width << ' ' << image_height << "\n255\n";
  for (int y=image_height-1; y>=0; y--) {
    for (int x=0; x<image_width; x++) {
      auto p = image[image_width * y + x];
      ofs << int(p.r) << ' ' << int(p.g) << ' ' << int(p.b) << '\n';
    }
  }
}
