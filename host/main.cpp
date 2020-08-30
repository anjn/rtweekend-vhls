#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "ocl_common.h"
#include "object.hpp"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
    return EXIT_FAILURE;
  }

  auto xclbin = argv[1];

  auto device = find_device();
  auto context = cl::Context(device);
  auto program = create_program(device, context, xclbin);
  auto q = cl::CommandQueue(context, device, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE);
  //auto q = cl::CommandQueue(context, device);

  auto rt_kernel = cl::Kernel(program, "rt");

  //const int samples_per_pixel = 2;
  //const int image_width = 64;
  //const int image_height = 36;

  const int samples_per_pixel = 16*16;
  const int image_width = 1920;
  const int image_height = 1080;

  host_buffer<float> host_image(image_width * image_height * 4);
  cl::Memory device_image = make_device_buffer(context, CL_MEM_READ_WRITE, host_image);

  host_buffer<object> host_objects;

  host_objects.push_back(make_sphere({ 0,      0, -1  }, 0.5, make_metal({1.0, 1.0, 1.0}, 0.0)));
  host_objects.push_back(make_sphere({ 0, -100.5, -1  }, 100, make_lambertian({0.8, 0.8, 0.3})));
  host_objects.push_back(make_sphere({-1,      0, -1  }, 0.4, make_metal({1.0, 0.8, 0.8}, 0.1)));
  host_objects.push_back(make_sphere({-0.5,   -0.4, -0.5}, 0.1, make_lambertian({0.3, 0.8, 0.3})));
  host_objects.push_back(make_sphere({ 1,      0, -1  }, 0.3, make_lambertian({0.3, 0.3, 0.8})));
  host_objects.push_back(make_sphere({ 0.4,   -0.4, -0.7  }, 0.1, make_lambertian({0.9, 0.8, 0.1})));

  cl::Memory device_objects = make_device_buffer(context, CL_MEM_READ_ONLY, host_objects);

  std::vector<cl::Event> events;

  auto render = [&](int sx, int sy, int ex, int ey, float ratio, std::vector<cl::Event> e_depends) {
    int arg = 0;
    rt_kernel.setArg(arg++, image_width);
    rt_kernel.setArg(arg++, image_height);
    rt_kernel.setArg(arg++, sx);
    rt_kernel.setArg(arg++, sy);
    rt_kernel.setArg(arg++, ex);
    rt_kernel.setArg(arg++, ey);
    rt_kernel.setArg(arg++, samples_per_pixel);
    rt_kernel.setArg(arg++, ratio);
    rt_kernel.setArg(arg++, int(host_objects.size()));
    rt_kernel.setArg(arg++, device_objects);
    rt_kernel.setArg(arg++, device_image);

    cl::Event e;
    q.enqueueTask(rt_kernel, &e_depends, &e);
    events.push_back(e);
  };

  cl::Event e_mem_wr;
  q.enqueueMigrateMemObjects({device_image, device_objects}, 0, nullptr, &e_mem_wr);

  for (int i=0; i<1; i++) {
    int sx = image_width / 1 * (i + 0);
    int ex = image_width / 1 * (i + 1);
    render(sx, 0, ex, image_height, 1.0f, {e_mem_wr});
  }
  //render(0, 0, image_width, image_height, 1.0f, {e_mem_wr});
  //render(0, 0, image_width, image_height, 0.5f, {e_mem_wr});

  q.enqueueMigrateMemObjects({device_image}, CL_MIGRATE_MEM_OBJECT_HOST, &events);
  q.finish();

  // Output in ppm format
  std::ofstream ofs("out.ppm");
  ofs << "P3\n" << image_width << ' ' << image_height << "\n255\n";
  for (int y=image_height-1; y>=0; y--) {
    for (int x=0; x<image_width; x++) {
      auto rgb = &host_image[(image_width * y + x) * 4];
      int ir = std::clamp(int(rgb[2]), 0, 255);
      int ig = std::clamp(int(rgb[1]), 0, 255);
      int ib = std::clamp(int(rgb[0]), 0, 255);
      ofs << ir << ' ' << ig << ' ' << ib << '\n';
    }
  }

  //device_image = cl::Memory();
  //rt_kernel = cl::Kernel();
  //q = cl::CommandQueue();
  //program = cl::Program();
  //context = cl::Context();
  //device = cl::Device();
}
