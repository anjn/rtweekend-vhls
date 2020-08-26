#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "ocl_common.h"

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

  const int samples_per_pixel = 8;
  const int image_width = 192;
  const int image_height = 108;

  host_buffer<uint32_t> host_image(image_width * image_height);
  cl::Memory device_image = make_device_buffer(context, CL_MEM_WRITE_ONLY, host_image);

  std::vector<cl::Event> events;

  auto render = [&](int sx, int sy, int ex, int ey) {
    int arg = 0;
    rt_kernel.setArg(arg++, image_width);
    rt_kernel.setArg(arg++, image_height);
    rt_kernel.setArg(arg++, sx);
    rt_kernel.setArg(arg++, sy);
    rt_kernel.setArg(arg++, ex);
    rt_kernel.setArg(arg++, ey);
    rt_kernel.setArg(arg++, samples_per_pixel);
    rt_kernel.setArg(arg++, device_image);

    cl::Event e;
    q.enqueueTask(rt_kernel, nullptr, &e);
    events.push_back(e);
  };

  //render(0, 0, image_width/2, image_height/2);
  //render(image_width/2, image_height/2, image_width, image_height);
  render(0, 0, image_width, image_height);

  q.enqueueMigrateMemObjects({device_image}, CL_MIGRATE_MEM_OBJECT_HOST, &events);
  q.finish();

  // Output in ppm format
  std::ofstream ofs("out.ppm");
  ofs << "P3\n" << image_width << ' ' << image_height << "\n255\n";
  for (int y=image_height-1; y>=0; y--) {
    for (int x=0; x<image_width; x++) {
      auto rgb = host_image[image_width * y + x];
      int ir = (rgb >>  0) & 0xff;
      int ig = (rgb >>  8) & 0xff;
      int ib = (rgb >> 16) & 0xff;
      ofs << ir << ' ' << ig << ' ' << ib << '\n';
    }
  }

  device_image = cl::Memory();
  rt_kernel = cl::Kernel();
  q = cl::CommandQueue();
  program = cl::Program();
  context = cl::Context();
  device = cl::Device();
}
