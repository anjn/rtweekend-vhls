#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "ocl_common.h"
#include "object.hpp"
#include "render_info.hpp"

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

  const int samples_per_pixel = 2;
  const int image_width = 64;
  const int image_height = 36;

  //const int samples_per_pixel = 256;
  //const int image_width = 1920;
  //const int image_height = 1080;

  host_buffer<float> host_image(image_width * image_height * 4);
  cl::Memory device_image = make_device_buffer(context, CL_MEM_READ_WRITE, host_image);

  host_buffer<object> host_objects;
  host_buffer<material> host_materials;

  host_objects.push_back(make_sphere({ 0.0,    0.0, -1.2}, 0.5)); host_materials.push_back(make_metal     ({1.0, 1.0, 1.0}, 0.0));
  host_objects.push_back(make_sphere({ 0.0, -100.5, -1.0}, 100)); host_materials.push_back(make_lambertian({0.8, 0.8, 0.3}     ));
  //host_objects.push_back(make_sphere({-1.0,   -0.1, -1.0}, 0.4)); host_materials.push_back(make_metal     ({1.0, 0.8, 0.8}, 0.1));
  //host_objects.push_back(make_sphere({-0.5,   -0.4, -0.5}, 0.1)); host_materials.push_back(make_lambertian({0.3, 0.8, 0.3}     ));
  //host_objects.push_back(make_sphere({ 1.0,   -0.2, -1.0}, 0.3)); host_materials.push_back(make_lambertian({0.3, 0.3, 0.8}     ));
  //host_objects.push_back(make_sphere({ 0.4,   -0.4, -0.7}, 0.1)); host_materials.push_back(make_lambertian({0.9, 0.8, 0.1}     ));
  //host_objects.push_back(make_sphere({ 0.1,   -0.4, -0.6}, 0.1)); host_materials.push_back(make_lambertian({0.9, 0.1, 0.1}     ));

  cl::Memory device_objects = make_device_buffer(context, CL_MEM_READ_ONLY, host_objects);
  cl::Memory device_materials = make_device_buffer(context, CL_MEM_READ_ONLY, host_materials);

  render_info p;
  p.image_w = image_width;
  p.image_h = image_height;
  p.start_x = 0;
  p.start_y = 0;
  p.end_x = image_width;
  p.end_y = image_height;
  p.samples_per_pixel = samples_per_pixel;
  p.output_ratio = 1.0f;
  p.num_objects = host_objects.size();

  p.pre_calculation();

  host_buffer<uint32_t> host_render;
  for (auto& v: to_array<uint32_t>(p)) host_render.push_back(v);
  cl::Memory device_render = make_device_buffer(context, CL_MEM_READ_ONLY, host_render);

  std::vector<cl::Event> events;

  auto render = [&](int sx, int sy, int ex, int ey, float ratio, std::vector<cl::Event> e_depends) {
    int arg = 0;
    rt_kernel.setArg(arg++, device_render);
    rt_kernel.setArg(arg++, device_objects);
    rt_kernel.setArg(arg++, device_materials);
    rt_kernel.setArg(arg++, device_image);

    cl::Event e;
    q.enqueueTask(rt_kernel, &e_depends, &e);
    events.push_back(e);
  };

  cl::Event e_mem_wr;
  q.enqueueMigrateMemObjects({device_render, device_objects, device_materials, device_image}, 0, nullptr, &e_mem_wr);

  for (int i=0; i<1; i++) {
    int sx = image_width / 1 * (i + 0);
    int ex = image_width / 1 * (i + 1);
    render(sx, 0, ex, image_height, 1.0f, {e_mem_wr});
  }

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
