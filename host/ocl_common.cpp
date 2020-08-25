#include <algorithm>
#include <iostream>
#include <fstream>
#include "ocl_common.h"

cl::Platform find_platform()
{
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  auto it = std::find_if(platforms.begin(), platforms.end(),
      [](cl::Platform& p) { return p.getInfo<CL_PLATFORM_NAME>() == "Xilinx"; });
  if (it == platforms.end())
    throw std::runtime_error("Error: Failed to find Xilinx platform");
  return *it;
}

cl::Device find_device(cl::Platform platform)
{
  std::vector<cl::Device> devices;
  platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
  std::cerr << devices[0].getInfo<CL_DEVICE_NAME>() << std::endl;
  return devices[0];
}

cl::Device find_device()
{
  return find_device(find_platform());
}

cl::Program create_program(cl::Device device, cl::Context context, const std::string& xclbin)
{
  // Check file size
  std::ifstream bin_file(xclbin.c_str(), std::ifstream::binary);
  bin_file.seekg (0, bin_file.end);
  unsigned nb = bin_file.tellg();
  bin_file.seekg (0, bin_file.beg);

  // Load to bufer
  char* buf = new char[nb];
  bin_file.read(buf, nb);

  cl::Program::Binaries binaries { {buf, nb} };
  return cl::Program(context, {device}, binaries);
}

uint64_t get_duration_ns(const cl::Event& event) {
  uint64_t start, end;
  event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_START, &start);
  event.getProfilingInfo<uint64_t>(CL_PROFILING_COMMAND_END, &end);
  return end - start;
}
