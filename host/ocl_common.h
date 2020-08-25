#include <string>
#include <vector>

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl2.hpp>
#include <CL/cl_ext.h>
#include <CL/cl_ext_xilinx.h>

template <typename T>
struct aligned_allocator
{
  using value_type = T;
  T* allocate(std::size_t num)
  {
    void* ptr = nullptr;
    if (posix_memalign(&ptr,4096,num*sizeof(T)))
      throw std::bad_alloc();
    return reinterpret_cast<T*>(ptr);
  }
  void deallocate(T* p, std::size_t num)
  {
    free(p);
  }
};


template<typename T>
using host_buffer = std::vector<T, aligned_allocator<T>>;

template<typename T, typename A>
cl::Buffer make_device_buffer(const cl::Context& context, cl_mem_flags flags, const std::vector<T,A>& host_buffer, int ddr_bank = -1)
{
  size_t size = sizeof(T) * host_buffer.size();
  cl_int err;
  cl::Buffer buffer;
  if (ddr_bank >= 0) {
    cl_mem_ext_ptr_t ext;
    ext.flags = XCL_MEM_TOPOLOGY | ddr_bank;
    ext.obj = (void*) host_buffer.data();
    ext.param = nullptr;
    buffer = cl::Buffer(context, flags | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR, size, &ext, &err);
  } else {
    buffer = cl::Buffer(context, flags | CL_MEM_USE_HOST_PTR, size, (void*) host_buffer.data(), &err);
  }
  if (err != CL_SUCCESS) {
    printf("Failed to allocate device buffer!\n");
    throw std::bad_alloc();
  }
  return buffer;
}

cl::Platform find_platform();
cl::Device find_device(cl::Platform platform);
cl::Device find_device();
cl::Program create_program(cl::Device device, cl::Context context, const std::string& xclbin);

uint64_t get_duration_ns(const cl::Event& event);
