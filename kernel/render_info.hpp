#pragma once
#include <cstdint>
#include <vector>

struct render_info
{
  uint16_t image_w;
  uint16_t image_h;
  uint16_t start_x;
  uint16_t start_y;
  uint16_t end_x;
  uint16_t end_y;
  uint16_t samples_per_pixel;
  float    output_ratio;
  uint16_t num_objects;

  // pre calc
  uint16_t render_w;
  uint16_t render_h;
  uint32_t num_pixels;
  uint32_t num_rays;

  void pre_calculation() {
    render_w = end_x - start_x;
    render_h = end_y - start_y;
    num_pixels = render_w * render_h;
    num_rays = num_pixels * samples_per_pixel;
  }

  static const int num_elements = 13;

  template<typename Archiver>
  void archive(Archiver a) {
#pragma HLS INLINE
    a( 0, image_w);
    a( 1, image_h);
    a( 2, start_x);
    a( 3, start_y);
    a( 4, end_x);
    a( 5, end_y);
    a( 6, samples_per_pixel);
    a( 7, output_ratio);
    a( 8, num_objects);
    a( 9, render_w);
    a(10, render_h);
    a(11, num_pixels);
    a(12, num_rays);
  }
};

template<typename T>
struct serializer
{
  T* arr;

  serializer(T* a) : arr(a) {}

  template<typename U, std::enable_if_t<std::is_integral<U>::value, std::nullptr_t> = nullptr>
  void operator()(int a, const U& v) {
    arr[a] = v;
  }

  template<typename U, std::enable_if_t<std::is_floating_point<U>::value, std::nullptr_t> = nullptr>
  void operator()(int a, const U& v) {
    union { U f; T i; } conv;
    conv.f = v;
    arr[a] = conv.i;
  }
};

template<typename T>
struct deserializer
{
  const T* arr;

  deserializer(const T* a) : arr(a) {
#pragma HLS INLINE
  }

  template<typename U, std::enable_if_t<std::is_integral<U>::value, std::nullptr_t> = nullptr>
  void operator()(int a, U& v) {
#pragma HLS INLINE
    static_assert(sizeof(T) >= sizeof(U));
    v = arr[a];
  }

  template<typename U, std::enable_if_t<std::is_floating_point<U>::value, std::nullptr_t> = nullptr>
  void operator()(int a, U& v) {
#pragma HLS INLINE
    static_assert(sizeof(T) >= sizeof(U));
    union { U f; T i; } conv;
    conv.i = arr[a];
    v = conv.f;
  }
};

template<typename T, typename A, int Size = sizeof(T) / sizeof(A)>
auto from_array(const A* arr)
{
#pragma HLS INLINE
  static_assert(sizeof(T) % sizeof(A) == 0);
  T obj;
//  for (int i=0; i<Size; i++) {
//#pragma HLS PIPELINE
//    reinterpret_cast<A*>(&obj)[i] = arr[i];
//  }
  obj.archive(deserializer<A>(arr));
  return obj;
}

template<typename A, typename T>
auto to_array(T& obj)
{
  std::vector<A> arr(T::num_elements);
  static_assert(sizeof(T) % sizeof(A) == 0);
//  for (int i=0; i<Size; i++) {
//    arr[i] = reinterpret_cast<const A*>(&obj)[i];
//  }
  obj.archive(serializer<A>(arr.data()));
  return arr;
}
