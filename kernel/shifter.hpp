#pragma once
#include <array>

template<typename T, int N>
struct shifter
{
  struct item {
    T data;
    bool valid;
  };

  std::array<item, N> arr;

  shifter() {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=arr complete
    for (int i=0; i<N; i++) {
      arr[i] = {T(), false};
    }
  }

  item shift(const T& v) {
    auto s = arr[0];
    for (int i=0; i<N-1; i++) {
      arr[i] = arr[i+1];
    }
    arr[N-1] = {v, true};
    return s;
  }

  item shift() {
    auto s = arr[0];
    for (int i=0; i<N-1; i++) {
      arr[i] = arr[i+1];
    }
    arr[N-1] = {T(), false};
    return s;
  }
};

