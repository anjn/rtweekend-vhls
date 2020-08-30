#pragma once

struct xorshift32
{
  static constexpr uint32_t min = 1;
  static constexpr uint32_t max = ~0;
  static constexpr float max_float = max;
  static constexpr float max_float_inv = 1.0f / max;

  uint32_t y;

  xorshift32(): y(2463534242) {}

  xorshift32(uint32_t seed): y(seed) {}

  uint32_t next() {
#pragma HLS INLINE
    y = y ^ (y << 13);
    y = y ^ (y >> 17);
    y = y ^ (y << 5);
    return y;
  }
};

template<typename T>
struct random_norm
{
  template<typename Random>
  static T gen(Random& r) {
    return T(r.next() - Random::min) / (Random::max - Random::min + 1);
  }

  template<typename Random>
  static T gen(Random& r, T min, T max) {
    return min + (max - min) * gen(r);
  }
};

template<int W, int I>
struct random_norm<ap_fixed<W,I>>
{
  template<typename Random>
  static ap_fixed<W,I> gen(Random& r) {
    ap_fixed<W,I> v = 0;
    v.range(W-I-1, 0) = r.next();
    return v;
  }
};
