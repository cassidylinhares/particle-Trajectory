#pragma once

#include <cstdio>

template <typename T>
struct VecR3 {
  T x{0};
  T y{0};
  T z{0};
};

template <typename T>
auto operator+(VecR3<T> a, VecR3<T> b) {
  return VecR3<T>{a.x + b.x, a.y + b.y, a.z + b.z};
}

template <typename T>
auto operator*(VecR3<T> a, T c) {
  return VecR3<T>{c * a.x, c * a.y, c*a.z};
}

template <typename T>
auto operator*(T c, VecR3<T> a) {
  return a * c;
}

template <typename T>
auto operator/(VecR3<T> a, T c) {
  return VecR3<T>{a.x / c, a.y / c, a.z/c};
}

struct TState {
  double t{0};
  VecR3<double> position{0, 0, 0};
  VecR3<double> velocity{0, 0, 0};
};

void print_tstate(TState s) {
  printf("%f %f %f %f %f %f %f\n", s.t, s.position.x, s.position.y, s.position.z, s.velocity.x,
         s.velocity.y, s.velocity.z);
}
