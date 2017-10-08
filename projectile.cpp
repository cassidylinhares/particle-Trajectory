#include "projectile.hpp"
#include <cmath>

const double m = 1;
const double q = 1;
const double Eo = 1;
const double w = 2;
const double _k = 0.5;
const double dt = 0.001;

auto force(TState s) { 
	double x = cos(_k*s.position.z - w*s.t);
	double y = sin(_k*s.position.z - w*s.t);
	double z = 0;
	return q * Eo/sqrt(2) * VecR3<double>{x, y, z}; 
}

auto verlet_step(TState s, VecR3<double> accel) {
  TState next;
  next.t = s.t + dt;
  next.position = s.position + (s.velocity * dt) + (accel * dt*dt* 0.5);
  next.velocity = s.velocity + (accel + (force(next) / m))*dt*0.5;
  return next;
}

void n_steps(unsigned n, TState state0) {
  print_tstate(state0);
  if (n == 0)
    return;
  else {
    auto state = state0;
    for (unsigned k = 0; k < n; ++k) {
      state = verlet_step(state, force(state) / m);
      print_tstate(state);
    }
  }
}

int main() {
  n_steps(7000, TState{0.0, {-Eo/(w*w*sqrt(2)), 0,0}, {0, Eo/(w*sqrt(2)),-0.25}});
  return 0;
}
