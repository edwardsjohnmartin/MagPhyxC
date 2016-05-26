#ifndef __EVENT_H__
#define __EVENT_H__

#include "./Dipole.h"

class Event {
 public:
  Event(const std::string& filename, const Dipole& d) : _n(1), _d(d) {
    _file = fopen(filename.c_str(), "w");
    fprintf(_file, "n, event_type, t, r, theta, phi,"
            "pr, ptheta, pphi, beta, E, dE\n");
  }

  ~Event() {
    fclose(_file);
  }

  int get_n() const { return _n; }

  bool log(const Dipole& new_d, const double t) {
    bool fired = false;
    // Log zero crossings
    if (isZeroCrossing(_d.get_theta(), new_d.get_theta())) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return d.get_theta();});
      event("theta = 0", logDipole, t);
      fired = true;
    }
    if (isZeroCrossing(_d.get_phi(), new_d.get_phi())) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return d.get_phi();});
      event("phi = 0", logDipole, t);
      fired = true;
    }
    if (isZeroCrossing(Physics::get_beta(_d), Physics::get_beta(new_d))) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return Physics::get_beta(d);});
      event("beta = 0", logDipole, t);
      fired = true;
    }
    if (isNegativeZeroCrossing(_d.get_pr(), new_d.get_pr())) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return d.get_pr();});
      event("pr = 0", logDipole, t);
      fired = true;
    }
    if (isZeroCrossing(_d.get_ptheta(), new_d.get_ptheta())) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return d.get_ptheta();});
      event("ptheta = 0", logDipole, t);
      fired = true;
    }
    if (isZeroCrossing(_d.get_pphi(), new_d.get_pphi())) {
      const Dipole logDipole = Dipole::interpolateZeroCrossing(
          _d, new_d, [](const Dipole& d) {return d.get_pphi();});
      event("pphi = 0", logDipole, t);
      fired = true;
    }
    _d = new_d;
    return fired;
  }

  void logCollision(const Dipole& new_d, const double t) {
    event("collision", new_d, t);
    _d = new_d;
  }

 private:
  void event(const std::string& name, const Dipole& d, const double t) {
    fprintf(_file, "%d,%s,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%.2e\n",
            _n, name.c_str(), t, d.get_r(),
            Physics::rad2deg(d.get_theta()), Physics::rad2deg(d.get_phi()),
            d.get_pr(), d.get_ptheta(), d.get_pphi(), Physics::get_beta(d),
            d.get_E(), d.get_dE());
    _n++;
  }

  static int sign(const double d) {
    const double EPSILON = 0.000000000001;
    return (d < -EPSILON) ? -1 : (d > EPSILON) ? 1 : 0;
  }

  static bool isZeroCrossing(const double a, const double b) {
    // Return true if we cross zero or if we move from non-zero to zero.
    if (sign(a) == 0) return false;
    return (sign(a) == -sign(b)) || (sign(b) == 0);
  }

  static bool isNegativeZeroCrossing(const double a, const double b) {
    // Return true if we cross zero from positive to negative or from
    // zero to negative.
    return (sign(a) > 0 && sign(b) <= 0);
  }

 private:
  FILE* _file;
  int _n;
  Dipole _d;
};

#endif
