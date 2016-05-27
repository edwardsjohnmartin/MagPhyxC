#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "./Dipole.h"
#include "./Physics.h"
#include "./Event.h"
#include "./Options.h"
#include "./Stepper.h"

using namespace std;

const Options::Dynamics default_dynamics = Options::BOUNCING;
const double default_n = 1e5;
const double default_h = 1e-2;
const double default_eps = 1e-10;

// Global Options object. Declared in Options.h.
Options o(default_n, default_h, default_eps, default_dynamics);

Dipole doSimulation(const Dipole& freeDipole, Event& event,
                       double h,
                       const int numEvents);

void printUsage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "SYNOPSIS\n");
  fprintf(stderr, "\tmagphyx [OPTIONS] (-i conditions | -f filename)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "DESCRIPTION\n");
  fprintf(stderr, 
          "\tmagphyx runs a magnet simulation given initial conditions\n"
          "\tspecified by either the -i or -f option. Events are output to\n"
          "\tevents.csv.\n"
          );
  fprintf(stderr, "\n");
  fprintf(stderr, "OPTIONS\n");
  fprintf(stderr, "\t-i r theta phi pr ptheta pphi\n");
  fprintf(stderr, "\t\tInitial conditions.\n");
  fprintf(stderr, "\t-f filename.csv\n");
  fprintf(stderr, 
          "\t\tInitial conditions are found in the second line of\n"
          "\t\tfilename.csv, which is in the same format as what is\n"
          "\t\texported from MagPhyx web version.\n"
          );
  fprintf(stderr, "\t-d (bouncing | rolling)\n");
  fprintf(stderr, "\t\tDynamics type. Default = bouncing.\n");
  fprintf(stderr, "\t-n numEvents\n");
  fprintf(stderr, 
          "\t\tExecutes the simulation until numEvents events occur.\n"
          "\t\tDefault = 1e5.\n");
  fprintf(stderr, "\t-h h\n");
  fprintf(stderr, "\t\tInitial step size. Default = 1e-2.\n");
  fprintf(stderr, "\t-e eps\n");
  fprintf(stderr, "\t\tError per step allowed. Note that this is error in\n"
          "\t\tterms of Runge-Kutta. The error in total energy will be\n"
          "\t\tsimilar to, but not bound by, this value. Default = 1e-10.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "EXAMPLES\n");
  fprintf(stderr, "\tmagphyx -n 1e6 -d bouncing -f init.csv\n");
  fprintf(stderr, "\tmagphyx -d rolling -i 1 3 -18.78982612 0 0 0\n");
  fprintf(stderr, "\n");
}

int main(int argc, char** argv) {
  int i = 1;
  bool stop = false;
  while (i < argc && !stop) {
    stop = true;
    if (o.ProcessArg(i, argv)) {
      stop = false;
    }
  }
  if (!o.initialized) {
    printUsage();
    return 1;
  }


  // Dipole freeDipole = initDipole(initFilename);
  Dipole freeDipole = o.dipole;
  // cout << freeDipole << endl;
  Event event("events.csv", freeDipole);
  freeDipole = doSimulation(freeDipole, event, o.h, o.numEvents);
  // cout << freeDipole << endl;
}

void printStateHeader() {
  printf("\n");
  printf(" %9s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
         "t", "h", "r", "theta", "phi", "pr", "ptheta", "pphi", "E", "dE");
  printf("-------------------------------------------"
         "-------------------------------------------"
         "--------------------------------------------\n");
}

void printState(const double t, const double h, const Dipole& d,
                const char prefix = ' ') {
  printf("%c%9.0lf %12.6g %12.5g %12.6g %12.6g %12.6g"
         "%12.6g %12.6g %12.6g %12.6g\n",
         prefix, t, h,
         d.get_r(), Physics::rad2deg(d.get_theta()),
         Physics::rad2deg(d.get_phi()),
         d.get_pr(), d.get_ptheta(), d.get_pphi(),
         d.get_E(), d.get_dE());
}

void printProgress(const int n, const Dipole& d, const bool fired) {
  if (!o.interactive && n % 1000 == 0 && fired) {
    printf("\r");
    printf("Num events = %-7d     dE = %-12e     ", n, d.get_dE());
    fflush(stdout);
  }
}

Dipole doSimulation(const Dipole& freeDipole, Event& event,
                       const double h_,
                       const int numEvents) {
  Stepper stepper(freeDipole, h_, o.eps);

  double t = 0.0;

  int n = 0;

  printf("\n");
  if (o.interactive) {
    printStateHeader();
    printState(0, h_, freeDipole);
  }
  printProgress(0, freeDipole, true);
  for (int i = 0; event.get_n() < numEvents; ++i) {
    try {
      stepper.step();
    } catch (logic_error& e) {
      printState(stepper.t, stepper.h, stepper.d);
      throw e;
    }

    if (o.interactive) {
      printState(stepper.t, stepper.h, stepper.d);
      cin.get();
    }

    // Keep theta and phi in the range [-180, 180]
    stepper.d.set_theta(Physics::rotate(stepper.d.get_theta()));
    stepper.d.set_phi(Physics::rotate(stepper.d.get_phi()));
    if (stepper.d.get_r() < 1) {
      // Handle collision. Iterate until we get close enough to reflect.
      stepper.undo();
      while (stepper.d.get_r() > 1.0000000000001) {
        stepper.stepHalf();
        if (stepper.d.get_r() < 1) {
          stepper.undo();
        } else {
          event.log(stepper.d, stepper.t);
          ++n;
        }
      }

      event.logCollision(stepper.d, stepper.t);
      printProgress(event.get_n(), stepper.d, true);
      // Specular reflection
      stepper.d.set_pr(-stepper.d.get_pr());

      stepper.reset();
    } else {
      const bool fired = event.log(stepper.d, stepper.t);
      printProgress(event.get_n(), stepper.d, fired);
      ++n;
    }
  }

  printf("\n");
  // printStateHeader();
  // printState(0, h_, freeDipole);
  // printState(stepper.t, stepper.h, stepper.d);

  printf("\n");
  printf("Results output to events.csv.\n");
  printf("\n");

  return freeDipole;
}
