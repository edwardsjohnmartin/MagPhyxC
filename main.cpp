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

const double default_n = 1e5;
const double default_h = 1e-2;
const double default_eps = 1e-10;
Options o(default_n, default_h, default_eps);

vector<string> split(string str, char delimiter);
Dipole initDipole(const string& filename);
Dipole doSimulation(const Dipole& freeDipole, Event& event,
                       double h,
                       const int numEvents);

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "\n");
    fprintf(stderr, "SYNOPSIS\n");
    fprintf(stderr, "\tmagphyx [OPTIONS] filename.csv\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DESCRIPTION\n");
    fprintf(stderr, "\tmagphyx runs a magnet simulation given initial\n"
            "\tconditions specified in the second line of filename.csv.\n"
            "\tfilename.csv is the same format as what is exported from\n"
            "\tMagPhyx web version. Events are output to events.csv.\n"
            );
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS\n");
    fprintf(stderr, "\t-n numEvents\n");
    fprintf(stderr, "\t\tExecutes the simulation until numEvents events "
            "occur. Default = 1e5.\n");
    fprintf(stderr, "\t-h h\n");
    fprintf(stderr, "\t\tInitial step size. Default = 1e-2.\n");
    fprintf(stderr, "\t-e eps\n");
    fprintf(stderr, "\t\tError per step allowed. Note that this is error in\n"
            "\t\tterms of Runge-Kutta. The error in total energy will be\n"
            "\t\tsimilar to, but not bound by, this value. Default = 1e-10.\n");
    fprintf(stderr, "\n");
    return 1;
  }

  int i = 1;
  bool stop = false;
  while (i < argc && !stop) {
    stop = true;
    if (o.ProcessArg(i, argv)) {
      stop = false;
    }
  }
  string initFilename = argv[i];

  Dipole freeDipole = initDipole(initFilename);
  // cout << freeDipole << endl;
  Event event("events.csv", freeDipole);
  freeDipole = doSimulation(freeDipole, event, o.h, o.numEvents);
  // cout << freeDipole << endl;
}

Dipole initDipole(const string& filename) {
  ifstream in(filename);
  string line;
  vector<string> tokens;

  // skip header line
  getline(in, line, '\r');
  // cout << line << endl;
  tokens = split(line, ',');
  if (tokens[0].compare("n") != 0) {
    throw logic_error("Illegal file");
  }
  
  // first data line
  getline(in, line, '\r');
  // cout << line << endl;
  tokens = split(line, ',');
  int i = 0;
  // cout << tokens[i] << endl;
  int n = stoi(tokens[i++]);
  string event_type = tokens[i++];
  double t = stod(tokens[i++]);
  double r = stod(tokens[i++]);
  double theta = Physics::deg2rad(stod(tokens[i++]));
  double phi = Physics::deg2rad(stod(tokens[i++]));
  double pr = stod(tokens[i++]);
  double ptheta = stod(tokens[i++]);
  double pphi = stod(tokens[i++]);
  double beta = stod(tokens[i++]);
  double E = stod(tokens[i++]);
  double dE = stod(tokens[i++]);

  return Dipole(r, theta, phi, pr, ptheta, pphi);
}

// You could also take an existing vector as a parameter.
vector<string> split(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok;
  
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }
  
  return internal;
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
  if (n % 1000 == 0 && fired) {
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
  printProgress(0, freeDipole, true);
  for (int i = 0; event.get_n() < numEvents; ++i) {
    try {
      stepper.step();
    } catch (logic_error& e) {
      printState(stepper.t, stepper.h, stepper.d);
      throw e;
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
