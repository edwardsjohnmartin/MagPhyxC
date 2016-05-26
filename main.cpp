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
Dipole updatePositions(const Dipole& freeDipole, Event& event,
                       double h,
                       const int numEvents);

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: MagPhyx [-n numSteps] [-h h] [-e eps] filename.csv" << endl;
    cerr << "  -n length of the simulation, in the form of number of events. "
         << "Default = 1e5" << endl;
    cerr << "  -h initial step size. Default = 1e-2" << endl;
    cerr << "  -e error per step allowed. Default = 1e-10" << endl;
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
  freeDipole = updatePositions(freeDipole, event, o.h, o.numEvents);
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

void printStatusHeader() {
  // printf(" %5s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
  //        "i", "h", "r", "theta", "phi", "pr", "ptheta", "pphi", "E", "dE");
  printf("\n");
  printf(" %5s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
         "t", "h", "r", "theta", "phi", "pr", "ptheta", "pphi", "E", "dE");
  printf("-----------------------------------------"
         "-----------------------------------------"
         "------------------------------------------\n");
}

void printStatus(const int i, const double t, const double h, const Dipole& d,
                 const int numSteps, const bool forcePrint = false) {
  const int n = (int)ceil(numSteps/20.0);
  if (i % n == 0 || forcePrint) {
    // printf("%s%5d %12.6g %12.5g %12.6g %12.6g %12.6g"
    //        "%12.6g %12.6g %12.6g %12.6g\n",
    printf("%s%5.0lf %12.6g %12.5g %12.6g %12.6g %12.6g"
           "%12.6g %12.6g %12.6g %12.6g\n",
        // (forcePrint?"*":" "), i, h,
        (forcePrint?"*":" "), t, h,
        d.get_r(), Physics::rad2deg(d.get_theta()),
        Physics::rad2deg(d.get_phi()),
        d.get_pr(), d.get_ptheta(), d.get_pphi(),
        d.get_E(), d.get_dE());
  }
}

void printStatus(const Event& event, const bool fired) {
  if (event.get_n() % 1000 == 0 && fired) {
    printf("\b\b\b\b\b\b\b%-7d", event.get_n());
    fflush(stdout);
  }
}

Dipole updatePositions(const Dipole& freeDipole, Event& event,
                       const double h_,
                       const int numEvents) {
  Stepper stepper(freeDipole, h_, o.eps);

  // Dipole d(freeDipole);
  double t = 0.0;
  // double h = o.h;

  int n = 0;
  // printStatusHeader();
  // printStatus(0, 0, stepper.h, stepper.d, numEvents);

  printf("\n# events = %-7d", 0);
  for (int i = 0; event.get_n() < numEvents; ++i) {
  // for (int i = 0; i < 10000; ++i) {
    try {
      stepper.step();
    } catch (logic_error& e) {
      printStatus(i, stepper.t, stepper.h, stepper.d, numEvents);
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
      printStatus(event, true);
      // Specular reflection
      stepper.d.set_pr(-stepper.d.get_pr());

      stepper.reset();
    } else {
      const bool fired = event.log(stepper.d, stepper.t);
      // printStatus(i, t, h, d, numEvents);
      printStatus(event, fired);
      ++n;
    }
  }

  printStatusHeader();
  printStatus(0, 0, h_, freeDipole, numEvents);
  printStatus(0, 0, stepper.h, stepper.d, numEvents);

  printf("\n");
  printf("Num events = %d\n", event.get_n());
  // printf("Num steps = %d\n", n);

  return freeDipole;
}
