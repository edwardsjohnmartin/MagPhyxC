/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include <vector>
#include <set>
#include <map>

// #include "./opencl/defs.h"

//------------------------------------------------------------------------------
// Options
//------------------------------------------------------------------------------
struct Options {
 public:
  std::vector<std::string> filenames;
  int numEvents;
  double h;
  double eps;
  bool bool_option;
  std::map<std::string, std::string> key2value;

 public:
  Options(const int numEvents_, const double h_, const double eps_)
      : numEvents(numEvents_), h(h_), eps(eps_), bool_option(true) {
    ReadOptionsFile();
  }

  bool ProcessArg(int& i, char** argv);

  std::string Value(
    const std::string& key, const std::string& default_value) const;
  bool BoolValue(const std::string& key, const bool default_value) const;
  int IntValue(const std::string& key, const int default_value) const;

  // The options file is key/value pairs, such as
  //   TEST_AMBIGUOUS_GPU 1
  //   DISPLAY_SOMETHING 0
  //   IMPORTANT_MATRIX 1 0 0 0 1 0 0 0 1
  void ReadOptionsFile();

};

// Global Options object
extern Options options;

#endif