/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2016 John Martin Edwards            **
 ** Idaho State University                            **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwajohn@isu.edu                               **
 ** or visit                                          **
 **    http://www2.cose.isu.edu/~edwajohn/            **
 *******************************************************/

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#if defined (WIN32)
	#include <functional>
#endif

#include "./Options.h"

using namespace std;

// Global Options object
// Options options;

bool Options::ProcessArg(int& i, char** argv) {
  Options& o = *this;
  int orig_i = i;
  if (strcmp(argv[i], "-n") == 0) {
    ++i;
    o.numEvents = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "-h") == 0) {
    ++i;
    o.h = atof(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "-e") == 0) {
    ++i;
    o.eps = atof(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "-b") == 0) {
    ++i;
    o.bool_option = false;
  }
  return i != orig_i;
}
// trim from start
string& ltrim(string &s) {
  s.erase(s.begin(), find_if(s.begin(),
                             s.end(), not1(ptr_fun<int, int>(isspace))));
  return s;
}

// trim from end
string& rtrim(string &s) {
  s.erase(find_if(s.rbegin(), s.rend(),
                  not1(ptr_fun<int, int>(isspace))).base(), s.end());
  return s;
}

// trim from both ends
string& trim(string &s) {
  return ltrim(rtrim(s));
}

void Options::ReadOptionsFile() {
  ifstream in("gvd.config");
  if (!in) return;
  while (!in.eof()) {
    string key;
    in >> key;
    string value;
    getline(in, value);
    
    if (!key.empty() && key[0] != '#')
      key2value[key] = trim(value);
  }
  in.close();
}

string Options::Value(
    const string& key, const string& default_value) const {
  if (key2value.find(key) == key2value.end())
    return default_value;
  const string value = key2value.find(key)->second;
  return value;
}

bool Options::BoolValue(
    const string& key, const bool default_value) const {
  const string value = Value(key, default_value?"true":"false");
  if (value == "0" || value == "false" || value == "False" || value == "FALSE")
    return false;
  return true;
}

int Options::IntValue(
    const string& key, const int default_value) const {
  stringstream ss;
  ss << default_value;
  const string value = Value(key, ss.str());
  return atoi(value.c_str());
}

