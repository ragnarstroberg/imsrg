#ifndef IMSRGProfiler_h
#define IMSRGProfiler_h 1

#include <unistd.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <map>

/// Profiling class with all static data members.
/// This is for keeping track of timing and memory usage, etc.
/// IMSRGProfiler methods should probably not be called inside parallel blocks.

using namespace std;

class IMSRGProfiler
{
 public:
  // timer and counter are declaired as static so that there's only one copy of each of them
  static map<string, double> timer; ///< For keeping timing information for various method calls
  static map<string, int> counter;
  static float start_time;

  IMSRGProfiler();
  map<string,size_t> CheckMem() const; // Kbytes  RSS  DIRTY
  map<string,float> GetTimes() const; // real  user  sys
  void PrintTimes() const;
  void PrintCounters() const;
  void PrintMemory() const;
  void PrintAll() const;
  size_t MaxMemUsage() const;
};

#endif
