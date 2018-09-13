#ifndef IMSRGProfiler_h
#define IMSRGProfiler_h 1

//#include <unistd.h>
//#include <sstream>
//#include <stdio.h>
//#include <iostream>
//#include <iomanip>
#include <map>
#include <string>

/// Profiling class with all static data members.
/// This is for keeping track of timing and memory usage, etc.
/// IMSRGProfiler methods should probably not be called inside parallel blocks.

//using namespace std;

class IMSRGProfiler
{
 public:
  // timer and counter are declaired as static so that there's only one copy of each of them
  static std::map<std::string, double> timer; ///< For keeping timing information for various method calls
  static std::map<std::string, int> counter;
  static float start_time;

  IMSRGProfiler();
  std::map<std::string,size_t> CheckMem(); // Kbytes  RSS  DIRTY
  std::map<std::string,float> GetTimes(); // real  user  sys
  void PrintTimes();
  void PrintCounters();
  void PrintMemory();
  void PrintAll();
  size_t MaxMemUsage();
};

#endif
