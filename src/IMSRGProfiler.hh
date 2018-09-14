///////////////////////////////////////////////////////////////////////////////////
//    IMSRGProfiler.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////


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
