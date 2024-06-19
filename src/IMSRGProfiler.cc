#include "IMSRGProfiler.hh"
#include <algorithm>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


std::map<std::string, double> IMSRGProfiler::timer;
std::map<std::string, int> IMSRGProfiler::counter;
double IMSRGProfiler::start_time = -1;

IMSRGProfiler::IMSRGProfiler()
{
  if (start_time < 0)
  {
    start_time = omp_get_wtime();
    counter["N_Threads"] = omp_get_max_threads();
  }
}
/// Check how much memory is being used.
///
std::map<std::string,size_t> IMSRGProfiler::CheckMem()
{
  std::map<std::string,size_t> s;
  s["Kbytes"]=0;s["RSS"]=0;s["DIRTY"]=0;

  std::ifstream ifs("/proc/self/statm"); // this is a pseudo file that provides information on the current process (linux only)
  if ( ifs.good() )
  {
    size_t npages;
    size_t page_size = sysconf(_SC_PAGESIZE);
    ifs >> npages >> npages;
    s["RSS"] = npages * page_size /1024; // divide by 1024 to convert to kilobytes
  }

/*
//#if defined(CHECKMEM) && ! defined(__APPLE__)
//  sprintf(commandstr,"pmap -x %d | tail -1",getpid()); // TODO make this more portable. On OSX, use vmmap. no idea for Windows...
  std::ostringstream commandstr;
  commandstr << "pmap -x " << getpid() << " | tail -1";
//#ifndef __APPLE__
//  FILE* output = popen(commandstr,"r");
  FILE* output = popen(commandstr.str().c_str(),"r");
  if (output==NULL or fgets(outbuf,500,output) == NULL)
    std::cout << " <<< IMSRGProfiler::CheckMem():  Problem reading output of pmap (pid = " << getpid() << ")" << std::endl;
  else
    std::istringstream(outbuf) >> buf >> buf >> s["Kbytes"] >> s["RSS"] >> s["DIRTY"];
//    std::istringstream(outbuf) >> commandstr >> buf >> s["Kbytes"] >> s["RSS"] >> s["DIRTY"];
*/

//#endif
  return s;
}

size_t IMSRGProfiler::MaxMemUsage()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  return (size_t) ru.ru_maxrss;
}

std::map<std::string,float> IMSRGProfiler::GetTimes()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  std::map<std::string,float> times;
  times["user"] = ru.ru_utime.tv_sec + 1e-6*ru.ru_utime.tv_usec;
  times["system"] = ru.ru_stime.tv_sec + 1e-6*ru.ru_stime.tv_usec;
  times["real"] = omp_get_wtime() - start_time;
  return times;
}

void IMSRGProfiler::PrintTimes()
{
   auto time_tot = GetTimes();
   
   std::cout << "====================== TIMES (s) ====================" << std::endl;
   std::cout.setf(std::ios::fixed);
   size_t max_func_length = 40;
   for ( auto it : timer ) {
    max_func_length = std::max( it.first.size(), max_func_length);
   }
   max_func_length += 5;
   for ( auto it : timer )
   {
     int nfill = (int) (20 * std::min( 1.0, it.second / time_tot["real"]));
     std::cout << std::setw(max_func_length) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(5) << std::right << it.second;
     std::cout << " (" << std::setw(4) << std::setprecision(1) << 100*it.second / time_tot["real"] << "%) |";
     for (int ifill=0; ifill<nfill; ifill++) std::cout << "*";
     for (int ifill=nfill; ifill<20; ifill++) std::cout << " ";
     std::cout << "|";
     std::cout  << std::endl;
     
   }
//   for (auto it : GetTimes())
   for (auto it : time_tot)
     std::cout << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(5) << std::right << it.second  << std::endl;
}

void IMSRGProfiler::PrintCounters()
{
   std::cout << "===================== COUNTERS =====================" << std::endl;
   std::cout.setf(std::ios::fixed);
   for ( auto it : counter )
     std::cout << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(0) << std::right << it.second  << std::endl;
}

void IMSRGProfiler::PrintMemory()
{
   std::cout << "===================== MEMORY (MB) ==================" << std::endl;
   for (auto it : CheckMem())
     std::cout << std::fixed << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(3) << std::right << it.second/1024. << std::endl;

   std::cout << std::fixed << std::setw(40) << std::left << "Max Used:  " << std::setw(12) << std::setprecision(3) << std::right << MaxMemUsage()/1024.  << std::endl;
}

void IMSRGProfiler::PrintAll()
{
  PrintCounters();
  PrintTimes();
  PrintMemory();
}


void IMSRGProfiler::Clear()
{
  timer.clear();
  counter.clear();
  start_time = omp_get_wtime();
}

