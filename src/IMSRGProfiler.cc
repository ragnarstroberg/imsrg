#include "IMSRGProfiler.hh"
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/resource.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <sstream>


std::map<std::string, double> IMSRGProfiler::timer;
std::map<std::string, int> IMSRGProfiler::counter;
float IMSRGProfiler::start_time = -1;

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
  char commandstr[100],outbuf[500],buf[100];
  std::map<std::string,size_t> s;
  sprintf(commandstr,"pmap -x %d | tail -1",getpid()); // TODO make this more portable. On OSX, use vmmap. no idea for Windows...
#ifndef __APPLE__
  FILE* output = popen(commandstr,"r");
  if (output==NULL or fgets(outbuf,500,output) == NULL)
    std::cout << " <<< IMSRGProfiler::CheckMem():  Problem reading output of pmap (pid = " << getpid() << ")" << std::endl;
  else
    std::istringstream(outbuf) >> commandstr >> buf >> s["Kbytes"] >> s["RSS"] >> s["DIRTY"];
#else
  s["Kbytes"]=0;s["RSS"]=0;s["DIRTY"]=0;
#endif
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
   for ( auto it : timer )
   {
     int nfill = (int) (20 * it.second / time_tot["real"]);
     std::cout << std::setw(40) << std::left << it.first + ":  " << std::setw(12) << std::setprecision(5) << std::right << it.second;
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


