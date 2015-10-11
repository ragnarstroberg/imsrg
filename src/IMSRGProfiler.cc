#include "IMSRGProfiler.hh"
#include <sys/time.h>
#include <sys/resource.h>


map<string, double> IMSRGProfiler::timer;
map<string, int> IMSRGProfiler::counter;


/// Check how much memory is being used.
///
map<string,size_t> IMSRGProfiler::CheckMem()
{
  char cmdstring[100],outbuf[500],buf[100];
  sprintf(cmdstring,"pmap -x %d | tail -1",getpid()); // TODO make this more portable. On OSX, use vmmap. no idea for Windows...
  FILE* output = popen(cmdstring,"r");
  fgets(outbuf,500,output);
  map<string,size_t> s;
  istringstream(outbuf) >> cmdstring >> buf >> s["Kbytes"] >> s["RSS"] >> s["DIRTY"];
  return s;
}

size_t IMSRGProfiler::MaxMemUsage()
{
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  return (size_t) ru.ru_maxrss;
}

void IMSRGProfiler::PrintTimes()
{
   cout << "====================== TIMES =======================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : timer )
   {
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(5) << std::right << it.second  << endl;
   }
}

void IMSRGProfiler::PrintCounters()
{
   cout << "===================== COUNTERS =====================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : counter )
   {
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(0) << std::right << it.second  << endl;
   }
}

void IMSRGProfiler::PrintAll()
{
  PrintCounters();
  PrintTimes();
}


