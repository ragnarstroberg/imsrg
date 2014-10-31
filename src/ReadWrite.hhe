
#ifndef ReadWrite_h
#define ReadWrite_h 1

#include <map>
#include <string>
#include "Operator.hh"


class ReadWrite
{
 public:
   void ReadSettingsFile(  char* filename);
   ModelSpace ReadModelSpace( const char* filename);
   void ReadBareTBME( const char* filename, Operator& Hbare);
   void CalculateKineticEnergy(Operator *Hbare);
   void WriteOneBody(Operator &, const char*);
   void WriteTwoBody(Operator &, const char*);
   //std::map<char*, char*> InputParameters;
   std::map<string,string> InputParameters;

};



#endif

