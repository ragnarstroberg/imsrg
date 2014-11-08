
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
   void WriteValenceOneBody(Operator &, const char*);
   void WriteValenceTwoBody(Operator &, const char*);
   void WriteNuShellX_sps( Operator& op, const char* filename);
   void WriteNuShellX_int( Operator& op, const char* filename);
   void WriteAntoine_int( Operator& op, const char* filename); // <- not implemented yet...

   std::map<string,string> InputParameters;
};



#endif

