
#ifndef ReadWrite_h
#define ReadWrite_h 1

#include "Operator.hh"


class ReadWrite
{
 public:
   ModelSpace ReadModelSpace( char* filename);
   void ReadBareInteraction( char* filename, Operator *Hbare);
   void ReadBareTBME( char* filename, Operator& Hbare);
   void CalculateKineticEnergy(Operator *Hbare);
   void WriteOneBody(Operator &, const char*);
   void WriteTwoBody(Operator &, const char*);

};



#endif

