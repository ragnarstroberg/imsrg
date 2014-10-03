
#ifndef ReadWrite_h
#define ReadWrite_h 1

#include "ModelSpace.hh"
#include "Operator.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "string.h"

class ReadWrite
{
 public:
   //void ReadModelSpace( char* filename, ModelSpace * modelspace);
   ModelSpace ReadModelSpace( char* filename);
   //void ReadBareInteraction( char* filename, Operator *Hbare);
   void ReadBareInteraction( char* filename, Operator *Hbare);
//   void ReadBareSPE( char* filename, Operator *Hbare);
   //void ReadBareTBME( char* filename, Operator *Hbare);
   void ReadBareTBME( char* filename, Operator& Hbare);
   void CalculateKineticEnergy(Operator *Hbare);

};



#endif

