
#ifndef ReadWrite_h
#define ReadWrite_h 1

#include <map>
#include <string>
#include "Operator.hh"

// ToDo: Add WriteOperator/ ReadOperator methods, which stores an operator in
// a convenient method for reading back in another calculation

class ReadWrite
{
 public:
   ReadWrite();
   void ReadSettingsFile(  string filename);
   ModelSpace ReadModelSpace( string filename);
   void ReadBareTBME( string filename, Operator& Hbare);
   void ReadBareTBME_Jason( string filename, Operator& Hbare);
   void ReadBareTBME_Navratil( string filename, Operator& Hbare);
   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax=-1);
   vector<int> MapDarmstadtOrbits( ModelSpace* modelspace, int& Emax);
   void WriteOneBody(Operator&, string);
   void WriteTwoBody(Operator&, string);
   void WriteValenceOneBody(Operator&, string);
   void WriteValenceTwoBody(Operator&, string);
   void WriteNuShellX_sps( Operator& op, string filename);
   void WriteNuShellX_int( Operator& op, string filename);
   void WriteAntoine_int( Operator& op, string filename); // <- not implemented yet...
   void WriteOperator(Operator& op, string filename);
   void ReadOperator(Operator& op, string filename);

   std::map<string,string> InputParameters;

   bool InGoodState(){return goodstate;};
   bool doCoM_corr;
   void SetCoMCorr(bool b){doCoM_corr = b;};
   bool goodstate;

};



#endif

