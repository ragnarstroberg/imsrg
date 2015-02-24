
#ifndef ReadWrite_h
#define ReadWrite_h 1

#include <map>
#include <string>
#include "Operator.hh"
//#include "Operator3.hh"

using namespace std;
class ReadWrite
{
 public:
   ReadWrite();
   void ReadSettingsFile(  string filename);
   ModelSpace ReadModelSpace( string filename);
   void ReadBareTBME( string filename, Operator& Hbare);
   void ReadBareTBME_Jason( string filename, Operator& Hbare);
//   template <class OPERATOR> void ReadBareTBME_Navratil( string filename, OPERATOR Hbare);
//template <class OPERATOR> void ReadWrite::ReadBareTBME_Navratil( string, Operator&);
//template <class OPERATOR> void ReadWrite::ReadBareTBME_Navratil( string, Operator3&);
//   void ReadBareTBME_Navratil( string filename, Operator3& Hbare);
   void ReadBareTBME_Navratil( string filename, Operator& Hbare);
//   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax=-1);
//   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare);
//   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax);
   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax=-1, int lmax=-1);
//   void Read_Darmstadt_3body( string filename, Operator3& Hbare, int Emax=-1);
//   void Read_Darmstadt_3body( string filename, Operator& Hbare, int Emax=-1);
   void Read_Darmstadt_3body( string filename, Operator& Hbare, int E1max, int E2max, int E3max);
//   vector<int> MapDarmstadtOrbits( ModelSpace* modelspace, int& Emax);
   void WriteOneBody(Operator&, string);
   void WriteTwoBody(Operator&, string);
   void WriteValenceOneBody(Operator&, string);
   void WriteValenceTwoBody(Operator&, string);
   void WriteNuShellX_sps( Operator& op, string filename);
   void WriteNuShellX_int( Operator& op, string filename);
   void WriteAntoine_int( Operator& op, string filename); // <- not implemented yet...
   void WriteOperator(Operator& op, string filename);
   void ReadOperator(Operator& op, string filename); // Not implemented yet...

   std::map<string,string> InputParameters;

   bool InGoodState(){return goodstate;};
   bool doCoM_corr;
   void SetCoMCorr(bool b){doCoM_corr = b;cout <<"Setting com_corr to "<< b << endl;};
   bool goodstate;

};

//template void ReadWrite::ReadBareTBME_Navratil<Operator&>( string, Operator&);
//template void ReadWrite::ReadBareTBME_Navratil<Operator3&>( string, Operator3&);


#endif

