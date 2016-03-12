
#ifndef ReadWrite_h
#define ReadWrite_h 1
#define BOOST_IOSTREAMS_NO_LIB 1

#include <map>
#include <string>
#include "Operator.hh"

using namespace std;


class ReadWrite
{
 public:
   ~ReadWrite();
   ReadWrite();
   void ReadSettingsFile(  string filename);
   void ReadTBME_Oslo( string filename, Operator& Hbare);
   void ReadTBME_OakRidge( string filename, Operator& Hbare);
   void ReadBareTBME_Jason( string filename, Operator& Hbare);
   void ReadBareTBME_Navratil( string filename, Operator& Hbare);
   void ReadBareTBME_Navratil_from_stream( istream& infile, Operator& Hbare);
   void ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int E1max, int E2max, int lmax);
   template<class T> void ReadBareTBME_Darmstadt_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int lmax);
   void Read_Darmstadt_3body( string filename, Operator& Hbare, int E1max, int E2max, int E3max);
   template<class T>void Read_Darmstadt_3body_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int E3max);
   void GetHDF5Basis( ModelSpace* modelspace, string filename, vector<array<int,5>>& Basis );
   void Read3bodyHDF5( string filename, Operator& op);
   void Read3bodyHDF5_new( string filename, Operator& op);
   void Write_me2j( string filename, Operator& op, int emax, int e2max, int lmax);
   void Write_me3j( string filename, Operator& op, int E1max, int E2max, int E3max);
   void WriteTBME_Navratil( string filename, Operator& Hbare);
   void WriteNuShellX_sps( Operator& op, string filename);
   void WriteNuShellX_int( Operator& op, string filename);
   void WriteNuShellX_op( Operator& op, string filename);
   void WriteNuShellX_intfile( Operator& op, string filename, string mode); 
   void WriteAntoine_int( Operator& op, string filename); // <- not implemented yet...
   void WriteOperator(Operator& op, string filename);
   void WriteOperatorHuman(Operator& op, string filename);
   void ReadOperator(Operator& op, string filename); 
   void CompareOperators(Operator& op1, Operator& op2, string filename);
   void ReadOneBody_Takayuki(string filename, Operator& Hbare);
   void ReadTwoBody_Takayuki(string filename, Operator& Hbare);
   void WriteOneBody_Takayuki(string filename, Operator& Hbare);
   void WriteTwoBody_Takayuki(string filename, Operator& Hbare);
   void WriteTensorOneBody(string filename, Operator& H, string opname);
   void WriteTensorTwoBody(string filename, Operator& H, string opname);
   void WriteOneBody_Oslo(string filename, Operator& Hbare);
   void WriteTwoBody_Oslo(string filename, Operator& Hbare);
   void ReadTwoBodyEngel(string filename, Operator& Op);
   void ReadTwoBodyEngel_from_stream(istream& infile, Operator& Op);
   void SetLECs(double c1, double c3, double c4, double cD, double cE);
   array<double,5> GetLECs(){return LECs;};
   void SetLECs_preset(string);
   void SetCoMCorr(bool b){doCoM_corr = b;cout <<"Setting com_corr to "<< b << endl;};
   void SetScratchDir( string d){scratch_dir = d;};
   string GetScratchDir(){return scratch_dir;};
   int GetAref(){return Aref;};
   int GetZref(){return Zref;};
   void SetAref(int a){Aref = a;};
   void SetZref(int z){Zref = z;};

   // Fields

   std::map<string,string> InputParameters;

   bool InGoodState(){return goodstate;};
   bool doCoM_corr;
   bool goodstate;
   array<double,5> LECs;
   string scratch_dir;
   string File2N;
   string File3N;
   int Aref;
   int Zref;   


};



/// Wrapper class so we can treat a vector of floats like a stream, using the extraction operator >>.
/// This is used for the binary version of ReadWrite::Read_Darmstadt_3body_from_stream().
class VectorStream 
{
 public:
  VectorStream(vector<float>& v) : vec(v), i(0) {};
//  VectorStream(vector<double>& v) : vec(v), i(0) {};
  VectorStream& operator>>(float& x) { x = vec[i++]; return (VectorStream&)(*this);}
//  VectorStream& operator>>(double& x) { x = vec[i++]; return (VectorStream&)(*this);}
  bool good(){ return i<vec.size(); };
  void getline(char[], int) {}; // Don't do nuthin'.
 private:
  vector<float>& vec;
//  vector<double>& vec;
  long long unsigned int i;
};

#endif

