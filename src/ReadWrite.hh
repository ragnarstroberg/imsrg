///////////////////////////////////////////////////////////////////////////////////
//    ReadWrite.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef ReadWrite_h
#define ReadWrite_h 1
#define BOOST_IOSTREAMS_NO_LIB 1

#include <map>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include "ModelSpace.hh"
#include "Operator.hh"
#include "ThreeBodyME.hh"
//#include "DaggerOperator.hh"
#include "Jacobi3BME.hh"


class ReadWrite
{

 public:
   ~ReadWrite();
   ReadWrite();
   void ReadTBME_Oslo( std::string filename, Operator& Hbare);
   void ReadTBME_OakRidge( std::string spname, std::string tbmename, Operator& Hbare, std::string format);
   void ReadBareTBME_Navratil( std::string filename, Operator& Hbare);
   void ReadBareTBME_Navratil_from_stream( std::istream& infile, Operator& Hbare);
   void ReadBareTBME_Darmstadt( std::string filename, Operator& Hbare, int E1max, int E2max, int lmax);
   template<class T> void ReadBareTBME_Darmstadt_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int lmax);
   void Read_Darmstadt_3body( std::string filename, Operator& Hbare, int E1max, int E2max, int E3max);
   size_t Count_Darmstadt_3body_to_read( Operator& Hbare, int E1max, int E2max, int E3max, std::vector<int>& orbits_remap, std::vector<size_t>& nread_list);
   template<class T>void Read_Darmstadt_3body_from_stream( T & infile, Operator& Hbare, int E1max, int E2max, int E3max);
//   void Store_Darmstadt_3body( std::vector<float>& ThreeBME, Operator& Hbare, int E1max, int E2max, int E3max);
   void Store_Darmstadt_3body( const std::vector<float>& ThreeBME, const std::vector<size_t>& nread_list, const std::vector<int>& orbits_remap, Operator& Hbare, int E1max, int E2max, int E3max);
#ifndef NO_HDF5
   void GetHDF5Basis( ModelSpace* modelspace, std::string filename, std::vector<std::array<int,5>>& Basis );
   void Read3bodyHDF5( std::string filename, Operator& op);
   void Read3bodyHDF5_new( std::string filename, Operator& op);
#endif
   size_t Jacobi2b_Channel_Hash(int S, int T, int Tz, int J);
   void Jacobi2b_Channel_UnHash(size_t key, int& S, int& T, int& Tz, int& J);
   void ReadDarmstadt_2bodyRel( std::string filename, Operator& Op );
   void Read2bCurrent_Navratil( std::string filename, Operator& Op);
   void Write_me2j( std::string filename, Operator& op, int emax, int e2max, int lmax);
   void Write_me3j( std::string filename, Operator& op, int E1max, int E2max, int E3max);
   void WriteTBME_Navratil( std::string filename, Operator& Hbare);
   void WriteNuShellX_sps( Operator& op, std::string filename);
   void WriteNuShellX_int( Operator& op, std::string filename);
   void WriteNuShellX_op( Operator& op, std::string filename);
   void ReadNuShellX_int( Operator& op, std::string filename);
   void ReadNuShellX_int_iso( Operator& op, std::string filename);
   void ReadNuShellX_sp( ModelSpace& ms, std::string filename);
   void WriteNuShellX_intfile( Operator& op, std::string filename, std::string mode);
   void WriteAntoine_int( Operator& op, std::string filename);
   void WriteAntoine_input( Operator& op, std::string filename, int A, int Z);
   void WriteOperator(Operator& op, std::string filename);
   void WriteOperatorHuman(Operator& op, std::string filename);
   void ReadOperator(Operator& op, std::string filename);
   void ReadOperatorHuman(Operator& op, std::string filename);
   void CompareOperators(Operator& op1, Operator& op2, std::string filename);
   void WriteTensorOneBody(std::string filename, Operator& H, std::string opname);
   void WriteTensorTwoBody(std::string filename, Operator& H, std::string opname);
   void WriteDaggerOperator( Operator& op, std::string filename, std::string opname="");

//   void WriteValence3body( ThreeBodyMEpn& threeBME, std::string filename );
   void WriteValence3body( ThreeBodyME& threeBME, std::string filename );

   void ReadBareTBME_Jason( std::string filename, Operator& Hbare);
   void ReadTensorOperator_Nathan( std::string filename1b, std::string filename2b, Operator& op);
   void ReadOperator_Nathan( std::string filename1b, std::string filename2b, Operator& op);
   void ReadOneBody_Takayuki(std::string filename, Operator& Hbare);
   void ReadTwoBody_Takayuki(std::string filename, Operator& Hbare);
   void WriteOneBody_Takayuki(std::string filename, Operator& Hbare);
   void WriteTwoBody_Takayuki(std::string filename, Operator& Hbare);
   void WriteOneBody_Simple(std::string filename, Operator& Hbare);
   void WriteOneBody_Oslo(std::string filename, Operator& Hbare);
   void WriteTwoBody_Oslo(std::string filename, Operator& Hbare);
   void ReadTwoBodyEngel(std::string filename, Operator& Op);
   void ReadTwoBodyEngel_from_stream(std::istream& infile, Operator& Op);
   void ReadRelCMOpFromJavier( std::string statefile, std::string MEfile, Operator& Op);
   void SetLECs(double c1, double c3, double c4, double cD, double cE);
   std::array<double,5> GetLECs(){return LECs;};
   void SetLECs_preset(std::string);
   void SetCoMCorr(bool b){doCoM_corr = b;std::cout <<"Setting com_corr to "<< b << std::endl;};
   void SetScratchDir( std::string d){scratch_dir = d; std::cout << "Setting scratch dir to " << d << std::endl;};
   std::string GetScratchDir(){return scratch_dir;};
   int GetAref(){return Aref;};
   int GetZref(){return Zref;};
   void SetAref(int a){Aref = a;};
   void SetZref(int z){Zref = z;};
   void Set3NFormat( std::string fmt ){format3N=fmt;};

   void ReadJacobi3NFiles( Jacobi3BME& jacobi3bme, std::string poi_name, std::string eig_name, std::string v3int_name );

   // added by T.Miyagi
   void ReadTokyo(std::string filename, Operator& op, std::string fmt);
   void ReadTokyo(std::string filename, Operator& op);
   void WriteTokyo(Operator& op, std::string filename, std::string mode);
   void WriteTokyoFull(Operator& op, std::string filename); // only for Hamiltonian
   void WriteTensorTokyo(std::string filename, Operator& op);
   Operator ReadOperator2b_Miyagi(std::string, ModelSpace&); // general operator me2j-like format
   void skip_comments(std::ifstream&);



   // added by A.Belley
//   void WriteOmega(std::string filename, std::string scratch, int size);
   void CopyFile(std::string file1, std::string file2);

   // Fields

//   std::map<std::string,std::string> InputParameters; // I believe this is very very deprecated

   bool InGoodState(){return goodstate;};
   bool doCoM_corr;
   bool goodstate;
   std::array<double,5> LECs;
   std::string scratch_dir;
   std::string File2N;
   std::string File3N;
   std::string format3N;
   int Aref;
   int Zref;


};



/// Wrapper class so we can treat a std::vector of floats like a stream, using the extraction operator >>.
/// This is used for the binary version of ReadWrite::Read_Darmstadt_3body_from_stream().
//template <typename T>
class VectorStream
{
 public:
//  VectorStream(std::vector<T>& v) : vec(v), i(0) {};
  VectorStream(std::vector<float>& v) : vec(v), i(0) {};
//  VectorStream(std::vector<double>& v) : vec(v), i(0) {};
//  VectorStream& operator>>(T& x) { x = vec[i++]; return (VectorStream&)(*this);}
  VectorStream& operator>>(float& x) { x = vec[i++]; return (VectorStream&)(*this);}
//  VectorStream& operator>>(double& x) { x = vec[i++]; return (VectorStream&)(*this);}
  bool good(){ return i<vec.size(); };
  void getline(char[], int) {}; // Don't do nuthin'.
  void read(char* buf, size_t len) {memcpy((void*)buf, (const void*)&vec[i], len);}; // Totally untested...
 private:
//  std::vector<T>& vec;
  std::vector<float>& vec;
//  std::vector<double>& vec;
  long long unsigned int i;
};




#endif

