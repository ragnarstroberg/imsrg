#include "ReadWrite.hh"
#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <ctime>
#include <unordered_map>
#include "omp.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#ifndef NO_HDF5
#include "H5Cpp.h"
#endif

#define LINESIZE 496
//#define HEADERSIZE 500
#define HEADERSIZE 255
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif
#ifndef HBARC
   #define HBARC 197.3269718 // hc in MeV * fm
#endif
#ifndef M_NUCLEON
   #define M_NUCLEON 938.9185
#endif

using namespace std;
#ifndef NO_HDF5
using namespace H5;
#endif

ReadWrite::~ReadWrite()
{
//  std::cout << "In ReadWrite destructor" << std::endl;
}

ReadWrite::ReadWrite()
: doCoM_corr(false), goodstate(true),LECs({-0.81,-3.20,5.40,1.271,-0.131}),File2N("none"),File3N("none"),format3N("me3j"),Aref(0),Zref(0) // default to the EM2.0_2.0 LECs
{
}

// This is old and deprecated.
void ReadWrite::ReadSettingsFile( std::string filename)
{
   char line[LINESIZE];
   std::string lstr;
   std::ifstream fin;
   fin.open(filename);
   if (! fin.good() )
   {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "!!!!  Trouble reading Settings file " << filename << "!!!!!!" << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      goodstate = false;
   }
   std::cout << "Reading settings file " << filename << std::endl;
   while (fin.getline(line,LINESIZE))
   {
      lstr = std::string(line);
      lstr = lstr.substr(0, lstr.find_first_of("#"));
      int colon = lstr.find_first_of(":");
      std::string param = lstr.substr(0, colon);
      if ( param.size() <1) continue;
      std::string value = lstr.substr(colon+1, lstr.length()-colon);
      value.erase(0,value.find_first_not_of(" \t\r\n"));
      value.erase(value.find_last_not_of(" \t\r\n")+1);
      if (value.size() < 1) continue;
      InputParameters[param] = value;
      std::cout << "parameter: [" << param << "] = [" << value << "]" << std::endl;
   }

}




/// Read two-body matrix elements from an Oslo-formatted file
void ReadWrite::ReadTBME_Oslo( std::string filename, Operator& Hbare)
{

  std::ifstream infile;
  char line[LINESIZE];
  int Tz,Par,J2,a,b,c,d;
  double fbuf[3];
  double tbme;
  int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
  File2N = filename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }

  infile.getline(line,LINESIZE);

  while (!strstr(line,"<ab|V|cd>") && !infile.eof()) // Skip lines until we see the header
  {
     infile.getline(line,LINESIZE);
  }

  // read the file one line at a time
  while ( infile >> Tz >> Par >> J2 >> a >> b >> c >> d >> tbme >> fbuf[0] >> fbuf[1] >> fbuf[2] )
  {
     // if the matrix element is outside the model space, ignore it.
     if (a>norbits or b>norbits or c>norbits or d>norbits) continue;
     a--; b--; c--; d--; // Fortran -> C  ==> 1 -> 0

     double com_corr = fbuf[2] * Hbare.GetModelSpace()->GetHbarOmega() / Hbare.GetModelSpace()->GetTargetMass();  

// NORMALIZATION: Read in normalized, antisymmetrized TBME's

     if (doCoM_corr)  tbme-=com_corr;

     Hbare.TwoBody.SetTBME(J2/2,Par,Tz,a,b,c,d, tbme ); // Don't do COM correction,

  }

  return;
}

/// Read two-body matrix elements from an Oslo-formatted file, as obtained from Gaute Hagen
void ReadWrite::ReadTBME_OakRidge( std::string spname, std::string tbmename, Operator& Hbare, std::string tbme_format )
{
  ModelSpace * modelspace = Hbare.GetModelSpace();
  std::ifstream spfile(spname);
  if (!spfile.good())
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << spname << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  std::cout << "reading Oak Ridge " << tbme_format << " format" << std::endl;

  File2N = tbmename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();
  double hw_file,spe,dummy;
  spfile >> hw_file;
  int index,n,l,j2,tz2;
  std::vector<index_t> orbit_remap;
  while( spfile >> index >> n >> l >> j2 >> tz2 >> spe >> dummy)
  {
    index_t orbit_index = modelspace->GetOrbitIndex(n,l,j2,tz2); // check which isospin convention is used
    if ((index_t)index > orbit_remap.size() ) orbit_remap.resize(index,-1);
//    std::cout << index << " " << n << " " << l << " " << j2 << " " << tz2 << " " << spe << std::endl;
    orbit_remap[index-1] = orbit_index;
  }
//  std::cout << "done with spe" << std::endl;

  spfile.close();

  auto openmode = std::ios::in;
  if (tbme_format=="binary") openmode |= std::ios::binary;
  std::ifstream tbmefile(tbmename, openmode);
//  std::ifstream tbmefile = (tbme_format=="binary") ? std::ifstream(tbmename, std::ios::in | std::ios::binary ) :  std::ifstream(tbmename );
//  if (tbme_format == "binary")
//  {
//    tbmefile = std::ifstream(tbmename, std::ios::in | std::ios::binary );
//  }
//  else
//  {
//    tbmefile = std::ifstream(tbmename );
//  }
//  std::ifstream tbmefile(tbmename, std::ios::in | std::ios::binary );
  if (!tbmefile.good())
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << tbmename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }

  // read 7 integers and 3 doubles
  while( tbmefile.good() )
  {
    int Tz,P,J2,aa,bb,cc,dd;
    double g1;//,g2,g3; // not using g2,g3 right now
    std::vector<int32_t> vint(7);
    std::vector<double> vdouble(3);
    if (tbme_format == "binary")
    {
      tbmefile.read( reinterpret_cast<char*>(vint.data()), vint.size()*sizeof(int32_t) );
      if (!tbmefile.good()) return;
      tbmefile.read( reinterpret_cast<char*>(vdouble.data()), vdouble.size()*sizeof(int64_t) );
      if (!tbmefile.good()) return;
    }
    else
    {
      tbmefile >> vint[0] >> vint[1] >> vint[2] >> vint[3] >> vint[4] >> vint[5] >> vint[6] >> vdouble[0] >> vdouble[1] >> vdouble[2];
    }
    Tz = vint[0];
    P  = vint[1];
    J2 = vint[2];
    aa = vint[3];
    bb = vint[4];
    cc = vint[5];
    dd = vint[6];

    g1 = vdouble[0];
//    g2 = vdouble[1];
//    g3 = vdouble[2];


//    double fnorm = 1.0;
//    if (Tz != 0)
//    {
//      if (aa==bb) fnorm /= sqrt(2.);
//      if (cc==dd) fnorm /= sqrt(2.);
//    }
    int a = orbit_remap[aa-1];
    int b = orbit_remap[bb-1];
    int c = orbit_remap[cc-1];
    int d = orbit_remap[dd-1];
//    double tbme = g1 * fnorm;
    double tbme = g1;
    if (a>=modelspace->GetNumberOrbits()) continue;
    if (b>=modelspace->GetNumberOrbits()) continue;
    if (c>=modelspace->GetNumberOrbits()) continue;
    if (d>=modelspace->GetNumberOrbits()) continue;
    Orbit& oa = modelspace->GetOrbit(a);
    Orbit& ob = modelspace->GetOrbit(b);
    Orbit& oc = modelspace->GetOrbit(c);
    Orbit& od = modelspace->GetOrbit(d);
    if ( 2*(oa.n+ob.n)+oa.l+ob.l > modelspace->GetE2max() ) continue;
    if ( 2*(oc.n+od.n)+oc.l+od.l > modelspace->GetE2max() ) continue;



    Hbare.TwoBody.SetTBME(J2/2,P,Tz,a,b,c,d, tbme ); 
  }


}



/// Read two-body matrix elements from an Oslo-formatted file
void ReadWrite::WriteTwoBody_Oslo( std::string filename, Operator& Op)
{

  std::ofstream outfile(filename);
  ModelSpace* modelspace = Op.GetModelSpace();
  int wint = 8;
  int wdouble = 12;
  int dprec = 6;
  outfile << fixed;
  outfile << setw(wint) << setprecision(dprec);

  if ( !outfile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble writing file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }

  outfile << "     ====> Interaction part" << std::endl;
  outfile << "Nucleon-Nucleon interaction model: n3lo" << std::endl;
  outfile << "Type of calculation: vlowk" << std::endl;
  outfile << "Number and value of starting energies:   1" << std::endl;
  outfile << "  0.000000E+00" << std::endl;
  outfile << "Total number of twobody matx elements:         ";
  std::streampos CountLocation = outfile.tellp();
  outfile << "0                                                 " << std::endl;
//467032          78000         311032          78000" << std::endl;
  outfile << "Matrix elements with the following legend, NOTE no hw/A for Hcom, pipj and rirj" << std::endl;
  outfile << "      Tz      Parity    J2        a        b        c        d     <ab|V|cd>    <ab|0|cd>    <ab|0|cd>    <ab|0|cd>" << std::endl;

  int linecounter = 0;
  array<int,3> Tzcounter = {0,0,0};
  for ( auto& itmat : Op.TwoBody.MatEl )
  {
    int ch = itmat.first[0]; // assume ch_bra == ch_ket
    auto& matrix = itmat.second;
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    int J2 = tbc.J*2;
    int Tz = tbc.Tz;
    int parity = tbc.parity;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
      Ket& bra = tbc.GetKet(ibra);
      for (int iket=ibra;iket<nkets;++iket)
      {
        Ket& ket = tbc.GetKet(iket);
        double tbme = matrix(ibra,iket);
        outfile << setw(wint) << Tz << " " << setw(wint) <<parity << " " << setw(wint) <<J2 << " " << setw(wint) << bra.p+1 << " " << setw(wint) <<bra.q+1 << " " << setw(wint) <<ket.p+1 << " " << setw(wint) <<ket.q+1 <<  " " << setw(wdouble) << setprecision(dprec) <<tbme << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << std::endl;
        ++linecounter;
        ++Tzcounter[Tz+1];
      }
    }

  }
  outfile.seekp( CountLocation );
  outfile << linecounter << "    " << Tzcounter[0] << "    " << Tzcounter[1] << "    " << Tzcounter[2];

  return;
}


//void ReadWrite::WriteOneBody_Oslo( std::string filename, Operator& Op)
void ReadWrite::WriteOneBody_Simple( std::string filename, Operator& Op)
{
  std::ofstream outfile(filename);
  ModelSpace* modelspace = Op.GetModelSpace();
  int wint = 3;
  int wdouble = 10;
  int dprec = 6;
  outfile << fixed;
  outfile << setw(wint) << setprecision(dprec);
  if ( !outfile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble writing file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  outfile << "  a     b      <a|V|b> " << std::endl;
  int norb = modelspace->GetNumberOrbits();
  for ( int i=0;i<norb;++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    for (int j : Op.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      outfile << setw(wint) << i << "   " << setw(wint) << j << "   " << setw(wdouble) << Op.OneBody(i,j) << std::endl;
    }
  }

}


void ReadWrite::WriteOneBody_Oslo( std::string filename, Operator& Op)
{
  std::ofstream outfile(filename);
  ModelSpace* modelspace = Op.GetModelSpace();

  outfile << "   ----> Oscillator parameters, Model space and single-particle data" << std::endl;
  outfile << "Mass number A of chosen nucleus (important for CoM corrections):          " << modelspace->GetTargetMass() << std::endl;
  outfile << "Oscillator length and energy: " << setw(12) << setprecision(6) << scientific 
          << HBARC / sqrt( M_NUCLEON * modelspace->GetHbarOmega() )
          << "  " << setw(12) << setprecision(6) << scientific <<  modelspace->GetHbarOmega() << std::endl;
  outfile << " Min and max value of partial wave ang. mom    0   " << modelspace->GetEmax() << std::endl;
  outfile << " Max value of relative orb mom or cm orb mom,  l or L=  " << modelspace->GetEmax() << std::endl;
  outfile << " Max value of relative n:  "  << modelspace->GetEmax()/2 << std::endl;
  outfile << " Max value of 2*n + l+ cm 2*N +L for large space:  " << modelspace->GetEmax() << std::endl;
  outfile << " Max value of 2*n + l+ cm 2*N +L for model space:  " << modelspace->GetEmax() << std::endl;
  outfile << " Total number of single-particle orbits  " << modelspace->GetNumberOrbits()  << std::endl;
  outfile << "Legend:         n     l     2j   tz    2n+l  HO-energy     evalence     particle/hole  inside/outside" << std::endl;

  for (int i=0;i<modelspace->GetNumberOrbits(); ++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    std::string ph = oi.cvq > 0 ? "particle" : "hole";
    std::string io = oi.cvq == 1 ? "inside" : "outside";
    int e = 2*oi.n+oi.l;
//    double hw = modelspace->GetHbarOmega();
    double spe = Op.OneBody(i,i);
    char line[512];
/// Switching order here to make EKK work with the MBPT code
//    sprintf(line,"Number: %3d %5d %5d %5d %5d %5d    %13.6e  %13.6e  %8s  %8s\n", i+1, oi.n, oi.l, oi.j2, oi.tz2, e, (e+1.5)*hw, spe, ph.c_str(), io.c_str());
    sprintf(line,"Number: %3d %5d %5d %5d %5d %5d    %13.6e  %13.6e  %8s  %8s\n", i+1, oi.n, oi.l, oi.j2, oi.tz2, e, spe, 0.0, ph.c_str(), io.c_str());
    outfile << line;
  }

}

void ReadWrite::ReadBareTBME_Jason( std::string filename, Operator& Hbare)
{

  std::ifstream infile;
  char line[LINESIZE];
  int Tz2,Par,J2,a,b,c,d;
  int na,nb,nc,nd;
  int la,lb,lc,ld;
  int ja,jb,jc,jd;
//  double fbuf[3];
  double tbme;
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norbits = modelspace->GetNumberOrbits();
  File2N = filename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }

//  infile.getline(line,LINESIZE);

//  while (!strstr(line,"<ab|V|cd>") && !infile.eof()) // Skip lines until we see the header
  //for (int& i : modelspace->holes ) // skip spe's at the top
  for (int i=0;i<6;++i ) // skip spe's at the top
  {
     infile.getline(line,LINESIZE);
     std::cout << "skipping line: " << line << std::endl;
  }

  // read the file one line at a time
  //while ( infile >> Tz >> Par >> J2 >> a >> b >> c >> d >> tbme >> fbuf[0] >> fbuf[1] >> fbuf[2] )
  while ( infile >> na >> la >> ja >> nb >> lb >> jb >> nc >> lc >> jc >> nd >> ld >> jd >> Tz2 >> J2 >> tbme )
  {
     int tza = Tz2<0 ? -1 : 1;
     int tzb = Tz2<2 ? -1 : 1;
     int tzc = Tz2<0 ? -1 : 1;
     int tzd = Tz2<2 ? -1 : 1;
     a = modelspace->Index1(na,la,ja,tza);
     b = modelspace->Index1(nb,lb,jb,tzb);
     c = modelspace->Index1(nc,lc,jc,tzc);
     d = modelspace->Index1(nd,ld,jd,tzd);
     // if the matrix element is outside the model space, ignore it.
     if (a>=norbits or b>=norbits or c>=norbits or d>=norbits) continue;
     Par = (la+lb)%2;
//     Hbare.SetTBME(J2/2,Par,Tz2/2,a,b,c,d, tbme ); 
     Hbare.TwoBody.SetTBME(J2/2,Par,Tz2/2,a,b,c,d, tbme ); 

  }

  return;
}


/// Read two body matrix elements from file formatted for Petr Navratil's
/// no-core shell model code.
/// These are given as normalized, JT coupled matix elements.
/// Matrix elements corresponding to orbits outside the modelspace are ignored.
void ReadWrite::ReadBareTBME_Navratil( std::string filename, Operator& Hbare)
{
  std::ifstream infile(filename);  
  if ( !infile.good() )
  {
    cerr << "************************************" << std::endl
         << "**    Trouble reading file  !!!   **" << filename << std::endl
         << "************************************" << std::endl;
    goodstate = false;
    return;
  }
  infile.close();
  File2N = filename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();

  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    std::ifstream infile(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadBareTBME_Navratil_from_stream(zipstream, Hbare);
  }
  else
  {
    std::ifstream infile(filename);
    ReadBareTBME_Navratil_from_stream(infile, Hbare);
  }
}

void ReadWrite::ReadBareTBME_Navratil_from_stream( std::istream& infile, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
//  int emax = modelspace->GetEmax();
  int norb = modelspace->GetNumberOrbits();
//  int nljmax = norb/2;
//  int herm = Hbare.IsHermitian() ? 1 : -1 ;
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2;
     orbits_remap[nlj] = i;
  }

//  std::ifstream infile;
//  infile.open(filename);


  // first line contains some information
//  int total_number_elements;
//  int N1max, N12max;
//  double hw, srg_lambda;
  double hw;
  int A = Hbare.GetModelSpace()->GetTargetMass();
  char header[500];
  infile.getline(header,500);
  hw = Hbare.GetModelSpace()->GetHbarOmega();
//  infile >> total_number_elements >> N1max >> N12max >> hw >> srg_lambda;

  int J,T;
  int nlj1,nlj2,nlj3,nlj4;
  double trel, h_ho_rel, vcoul, vpn, vpp, vnn;
  // Read the TBME
  while( infile >> nlj1 >> nlj2 >> nlj3 >> nlj4 >> J >> T >> trel >> h_ho_rel >> vcoul >> vpn >> vpp >> vnn )
  {
    --nlj1;--nlj2;--nlj3;--nlj4;  // Fortran -> C indexing
    auto it_a = orbits_remap.find(nlj1);
    auto it_b = orbits_remap.find(nlj2);
    auto it_c = orbits_remap.find(nlj3);
    auto it_d = orbits_remap.find(nlj4);
    if (it_a==orbits_remap.end() or it_b==orbits_remap.end() or it_c==orbits_remap.end() or it_d==orbits_remap.end() ) continue;
    int a = orbits_remap[nlj1];
    int b = orbits_remap[nlj2];
    int c = orbits_remap[nlj3];
    int d = orbits_remap[nlj4];

  // from Petr: The Trel and the H_HO_rel must be scaled by (2/A) hbar Omega
    if (doCoM_corr)
    {
       double CoMcorr = trel * 2 * hw / A;
       vpn += CoMcorr;
       vpp += CoMcorr;
       vnn += CoMcorr;
    }
    
    Orbit oa = modelspace->GetOrbit(a);
    Orbit ob = modelspace->GetOrbit(b);
    int parity = (oa.l + ob.l) % 2;

    if (T==1)
    {
      Hbare.TwoBody.SetTBME(J,parity,-1,a,b,c,d,vpp);
      Hbare.TwoBody.SetTBME(J,parity,1,a+1,b+1,c+1,d+1,vnn);
    }

    if (std::abs(vpn)>1e-6)
    {
      Hbare.TwoBody.Set_pn_TBME_from_iso(J,T,0,a,b,c,d,vpn);
    }

  }
  
}



void ReadWrite::WriteTBME_Navratil( std::string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2;
//     orbits_remap[nlj] = i;
     orbits_remap[i]   = nlj+1;
//     orbits_remap[i+1] = nlj+1;
  }

  std::ofstream outfile(filename);

  double hw = modelspace->GetHbarOmega();
  hw = Hbare.GetModelSpace()->GetHbarOmega();
  double srg_lambda = 0;
  outfile << 0 << "    " << modelspace->GetEmax() << "    " << 2*modelspace->GetEmax() << "   " << hw << "     " << srg_lambda << std::endl;

  outfile << setiosflags(std::ios::fixed);

  double trel=0, h_ho_rel=0, vcoul=0;

  for (int a=0; a<norb; a+=2)
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b=0; b<=a; b+=2)
    {
      Orbit& ob = modelspace->GetOrbit(b);
      int jab_min = std::abs(oa.j2-ob.j2)/2;
      int jab_max = (oa.j2+ob.j2)/2;
      for (int c=0; c<=a; c+=2)
      {
        Orbit& oc = modelspace->GetOrbit(c);
        for (int d=0; d<=c; d+=2)
        {
          Orbit& od = modelspace->GetOrbit(d);
          if ( (oa.l + ob.l + oc.l + od.l)%2 > 0) continue;
          int jcd_min = std::abs(oc.j2-od.j2)/2;
          int jcd_max = (oc.j2+od.j2)/2;
          int jmin = std::max(jab_min,jcd_min);
          int jmax = std::min(jab_max,jcd_max);
          for (int J=jmin; J<=jmax; ++J)
          {
//            double vpp =  Hbare.TwoBody.Get_iso_TBME_from_pn(J, 1, -1, a, b, c, d);
//            double vnn =  Hbare.TwoBody.Get_iso_TBME_from_pn(J, 1,  1, a, b, c, d);
            double vpp =  Hbare.TwoBody.GetTBME_J(J, a, b, c, d);
            double vnn =  Hbare.TwoBody.GetTBME_J(J, a+1, b+1, c+1, d+1);
            double v10 =  Hbare.TwoBody.Get_iso_TBME_from_pn(J, 1,  0, a, b, c, d);
            double v00 =  Hbare.TwoBody.Get_iso_TBME_from_pn(J, 0,  0, a, b, c, d);
            if (a==b)
            {
              vpp /= SQRT2;
              vnn /= SQRT2;
            }
            if (c==d)
            {
              vpp /= SQRT2;
              vnn /= SQRT2;
            }
            if (std::abs(vpp)>1e-7 or std::abs(vnn)>1e-7 or std::abs(v10)>1e-7)
            {
            outfile << setw(3) << orbits_remap.at(a) << " "
                    << setw(3) << orbits_remap.at(b) << " "
                    << setw(3) << orbits_remap.at(c) << " "
                    << setw(3) << orbits_remap.at(d) << " "
                    << setw(3) << J << "    " << 1 << " "
                    << setw(10) << setprecision(6)
                    << trel << " "  << setw(10) << setprecision(6)    << h_ho_rel << " "  << setw(10) << setprecision(6)  << vcoul << " "
                    << setw(10) << setprecision(6)
                    << v10 << " "
                    << setw(10) << setprecision(6)
                    << vpp << " "
                    << setw(10) << setprecision(6)
                    << vnn << std::endl;
            }
            if (std::abs(v00)>1e-7)
            {
            outfile << setw(3) << orbits_remap.at(a) << " "
                    << setw(3) << orbits_remap.at(b) << " "
                    << setw(3) << orbits_remap.at(c) << " "
                    << setw(3) << orbits_remap.at(d) << " "
                    << setw(3) << J << "    " << 0 << " "
                    << setw(10) << setprecision(6)
                    << trel << " "  << setw(10) << setprecision(6)    << h_ho_rel << " "  << setw(10) << setprecision(6)  << vcoul << " "
                    << setw(10) << setprecision(6)
                    << v00 << " "
                    << setw(10) << setprecision(6)
                    << 0.0 << " "
                    << setw(10) << setprecision(6)
                    << 0.0 << std::endl;
            }
          }
        }
      }
    }
  }
}








/// Decide if the file is gzipped or ascii, create a stream, then call ReadBareTBME_Darmstadt_from_stream().
void ReadWrite::ReadBareTBME_Darmstadt( std::string filename, Operator& Hbare, int emax, int Emax, int lmax)
{

  File2N = filename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();
  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    std::ifstream infile(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadBareTBME_Darmstadt_from_stream(zipstream, Hbare,  emax, Emax, lmax);
  }
  else if (filename.substr( filename.find_last_of(".")) == ".bin")
  {
    std::ifstream infile(filename, std::ios::binary);    
    if ( !infile.good() )
    {
      cerr << "problem opening " << filename << ". Exiting." << std::endl;
      return ;
    }


    infile.seekg(0, infile.end);
    int n_elem = infile.tellg();
    infile.seekg(0, infile.beg);
    n_elem -= infile.tellg();
    n_elem -= HEADERSIZE;
    n_elem /= sizeof(float);
    
    char header[HEADERSIZE];
//    infile.read((char*)&n_elem,sizeof(int));
    infile.read(header,HEADERSIZE);
    
//    std::vector<double> v(n_elem);
//    infile.read((char*)&v[0], n_elem*sizeof(double));
    std::vector<float> v(n_elem);
    infile.read((char*)&v[0], n_elem*sizeof(float));
    infile.close();
    VectorStream vecstream(v);
    std::cout << "n_elem = " << n_elem << std::endl;
    ReadBareTBME_Darmstadt_from_stream(vecstream, Hbare,  emax, Emax, lmax);
  }
  else
  {
    std::ifstream infile(filename);
    ReadBareTBME_Darmstadt_from_stream(infile, Hbare,  emax, Emax, lmax);
  }
}


/// Decide the file format from the extension -- .me3j (Darmstadt group format, human-readable), .gz (gzipped me3j, less storage),
/// .bin (me3j converted to binary, faster to read), .h5 (HDF5 format). Default is to assume .me3j.
/// For the first three, the file is converted to a stream and sent to ReadDarmstadt_3body_from_stream().
/// For the HDF5 format, a separate function is called: Read3bodyHDF5().
void ReadWrite::Read_Darmstadt_3body( std::string filename, Operator& Hbare, int E1max, int E2max, int E3max)
{

  double start_time = omp_get_wtime();
  std::string extension = filename.substr( filename.find_last_of("."));
  File3N = filename;
  Aref = Hbare.GetModelSpace()->GetAref();
  Zref = Hbare.GetModelSpace()->GetZref();

  if (extension == ".me3j")
  {
    std::ifstream infile(filename);
    Read_Darmstadt_3body_from_stream(infile, Hbare,  E1max, E2max, E3max);
  }
  else if ( extension == ".gz")
  {
    std::ifstream infile(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    Read_Darmstadt_3body_from_stream(zipstream, Hbare,  E1max, E2max, E3max);
  }
  else if (extension == ".bin")
  {
    std::ifstream infile(filename, std::ios::binary);    
    if ( !infile.good() )
    {
      cerr << "problem opening " << filename << ". Exiting." << std::endl;
      return ;
    }
    

    infile.seekg(0, infile.end);
//    int n_elem = infile.tellg();
    size_t n_elem = infile.tellg();
    infile.seekg(0, infile.beg);
    n_elem -= infile.tellg();
    n_elem -= HEADERSIZE;
    
    char header[HEADERSIZE];

    n_elem /= sizeof(float);

//    infile.read((char*)&n_elem,sizeof(int));
    infile.read(header,HEADERSIZE);
    
//    std::vector<float> v(n_elem);
//    infile.read((char*)&v[0], n_elem*sizeof(float));
    std::vector<float> v(n_elem);
    infile.read((char*)&v[0], n_elem*sizeof(float));
    infile.close();
    VectorStream vecstream(v);
    std::cout << "n_elem = " << n_elem <<  std::endl;
    Read_Darmstadt_3body_from_stream(vecstream, Hbare,  E1max, E2max, E3max);
  }
  else if (extension == ".h5")
  {
#ifndef NO_HDF5
    Read3bodyHDF5_new( filename, Hbare );
#else
    std::cout << "!!!!! ERROR: COMPILED WITH FLAG NO_HDF5 -> Can't read .h5 files !!!!!!!!!!!" << std::endl;
    exit(EXIT_FAILURE);
#endif
  }
  else
  {
    std::cout << "assuming " << filename << " is of me3j format ... " << std::endl;
    std::ifstream infile(filename);
    Read_Darmstadt_3body_from_stream(infile, Hbare,  E1max, E2max, E3max);
  }

  Hbare.profiler.timer["Read_3body_file"] += omp_get_wtime() - start_time;
}




/// Read TBME's from a file formatted by the Darmstadt group.
/// The file contains just the matrix elements, and the corresponding quantum numbers
/// are inferred. This means that the model space of the file must also be specified.
/// emax refers to the maximum single-particle oscillator shell. Emax refers to the
/// maximum of the sum of two single particles. lmax refers to the maximum single-particle \f$ \ell \f$.
/// If Emax is not specified, it should be 2 \f$\times\f$ emax. If lmax is not specified, it should be emax.
/// Also note that the matrix elements are given in un-normalized form, i.e. they are just
/// the M scheme matrix elements multiplied by Clebsch-Gordan coefficients for JT coupling.
//void ReadWrite::ReadBareTBME_Darmstadt_from_stream( istream& infile, Operator& Hbare, int emax, int Emax, int lmax)
template<class T>
void ReadWrite::ReadBareTBME_Darmstadt_from_stream( T& infile, Operator& Hbare, int emax, int Emax, int lmax)
{
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::vector<int> orbits_remap;

  if (emax < 0)  emax = modelspace->Emax;
  if (lmax < 0)  lmax = emax;

  for (int e=0; e<=min(emax,modelspace->Emax); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size()-1;

//  double tbme_pp,tbme_nn,tbme_10,tbme_00;
  float tbme_pp,tbme_nn,tbme_10,tbme_00;
  // skip the first line
  char line[LINESIZE];
  infile.getline(line,LINESIZE);

  for(int nlj1=0; nlj1<=nljmax; ++nlj1)
  {
    int a =  orbits_remap[nlj1];
    Orbit & o1 = modelspace->GetOrbit(a);
    int e1 = 2*o1.n + o1.l;
    if (e1 > modelspace->Emax) break;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_remap[nlj2];
      Orbit & o2 = modelspace->GetOrbit(b);
      int e2 = 2*o2.n + o2.l;
      if (e1+e2 > Emax) break;
      int parity = (o1.l + o2.l) % 2;

      for(int nlj3=0; nlj3<=nlj1; ++nlj3)
      {
        int c =  orbits_remap[nlj3];
        Orbit & o3 = modelspace->GetOrbit(c);
        int e3 = 2*o3.n + o3.l;

        for(int nlj4=0; nlj4<=(nlj3==nlj1 ? nlj2 : nlj3); ++nlj4)
        {
          int d =  orbits_remap[nlj4];
          Orbit & o4 = modelspace->GetOrbit(d);
          int e4 = 2*o4.n + o4.l;
          if (e3+e4 > Emax) break;
          if ( (o1.l + o2.l + o3.l + o4.l)%2 != 0) continue;
          int Jmin = std::max( std::abs(o1.j2 - o2.j2), std::abs(o3.j2 - o4.j2) )/2;
          int Jmax = std::min (o1.j2 + o2.j2, o3.j2+o4.j2)/2;
          if (Jmin > Jmax) continue;
          for (int J=Jmin; J<=Jmax; ++J)
          {

             // File is read here.
             // Matrix elements are written in the file with (T,Tz) = (0,0) (1,1) (1,0) (1,-1)
             infile >> tbme_00 >> tbme_nn >> tbme_10 >> tbme_pp;

             if (a>=norb or b>=norb or c>=norb or d>=norb) continue;

             // Normalization. The TBMEs are read in un-normalized.
             double norm_factor = 1;
             if (a==b)  norm_factor /= SQRT2;
             if (c==d)  norm_factor /= SQRT2;

//             std::cout << a << " " << b << " " << c << " " << d << "   " << J << "   "
//                  << fixed << setprecision(7) << setw(11) << tbme_00 << " " << tbme_nn << " " << tbme_10 << " " << tbme_pp << std::endl;

             if (norm_factor>0.9 or J%2==0)
             {
                Hbare.TwoBody.SetTBME(J,parity,-1,a,b,c,d,tbme_pp*norm_factor);
                Hbare.TwoBody.SetTBME(J,parity,1,a+1,b+1,c+1,d+1,tbme_nn*norm_factor);
                Hbare.TwoBody.Set_pn_TBME_from_iso(J,1,0,a,b,c,d,tbme_10*norm_factor);
             }
             if (norm_factor>0.9 or J%2!=0)
             { 
                Hbare.TwoBody.Set_pn_TBME_from_iso(J,0,0,a,b,c,d,tbme_00*norm_factor);
             }

          }
        }
      }
    }
  }

}


///// Read me3j format three-body matrix elements. Pass in E1max, E2max, E3max for the file, so that it can be properly interpreted.
///// The modelspace truncation doesn't need to coincide with the file truncation. For example, you could have an emax=10 modelspace
///// and read from an emax=14 file, and the matrix elements with emax>10 would be ignored.
//template <class T>
//void ReadWrite::Read_Darmstadt_3body_from_stream( T& infile, Operator& Hbare, int E1max, int E2max, int E3max)
//{
//  if ( !infile.good() )
//  {
//     cerr << "************************************" << std::endl
//          << "**    Trouble reading file  !!!   **" << std::endl
//          << "************************************" << std::endl;
//     goodstate = false;
//     return;
//  }
//  if (Hbare.particle_rank < 3)
//  {
//    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
//    cerr << " Oops. Looks like we're trying to read 3body matrix elements to a " << Hbare.particle_rank << "-body operator. For shame..." << std::endl;
//    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
//    goodstate = false;
//    return;
//  }
//  ModelSpace * modelspace = Hbare.GetModelSpace();
//  int e1max = modelspace->GetEmax();
//  int e2max = modelspace->GetE2max(); // not used yet
//  int e3max = modelspace->GetE3max();
//  int lmax3 = modelspace->GetLmax3();
//  std::cout << "Reading 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << std::endl;
//
//  std::vector<int> orbits_remap(0);
//  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.
//
//  for (int e=0; e<=min(E1max,e1max); ++e)
//  {
//    int lmin = e%2;
//    for (int l=lmin; l<=min(e,lmax); l+=2)
//    {
//      int n = (e-l)/2;
//      int twojMin = std::abs(2*l-1);
//      int twojMax = 2*l+1;
//      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
//      {
//         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
//      }
//    }
//  }
//  int nljmax = orbits_remap.size();
//
//
//  // skip the first line
//  char line[LINESIZE];
//  infile.getline(line,LINESIZE);
//
////  std::cout << "In read3body nthreads = " << omp_get_num_threads() << std::endl;
////  omp_set_num_threads(min(2,omp_get_max_threads())); // it's not clear that this is actually helping...
//
//  // begin giant nested loops
//  int nread = 0;
//  int nkept = 0;
//  for(int nlj1=0; nlj1<nljmax; ++nlj1)
//  {
//    int a =  orbits_remap[nlj1];
//    Orbit & oa = modelspace->GetOrbit(a);
//    int ea = 2*oa.n + oa.l;
//    if (ea > E1max) break;
//    if (ea > e1max) break;
//    if (ea > e3max) break;
////    std::cout << setw(5) << setprecision(2) << nlj1*(nlj1+1.)/(nljmax*(nljmax+1))*100 << " % done" << '\r';
////    std::cout.flush();
//
//    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
//    {
//      int b =  orbits_remap[nlj2];
//      Orbit & ob = modelspace->GetOrbit(b);
//      int eb = 2*ob.n + ob.l;
//      if ( (ea+eb) > E2max) break;
//
//      for(int nlj3=0; nlj3<=nlj2; ++nlj3)
//      {
//        int c =  orbits_remap[nlj3];
//        Orbit & oc = modelspace->GetOrbit(c);
//        int ec = 2*oc.n + oc.l;
//        if ( (ea+eb+ec) > E3max) break;
//
//        // Get J limits for bra <abc|
//        int JabMax  = (oa.j2 + ob.j2)/2;
//        int JabMin  = std::abs(oa.j2 - ob.j2)/2;
//
//        int twoJCMindownbra;
//        if (std::abs(oa.j2 - ob.j2) >oc.j2)
//           twoJCMindownbra = std::abs(oa.j2 - ob.j2)-oc.j2;
//        else if (oc.j2 < (oa.j2+ob.j2) )
//           twoJCMindownbra = 1;
//        else
//           twoJCMindownbra = oc.j2 - oa.j2 - ob.j2;
//        int twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2;
//
//
//        // now loop over possible ket orbits
//        for(int nnlj1=0; nnlj1<=nlj1; ++nnlj1)
//        {
//          int d =  orbits_remap[nnlj1];
//          Orbit & od = modelspace->GetOrbit(d);
//          int ed = 2*od.n + od.l;
//
//          for(int nnlj2=0; nnlj2 <= ((nlj1 == nnlj1) ? nlj2 : nnlj1); ++nnlj2)
//          {
//            int e =  orbits_remap[nnlj2];
//            Orbit & oe = modelspace->GetOrbit(e);
//            int ee = 2*oe.n + oe.l;
//
//            int nnlj3Max = (nlj1 == nnlj1 and nlj2 == nnlj2) ? nlj3 : nnlj2;
//            for(int nnlj3=0; nnlj3 <= nnlj3Max; ++nnlj3)
//            {
//              int f =  orbits_remap[nnlj3];
//              Orbit & of = modelspace->GetOrbit(f);
//              int ef = 2*of.n + of.l;
//              if ( (ed+ee+ef) > E3max) break;
//              // check parity
//              if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2 !=0 ) continue;
//
//              // Get J limits for ket |def>
//              int JJabMax = (od.j2 + oe.j2)/2;
//              int JJabMin = std::abs(od.j2 - oe.j2)/2;
//
//              int twoJCMindownket;
//              if ( std::abs(od.j2 - oe.j2) > of.j2 )
//                 twoJCMindownket = std::abs(od.j2 - oe.j2) - of.j2;
//              else if ( of.j2 < (od.j2+oe.j2) )
//                 twoJCMindownket = 1;
//              else
//                 twoJCMindownket = of.j2 - od.j2 - oe.j2;
//
//              int twoJCMaxupket = od.j2 + oe.j2 + of.j2;
//
//              int twoJCMindown = max(twoJCMindownbra, twoJCMindownket);
//              int twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket);
//              if (twoJCMindown > twoJCMaxup) continue;
//
//              //inner loops
//              for(int Jab = JabMin; Jab <= JabMax; Jab++)
//              {
//               for(int JJab = JJabMin; JJab <= JJabMax; JJab++)
//               {
//                //summation bounds for twoJC
//                int twoJCMin = max( std::abs(2*Jab - oc.j2), std::abs(2*JJab - of.j2));
//                int twoJCMax = min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
//       
//                // read all the ME for this range of J,T into block
//                if (twoJCMin>twoJCMax) continue;
//                size_t blocksize = ((twoJCMax-twoJCMin)/2+1)*5;
////                std::cout << "constructing block of size " << blocksize << "  =5* ((" << twoJCMax << " - " << twoJCMin << ")/2+1)" << std::endl;
//                std::vector<float> block(blocksize,0);
//                for (size_t iblock=0;iblock<blocksize;iblock++) infile >> block[iblock];
//                nread += blocksize;
////                std::cout << "Done making block" << std::endl;
//
////                for(int twoJC = twoJCMin; twoJC <= twoJCMax; twoJC += 2)
//                // note that the maximum number of threads is set at the beginning of the nested loops
//                double t_start_loop = omp_get_wtime();
////                #pragma omp parallel for schedule(dynamic,1) // parallelize in the J loop because they can't interfere with each other
//                for(int JTind = 0; JTind <= (twoJCMax-twoJCMin)+1; JTind++)
//                {
// //                std::cout << "   Read 3body threadnum = " << omp_get_thread_num() << std::endl;
//                 int twoJC = twoJCMin + (JTind/2)*2;
//                 int twoT = 1+(JTind%2)*2;
//                 for(int tab = 0; tab <= 1; tab++) // the total isospin loop can be replaced by i+=5
//                 {
//                  for(int ttab = 0; ttab <= 1; ttab++)
//                  {
//                   //summation bounds
//                   if ( twoT > min( 2*tab+1, 2*ttab+1) ) continue;
////                   int twoTMin = 1; // twoTMin can just be used as 1
////                   int twoTMax = min( 2*tab +1, 2*ttab +1);
//       
////                   for(int twoT = twoTMin; twoT <= twoTMax; twoT += 2)
////                   {
//////                    double V;
////                    float V = 0;
////                    infile >> V;
//                    float V = block[5*(twoJC-twoJCMin)/2+2*tab+ttab+(twoT-1)/2];
////                    ++nread;
//                    bool autozero = false;
//                    if (oa.l>lmax3 or ob.l>lmax3 or oc.l>lmax3 or od.l>lmax3 or oe.l>lmax3 or of.l>lmax3) V=0;
//
//                    
//
//                    if ( a==b and (tab+Jab)%2==0 ) autozero = true;
//                    if ( d==e and (ttab+JJab)%2==0 ) autozero = true;
//                    if ( a==b and a==c and twoT==3 and oa.j2<3 ) autozero = true;
//                    if ( d==e and d==f and twoT==3 and od.j2<3 ) autozero = true;
//
//                       if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
//                          and (ea+eb+ec<=e3max) and (ed+ee+ef<=e3max) )
//                       {
//                         #pragma omp atomic
//                         ++nkept;
//                       }
//
//                    if (not autozero and std::abs(V)>1e-5)
//                    {
////                       double V0 = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
////                       V0 = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
//                       if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
//                          and (ea+eb+ec<=e3max) and (ed+ee+ef<=e3max) )
//                       {
////                        std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << Jab << " " << JJab << " " << twoJC << " " << tab << " " << ttab << " " << twoT << " " << V << std::endl;
//                        Hbare.ThreeBody.SetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f, V);
//                       }
//                    }
//
//                    if (autozero)
//                    {
//                       if (std::abs(V) > 1e-6 and ea<=e1max and eb<=e1max and ec<=e1max)
//                       {
//                          std::cout << " <-------- AAAAHHHH!!!!!!!! Reading 3body file and this should be zero, but it's " << V << std::endl;
//                          goodstate = false;
//                       }
//                    }
//       
////                   }//twoT
//                  }//ttab
//                 }//tab
//                }//twoJ
//                modelspace->profiler.timer["Read_3body_inner_loop"] += omp_get_wtime() - t_start_loop;
//                if (not goodstate or not infile.good()) return;
//               }//JJab
//              }//Jab
//
//            }
//          }
//        }
//      }
//    }
//  }
//  std::cout << "Read in " << nread << " floating point numbers (" << nread * sizeof(float)/1024./1024./1024. << " GB)" << std::endl;
//  std::cout << "Stored " << nkept << " floating point numbers (" << nkept * sizeof(float)/1024./1024./1024. << " GB)" << std::endl;
//
//}





/// Read me3j format three-body matrix elements. Pass in E1max, E2max, E3max for the file, so that it can be properly interpreted.
/// The modelspace truncation doesn't need to coincide with the file truncation. For example, you could have an emax=10 modelspace
/// and read from an emax=14 file, and the matrix elements with emax>10 would be ignored.
size_t ReadWrite::Count_Darmstadt_3body_to_read( Operator& Hbare, int E1max, int E2max, int E3max, std::vector<int>& orbits_remap, std::vector<size_t>& nread_list)
{
  double t_start = omp_get_wtime();
//  if ( !infile.good() )
//  {
//     cerr << "************************************" << std::endl
//          << "**    Trouble reading file  !!!   **" << std::endl
//          << "************************************" << std::endl;
//     goodstate = false;
//     return;
//  }
//  if (Hbare.particle_rank < 3)
//  {
//    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
//    cerr << " Oops. Looks like we're trying to read 3body matrix elements to a " << Hbare.particle_rank << "-body operator. For shame..." << std::endl;
//    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
//    goodstate = false;
//    return;
//  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int e1max = modelspace->GetEmax();
  int e2max = modelspace->GetE2max(); // not used yet
  int e3max = modelspace->GetE3max();
//  int lmax3 = modelspace->GetLmax3();
  std::cout << "Reading 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << std::endl;

//  std::vector<int> orbits_remap(0);
  orbits_remap.clear();
//  nread_list.clear();
  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.

  for (int e=0; e<=min(E1max,e1max); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size();


  // skip the first line

  // begin giant nested loops
  size_t nread = 0;
//  int nkept = 0;

//  std::vector<size_t> nread_list;

  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
    nread_list.push_back(nread);
    int a =  orbits_remap[nlj1];
    Orbit & oa = modelspace->GetOrbit(a);
    int ea = 2*oa.n + oa.l;
    if (ea > E1max) break;
    if (ea > e1max) break;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_remap[nlj2];
      Orbit & ob = modelspace->GetOrbit(b);
      int eb = 2*ob.n + ob.l;
//      if ( ((size_t) (nlj1)*(nlj1+1)/2 + nlj2) != nread_list.size())
//      {
//        std::cout << " woops, snafu in constructing nread_list.   nlj1 = " << nlj1 << " nlj2 = " << nlj2 << "  nread_list.size() = " << nread_list.size() << "   nread = " << nread << std::endl;
//      }
//      nread_list.push_back(nread);
      if ( (ea+eb) > E2max) break;

      for(int nlj3=0; nlj3<=nlj2; ++nlj3)
      {
        int c =  orbits_remap[nlj3];
        Orbit & oc = modelspace->GetOrbit(c);
        int ec = 2*oc.n + oc.l;
        if ( (ea+eb+ec) > E3max) break;

        // Get J limits for bra <abc|
        int JabMax  = (oa.j2 + ob.j2)/2;
        int JabMin  = std::abs(oa.j2 - ob.j2)/2;

        int twoJCMindownbra;
        if (std::abs(oa.j2 - ob.j2) >oc.j2)
           twoJCMindownbra = std::abs(oa.j2 - ob.j2)-oc.j2;
        else if (oc.j2 < (oa.j2+ob.j2) )
           twoJCMindownbra = 1;
        else
           twoJCMindownbra = oc.j2 - oa.j2 - ob.j2;
        int twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2;


        // now loop over possible ket orbits
        for(int nnlj1=0; nnlj1<=nlj1; ++nnlj1)
        {
          int d =  orbits_remap[nnlj1];
          Orbit & od = modelspace->GetOrbit(d);
          int ed = 2*od.n + od.l;

          for(int nnlj2=0; nnlj2 <= ((nlj1 == nnlj1) ? nlj2 : nnlj1); ++nnlj2)
          {
            int e =  orbits_remap[nnlj2];
            Orbit & oe = modelspace->GetOrbit(e);
            int ee = 2*oe.n + oe.l;

            int nnlj3Max = (nlj1 == nnlj1 and nlj2 == nnlj2) ? nlj3 : nnlj2;
            for(int nnlj3=0; nnlj3 <= nnlj3Max; ++nnlj3)
            {
              int f =  orbits_remap[nnlj3];
              Orbit & of = modelspace->GetOrbit(f);
              int ef = 2*of.n + of.l;
              if ( (ed+ee+ef) > E3max) break;
              // check parity
              if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2 !=0 ) continue;

              // Get J limits for ket |def>
              int JJabMax = (od.j2 + oe.j2)/2;
              int JJabMin = std::abs(od.j2 - oe.j2)/2;

              int twoJCMindownket;
              if ( std::abs(od.j2 - oe.j2) > of.j2 )
                 twoJCMindownket = std::abs(od.j2 - oe.j2) - of.j2;
              else if ( of.j2 < (od.j2+oe.j2) )
                 twoJCMindownket = 1;
              else
                 twoJCMindownket = of.j2 - od.j2 - oe.j2;

              int twoJCMaxupket = od.j2 + oe.j2 + of.j2;

              int twoJCMindown = max(twoJCMindownbra, twoJCMindownket);
              int twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket);
              if (twoJCMindown > twoJCMaxup) continue;

              //inner loops
              for(int Jab = JabMin; Jab <= JabMax; Jab++)
              {
               for(int JJab = JJabMin; JJab <= JJabMax; JJab++)
               {
                //summation bounds for twoJC
                int twoJCMin = std::max( std::abs(2*Jab - oc.j2), std::abs(2*JJab - of.j2));
                int twoJCMax = std::min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
       
                // read all the ME for this range of J,T into block
                if (twoJCMin>twoJCMax) continue;
                size_t blocksize = ((twoJCMax-twoJCMin)/2+1)*5;
                nread += blocksize;


//                if (not goodstate or not infile.good()) return;
               }//JJab
              }//Jab

            }
          }
        }
      }
    }
  }
  
  modelspace->profiler.timer["Count_3BME"] += omp_get_wtime() - t_start;
  return nread;
}

template <class T>
void ReadWrite::Read_Darmstadt_3body_from_stream( T& infile, Operator& Hbare, int E1max, int E2max, int E3max)
{

  double t_start = omp_get_wtime();
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  if (Hbare.particle_rank < 3)
  {
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
    cerr << " Oops. Looks like we're trying to read 3body matrix elements to a " << Hbare.particle_rank << "-body operator. For shame..." << std::endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
    goodstate = false;
    return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  std::vector<int> orbits_remap;
  std::vector<size_t> nread_list;
  size_t nread = Count_Darmstadt_3body_to_read( Hbare, E1max, E2max, E3max, orbits_remap, nread_list);




  std::vector<float> ThreeBME(nread,0.);

//  #define BUFFSIZE3N 1024*1000
  if (format3N == "me3j")
  {
    char line[LINESIZE];
    infile.getline(line,LINESIZE);  // read the header
//    char buff[BUFFSIZE3N];
//    size_t read_so_far = 0;
//    while (read_so_far < nread and infile.good())
//    {
//      infile.read(buff, BUFFSIZE3N);
//   //.... other stuff...   
//    }
    for (size_t i=0;i<nread;++i) infile >> ThreeBME[i];
  }
  else if (format3N == "navratil" or format3N == "Navratil")
  {
    uint32_t delimiter; // This is machine-dependent. This seems to work on the local cluster...
    float v;
    for (size_t i=0;i<nread;++i)
    {
       infile.read((char*)&delimiter, sizeof(delimiter));
       infile.read((char*)&v,         sizeof(v));
       infile.read((char*)&delimiter, sizeof(delimiter));
       ThreeBME[i] = v;
    }
  }

  modelspace->profiler.timer["Read_3BME"] += omp_get_wtime() - t_start;
  std::cout << "Read in " << nread << " floating point numbers (" << nread * sizeof(float)/1024./1024./1024. << " GB)" << std::endl;
  Store_Darmstadt_3body( ThreeBME, nread_list, orbits_remap, Hbare, E1max, E2max, E3max);

}



/// Read me3j format three-body matrix elements. Pass in E1max, E2max, E3max for the file, so that it can be properly interpreted.
/// The modelspace truncation doesn't need to coincide with the file truncation. For example, you could have an emax=10 modelspace
/// and read from an emax=14 file, and the matrix elements with emax>10 would be ignored.
void ReadWrite::Store_Darmstadt_3body( const std::vector<float>& ThreeBME, const std::vector<size_t>& nread_list, const std::vector<int>& orbits_remap, Operator& Hbare, int E1max, int E2max, int E3max)
{

  double t_start = omp_get_wtime();
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int e1max = modelspace->GetEmax();
//  int e2max = modelspace->GetE2max(); // not used yet
  int e3max = modelspace->GetE3max();
  int lmax3 = modelspace->GetLmax3();
//  std::cout << "Reading 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << std::endl;

//  std::vector<int> orbits_remap(0);
//  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.
//
//  for (int e=0; e<=min(E1max,e1max); ++e)
//  {
//    int lmin = e%2;
//    for (int l=lmin; l<=min(e,lmax); l+=2)
//    {
//      int n = (e-l)/2;
//      int twojMin = std::abs(2*l-1);
//      int twojMax = 2*l+1;
//      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
//      {
//         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
//      }
//    }
//  }
  int nljmax = orbits_remap.size();


  std::cout << "begin storing 3N matrix elements" << std::endl;
  // begin giant nested loops
  size_t nkept = 0;
  modelspace->PreCalculateSixJ(); // Get all the sixJ so we don't have to worry about threading issues
  // combine the first two loops into one to better utilize more threads
//  #pragma omp parallel for schedule(dynamic,1) reduction(+ : nkept)  
//  for (int index12=0; index12< nljmax*(nljmax+1)/2; ++index12)
//  {
  #pragma omp parallel for schedule(dynamic,1) reduction(+ : nkept)  
  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
//    int nlj1 = int( (sqrt(8*index12+1)-1)/2);
//    int nlj2 = index12 - nlj1*(nlj1+1)/2;
//    size_t nread = nread_list[index12];
    size_t nread = nread_list[nlj1];
    int a =  orbits_remap[nlj1];
    Orbit & oa = modelspace->GetOrbit(a);
    int ea = 2*oa.n + oa.l;
//    if (ea > E1max) break;
//    if (ea > e1max) break;
//    if (ea > e3max) break;
    if (ea > E1max) continue;
    if (ea > e1max) continue;
    if (ea > e3max) continue;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_remap[nlj2];
      Orbit & ob = modelspace->GetOrbit(b);
      int eb = 2*ob.n + ob.l;
      if ( (ea+eb) > E2max) break;
//      if ( (ea+eb) > E2max) continue;

      for(int nlj3=0; nlj3<=nlj2; ++nlj3)
      {
        int c =  orbits_remap[nlj3];
        Orbit & oc = modelspace->GetOrbit(c);
        int ec = 2*oc.n + oc.l;
        if ( (ea+eb+ec) > E3max) break;

        // Get J limits for bra <abc|
        int JabMax  = (oa.j2 + ob.j2)/2;
        int JabMin  = std::abs(oa.j2 - ob.j2)/2;

        int twoJCMindownbra;
        if (std::abs(oa.j2 - ob.j2) >oc.j2)
           twoJCMindownbra = std::abs(oa.j2 - ob.j2)-oc.j2;
        else if (oc.j2 < (oa.j2+ob.j2) )
           twoJCMindownbra = 1;
        else
           twoJCMindownbra = oc.j2 - oa.j2 - ob.j2;
        int twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2;


        // now loop over possible ket orbits
        for(int nnlj1=0; nnlj1<=nlj1; ++nnlj1)
        {
          int d =  orbits_remap[nnlj1];
          Orbit & od = modelspace->GetOrbit(d);
          int ed = 2*od.n + od.l;

          for(int nnlj2=0; nnlj2 <= ((nlj1 == nnlj1) ? nlj2 : nnlj1); ++nnlj2)
          {
            int e =  orbits_remap[nnlj2];
            Orbit & oe = modelspace->GetOrbit(e);
            int ee = 2*oe.n + oe.l;

            int nnlj3Max = (nlj1 == nnlj1 and nlj2 == nnlj2) ? nlj3 : nnlj2;
            for(int nnlj3=0; nnlj3 <= nnlj3Max; ++nnlj3)
            {
              int f =  orbits_remap[nnlj3];
              Orbit & of = modelspace->GetOrbit(f);
              int ef = 2*of.n + of.l;
              if ( (ed+ee+ef) > E3max) break;
              // check parity
              if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2 !=0 ) continue;

              // Get J limits for ket |def>
              int JJabMax = (od.j2 + oe.j2)/2;
              int JJabMin = std::abs(od.j2 - oe.j2)/2;

              int twoJCMindownket;
              if ( std::abs(od.j2 - oe.j2) > of.j2 )
                 twoJCMindownket = std::abs(od.j2 - oe.j2) - of.j2;
              else if ( of.j2 < (od.j2+oe.j2) )
                 twoJCMindownket = 1;
              else
                 twoJCMindownket = of.j2 - od.j2 - oe.j2;

              int twoJCMaxupket = od.j2 + oe.j2 + of.j2;

              int twoJCMindown = max(twoJCMindownbra, twoJCMindownket);
              int twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket);
              if (twoJCMindown > twoJCMaxup) continue;

              //inner loops
              for(int Jab = JabMin; Jab <= JabMax; Jab++)
              {
               for(int JJab = JJabMin; JJab <= JJabMax; JJab++)
               {
                //summation bounds for twoJC
                int twoJCMin = std::max( std::abs(2*Jab - oc.j2), std::abs(2*JJab - of.j2));
                int twoJCMax = std::min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
       
                // read all the ME for this range of J,T into block
                if (twoJCMin>twoJCMax) continue;
                size_t blocksize = ((twoJCMax-twoJCMin)/2+1)*5;

                for(int JTind = 0; JTind <= (twoJCMax-twoJCMin)+1; JTind++)
                {
                 int twoJC = twoJCMin + (JTind/2)*2;
                 int twoT = 1+(JTind%2)*2;
                 for(int tab = 0; tab <= 1; tab++) // the total isospin loop can be replaced by i+=5
                 {
                  for(int ttab = 0; ttab <= 1; ttab++)
                  {
                   //summation bounds
                   if ( twoT > min( 2*tab+1, 2*ttab+1) ) continue;
//                   int twoTMin = 1; // twoTMin can just be used as 1
//                   int twoTMax = min( 2*tab +1, 2*ttab +1);
       
                    size_t index_ab = 5*(twoJC-twoJCMin)/2+2*tab+ttab+(twoT-1)/2;
                    if (nread+index_ab >=ThreeBME.size())
                    {
                      std::cout << "OH NO!!! trying to access element " << nread << "+" << index_ab << " = " << nread+index_ab << "  which is >= "<< ThreeBME.size() << std::endl;
                    }
                    float V;
                    V = ThreeBME[nread + index_ab ];
                    bool autozero = false;
                    if (oa.l>lmax3 or ob.l>lmax3 or oc.l>lmax3 or od.l>lmax3 or oe.l>lmax3 or of.l>lmax3) V=0;


                    if ( ( a==b and (tab+Jab)%2==0 )
                      or ( d==e and (ttab+JJab)%2==0 )
                      or ( a==b and a==c and twoT==3 and oa.j2<3 )
                      or ( d==e and d==f and twoT==3 and od.j2<3 )) autozero = true;

                    if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
                          and (ea+eb+ec<=e3max) and (ed+ee+ef<=e3max) )
                    {
                      ++nkept;
                      if ( std::abs(V)>1e-6 )
                      {
                        if (not autozero )
                        {
//                            std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << Jab << " " << JJab << " " << twoJC << " " << tab << " " << ttab << " " << twoT << " " << V << std::endl;
                            Hbare.ThreeBody.SetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f, V);
                        }
                        else if (autozero)
                        {
                            printf(" <--------- AAAHHHH!!!!!! Reading 3body file. <%d %d %d  %d %d |V| %d %d %d  %d %d>_(%d %d) should be zero but its %f.  nread = %lu index_ab = %lu\n",a,b,c,Jab,tab,d,e,f,JJab,ttab,twoJC,twoT,V,nread,index_ab);
//                              std::cout << " <-------- AAAAHHHH!!!!!!!! Reading 3body file and this should be zero, but it's " << V << std::endl;
//                            std::cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << Jab << " " << JJab << " " << twoJC << " " << tab << " " << ttab << " " << twoT << " " << V << std::endl;
//                            std::cout << "nread = " << nread << "  index_ab = " << index_ab << std::endl;
                              goodstate = false;
                        }
                      }
                    }
       
                  }//ttab
                 }//tab
                }//twoJ
                nread += blocksize;
               }//JJab
              }//Jab

            } //nnlj3
          } //nnlj2
        } //nnlj1
      } //nlj3
    } //nlj2
  } //nlj1
  
  std::cout << "Stored " << nkept << " floating point numbers (" << nkept * sizeof(float)/1024./1024./1024. << " GB)" << std::endl;

  modelspace->profiler.timer["Store_3BME"] += omp_get_wtime() - t_start;
}






#ifndef NO_HDF5
/// Read three-body basis from HDF5 formatted file. This routine was ported to C++ from
/// a C routine by Heiko Hergert, with as little modification as possible.
void ReadWrite::GetHDF5Basis( ModelSpace* modelspace, std::string filename, std::vector<array<int,5>>& Basis)
{
  H5File file(filename, H5F_ACC_RDONLY);
  // The parameter alpha enumerates the different 3body states |abc> coupled to J12 and J (no isospin)
  DataSet basis = file.openDataSet("alphas");
  DataSpace basis_dspace = basis.getSpace();

  int nDim = basis_dspace.getSimpleExtentNdims();
  hsize_t iDim[6];
  int status = basis_dspace.getSimpleExtentDims(iDim,NULL);
  if (status != nDim)
  {
     cerr << "Error: Failed to read dataset dimensions!" << std::endl;
     return;
  }
  
  int alpha_max = iDim[0]; // alpha_max is the largest dimension
  for (int i=0;i<nDim;++i)
    alpha_max = max(alpha_max, int(iDim[i]));

  // Generate a 2d buffer in contiguous memory
  int** dbuf = new int*[iDim[0]];
  dbuf[0] = new int[iDim[0]*iDim[1]];
  for (hsize_t i=1;i<iDim[0];++i)
  {
    dbuf[i] = dbuf[i-1] + iDim[1];
  }

  basis.read(&dbuf[0][0], PredType::NATIVE_INT);

  Basis.resize(alpha_max+1);

  for( int alpha=1; alpha<=alpha_max; ++alpha)
  {
    int alphap = dbuf[alpha-1][0];
    int n1     = dbuf[alpha-1][1];
    int l1     = dbuf[alpha-1][2];
    int twoj1  = dbuf[alpha-1][3];
    int n2     = dbuf[alpha-1][4];
    int l2     = dbuf[alpha-1][5];
    int twoj2  = dbuf[alpha-1][6];
    int n3     = dbuf[alpha-1][7];
    int l3     = dbuf[alpha-1][8];
    int twoj3  = dbuf[alpha-1][9];
    int J12    = dbuf[alpha-1][10];
    int twoJ   = dbuf[alpha-1][11];

    int o1 = modelspace->GetOrbitIndex(n1,l1,twoj1,-1);
    int o2 = modelspace->GetOrbitIndex(n2,l2,twoj2,-1);
    int o3 = modelspace->GetOrbitIndex(n3,l3,twoj3,-1);
    if (alpha != alphap)
    {
      cerr << "Error. alpha != alphap " << std::endl;
      return;
    }

    Basis[alpha] = {o1,o2,o3,J12,twoJ};
    
  }

  delete[] dbuf[0];
  delete[] dbuf;

}
#endif

#ifndef NO_HDF5
/// Read three-body matrix elements from HDF5 formatted file. This routine was ported to C++ from
/// a C routine by Heiko Hergert, with as little modification as possible.
void ReadWrite::Read3bodyHDF5( std::string filename,Operator& op )
{

  const int SLABSIZE = 10000000;
  File3N = filename;
  Aref = op.GetModelSpace()->GetAref();
  Zref = op.GetModelSpace()->GetZref();

  ModelSpace* modelspace = op.GetModelSpace();
  std::vector<array<int,5>> Basis;
  GetHDF5Basis(modelspace, filename, Basis);


  H5File file(filename, H5F_ACC_RDONLY);
  DataSet label = file.openDataSet("vtnf_labels");
  DataSpace label_dspace = label.getSpace();
  DataSet value = file.openDataSet("vtnf");
  DataSpace value_dspace = value.getSpace();

  int label_nDim = label_dspace.getSimpleExtentNdims();
  if (label_nDim != 2)
  {
    cerr << "Error. Expected label_nDim==2, but got << label_nDim." << std::endl;
    return;
  }
  hsize_t label_maxDim[2];
  int label_status = label_dspace.getSimpleExtentDims(label_maxDim,NULL);
  if (label_status != label_nDim)
  {
    cerr << "Error. failed to read dataset dimensions for label." << std::endl;
    return;
  }
  
  hsize_t label_curDim[2];
  label_curDim[0] = min(SLABSIZE, int(label_maxDim[0]));
  label_curDim[1] = 7 ;
  
  DataSpace label_buf_dspace(2,label_curDim);

  // Generate a 2d buffer in contiguous memory
  int **label_buf = new int*[label_curDim[0]];
  label_buf[0] = new int[label_curDim[0] * label_curDim[1]];
  for (hsize_t i=1; i<label_curDim[0]; ++i)
  {
    label_buf[i] = label_buf[i-1] + label_curDim[1];  
  }
  
  int value_nDim = value_dspace.getSimpleExtentNdims();
  if (value_nDim != 2)
  {
    cerr << "Error. Expected value_nDim==2, but got << value_nDim." << std::endl;
    return;
  }
  hsize_t value_maxDim[2];
  int value_status = value_dspace.getSimpleExtentDims(value_maxDim,NULL);
  if (value_status != value_nDim)
  {
    cerr << "Error. failed to read dataset dimensions for value." << std::endl;
    return;
  }
  
  hsize_t value_curDim[2];
  value_curDim[0] = min(SLABSIZE, int(value_maxDim[0]));
  value_curDim[1] = 1 ;

  DataSpace value_buf_dspace(1,value_curDim);
  
  // Generate a 1d buffer in contiguous memory, also known as an array...
  double *value_buf = new double[value_curDim[0]];

  // break the file into slstd::abs for reading
  int nSlstd::abs = (int)((double)value_maxDim[0]/((double)SLABSIZE)) + 1;

  hsize_t stride[2] = {1,1};
  hsize_t count[2] = {1,1};

  // loop through the slstd::abs
  for ( int n=0; n<nSlstd::abs; ++n)
  {
    hsize_t start[2] = { n*value_curDim[0], 0};
    hsize_t label_block[2];
    hsize_t value_block[2];
    if (n==nSlstd::abs-1)
    {
      label_block[0] = label_maxDim[0]-(nSlstd::abs-1)*SLABSIZE;
      label_block[1] = label_maxDim[1];
      value_block[0] = value_maxDim[0]-(nSlstd::abs-1)*SLABSIZE;
      value_block[1] = value_maxDim[1];

      // Not clear exactly why this needs to be done.
      label_buf_dspace.close();
      value_buf_dspace.close();
      label_buf_dspace = DataSpace(2,label_block);
      value_buf_dspace = DataSpace(2,value_block);

    }
    else
    {
      label_block[0] = label_curDim[0];
      label_block[1] = label_curDim[1];
      value_block[0] = value_curDim[0];
      value_block[1] = value_curDim[1];
    }

    label_dspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, label_block);
    value_dspace.selectHyperslab( H5S_SELECT_SET, count, start, stride, value_block);

    // Read the label data into label_buf, and matrix elements into value_buf
    label.read( &label_buf[0][0], PredType::NATIVE_INT, label_buf_dspace, label_dspace );
    value.read( &value_buf[0], PredType::NATIVE_DOUBLE, value_buf_dspace, value_dspace );


    for (hsize_t i=0; i<value_block[0]; ++i)
    {
       int  alpha  = label_buf[i][0];
       int  T12    = label_buf[i][1]/2;
       int  alphap = label_buf[i][2];
       int  TT12   = label_buf[i][3]/2;
       int  twoT   = label_buf[i][4];
       int  twoJ   = label_buf[i][5];
       int  Pi     = label_buf[i][6];

       if (alpha<alphap) continue;

       double me   = value_buf[i];
       me *= HBARC;
//       if (alpha != alphap) me *=0.5;

       int a    = Basis[alpha][0];
       int b    = Basis[alpha][1];
       int c    = Basis[alpha][2];
       int J12  = Basis[alpha][3];
       int J2   = Basis[alpha][4];

       int d    = Basis[alphap][0];
       int e    = Basis[alphap][1];
       int f    = Basis[alphap][2];
       int JJ12 = Basis[alphap][3];
       int J2p  = Basis[alphap][4];

       int norb = modelspace->GetNumberOrbits();
       if (a>=norb or b>=norb or c>=norb or d>=norb or e>=norb or f>=norb) continue;

       Orbit& oa = modelspace->GetOrbit(a);
       Orbit& ob = modelspace->GetOrbit(b);
       Orbit& oc = modelspace->GetOrbit(c);
       Orbit& od = modelspace->GetOrbit(d);
       Orbit& oe = modelspace->GetOrbit(e);
       Orbit& of = modelspace->GetOrbit(f);
       int parity_abc = ( oa.l+ob.l+oc.l )%2;
       int parity_def = ( od.l+oe.l+of.l )%2;
       if (parity_abc != Pi or parity_def != Pi)
       {
         cerr << "Error. Mismatching parity !  < "  << parity_abc << " " << parity_def << " " << Pi << "    " << std::endl;
       }
       if (J2 != twoJ or J2p != twoJ)
       {
         cerr << "Error. Mismatching total J! " << J2 << " " << J2p << " " << twoJ << "   alphas = " << alpha << ", " << alphap << std::endl;
       }

       me *= 0.5; // According to Heiko, this shouldn't be here. But comparing matrix elements with Petr's interaction suggests otherwise.
       me *= modelspace->phase(oa.n+ob.n+oc.n+od.n+oe.n+of.n); // shamelessly copying Heiko. I don't understand this.

       op.ThreeBody.SetME(J12,JJ12,twoJ,T12,TT12,twoT,a,b,c,d,e,f, me);
       if (a==d and b==e and c==f and ( J12!=JJ12 ) )
          op.ThreeBody.SetME(JJ12,J12,twoJ,TT12,T12,twoT,a,b,c,d,e,f, me);

    } //loop over matrix elements
  } // loop over slstd::abs
  delete[] label_buf[0];
  delete[] label_buf;
  delete[] value_buf;
  std::cout << "Writing me3j file..." << std::endl;
  Write_me3j(filename + "_to_me3j", op, 2, 24, 12);
  std::cout << "done" << std::endl;
}
#endif





#ifndef NO_HDF5
// THIS ONE SEEMS TO WORK, SO FAR

void ReadWrite::Read3bodyHDF5_new( std::string filename,Operator& op )
{

  ModelSpace* modelspace = op.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  File3N = filename;
  Aref = op.GetModelSpace()->GetAref();
  Zref = op.GetModelSpace()->GetZref();

  int t12p_list[5] = {0,0,1,1,1};
  int t12_list[5]  = {0,1,0,1,1};
  int twoT_list[5] = {1,1,1,1,3};

  H5File file(filename, H5F_ACC_RDONLY);

  DataSet basis = file.openDataSet("alphas");
  DataSpace basis_dspace = basis.getSpace();
  hsize_t iDim_basis[6];
  basis_dspace.getSimpleExtentDims(iDim_basis,NULL);


  // Generate a 2d buffer in contiguous memory
  //TODO: Do this with std::vectors so we don't have to worry about new and delete[]
  int** dbuf = new int*[iDim_basis[0]];
  dbuf[0] = new int[iDim_basis[0]*iDim_basis[1]];
  for (hsize_t i=1;i<iDim_basis[0];++i)
  {
    dbuf[i] = dbuf[i-1] + iDim_basis[1];
  }

  basis.read(&dbuf[0][0], PredType::NATIVE_INT);


  DataSet value = file.openDataSet("vtnf");
  DataSpace value_dspace = value.getSpace();
  hsize_t value_maxDim[2];
  value_dspace.getSimpleExtentDims(value_maxDim,NULL);

  hsize_t value_curDim[2];
  value_curDim[0] = value_maxDim[0];
  value_curDim[1] = 5 ; // c1 c3 c4 cD cE
  
  // This needs to be a 2d array now
  float **value_buf = new float*[value_maxDim[0]];
  value_buf[0] = new float[value_curDim[0] * value_curDim[1]];
  for (hsize_t i=1; i<value_curDim[0]; ++i)
  {
    value_buf[i] = value_buf[i-1] + value_curDim[1];  
  }

  value.read(&value_buf[0][0], PredType::NATIVE_FLOAT);

  int alpha_max = iDim_basis[0];

//  int i=-5; 
  long long i=-5; 
  for (int alphaspp=0;alphaspp<alpha_max;++alphaspp)
  {
    int lap = dbuf[alphaspp][2];
    int lbp = dbuf[alphaspp][5];
    int lcp = dbuf[alphaspp][8];

    int ap = modelspace->GetOrbitIndex(dbuf[alphaspp][1],lap,dbuf[alphaspp][3],-1);
    int bp = modelspace->GetOrbitIndex(dbuf[alphaspp][4],lbp,dbuf[alphaspp][6],-1);
    int cp = modelspace->GetOrbitIndex(dbuf[alphaspp][7],lcp,dbuf[alphaspp][9],-1);
    int j12p = dbuf[alphaspp][10];
    int jtotp = dbuf[alphaspp][11];
    if (ap > norb) break;

    for (int alphasp=alphaspp; alphasp<alpha_max;++alphasp)
    {
      int la = dbuf[alphasp][2];
      int lb = dbuf[alphasp][5];
      int lc = dbuf[alphasp][8];
      int a = modelspace->GetOrbitIndex(dbuf[alphasp][1],la,dbuf[alphasp][3],-1);
      int b = modelspace->GetOrbitIndex(dbuf[alphasp][4],lb,dbuf[alphasp][6],-1);
      int c = modelspace->GetOrbitIndex(dbuf[alphasp][7],lc,dbuf[alphasp][9],-1);
      int j12 = dbuf[alphasp][10];
      int jtot = dbuf[alphasp][11];
      if (jtot != jtotp or (lap+lbp+lcp+la+lb+lc)%2>0) continue; 
      i+=5;
      if (ap>=norb or bp>=norb or cp>=norb) continue;
      if (a>=norb or b>=norb or c>=norb) continue;
      
      for (hsize_t k_iso=0;k_iso<5;++k_iso)
      {
       int T12  = t12p_list[k_iso];
       int TT12 = t12_list[k_iso];
       int twoT = twoT_list[k_iso];
       float *me = value_buf[i+k_iso];
       float summed_me = 0;
       for (int ii=0;ii<5;++ii) summed_me += LECs[ii] * me[ii] ;
       summed_me *= HBARC;
       // Phase due to different conventions for HO wave functions.
       // Now obsolete -- Feb 2016
//       summed_me *= modelspace->phase(dbuf[alphasp][1]+dbuf[alphasp][4]+dbuf[alphasp][7]+dbuf[alphaspp][1]+dbuf[alphaspp][4]+dbuf[alphaspp][7]);

       if ( (ap==bp and (j12p+T12)%2 !=1) or ( a==b  and (j12+TT12)%2 !=1 ) )
       {
         if ( std::abs(summed_me)>1.0e-6  )
         {
           std::cout << "AAHH!!  by J+T symmetry should be zero!" << std::endl;
         }
       }
       else
       {
//        std::cout << "read <" << ap << " " << bp << " " << cp << " | V | " << a << " " << b << " " << c << "  (" << j12p << " " << j12 << " " << jtot << ")  ( " << T12 << " " << TT12 << " " << twoT << std::endl;
        op.ThreeBody.SetME(j12p,j12,jtot,T12,TT12,twoT,ap,bp,cp,a,b,c, summed_me);
        if (a==ap and b==bp and c==cp and j12 != j12p) // we're only looping through alphap > alphaspp, while I'm set up to read in all J,T possibilities for a given set of orbits
        {
          op.ThreeBody.SetME(j12,j12p,jtot,TT12,T12,twoT,ap,bp,cp,a,b,c, summed_me);
        }
       }
      }

    }
  }

  delete[] dbuf[0];
  delete[] dbuf;
  delete[] value_buf[0];
  delete[] value_buf;

}
#endif




void ReadWrite::ReadOperator_Nathan( std::string filename1b, std::string filename2b, Operator& op)
{
  std::ifstream infile(filename1b);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  index_t a,b,c,d,J;
  double me;
  int herm = op.IsHermitian() ? 1 : -1;
  index_t norb = op.GetModelSpace()->GetNumberOrbits();
  while ( infile >> a >> b >> me  )
  {
    if (a>=norb or b>=norb) continue;
    op.OneBody(a,b) = me;
    op.OneBody(b,a) = herm*me;
  }
  infile.close();

  infile.open(filename2b);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
//  char header[500];
//  infile.getline(header,500);
  while ( infile >> a >> b >> c >> d >> J >> me )
  {
    if (a>=norb or b>=norb) continue;
    if (c>=norb or d>=norb) continue;
//    if (a==b) me /= sqrt(2);
//    if (c==d) me /= sqrt(2);
    if (a==b) me /= SQRT2;
    if (c==d) me /= SQRT2;
    op.TwoBody.SetTBME_J(J,a,b,c,d,me);
  }
  infile.close();


}


void ReadWrite::ReadTensorOperator_Nathan( std::string filename1b, std::string filename2b, Operator& op)
{
  std::ifstream infile(filename1b);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening 1b file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  index_t a,b,c,d,J1,J2;
  double me;
  ModelSpace* modelspace = op.GetModelSpace();
  int herm = op.IsHermitian() ? 1 : -1;
  index_t norb = op.GetModelSpace()->GetNumberOrbits();
  while ( infile >> a >> b >> me  )
  {
    if (a>=norb or b>=norb) continue;
    int flipphase = herm*modelspace->phase( (modelspace->GetOrbit(a).j2 - modelspace->GetOrbit(b).j2)/2);
    op.OneBody(a,b) = me;
    op.OneBody(b,a) = flipphase*me;
  }
  infile.close();
  std::cout << "Done reading 1b file." << std::endl;

  infile.open(filename2b);
  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening 2b file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
//  char header[500];
//  infile.getline(header,500);
  while ( infile >> a >> b >> c >> d >> J1 >> J2 >> me )
  {
    if (a>=norb or b>=norb) continue;
    if (c>=norb or d>=norb) continue;
    op.TwoBody.SetTBME_J(J1,J2,a,b,c,d,me);
  }
  infile.close();
  std::cout << "Done reading 2b file." << std::endl;


}

/*
void ReadWrite::Read2bCurrent_Navratil( std::string filename, Operator& Op)
{
   std::cout << "Begin Read2bCurrent_Navratil" << std::endl;
//  std::ifstream infile(filename);
    std::ifstream infilegz(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream infile;
    infile.push(boost::iostreams::gzip_decompressor());
    infile.push(infilegz);


  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening file  !!!   **  " << filename<< std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Op.GetModelSpace();
  std::unordered_map<int,int> orbits_remap;
  int norb = modelspace->GetNumberOrbits();
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2;
     orbits_remap[nlj] = i;
  }

  std::cout << "Done remapping orbits" << std::endl;

  // locally-defined isospin clebsch gordan coefficent
  // hopefully faster than calling AngMom::CG
  auto isospinCG = [](int tz1, int tz2, int T12){
    if (std::abs(tz1+tz2) > 2*T12) return 0.0;
    if (tz1==tz2) return (double)T12;
    if (T12==1) return 1/sqrt(2.0);
    return tz1 / sqrt(2.0);
  };

  std::string strbuf;
  int J_op,T_op,ain,bin,cin,din,a,b,c,d,j12,t12,j34,t34,tza,tzb;
  double mat_el;

  for (int i=0;i<5;++i) infile >> strbuf; // Five useless lines...
  infile >> J_op >> T_op;
  for (int i=0;i<8;++i) infile >> strbuf; // Another eight useless lines...

  std::cout << "Start while loop" << std::endl;
  while( infile >> ain >> bin >> cin >> din >> j12 >> j34 >> t12 >> t34 >> tza >> tzb >> mat_el)
  {
    if (cin<ain or (cin==ain and din<bin)) continue;
    if ((ain==cin and bin==din) and (j34<j12 or (j34==j12 and t34<t12))) continue;
   ain--; bin--;cin--;din--; // Fortran to C
    if (    orbits_remap.find(ain)==orbits_remap.end() or orbits_remap.find(bin)==orbits_remap.end()
         or orbits_remap.find(cin)==orbits_remap.end() or orbits_remap.find(din)==orbits_remap.end() ) continue;
//   if (true)
//   if (ain==0 and bin==0 and cin==0 and din==0 and j12==1 and j34==0)
//   if (ain==1 and bin==1 and cin==1 and din==2 and j12==0 and j34==1)
   if (ain==1 and bin==2 and cin==1 and din==2 and j12==1 and j34==1)
   {
    std::cout << "< " << ain << " " << bin << " " << j12 << " " << t12 << " ||| Op J=" << J_op << " T=" <<T_op << " ||| "
                 << cin << " " << din << " " << j34 << " " << t34 << " > = " << mat_el << std::endl;
   }
    a = orbits_remap[ain];  // remapping gives proton index. neutron is a+1.
    b = orbits_remap[bin];
    c = orbits_remap[cin];
    d = orbits_remap[din];
//    std::cout << "    --> " << a << " " << b << " " << c << " " << d << std::endl;
  //  mat_el = < a b j12 t12 ||| OP (J,T) ||| c d j34 t34>    (doubly-reduced)
  //   (note that tza and tzb from the file are not meaningful here and should be ignored).
  //   also note that Petr uses the particle-physics convention tz|p> = +1|p>
  //   which is the opposite of what is used elsewhere in this code.

    for ( int tza : {-1,1} )
    {
     for ( int tzb : {-1,1} )
     {
      if (ain==bin and tzb<tza) continue;
      if (std::abs(tza+tzb)>2*t12) continue;
      double cg12 = isospinCG(tza,tzb,t12);
      if (ain==bin and (tza!=tzb)) cg12 *= sqrt(2.0);
      for ( int tzc : {-1,1} )
      {
       for ( int tzd : {-1,1} )
       {
        if (cin==din and tzd>tzc) continue;
        if (std::abs(tzc+tzd)>2*t34) continue;
        double Tz_op = (tza+tzb-tzc-tzd)/2;
        if (std::abs(Tz_op)!=T_op) continue;
        double cg34 = isospinCG(tzc,tzd,t34);
        if (cin==din and (tzc!=tzd)) cg34 *= sqrt(2.0);
        if ( (ain==cin) and (bin==din) and (j12==j34))
        {
          if (tzc<tza or (tzc==tza and tzd<tzb) ) continue;
        }
        double WignerEckartCoeff = AngMom::CG(t34,(tzc+tzd)/2,T_op,Tz_op,t12,(tza+tzb)/2) / sqrt(2*t12+1.);
        if (Tz_op < 0) WignerEckartCoeff *= -1; // If tau is a tensor operator, there needs to a relative minus sign for t+ and t-
        // Note that this strategy could come back to bite me if the file
        // contains redundant information. There are safer but more memory-intensive ways to do this...
//        if (true  )
//        if (ain==0 and bin==0 and cin==0 and din==0 and ((j12==1 and j34==0)or(j12==0 and j34==1)) )
//        if (ain==1 and bin==1 and cin==1 and din==2 and j12==0 and j34==1)
        if (ain==1 and bin==2 and cin==1 and din==2 and j12==1 and j34==1)
        {
        std::cout << "tza,tzb,t12, tzc,tzd,t34, cg12,cg34,WE        : " << tza << " " << tzb << " " << t12 << ",   " << tzc << " " << tzd << " " << t34
             << "   " << cg12 << " " << cg34 << " " << WignerEckartCoeff << " ( from " << AngMom::CG(t34,(tzc+tzd)/2,T_op,Tz_op,t12,(tza+tzb)/2) << " / " <<sqrt(2*t12+1.) << ") " << std::endl;
        std::cout << "     cg12 is < 1/2 " << tza << "/2 1/2 " << tzb << "/2 | " << t12 << " " << (tza+tzb)/2 << " > " << std::endl;
        std::cout << "     cg34 is < 1/2 " << tzc << "/2 1/2 " << tzd << "/2 | " << t34 << " " << (tzc+tzd)/2 << " > " << std::endl;
        std::cout << "     W-E Clebsch is < " << t34 << " " << (tzc+tzd)/2 << "  " << T_op << " " << Tz_op << " | " << t12 << " " << (tza+tzb)/2 << " > " << std::endl;
        std::cout << a+(1-tza)/2<<  "   " << b+(1-tzb)/2<<  "   " << c+(1-tzc)/2<<  "   " << d+(1-tzd)/2 << std::endl;
        std::cout << "       Adding " << mat_el*cg12*cg34*WignerEckartCoeff
             << " to " << Op.TwoBody.GetTBME_J_norm(j12, j34, a+(1-tza)/2, b+(1-tzb)/2, c+(1-tzc)/2, d+(1-tzd)/2  ) << std::endl;
        }
        Op.TwoBody.AddToTBME_J( j12, j34, a+(1-tza)/2, b+(1-tzb)/2, c+(1-tzc)/2, d+(1-tzd)/2, mat_el*cg12*cg34*WignerEckartCoeff);
//        if (ain==0 and bin==0 and cin==0 and din==0 and j12==1 and j34==0)
//        if (ain==0 and bin==0 and cin==0 and din==0 and ((j12==1 and j34==0)or(j12==0 and j34==1)) )
//        if (ain==1 and bin==1 and cin==1 and din==2 and j12==0 and j34==1)
        if (ain==1 and bin==2 and cin==1 and din==2 and j12==1 and j34==1)
        {
          std::cout << "Afterwards, it's " << Op.TwoBody.GetTBME_J_norm(j12, j34, a+(1-tza)/2, b+(1-tzb)/2, c+(1-tzc)/2, d+(1-tzd)/2  ) << std::endl;
          std::cout << "$$$$ <2 2 J0 || 2 5 J1> = " << Op.TwoBody.GetTBME_J_norm(0,1,2,2,2,5)
               << "     <5 2 J1 || 2 2 J0> = " << Op.TwoBody.GetTBME_J_norm(1,0,5,2,2,2) << std::endl;
        }
       }
      }
     }
    }
  }
}



*/


uint64_t Petr2BC_hash(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t jab, uint64_t jcd, uint64_t tab, uint64_t tcd )
{
  return    a + (b <<8 ) + (c << 16) + (d << 24) + (jab << 32) + (jcd << 40) + (tab << 48) + (tcd << 56); 
}


void ReadWrite::Read2bCurrent_Navratil( std::string filename, Operator& Op)
{
    std::ifstream infilegz(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream infile;
    infile.push(boost::iostreams::gzip_decompressor());
    infile.push(infilegz);


  if ( !infile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening file  !!!   **  " << filename<< std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Op.GetModelSpace();
  std::unordered_map<int,int> orbits_remap;
  int norb = modelspace->GetNumberOrbits();
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2;
     orbits_remap[nlj] = i;
  }


  std::string strbuf;
  int J_op,T_op,ain,bin,cin,din,a,b,c,d,j12,t12,j34,t34,tza,tzb;
  double mat_el;

  map<uint64_t, double> DoubleReducedME;

  using namespace AngMom;

  for (int i=0;i<5;++i) infile >> strbuf; // Five useless lines...
  infile >> J_op >> T_op;
  for (int i=0;i<8;++i) infile >> strbuf; // Another eight useless lines...

  while( infile >> ain >> bin >> cin >> din >> j12 >> j34 >> t12 >> t34 >> tza >> tzb >> mat_el)
  {
    ain--; bin--;cin--;din--; // Fortran to C
    if (    orbits_remap.find(ain)==orbits_remap.end() or orbits_remap.find(bin)==orbits_remap.end()
         or orbits_remap.find(cin)==orbits_remap.end() or orbits_remap.find(din)==orbits_remap.end() ) continue;
    a = orbits_remap[ain];  // remapping gives proton index. neutron is a+1.
    b = orbits_remap[bin];
    c = orbits_remap[cin];
    d = orbits_remap[din];
    if (a > b)
    {
      Orbit& oa = modelspace->GetOrbit(a);
      Orbit& ob = modelspace->GetOrbit(b);
      swap(a,b);
      mat_el *= modelspace->phase( (oa.j2 + ob.j2)/2 - j12 - t12);
    }
    if (c > d)
    {
      Orbit& oc = modelspace->GetOrbit(c);
      Orbit& od = modelspace->GetOrbit(d);
      swap(c,d);
      mat_el *= modelspace->phase( (oc.j2 + od.j2)/2 - j34 - t34);
    }
    
    uint64_t hashkey = Petr2BC_hash(a,b,c,d,j12,j34,t12,t34);
    DoubleReducedME[hashkey] = mat_el;
  }

  for ( auto& itmat : Op.TwoBody.MatEl )
  {

    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];
    auto& TBME = itmat.second;
    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int Jbra = tbc_bra.J;
    int Jket = tbc_ket.J;
    int Tzbra = -tbc_bra.Tz;
    int Tzket = -tbc_ket.Tz;

    int nbra = tbc_bra.GetNumberKets();
    int nket = tbc_ket.GetNumberKets();

    for (int ibra=0; ibra<nbra; ++ibra)
    {
      Ket& bra = tbc_bra.GetKet(ibra);
      int a = bra.p;
      int b = bra.q;
      a -= a%2;
      b -= b%2;
      int tza = -bra.op->tz2;
      int tzb = -bra.oq->tz2;
      for (int iket=0; iket<nket; ++iket)
      {
        Ket& ket = tbc_ket.GetKet(iket);
        int c = ket.p;
        int d = ket.q;
        c -= c%2;
        d -= d%2;
        int tzc = -ket.op->tz2;
        int tzd = -ket.oq->tz2;
        double tbme = 0;
        for (int Tbra = std::abs(Tzbra); Tbra<=1; Tbra++)
        {
          if (a==b and ((Tbra+Jbra%2)<1)) continue;
          double iso_clebsch_bra = CG(0.5,tza*0.5,0.5,tzb*0.5,Tbra,Tzbra);
          if (a==b and Tzbra==0) iso_clebsch_bra *= sqrt(2);
          for (int Tket = std::abs(Tzket); Tket<=1; Tket++)
          {
            double iso_clebsch_ket = CG(0.5,tzc*0.5,0.5,tzd*0.5,Tket,Tzket);
            if (c==d and Tzket==0) iso_clebsch_ket *= sqrt(2);
            double WignerEckart_factor = CG(Tket,Tzket,T_op,Tzbra-Tzket, Tbra,Tzbra) /sqrt(2*Tbra+1.) * (Tzbra-Tzket);
            uint64_t hash_key = Petr2BC_hash(a,b,c,d,Jbra,Jket,Tbra,Tket);
            tbme += DoubleReducedME[hash_key] * iso_clebsch_bra * iso_clebsch_ket * WignerEckart_factor;

          }
        }
        TBME(ibra,iket) = tbme;
      }
    }

  }

}


void ReadWrite::Write_me2j( std::string outfilename, Operator& Hbare, int emax, int Emax, int lmax)
{
  std::ofstream outfile(outfilename);
  if ( !outfile.good() )
  {
     cerr << "************************************" << std::endl
          << "**    Trouble opening file  !!!   **" << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  std::vector<int> orbits_remap;

  if (emax < 0)  emax = modelspace->GetEmax();
  if (lmax < 0)  lmax = emax;

  for (int e=0; e<=min(emax,modelspace->GetEmax()); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size()-1;

//  double tbme_pp,tbme_nn,tbme_10,tbme_00;
  float tbme_pp,tbme_nn,tbme_10,tbme_00;
  // skip the first line
  time_t time_now = chrono::system_clock::to_time_t(chrono::system_clock::now());
//  outfile << "    generated by IMSRG code on " << ctime(&time_now)<< std::endl;
  outfile << "    generated by IMSRG code on " << ctime(&time_now);
  int icount = 0;

  outfile << setiosflags(std::ios::fixed);
  std::cout << "Writing file " << outfilename << "  emax =  " << emax << "  e2max = " << Emax << "  lmax = " << lmax << "  nljmax = " << nljmax << std::endl;

  for(int nlj1=0; nlj1<=nljmax; ++nlj1)
  {
    int a =  orbits_remap[nlj1];
    Orbit & o1 = modelspace->GetOrbit(a);
    int e1 = 2*o1.n + o1.l;
    if (e1 > emax) break;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_remap[nlj2];
      Orbit & o2 = modelspace->GetOrbit(b);
      int e2 = 2*o2.n + o2.l;
      if (e1+e2 > Emax) break;
      int parity = (o1.l + o2.l) % 2;

      for(int nlj3=0; nlj3<=nlj1; ++nlj3)
      {
        int c =  orbits_remap[nlj3];
        Orbit & o3 = modelspace->GetOrbit(c);
        int e3 = 2*o3.n + o3.l;

        for(int nlj4=0; nlj4<=(nlj3==nlj1 ? nlj2 : nlj3); ++nlj4)
        {
          int d =  orbits_remap[nlj4];
          Orbit & o4 = modelspace->GetOrbit(d);
          int e4 = 2*o4.n + o4.l;
          if (e3+e4 > Emax) break;
          if ( (o1.l + o2.l + o3.l + o4.l)%2 != 0) continue;
          int Jmin = std::max( std::abs(o1.j2 - o2.j2), std::abs(o3.j2 - o4.j2) )/2;
          int Jmax = std::min (o1.j2 + o2.j2, o3.j2+o4.j2)/2;
          if (Jmin > Jmax) continue;
          for (int J=Jmin; J<=Jmax; ++J)
          {
             // me2j format is unnormalized
             double norm_factor = 1;
             if (a==b)  norm_factor *= SQRT2;
             if (c==d)  norm_factor *= SQRT2;

             // Matrix elements are written in the file with (T,Tz) = (0,0) (1,1) (1,0) (1,-1)
             tbme_pp = Hbare.TwoBody.GetTBME(J,parity,-1,a,b,c,d);        // unnormalized
             tbme_nn = Hbare.TwoBody.GetTBME(J,parity,1,a+1,b+1,c+1,d+1); // unnormalized
             tbme_10 = Hbare.TwoBody.Get_iso_TBME_from_pn(J,1,0,a,b,c,d); // normalized
             tbme_00 = Hbare.TwoBody.Get_iso_TBME_from_pn(J,0,0,a,b,c,d); // normalized

             
             outfile << setprecision(7) << setw(12) << tbme_00*norm_factor<< " "  ;
             if ((icount++)%10==9)
             {
               outfile << std::endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_nn<< " "  ;
             if ((icount++)%10==9)
             {
               outfile << std::endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_10*norm_factor << " "  ;
             if ((icount++)%10==9)
             {
               outfile << std::endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_pp << " "  ;
             if ((icount++)%10==9)
             {
               outfile << std::endl;
             }


          }
        }
      }
    }
  }
  if (icount%10 !=9) outfile << std::endl;

}






void ReadWrite::Write_me3j( std::string ofilename, Operator& Hbare, int E1max, int E2max, int E3max)
{
  std::ofstream outfile(ofilename);

  if (Hbare.particle_rank < 3)
  {
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
    cerr << " Oops. Looks like we're trying to write 3body matrix elements from a " << Hbare.particle_rank << "-body operator. For shame..." << std::endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << std::endl;
    goodstate = false;
    return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int e1max = modelspace->GetEmax();
  int e2max = modelspace->GetE2max(); // not used yet
  int e3max = modelspace->GetE3max();
  std::cout << "Writing 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << std::endl;

  std::vector<int> orbits_remap(0);
  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.

  for (int e=0; e<=min(E1max,e1max); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size();


  outfile << "(*** nsuite/me3b/v0.0.0 (Dec 22 2010) me3j-f2 ***)" << std::endl;
  outfile << setiosflags(std::ios::fixed);
  // begin giant nested loops
  int icount = 0;
  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
    int a =  orbits_remap[nlj1];
    Orbit & oa = modelspace->GetOrbit(a);
    int ea = 2*oa.n + oa.l;
    if (ea > E1max) break;
    if (ea > e1max) break;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_remap[nlj2];
      Orbit & ob = modelspace->GetOrbit(b);
      int eb = 2*ob.n + ob.l;
      if ( (ea+eb) > E2max) break;

      for(int nlj3=0; nlj3<=nlj2; ++nlj3)
      {
        int c =  orbits_remap[nlj3];
        Orbit & oc = modelspace->GetOrbit(c);
        int ec = 2*oc.n + oc.l;
        if ( (ea+eb+ec) > E3max) break;

        // Get J limits for bra <abc|
        int JabMax  = (oa.j2 + ob.j2)/2;
        int JabMin  = std::abs(oa.j2 - ob.j2)/2;

        int twoJCMindownbra;
        if (std::abs(oa.j2 - ob.j2) >oc.j2)
           twoJCMindownbra = std::abs(oa.j2 - ob.j2)-oc.j2;
        else if (oc.j2 < (oa.j2+ob.j2) )
           twoJCMindownbra = 1;
        else
           twoJCMindownbra = oc.j2 - oa.j2 - ob.j2;
        int twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2;


        // now loop over possible ket orbits
        for(int nnlj1=0; nnlj1<=nlj1; ++nnlj1)
        {
          int d =  orbits_remap[nnlj1];
          Orbit & od = modelspace->GetOrbit(d);
          int ed = 2*od.n + od.l;

          for(int nnlj2=0; nnlj2 <= ((nlj1 == nnlj1) ? nlj2 : nnlj1); ++nnlj2)
          {
            int e =  orbits_remap[nnlj2];
            Orbit & oe = modelspace->GetOrbit(e);
            int ee = 2*oe.n + oe.l;

            int nnlj3Max = (nlj1 == nnlj1 and nlj2 == nnlj2) ? nlj3 : nnlj2;
            for(int nnlj3=0; nnlj3 <= nnlj3Max; ++nnlj3)
            {
              int f =  orbits_remap[nnlj3];
              Orbit & of = modelspace->GetOrbit(f);
              int ef = 2*of.n + of.l;
              if ( (ed+ee+ef) > E3max) break;
              // check parity
              if ( (oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2 !=0 ) continue;

              // Get J limits for ket |def>
              int JJabMax = (od.j2 + oe.j2)/2;
              int JJabMin = std::abs(od.j2 - oe.j2)/2;

              int twoJCMindownket;
              if ( std::abs(od.j2 - oe.j2) > of.j2 )
                 twoJCMindownket = std::abs(od.j2 - oe.j2) - of.j2;
              else if ( of.j2 < (od.j2+oe.j2) )
                 twoJCMindownket = 1;
              else
                 twoJCMindownket = of.j2 - od.j2 - oe.j2;

              int twoJCMaxupket = od.j2 + oe.j2 + of.j2;

              int twoJCMindown = std::max(twoJCMindownbra, twoJCMindownket);
              int twoJCMaxup = std::min(twoJCMaxupbra, twoJCMaxupket);
              if (twoJCMindown > twoJCMaxup) continue;

              //inner loops
              for(int Jab = JabMin; Jab <= JabMax; Jab++)
              {
               for(int JJab = JJabMin; JJab <= JJabMax; JJab++)
               {
                //summation bounds for twoJC
                int twoJCMin = max( std::abs(2*Jab - oc.j2), std::abs(2*JJab - of.j2));
                int twoJCMax = min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
       
                for(int twoJC = twoJCMin; twoJC <= twoJCMax; twoJC += 2)
                {
//                 ++useless_counter;
//                 if (a<5 and b<5 and c<5 and d<5 and e<5 and f<5)
//                 {
//                 std::cout << "#" << useless_counter << "  " << a << "-" << b << "-" << c << "-" << "  " << d << "-" << e << "-" << f << "  "
//                      << Jab << " " << JJab << " " << twoJC << std::endl;
//                 }
                 for(int tab = 0; tab <= 1; tab++) // the total isospin loop can be replaced by i+=5
                 {
                  for(int ttab = 0; ttab <= 1; ttab++)
                  {
                   //summation bounds
                   int twoTMin = 1; // twoTMin can just be used as 1
                   int twoTMax = min( 2*tab +1, 2*ttab +1);
       
                   for(int twoT = twoTMin; twoT <= twoTMax; twoT += 2)
                   {
                    double V = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
                    if ((a==b and (Jab+tab)%2!=1) or (d==e and (JJab+ttab)%2!=1) )
                    {
                      if ( std::abs(V) > 1e-4 )  // There may be some numerical noise from using floats and 6Js at the level of 1e-6. Ignore that.
                      {
                         std::cout << "!!! Warning: <"
                              << a << " " << b << " " << c << " || " << d << " " << e << " " << f << "> ("
                              << Jab << " " << JJab << " " << twoJC << ") , (" << tab << " " << ttab << " " << twoT
                              <<") should be zero, but it's " << V << std::endl;
                      }
                      V = 0.0; // It should be zero, so set it to zero.
                    }
                    outfile << setprecision(7) << setw(12) << V << " "  ;
//                    if ((icount++)%10==9)
//                      outfile << std::endl;
                   }//twoT
                  }//ttab
                 }//tab
                 if (icount%10 == 5)
                      outfile << std::endl;
                 icount += 5;
//                if (not infile.good() ) break;
                }//twoJ
               }//JJab
       
              }//Jab


            }
          }
        }
      }
    }
  }
  if (icount%10 !=9) outfile << std::endl;

}







void ReadWrite::WriteNuShellX_op(Operator& op, std::string filename)
{
  WriteNuShellX_intfile(op,filename,"op");
}
void ReadWrite::WriteNuShellX_int(Operator& op, std::string filename)
{
  WriteNuShellX_intfile(op,filename,"int");
}

/// Write the valence space part of the interaction to a NuShellX *.int file.
/// Note that for operators other than the Hamiltonian
/// NuShellX assumes identical orbits for protons and neutrons,
/// so that the pnpn interaction should be equal to the pnnp interaction.
/// This is only approximately true for interactions generated with IMSRG,
/// so some averaging is required.
void ReadWrite::WriteNuShellX_intfile(Operator& op, std::string filename, std::string mode)
{
   std::ofstream intfile;
   intfile.open(filename, std::ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();

   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   int pdouble = 6; // precision for printing doubles

   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   std::vector<int> valence_protons(modelspace->valence.size());
   std::vector<int> valence_neutrons(modelspace->valence.size());
   auto it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->proton_orbits.begin(), modelspace->proton_orbits.end(),valence_protons.begin());
   valence_protons.resize(it-valence_protons.begin());
   it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->neutron_orbits.begin(), modelspace->neutron_orbits.end(),valence_neutrons.begin());
   valence_neutrons.resize(it-valence_neutrons.begin());

   // construct conversion maps from local orbit index to nushell index
   map<int,int> orb2nushell;
   map<int,int> nushell2orb;
   int counter = 1;
   // protons first
   for ( auto i : valence_protons ) orb2nushell[i] = counter++;
   // then neutrons
   for ( auto i : valence_neutrons ) orb2nushell[i] = counter++;
   for ( auto& it : orb2nushell) nushell2orb[it.second] = it.first;

   // Get A of the core
   int Acore=0;
   for (auto& i : modelspace->core)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 + 1;
   }
#ifdef BUILDVERSION
   intfile << "! shell model effective interaction generated by IMSRG version " << BUILDVERSION << std::endl;
#else
   intfile << "! shell model effective interaction generated by IMSRG" << std::endl;
#endif

   intfile << "! input 2N: " << File2N.substr( File2N.find_last_of("/\\")+1 ) << std::endl;
   intfile << "! input 3N: " << File3N.substr( File3N.find_last_of("/\\")+1 ) << std::endl;
   intfile << "! e1max: " << modelspace->GetEmax() << "  e2max: " << modelspace->GetE2max() << "   e3max: " << modelspace->GetE3max() << "   hw: " << modelspace->GetHbarOmega();
   intfile << "   Aref: " << Aref << "  Zref: " << Zref << "  A_for_kinetic_energy: " << modelspace->GetTargetMass() << std::endl;
   intfile << "! Zero body term: " << op.ZeroBody << std::endl;
   intfile << "! Index   n l j tz" << std::endl;

   for ( auto& it : nushell2orb )
   {
      Orbit& oi = modelspace->GetOrbit(it.second);
      intfile << "!  " << it.first << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << std::endl;
   }
   intfile << "!" << std::endl;
   intfile << "-999  ";
   for ( auto& it : nushell2orb )
   {
      intfile << op.OneBody(it.second,it.second) << "  ";
   }

   intfile << "  " << Acore << "  " << Acore+2 << "  0.00000 " << std::endl; // No mass dependence for now...

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      for (auto& ibra: tbc.GetKetIndex_vv() )
      {
         Ket &bra = tbc.GetKet(ibra);
         int a = bra.p;
         int b = bra.q;
         Orbit& oa = modelspace->GetOrbit(a);
         Orbit& ob = modelspace->GetOrbit(b);
         for (auto& iket: tbc.GetKetIndex_vv())
         {
            if (iket<ibra) continue;
            Ket &ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            Orbit& oc = modelspace->GetOrbit(c);
            Orbit& od = modelspace->GetOrbit(d);

            // don't pull a_ind and b_ind out of this loop, on account of the swap below.
            int a_ind = orb2nushell[a];
            int b_ind = orb2nushell[b];
            int c_ind = orb2nushell[c];
            int d_ind = orb2nushell[d];
            int T = std::abs(tbc.Tz);

            double tbme = op.TwoBody.GetTBME_norm(ch,a,b,c,d);
            // NuShellX quirk: even though it uses pn formalism, it requires that Vpnpn = Vnpnp,
            //  i.e. the spatial wf for protons and neutrons are assumed to be the same.
            // This seems to only be an issue for expectation values of operators,
            // so the mode "op" averages Vpnpn and Vnpnp to make them equal.
            if (mode=="op" and oa.tz2!=ob.tz2)
            {
               int aa = a - oa.tz2;
               int bb = b - ob.tz2;
               int cc = c - oc.tz2;
               int dd = d - od.tz2;
               tbme += op.TwoBody.GetTBME_norm(ch,aa,bb,cc,dd);
               tbme /= 2;
            }

            if ( std::abs(tbme) < 1e-6) continue;
            if (T==0)
            {
               if ( oa.j2 != ob.j2 or oa.l != ob.l or oa.n != ob.n ) tbme *= SQRT2; // pn TBMEs are unnormalized
               if ( oc.j2 != od.j2 or oc.l != od.l or oc.n != od.n ) tbme *= SQRT2; // pn TBMEs are unnormalized
               T = (tbc.J+1)%2;
            }
            // in NuShellX, the proton orbits must come first. This can be achieved by
            // ensuring that the bra and ket indices are in ascending order.
            if (a_ind > b_ind)
            {
               swap(a_ind,b_ind);
               tbme *= bra.Phase(tbc.J);
            }
            if (c_ind > d_ind)
            {
               swap(c_ind,d_ind);
               tbme *= ket.Phase(tbc.J);
            }
            if ((a_ind > c_ind) or (c_ind==a_ind and b_ind>d_ind) )
            {
              swap(a_ind,c_ind);
              swap(b_ind,d_ind);
            }

            intfile  << setw(wint) << a_ind << setw(wint) << b_ind << setw(wint) << c_ind << setw(wint) << d_ind << "    "
                     << setw(wint) << tbc.J << " " << setw(wint) << T     << "       "
                     << setw(wdouble) << setiosflags(std::ios::fixed) << setprecision(pdouble) << tbme
                     << std::endl;
         }
      }
   }
   intfile.close();
}





/// Write the valence model space to a NuShellX *.sp file.
void ReadWrite::WriteNuShellX_sps(Operator& op, std::string filename)
{
   std::ofstream spfile;
   spfile.open(filename, std::ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();

   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   std::vector<int> valence_protons(modelspace->valence.size());
   std::vector<int> valence_neutrons(modelspace->valence.size());
   auto it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->proton_orbits.begin(), modelspace->proton_orbits.end(),valence_protons.begin());
   valence_protons.resize(it-valence_protons.begin());
   it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->neutron_orbits.begin(), modelspace->neutron_orbits.end(),valence_neutrons.begin());
   valence_neutrons.resize(it-valence_neutrons.begin());

   
   // construct conversion maps from local orbit index to nushell index
   map<int,int> orb2nushell;
   map<int,int> nushell2orb;
   int counter = 1;
   // protons first
   for ( auto i : valence_protons ) orb2nushell[i] = counter++;
   // then neutrons
   for ( auto i : valence_neutrons ) orb2nushell[i] = counter++;
   for ( auto& it : orb2nushell) nushell2orb[it.second] = it.first;

   // Get A,Z of the core
   int Acore=0, Zcore=0;
   for (auto& i : modelspace->core)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 + 1;
      if (oi.tz2 < 0)
        Zcore += oi.j2+1;
   }

   // Write to file. In NuShellX convention, radial QN n starts at 1
   spfile << "! modelspace for IMSRG interaction" << std::endl;
   spfile << "pn" << std::endl; // proton-neutron formalism
   spfile << Acore << " " << Zcore << std::endl;
   spfile << modelspace->valence.size() << std::endl;
   spfile << "2 " << valence_protons.size() << " " << valence_neutrons.size() << std::endl;

   for (auto& it : nushell2orb )
   {
     Orbit& oi = modelspace->GetOrbit(it.second);
     spfile << it.first << " " << oi.n+1 << " " << oi.l << " " << oi.j2 << std::endl;
   }

   spfile << std::endl;
   spfile.close();

}


void ReadWrite::ReadNuShellX_int(Operator& op, std::string filename)
{
  ModelSpace* modelspace = op.GetModelSpace();
  std::ifstream intfile(filename);
  char line[500];
  while(intfile.getline(line,500))
  {
    if (line[0] != '!') break;
    std::string linestr(line);
    if ( linestr.find("Index") != std::string::npos ) break;
  }
  map<index_t,index_t> orbit_map;
  while(intfile.getline(line,500))
  {
    if (line[0] != '!') break;
    char exclamation_point;
    index_t indx;
    int n,l,j2,tz2;
    std::string j,tz;
    std::string linestr(line);
    if ( linestr.find("/") == std::string::npos ) break;
    std::istringstream(line) >> exclamation_point >> indx >> n >> l >> j >> tz;
    std::istringstream(j.substr(0,j.find("/")-j.front())) >> j2;
    std::istringstream(tz.substr(0,j.find("/")-tz.front())) >> tz2;
    orbit_map[indx] = modelspace->GetOrbitIndex(n,l,j2,tz2);
  }

  double dummy;
  intfile >> dummy; // read the -999 that doesn't mean anything
  for (auto& orb : orbit_map )
  {
    intfile >> dummy;
    op.OneBody(orb.second,orb.second) = dummy;
  }

  double acore_nushell, a0_nushell, exp_nushell;
  intfile >> acore_nushell >> a0_nushell >> exp_nushell;
  double scalefactor = pow( a0_nushell / modelspace->GetAref(), exp_nushell) ;
//  for (int i=0;i<3;++i) intfile >> dummy; // A-dependence parameters
//  std::cout << op.OneBody << std::endl;
  int a,b,c,d,J,Tprime;
  double V;
  while( intfile >> a >> b >> c >> d >> J >> Tprime >> V)
  {

    V *= scalefactor;
    Orbit& oa = modelspace->GetOrbit(orbit_map[a]);
    Orbit& ob = modelspace->GetOrbit(orbit_map[b]);
    Orbit& oc = modelspace->GetOrbit(orbit_map[c]);
    Orbit& od = modelspace->GetOrbit(orbit_map[d]);
    if (oa.tz2 != ob.tz2)
    {
       if ( (oa.j2 != ob.j2) or (oa.l != ob.l) or (oa.n != ob.n) ) V /= SQRT2; // pn TBMEs are unnormalized
       if ( (oc.j2 != od.j2) or (oc.l != od.l) or (oc.n != od.n) ) V /= SQRT2; // pn TBMEs are unnormalized
    }
    op.TwoBody.SetTBME_J(J,orbit_map[a],orbit_map[b],orbit_map[c],orbit_map[d],V);
  }

}


void ReadWrite::ReadNuShellX_int_iso(Operator& op, std::string filename)
{
  ModelSpace* modelspace = op.GetModelSpace();
  std::ifstream intfile(filename);
  char line[500];
  while(intfile.getline(line,500))
  {
    if (line[0] != '!') break;
    std::string linestr(line);
    if ( linestr.find("Index") != std::string::npos ) break;
  }
  map<index_t,index_t> orbit_map;
  while(intfile.getline(line,500))
  {
    if (line[0] != '!') break;
    char exclamation_point;
    index_t indx;
    int n,l,j2,tz2;
    std::string j,tz;
    std::string linestr(line);
    if ( linestr.find("/") == std::string::npos ) break;
    std::istringstream(line) >> exclamation_point >> indx >> n >> l >> j >> tz;
    std::istringstream(j.substr(0,j.find("/")-j.front())) >> j2;
    std::istringstream(tz.substr(0,j.find("/")-tz.front())) >> tz2;
    orbit_map[indx] = modelspace->GetOrbitIndex(n,l,j2,-1);
  }

  double dummy;
  intfile >> dummy; // read the -999 that doesn't mean anything
//  std::cout << " Reading SPE line" << std::endl;
  for (auto& orb : orbit_map )
  {
    intfile >> dummy;
    op.OneBody(orb.second,orb.second) = dummy;
    op.OneBody(orb.second+1,orb.second+1) = dummy;
//    std::cout << orb.second << " " << dummy << std::endl;
  }

  double acore_nushell, a0_nushell, exp_nushell;
  intfile >> acore_nushell >> a0_nushell >> exp_nushell;
  double scalefactor = pow( a0_nushell / modelspace->GetAref(), exp_nushell) ;
//  std::cout << "acore_nushell = " << acore_nushell << "   a0_nushell = " << a0_nushell << "  exp_nushell = " << exp_nushell << std::endl;
//  for (int i=0;i<3;++i) intfile >> dummy; // A-dependence parameters
//  std::cout << op.OneBody << std::endl;
  uint64_t a,b,c,d,J,T;
  double V;
  map<uint64_t,double> IsoTBME;
  while( intfile >> a >> b >> c >> d >> J >> T >> V)
  {
    V *= scalefactor;
    uint64_t hashkey = ( ((uint64_t) J) + ((uint64_t)T<<10) + (a<<20) + (b<<30) + (c<<40) + (d<<50)  );
//    std::cout << "setting " << a << " " << b << " " << c << " " << d << " " << J << " " << T << " ->  "<< hashkey << "   " << V << std::endl;
    IsoTBME[hashkey] = V;
  }

  for ( auto& a_indx : orbit_map )
  {
   Orbit& oa = modelspace->GetOrbit(a_indx.second);
   for ( auto& b_indx : orbit_map )
   {
    Orbit& ob = modelspace->GetOrbit(b_indx.second);
    for ( auto& c_indx : orbit_map )
    {
     Orbit& oc = modelspace->GetOrbit(c_indx.second);
     for ( auto& d_indx : orbit_map )
     {
      Orbit& od = modelspace->GetOrbit(d_indx.second);
      int Jmin = std::max( std::abs(oa.j2-ob.j2), std::abs(oc.j2-od.j2) )/2;
      int Jmax = std::min( oa.j2+ob.j2, oc.j2+od.j2 )/2;
      for (int J=Jmin; J<=Jmax; ++J)
      {
        uint64_t hashT1 = ( J + (1L<<10) + (a_indx.first<<20) + (b_indx.first<<30) + (c_indx.first<<40) + (d_indx.first<<50) );
        uint64_t hashT0 = ( J + (0L<<10) + (a_indx.first<<20) + (b_indx.first<<30) + (c_indx.first<<40) + (d_indx.first<<50) );
        double V1 = IsoTBME[hashT1];
        double V0 = IsoTBME[hashT0];
//        std::cout << J << " " <<  a_indx.first   << " " <<    b_indx.first << " " <<    c_indx.first << " " <<    d_indx.first << "     " << hashT1 << ": " <<    V1 << "   " << hashT0 << "  " << V0 << std::endl;
        double Vpp = V1;
        double Vpnpn = (V1 + V0) * 0.5;
        double Vpnnp = (V1 - V0) * 0.5;
        double Vnppn = (V1 - V0) * 0.5;
        double Vnpnp = (V1 + V0) * 0.5;

        if ( std::abs(Vpp)>1e-6 )
        {
//        std::cout << "Vnn: " << J << " " <<  a_indx.second   << " " <<    b_indx.second << " " <<    c_indx.second << " " <<    d_indx.second << " " <<    Vpp << std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second,   b_indx.second,   c_indx.second,   d_indx.second,   Vpp);
//        std::cout << "Vpp: " << J << " " <<  a_indx.second+1 << " " <<    b_indx.second+1 << " " <<  c_indx.second+1 << " " <<  d_indx.second+1 << " " <<  Vpp<< std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second+1, b_indx.second+1, c_indx.second+1, d_indx.second+1, Vpp);
        }
        if (std::abs(Vpnpn)>1e-6)
        {
//        std::cout << "Vpnpn: " << J << " " <<  a_indx.second+1 << " " <<    b_indx.second << " " <<    c_indx.second+1 << " " <<  d_indx.second << " " <<    Vpnpn<< std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second+1, b_indx.second,   c_indx.second+1, d_indx.second,   Vpnpn);
        }
        if (std::abs(Vpnnp)>1e-6)
        {
//        std::cout << "Vpnnp: " << J << " " <<  a_indx.second+1 << " " <<    b_indx.second << " " <<    c_indx.second << " " <<    d_indx.second+1 << " " <<  Vpnnp<< std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second+1, b_indx.second,   c_indx.second,   d_indx.second+1, Vpnnp);
        }
        if (std::abs(Vnppn)>1e-6)
        {
//        std::cout << "Vnppn: " << J << " " <<  a_indx.second   << " " <<    b_indx.second+1 << " " <<  c_indx.second+1 << " " <<  d_indx.second << " " <<    Vnppn<< std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second,   b_indx.second+1, c_indx.second+1, d_indx.second,   Vnppn);
        }
        if (std::abs(Vnpnp)>1e-6)
        {
//        std::cout << "Vnpnp: " << J << " " <<  a_indx.second   << " " <<    b_indx.second+1 << " " <<  c_indx.second << " " <<    d_indx.second+1 << " " <<  Vnpnp<< std::endl;
        op.TwoBody.SetTBME_J(J, a_indx.second,   b_indx.second+1, c_indx.second,   d_indx.second+1, Vnpnp);
        }
      }
     }
    }
   }
  }

}




/// This now appears to be working properly
void ReadWrite::WriteAntoine_int(Operator& op, std::string filename)
{
   std::ofstream intfile;
   intfile.open(filename, std::ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int nvalence_orbits = modelspace->valence.size();
   int Acore = modelspace->GetAref(); // Final interaction should have ref = core
   int Zcore = modelspace->GetAref(); // Final interaction should have ref = core
   int idens = 0; // no density dependence
   std::string title = "IMSRG INTERACTION";

   std::vector<int> nlj_labels; // list of antoine-style nlj labels
   map<index_t,index_t> orbit_map; // map from index number to nlj index
   std::vector<double> spe(nvalence_orbits,0.0); // one-body terms
   for (auto& v : modelspace->valence)
   {
     Orbit& ov = modelspace->GetOrbit(v);
//     int label = ov.n*1000 + ov.l*100 + ov.j2;
     int label = ((ov.n*2+(ov.tz2+1)/2 )%10)*1000 + ov.l*100 + ov.j2;
     if ( find(nlj_labels.begin(),nlj_labels.end(),label) == nlj_labels.end() )
        nlj_labels.push_back(label);
     orbit_map[v] = find(nlj_labels.begin(),nlj_labels.end(),label) - nlj_labels.begin();
     spe[orbit_map[v]] = op.OneBody(v,v);
   }

   intfile << title << std::endl;
   // 2 indicates pn treated separately
   intfile << "2 " << nvalence_orbits << " ";
   // write NLJ of the orbits
   for (auto& nlj : nlj_labels) intfile << nlj << " ";
   intfile << std::endl;
   // write SPEs
   for (auto& s : spe) intfile << " " << s;
   intfile << std::endl;
   for (auto& s : spe) intfile << " " << s;
   intfile << std::endl;
   // No density dependence
   intfile << idens << " " << Zcore << " " << Acore-Zcore << " " << 0.0 << " " << 0.0 << std::endl;

   // Start loop over TBMEs
   for (auto a : modelspace->valence)
   {
     Orbit& oa = modelspace->GetOrbit(a);
     int nlja = nlj_labels[orbit_map[a]];
     for (auto b : modelspace->valence)
     {
//       if (b<a) continue;
       Orbit& ob = modelspace->GetOrbit(b);
       if (oa.tz2 > ob.tz2) continue;
       if (oa.tz2==ob.tz2 and b<a) continue;
       int nljb = nlj_labels[orbit_map[b]];
       for (auto c : modelspace->valence)
       {
         if (c<a) continue;
         Orbit& oc = modelspace->GetOrbit(c);
         int nljc = nlj_labels[orbit_map[c]];
         for (auto d : modelspace->valence)
         {
           Orbit& od = modelspace->GetOrbit(d);
           if (oc.tz2 > od.tz2) continue;
           if (oc.tz2==od.tz2 and (d<c or (c==a and d<b))) continue;
           if ((oa.tz2+ob.tz2) != (oc.tz2+od.tz2) ) continue;
           if ( (oa.l+ob.l+oc.l+od.l)%2>0 ) continue;
//           std::cout << a << " " << b << " " << c << " " << d << std::endl;
           int nljd = nlj_labels[orbit_map[d]];
           int Tmin = std::abs(oa.tz2+ob.tz2) -1; // -1 means pn, 1 means pp or nn. T loop goes std::abs(Tmin) to Tmax
           int Tmax = 1; // always 1.
           int Jmin = std::max(std::abs(oa.j2-ob.j2),std::abs(oc.j2-od.j2))/2;
           int Jmax = std::min(oa.j2+ob.j2,oc.j2+od.j2)/2;
           if (Jmin<=Jmax)
           {
             intfile << setw(2) << Tmin << " " << setw(2) << Tmax << " " << nlja << " " << nljb << " " << nljc << " " << nljd << " " << setw(2) << Jmin << " " << setw(2) << Jmax << std::endl;
             for (int J=Jmin;J<=Jmax;++J)
             {
                intfile.setf(std::ios::fixed, std::ios::floatfield );
                intfile << " " << setw(10) << setprecision(6) <<  op.TwoBody.GetTBME_J_norm(J, a, b, c, d);
             }
             intfile << std::endl;
           }
         }
       }
     }
   }
   
}


/// Generate an input file for antoine, which may be edited afterwards if need be
void ReadWrite::WriteAntoine_input(Operator& op, std::string filename,int A, int Z)
{
   std::ofstream inputfile;
   inputfile.open(filename, std::ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
//   int nvalence_orbits = modelspace->valence.size();
   int Acore = modelspace->GetAref(); // Final interaction should have ref = core
   int Zcore = modelspace->GetZref(); // Final interaction should have ref = core
   int Aval = A-Acore;
   int Zval = Z-Zcore; 
   int MPRO = Aval%2; // M projection to use
   std::vector<int> JVAL = {MPRO,MPRO+2,MPRO+4,MPRO+6}; // 2*J values to calculate
   int N_for_each_J = 2; // number to calculate for each J. This could be different for each, if desired.

   std::vector<int> proton_shells;
   std::vector<int> neutron_shells;

   for (auto& v : modelspace->valence)
   {
      Orbit& ov = modelspace->GetOrbit(v);
      int nlj =  ((ov.n*2+(ov.tz2+1)/2 )%10)*1000 + ov.l*100 + ov.j2;
      if (ov.tz2 < 0)
        proton_shells.push_back(nlj);
      else
        neutron_shells.push_back(nlj);
   }


   int CAS = 4; // usual Lanczos calculation with random pivot
//   int KTEXT = 1; // will there be a name given?
   int KTEXT = 0; // will there be a name given?
   std::string TEXT = "IMSRG"; // name that goes somewhere?
   int IPRI = 0; // printing option. IPRI=1 prints more stuff

//   inputfile << CAS << " " << KTEXT << " " << IPRI << " " << TEXT << std::endl;
   inputfile << CAS << " " << KTEXT << " " << IPRI << std::endl;
   
   int Fil1 = 50; // file where we write std::vectors
   int NLEC1 = 0; // number of std::vectors read at the beginning
   int NCAL = JVAL.size(); // number of std::vectors to treat. 0 means "all"
   
   inputfile  << Fil1 << " " << NLEC1 << " " << NCAL << std::endl;

   // do protons first
   inputfile << Zval << " " << proton_shells.size();
   for (auto& p : proton_shells ) inputfile << " " << p;
   for (size_t i=0;i<proton_shells.size();++i ) inputfile << " " << 0; // weight factors for orbits
   inputfile << " " << 0 << std::endl; // SAUT = maximum weight factor for protons

   // do neutrons second
   inputfile << Aval-Zval << " " << proton_shells.size();
   for (auto& n : neutron_shells ) inputfile << " " << n;
   for (size_t i=0;i<neutron_shells.size();++i ) inputfile << " " << 0; // weight factors for orbits
   inputfile << " " << 0 << std::endl; // SAUT = maximum weight factor for neutrons

   int PARI = ((Zval%2) * (proton_shells[0]%1000)/100 + ((Aval-Zval)%2) * (neutron_shells[0]%1000)/100 )%2; // parity. 0=+ , 1=-. make a reasonable guess.
   int JUMP = 0; // maximum weight factor for protons + neutrons

   inputfile << " " << MPRO << " " << PARI << " " << JUMP << std::endl;

   int FLNUC = 90; // file number for hamiltonian file
   int COUL = 0; // 0 means no Coulomb. we'll do that ourselves, thank you.
   inputfile << FLNUC << " " << COUL << std::endl;

   for (auto jv : JVAL) inputfile << jv << " ";
   inputfile << std::endl;
   for (size_t i=0;i<JVAL.size();++i) inputfile << N_for_each_J << " ";
   inputfile << std::endl;

   int NLOOP = 200; // maximum number of lanczos iterations
   double ZFIT = 0.0005;  // test of convergence
   int ORTH = 0; // orthogonalize each state to all states in the std::vector. usual is 0=no
   inputfile << NLOOP << " "<< ZFIT << " " << ORTH << std::endl;

}



/// Write an operator to a plain-text file
void ReadWrite::WriteOperatorHuman(Operator& op, std::string filename)
{
   std::ofstream opfile;
   opfile.open(filename, std::ofstream::out);
   if (not opfile.good() )
   {
     std::cout << "Trouble opening " << filename << ". Aborting WriteOperator." << std::endl;
     return;
   }
   ModelSpace * modelspace = op.GetModelSpace();

   if (op.IsHermitian() )
   {
      opfile << "Hermitian" << std::endl;
   }
   else if (op.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << std::endl;
   }
   else
   {
      opfile << "Non-Hermitian" << std::endl;
   }
   opfile << op.GetJRank() << "  " << op.GetTRank() << "  " << op.GetParity() << std::endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op.ZeroBody << std::endl;

   opfile << "$OneBody:\t" << std::endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (std::abs(op.OneBody(i,j)) > 0)
            opfile << fixed << setw(3) << i << "\t" << fixed << setw(3) << j << "\t" << fixed << setw(18) << setprecision(12) << op.OneBody(i,j) << std::endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << std::endl;

   for ( auto& it : op.TwoBody.MatEl )
   {
      int chbra = it.first[0];
      int chket = it.first[1];
      int nbras = it.second.n_rows;
      int nkets = it.second.n_cols;
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(chbra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(chket);
      for (int ibra=0; ibra<nbras; ++ibra)
      {
        Ket& bra = tbc_bra.GetKet(ibra);
        for (int iket=0; iket<nkets; ++iket)
        {
          Ket& ket = tbc_ket.GetKet(iket);
           double tbme = it.second(ibra,iket);
//           if (bra.p == bra.q) tbme *= sqrt(2); // For comparison with Nathan CHANGE THIS
//           if (ket.p == ket.q) tbme *= sqrt(2); // For comparison with Nathan CHANGE THIS
           if ( std::abs(tbme) > 1e-7 )
           {
             opfile << setw(2) << tbc_bra.J << " " << setw(2) << tbc_bra.parity << " " << setw(3) << tbc_bra.Tz  << "    "
                    << setw(2) << tbc_ket.J << " " << setw(2) << tbc_ket.parity << " "  << setw(3) << tbc_ket.Tz  << "    "
                  << setw(3) << bra.p << " "  << setw(3) << bra.q  << " "  << setw(3) << ket.p << " "  << setw(3) << ket.q  << "   "
                  << fixed << setw(18) << setprecision(12) << tbme << std::endl;
           }
        }
      }
   }

   opfile.close();

}





/// Write an operator to a plain-text file
void ReadWrite::WriteOperator(Operator& op, std::string filename)
{
   std::ofstream opfile;
   opfile.open(filename, std::ofstream::out);
   if (not opfile.good() )
   {
     std::cout << "Trouble opening " << filename << ". Aborting WriteOperator." << std::endl;
     return;
   }
   ModelSpace * modelspace = op.GetModelSpace();

   if (op.IsHermitian() )
   {
      opfile << "Hermitian" << std::endl;
   }
   else if (op.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << std::endl;
   }
   else
   {
      opfile << "Non-Hermitian" << std::endl;
   }
   opfile << op.GetJRank() << "  " << op.GetTRank() << "  " << op.GetParity() << std::endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op.ZeroBody << std::endl;

   opfile << "$OneBody:\t" << std::endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (std::abs(op.OneBody(i,j)) > 0)
            opfile << i << "\t" << j << "\t" << setprecision(10) << op.OneBody(i,j) << std::endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << std::endl;

   for ( auto& it : op.TwoBody.MatEl )
   {
      int chbra = it.first[0];
      int chket = it.first[1];
      int nbras = it.second.n_rows;
      int nkets = it.second.n_cols;
      for (int ibra=0; ibra<nbras; ++ibra)
      {
        for (int iket=0; iket<nkets; ++iket)
        {
           double tbme = it.second(ibra,iket);
           if ( std::abs(tbme) > 1e-7 )
           {
             opfile << setw(4) << chbra << " " << setw(4) << chket << "   "
                  << setw(4) << ibra  << " " << setw(4) << iket  << "   "
                  << setw(10) << setprecision(6) << tbme << std::endl;
           }
        }
      }
   }

   opfile.close();

}

/// Read an operator from a plain-text file
void ReadWrite::ReadOperator(Operator &op, std::string filename)
{
   std::ifstream opfile;
   opfile.open(filename);
   if (not opfile.good() )
   {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
   }

   std::string tmpstr;
   int i,j,chbra,chket;
   double v;

   opfile >> tmpstr;
   if (tmpstr == "Hermitian")
   {
      op.SetHermitian();
   }
   else if (tmpstr == "Anti-Hermitian")
   {
      op.SetAntiHermitian();
   }
   else
   {
      op.SetNonHermitian();
   }
   int jrank,trank,parity;
   opfile >> jrank >> trank >> parity;

   opfile >> tmpstr >> v;
   op.ZeroBody = v;

   getline(opfile, tmpstr);
   getline(opfile, tmpstr);
   getline(opfile, tmpstr);
   while (tmpstr[0] != '$')
   {
      std::stringstream ss(tmpstr);
      ss >> i >> j >> v;
      op.OneBody(i,j) = v;
      if ( op.IsHermitian() )
         op.OneBody(j,i) = v;
      else if ( op.IsAntiHermitian() )
         op.OneBody(j,i) = -v;
      getline(opfile, tmpstr);
   }

  while(opfile >> chbra >> chket >> i >> j >> v)
  {
    op.TwoBody.SetTBME(chbra,chket,i,j,v);
  }

   opfile.close();
   
}

/// Read an operator from a plain-text file
void ReadWrite::ReadOperatorHuman(Operator &op, std::string filename)
{
   std::ifstream opfile;
   opfile.open(filename);
   if (not opfile.good() )
   {
     cerr << "************************************" << std::endl
          << "**    Trouble reading file  !!!   **" << filename << std::endl
          << "************************************" << std::endl;
     goodstate = false;
     return;
   }

   std::string tmpstr;
   int i,j,k,l,Jbra,Jket,pbra,pket,Tzbra,Tzket;
   double v;

   opfile >> tmpstr;
   if (tmpstr == "Hermitian")
   {
      op.SetHermitian();
   }
   else if (tmpstr == "Anti-Hermitian")
   {
      op.SetAntiHermitian();
   }
   else
   {
      op.SetNonHermitian();
   }
   int jrank,trank,parity;
   opfile >> jrank >> trank >> parity;

   opfile >> tmpstr >> v;
   op.ZeroBody = v;

   getline(opfile, tmpstr);
   getline(opfile, tmpstr);
   getline(opfile, tmpstr);
   while (tmpstr[0] != '$')
   {
      std::stringstream ss(tmpstr);
      ss >> i >> j >> v;
      op.OneBody(i,j) = v;
      if ( op.IsHermitian() )
         op.OneBody(j,i) = v;
      else if ( op.IsAntiHermitian() )
         op.OneBody(j,i) = -v;
      getline(opfile, tmpstr);
   }

//  while(opfile >> chbra >> chket >> i >> j >> v)
  while(opfile >> Jbra >> pbra >> Tzbra >> Jket >> pket >> Tzket >> i >> j >> k >> l >> v)
  {
    op.TwoBody.SetTBME(Jbra,pbra,Tzbra,Jket,pket,Tzket,i,j,k,l,v);
  }

   opfile.close();
   
}




/// Write an operator to a plain-text file
void ReadWrite::CompareOperators(Operator& op1, Operator& op2, std::string filename)
{
   std::ofstream opfile;
   opfile.open(filename, std::ofstream::out);
   if (not opfile.good() )
   {
     std::cout << "Trouble opening " << filename << ". Aborting WriteOperator." << std::endl;
     return;
   }
   ModelSpace * modelspace = op1.GetModelSpace();

   if (op1.IsHermitian() )
   {
      opfile << "Hermitian" << std::endl;
   }
   else if (op1.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << std::endl;
   }
   else
   {
      opfile << "Non-Hermitian" << std::endl;
   }
   opfile << op1.GetJRank() << "  " << op1.GetTRank() << "  " << op1.GetParity() << std::endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op1.ZeroBody << "   " << op2.ZeroBody << std::endl;

   opfile << "$OneBody:\t" << std::endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op1.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (std::abs(op1.OneBody(i,j)) > 0 or std::abs(op2.OneBody(i,j))>0 )
            opfile << i << "\t" << j << "\t" << setprecision(10) << op1.OneBody(i,j) << "   " << op2.OneBody(i,j) << std::endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << std::endl;

   for ( auto& it : op1.TwoBody.MatEl )
   {
      int chbra = it.first[0];
      int chket = it.first[1];
      int nbras = it.second.n_rows;
      int nkets = it.second.n_cols;
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(chbra);
      for (int ibra=0; ibra<nbras; ++ibra)
      {
        Ket& bra = tbc_bra.GetKet(ibra);
        for (int iket=0; iket<nkets; ++iket)
        {
          Ket& ket = tbc_bra.GetKet(iket);
           double tbme1 = it.second(ibra,iket);
           double tbme2 = op2.TwoBody.GetMatrix(chbra,chket)(ibra,iket);
           if ( std::abs(tbme1) > 1e-7 or std::abs(tbme2)>1e-7 )
           {
             opfile << setw(4) << tbc_bra.J << " " << tbc_bra.parity << " " << tbc_bra.Tz  << "    "
                  << setw(4) << bra.p << " " << bra.q  << " " << ket.p << " " << ket.q  << "   "
                  << setw(10) << setprecision(6) << tbme1 << "  "
                  << setw(10) << setprecision(6) << tbme2 << std::endl;
           }
        }
      }
   }

   opfile.close();

}



void ReadWrite::ReadOneBody_Takayuki(std::string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2 + 1;
     orbits_remap[nlj] = i;
  }


  std::ifstream infile(filename);
  int a,b,tza,tzb;
  double me;
  while( infile >> tza >> a >> tzb >> b >> me )
  {
    int aa = orbits_remap.at(a) + (tza+1)/2;
    int bb = orbits_remap.at(b) + (tzb+1)/2;
    Hbare.OneBody(aa,bb) = me;
    if (Hbare.IsHermitian())
      Hbare.OneBody(bb,aa) = me;
    else if (Hbare.IsAntiHermitian())
      Hbare.OneBody(bb,aa) = -me;
  }
}

void ReadWrite::ReadTwoBody_Takayuki(std::string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2 + 1;
     orbits_remap[nlj] = i;
  }


  std::string dummy;
  std::ifstream infile(filename);
  getline( infile, dummy );
  getline( infile, dummy );
  int a,b,c,d,tza,tzb,tzc,tzd,J;
  double me,pipj,rirj;
  while( infile >> tza >> a >> tzb >> b >> tzc >> c >> tzd >> d >> J >> me >> pipj >> rirj )
  {
    if ( orbits_remap.find(a) == orbits_remap.end() ) continue;
    if ( orbits_remap.find(b) == orbits_remap.end() ) continue;
    if ( orbits_remap.find(c) == orbits_remap.end() ) continue;
    if ( orbits_remap.find(d) == orbits_remap.end() ) continue;
    int aa = orbits_remap.at(a) + (tza+1)/2;
    int bb = orbits_remap.at(b) + (tzb+1)/2;
    int cc = orbits_remap.at(c) + (tzc+1)/2;
    int dd = orbits_remap.at(d) + (tzd+1)/2;
    if ( (aa==bb or cc==dd) and (J%2)>0 ) continue;
    if (std::abs(me)<1e-6) continue;
    Hbare.TwoBody.SetTBME_J(J,aa,bb,cc,dd,me);
  }
}


void ReadWrite::WriteOneBody_Takayuki(std::string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2 + 1;
     orbits_remap[i]   = nlj;
     orbits_remap[i+1] = nlj;
  }

  std::ofstream outfile(filename);
  outfile << setiosflags(std::ios::fixed);

  for (int a=0; a<norb; ++a)
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
    {
      Orbit& ob = modelspace->GetOrbit(b);
      double me = Hbare.OneBody(a,b);
      if (std::abs(me) > 1e-7)
      {
      outfile << setw(3) << oa.tz2 << " " << setw(3) << orbits_remap.at(a) << " "
              << setw(3) << ob.tz2 << " " << setw(3) << orbits_remap.at(b) << " " 
              << setw(12) << setprecision(8) <<  me << std::endl;
      }
    }
  }

}



void ReadWrite::WriteTwoBody_Takayuki(std::string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  std::unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + std::max(oi.l-1,0) + (oi.j2 - std::abs(2*oi.l-1))/2 + 1;
//     orbits_remap[nlj] = i;
     orbits_remap[i]   = nlj;
     orbits_remap[i+1] = nlj;
  }

  std::ofstream outfile(filename);
  outfile << setiosflags(std::ios::fixed);

  for ( auto& itmat : Hbare.TwoBody.MatEl )
  {
    int ch_ket = itmat.first[1];
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch_ket);
    int J = tbc.J;
    int nkets = tbc.GetNumberKets();
    for (int ibra=0;ibra<nkets;++ibra)
    {
      Ket& bra = tbc.GetKet(ibra);
      for (int iket=0; iket<=ibra; ++iket)
      {
        Ket& ket = tbc.GetKet(iket);
        double tbme = itmat.second(ibra,iket);
        if (std::abs(tbme)<1e-8) continue;
        outfile << setw(3) << bra.op->tz2 << " " << setw(3) << orbits_remap.at(bra.p) << " "
                << setw(3) << bra.oq->tz2 << " " << setw(3) << orbits_remap.at(bra.q) << " "
                << setw(3) << ket.op->tz2 << " " << setw(3) << orbits_remap.at(ket.p) << " "
                << setw(3) << ket.oq->tz2 << " " << setw(3) << orbits_remap.at(ket.q) << " "
                << setw(3) << J << setw(12) << setprecision(8) << tbme << std::endl;
      }
    }
  }


}


void ReadWrite::WriteTensorOneBody(std::string filename, Operator& Op, std::string opname)
{
   std::ofstream outfile(filename);
   ModelSpace * modelspace = Op.GetModelSpace();
   int wint = 4; // width for printing integers

   int wdouble = 16; // width for printing doubles
   int pdouble = 9; // precision for printing doubles


   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   std::vector<int> valence_protons(modelspace->valence.size());
   std::vector<int> valence_neutrons(modelspace->valence.size());
   auto it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->proton_orbits.begin(), modelspace->proton_orbits.end(),valence_protons.begin());
   valence_protons.resize(it-valence_protons.begin());
   it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->neutron_orbits.begin(), modelspace->neutron_orbits.end(),valence_neutrons.begin());
   valence_neutrons.resize(it-valence_neutrons.begin());
   // construct conversion maps from local orbit index to nushell index
   map<int,int> orb2nushell;
   map<int,int> nushell2orb;
   int counter = 1;
   // protons first
   for ( auto i : valence_protons ) orb2nushell[i] = counter++;
   // then neutrons
   for ( auto i : valence_neutrons ) orb2nushell[i] = counter++;
   for ( auto& it : orb2nushell) nushell2orb[it.second] = it.first;

   outfile << fixed << setprecision(pdouble);
   outfile << "!  One-body matrix elements for tensor operator: " << opname << "   generated with IM-SRG" << std::endl;
   outfile << "!  Rank_J :  " << Op.GetJRank() << std::endl;
   outfile << "!  Rank_T :  " << Op.GetTRank() << std::endl;
   outfile << "!  Parity :  " << showpos << 1-2*Op.GetParity() << noshowpos << std::endl;
   outfile << "!  Zero body term:  " << Op.ZeroBody << std::endl;
   outfile << "!  index   n   l   2j   2tz " << std::endl;

   for ( auto& it : nushell2orb )
   {
      Orbit& oi = modelspace->GetOrbit(it.second);
      outfile << "! " << setw(5) << it.first << "  " << oi.n << "  " << oi.l << " " << setw(2) << oi.j2 << " " << setw(2) << oi.tz2  << std::endl;
   }



   outfile << "!" << std::endl << "!  a   b   < a || Op || b > " << std::endl;

   for ( auto a : modelspace->valence )
   {
      int a_ind = orb2nushell[a];
      for ( auto b : modelspace->valence )
     {
        double me = Op.OneBody(a,b);
        if ( std::abs(me) < 1e-7 ) continue;
        int b_ind = orb2nushell[b];
        outfile << setw(wint) << a_ind << " " << setw(wint) << b_ind << " " << fixed << setw(wdouble) << setprecision(pdouble) <<  me << std::endl;
     }
   }
}

void ReadWrite::WriteTensorTwoBody(std::string filename, Operator& Op, std::string opname)
{
   std::ofstream outfile(filename);
   ModelSpace * modelspace = Op.GetModelSpace();
   int wint = 4; // width for printing integers

   int wdouble = 16; // width for printing doubles
   int pdouble = 9; // precision for printing doubles

   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   std::vector<int> valence_protons(modelspace->valence.size());
   std::vector<int> valence_neutrons(modelspace->valence.size());
   auto it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->proton_orbits.begin(), modelspace->proton_orbits.end(),valence_protons.begin());
   valence_protons.resize(it-valence_protons.begin());
   it = set_intersection(modelspace->valence.begin(), modelspace->valence.end(), modelspace->neutron_orbits.begin(), modelspace->neutron_orbits.end(),valence_neutrons.begin());
   valence_neutrons.resize(it-valence_neutrons.begin());
   // construct conversion maps from local orbit index to nushell index
   map<int,int> orb2nushell;
   map<int,int> nushell2orb;
   int counter = 1;
   // protons first
   for ( auto i : valence_protons ) orb2nushell[i] = counter++;
   // then neutrons
   for ( auto i : valence_neutrons ) orb2nushell[i] = counter++;
   for ( auto& it : orb2nushell) nushell2orb[it.second] = it.first;



   outfile << fixed << setprecision(pdouble);
   outfile << "!  Two-body matrix elements for tensor operator: " << opname << "   generated with IM-SRG" << std::endl;
   outfile << "!  Rank_J :  " << Op.GetJRank() << std::endl;
   outfile << "!  Rank_T :  " << Op.GetTRank() << std::endl;
   outfile << "!  Parity :  " << showpos << 1-2*Op.GetParity() << noshowpos << std::endl;
   outfile << "!  Zero body term:  " << Op.ZeroBody << std::endl;
   outfile << "!  index   n   l   2j   2tz " << std::endl;


   for ( auto& it : nushell2orb )
   {
      Orbit& oi = modelspace->GetOrbit(it.second);
      outfile << "! " << setw(5) << it.first << "  " << oi.n << "  " << oi.l << " " << setw(2) << oi.j2  << " " << setw(2) << oi.tz2  << std::endl;
   }



   outfile << "!  a    b    c    d     Jab  Jcd   <ab Jab || Op || cd Jcd>" << std::endl;
   for ( auto& itmat : Op.TwoBody.MatEl )
   {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(itmat.first[0]);
     TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(itmat.first[1]);
     auto& matrix = itmat.second;
     for (auto& ibra: tbc_bra.GetKetIndex_vv() )
     {
       Ket& bra = tbc_bra.GetKet(ibra);
       int a_ind = orb2nushell[bra.p];
       int b_ind = orb2nushell[bra.q];
//       for ( int iket=0; iket<nkets; ++iket)
       for (auto& iket: tbc_ket.GetKetIndex_vv() )
       {
         double me = matrix(ibra,iket);
         if (std::abs(me) < 1e-7) continue;
         Ket& ket = tbc_ket.GetKet(iket);
         int c_ind = orb2nushell[ket.p];
         int d_ind = orb2nushell[ket.q];
         outfile << setw(wint) << a_ind << " " << setw(wint) << b_ind << " " << setw(wint) << c_ind << " " << setw(wint) << d_ind << "   "
                 << setw(wint) << tbc_bra.J << " " << setw(wint) << tbc_ket.J << "   " << setw(wdouble) << setprecision(pdouble) << me << std::endl;
         
       }
     }

   }

}


void ReadWrite::ReadTwoBodyEngel(std::string filename, Operator& Op)
{
  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    std::ifstream infile(filename, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadTwoBodyEngel_from_stream(zipstream, Op);
  }
  else
  {
    std::ifstream infile(filename);
    ReadTwoBodyEngel_from_stream(infile, Op);
  }
}

// These are double beta decay matrix elements, so they are angular momentum scalars
// and have delta Tz=2
void ReadWrite::ReadTwoBodyEngel_from_stream( istream& infile, Operator& Op)
{
  int a,b,c,d,J;
  double tbme;
  ModelSpace* modelspace = Op.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  while( infile >> a >> b >> c >> d >> J >> tbme )
  {
    a = a-2;
    b = b-2;
    c = c-1;
    d = d-1;
    if (a >= norb) break;
    if (b >= norb or c>=norb or d>=norb) continue;
    Op.TwoBody.SetTBME_J(J,a,b,c,d,tbme);
  }
}
/*
void ReadWrite::ReadRelCMOpFromJavier( std::string statefilename, std::string MEfilename, Operator& Op)
{
  std::ifstream statefile(statefilename);
  std::ifstream MEfile(MEfilename);
  if (not MEfile.good() )
  {
    std::cout << "Trouble reading " << MEfilename << std::endl;
  }


  // First, read in the file which lists the state labelling.
  struct state_t {  int e12; int n; int N; int J; int S; int L; int lam; int LAM; int T; int Tz; };
  state_t tmp_state;            // temporary struct to read the data in.
  std::vector<state_t> statelist(1); // initialize with an empty state, since Javier starts indexing at 1 (boo Fortan).
  statefile.ignore(500,'\n');   // skip the first line
  int index;
  std::cout << "Read State File..." << std::endl;
  while( statefile >> index >> tmp_state.e12 >> tmp_state.n >> tmp_state.N >> tmp_state.J >> tmp_state.S >> tmp_state.L >> tmp_state.lam >> tmp_state.LAM >> tmp_state.T >> tmp_state.Tz )
  {
     statelist.push_back(tmp_state); // add it to the list
     std::cout << index << " " << tmp_state.e12 << " " << tmp_state.n << " " << tmp_state.N << " " << tmp_state.J << " " << tmp_state.S << " " << tmp_state.L << " " << tmp_state.lam << " " << tmp_state.LAM << " " << tmp_state.T << " " << tmp_state.Tz << std::endl;
  }
  statefile.close();


  // Second, read the file which contains the relative/cm matrix elements.
  int bra_index, ket_index;
  double MErel,MEcm;
  for (int i=0;i<6;i++) MEfile.ignore(500,'\n'); // skip header info, since I don't do anything with it (for now).
  std::cout << "Read ME File..." << std::endl;
  while( MEfile >> bra_index >> ket_index >> MErel >> MEcm )
  {
    state_t& bra = statelist[bra_index];
    state_t& ket = statelist[ket_index];
    std::cout << bra.n << " " <<  bra.lam << " " <<  bra.N << " " <<  bra.LAM << " " <<  bra.L << " " <<  bra.S << " " <<  bra.J << " " <<  bra.T << " " <<  bra.Tz << " " 
         << ket.n << " " <<  ket.lam << " " <<  ket.N << " " <<  ket.LAM << " " <<  ket.L << " " <<  ket.S << " " <<  ket.J << " " <<  ket.T << " " <<  ket.Tz << " " 
         << "(" << bra_index << "," << ket_index << ") " << MErel << " " <<  MEcm << std::endl;
    Op.TwoBody.AddToTBME_RelCM(bra.n, bra.lam, bra.N, bra.LAM, bra.L, bra.S, bra.J, bra.T, bra.Tz, 
                               ket.n, ket.lam, ket.N, ket.LAM, ket.L, ket.S, ket.J, ket.T, ket.Tz, MErel, MEcm);
  }

}
*/


struct javier_state_t {  int e12; int n; int N; int J; int S; int L; int lam; int LAM; int T; int Tz;
                  javier_state_t(){};
                  javier_state_t(int e12_in, int n_in, int N_in, int J_in, int S_in, int L_in, int lam_in, int LAM_in, int T_in, int Tz_in) :
                    e12(e12_in), n(n_in), N(N_in), J(J_in), S(S_in), L(L_in), lam(lam_in), LAM(LAM_in), T(T_in), Tz(Tz_in) {};
};

// annoyingly have to provide this as well
bool operator ==( const javier_state_t& st1, const javier_state_t& st2 )
{
  return (st1.e12==st2.e12 and st1.n==st2.n and st1.N==st2.N and st1.J==st2.J and st1.S==st2.S and st1.L==st2.L and st1.lam==st2.lam and st1.T==st2.T and st1.Tz==st2.Tz);
}

ostream& operator<< (ostream& stream, const javier_state_t& st)
{
  return stream << "e12=" << st.e12
         << ", n=" << st.n
         << ", N=" << st.N
         << ", J=" << st.J
         << ", S=" << st.S
         << ", L=" << st.L
         << ", lam=" << st.lam
         << ", LAM=" << st.LAM
         << ", T=" << st.T
         << ", Tz=" << st.Tz;
}

namespace std{
template<>
struct hash<javier_state_t>
{
  size_t operator()(const javier_state_t& st) const
  {
    index_t key = ((index_t)st.e12  <<32)^ // can be 0...32 -> 5 bits
                  ((index_t)st.n    <<28)^ // can be 0...16 -> 4 bits
                  ((index_t)st.N    <<24)^ // can be 0...16 -> 4 bits
                  ((index_t)st.J    <<19)^ // can be 0...32 -> 5 bits
                  ((index_t)st.S    <<18)^ // can be 0,1    -> 1 bit
                  ((index_t)st.L    <<13)^ // can be 0...32 -> 5 bits
                  ((index_t)st.lam  << 8)^ // can be 0...32 -> 5 bits
                  ((index_t)st.LAM  << 3)^ // can be 0...32 -> 5 bits
                  ((index_t)st.T    << 2)^ // can be 0,1    -> 1 bit
                  ((index_t)(st.Tz+1)  );  // can be 0,1,2  -> 2 bits
     return key;
  }
};
}

void ReadWrite::ReadRelCMOpFromJavier( std::string statefilename, std::string MEfilename, Operator& Op)
{
  std::ifstream statefile;
  std::ifstream MEfile;
  statefile.open(statefilename);
  MEfile.open(MEfilename);
  if (not MEfile.good() )
  {
    std::cout << "Trouble reading " << MEfilename << std::endl;
  }

  ModelSpace* modelspace = Op.GetModelSpace();

  std::cout << "Reading " << statefilename << std::endl;
  // First, read in the file which lists the state labelling.
  javier_state_t tmp_state;            // temporary struct to read the data in.
  std::vector<javier_state_t> statelist(1); // initialize with an empty state, since Javier starts indexing at 1 (boo Fortan).
  std::unordered_map<javier_state_t,index_t> statemap; // initialize with an empty state, since Javier starts indexing at 1 (boo Fortan).
  statefile.ignore(500,'\n');   // skip the first line
  index_t index;
  while( statefile >> index >> tmp_state.e12 >> tmp_state.n >> tmp_state.N >> tmp_state.J >> tmp_state.S >> tmp_state.L >> tmp_state.lam >> tmp_state.LAM >> tmp_state.T >> tmp_state.Tz )
  {
     tmp_state.Tz *=-1;
     statelist.push_back(tmp_state); // add it to the list
     statemap[tmp_state] = index; // add it to the map
  }
  statefile.close();


  std::cout << "Reading " << MEfilename << std::endl;
  // Second, read the file which contains the relative/cm matrix elements.
  arma::sp_mat matrel(statelist.size(), statelist.size());
  arma::sp_mat matcm(statelist.size(), statelist.size());
  int bra_index, ket_index;
  double MErel,MEcm;
  for (int i=0;i<6;i++) MEfile.ignore(500,'\n'); // skip header info, since I don't do anything with it (for now).
  std::cout << "Read ME File..." << std::endl;
  while( MEfile >> bra_index >> ket_index >> MErel >> MEcm )
  {
    matrel(bra_index,ket_index) = MErel;
    matcm(bra_index,ket_index) = MEcm;
  }

  std::cout << "Filling Op" << std::endl;
  // Finally, loop over all the lab frame matrix elements, and fill them.
  // This is the dumbest, most straightforward way I can think of doing it.
  for ( auto& itmat : Op.TwoBody.MatEl )
  {
    index_t ch_bra = itmat.first[0];
    index_t ch_ket = itmat.first[1];
    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int nkets_bra = tbc_bra.GetNumberKets();
    int nkets_ket = tbc_ket.GetNumberKets();
    int Jab = tbc_bra.J;
    int Jcd = tbc_ket.J;
    int Tzab = tbc_bra.Tz;
    int Tzcd = tbc_ket.Tz;
    for (int ibra=0; ibra<nkets_bra; ++ibra)
    {
      std::cout << "Getting Ket ibra=" << ibra << std::endl;
      Ket& bra = tbc_bra.GetKet(ibra);
      int na = bra.op->n;
      int la = bra.op->l;
      double ja = 0.5*bra.op->j2;
      double ta = 0.5*bra.op->tz2;
      int nb = bra.oq->n;
      int lb = bra.oq->l;
      double jb = 0.5*bra.oq->j2;
      double tb = 0.5*bra.oq->tz2;
      int Lab_min = std::max(std::abs(la-lb),Jab-1);
      int Lab_max = std::min(la+lb,Jab+1);
      int eab = 2*(na+nb)+la+lb;
      std::cout << "eab =  " << eab << std::endl;
      for (int iket=0;iket<nkets_ket; ++iket)
      {
        std::cout << "Getting Ket iket=" << iket << std::endl;
        Ket& ket = tbc_ket.GetKet(iket);
        int nc = ket.op->n;
        int lc = ket.op->l;
        double jc = 0.5*ket.op->j2;
        double tc = 0.5*ket.op->tz2;
        int nd = ket.oq->n;
        int ld = ket.oq->l;
        double jd = 0.5*ket.oq->j2;
        double td = 0.5*ket.oq->tz2;
        int Lcd_min = std::max(std::abs(lc-ld),Jcd-1);
        int Lcd_max = std::min(lc+ld,Jcd+1);
        int ecd = 2*(nc+nd)+lc+ld;
        std::cout << " ecd = " << ecd << std::endl;
        if (eab==0 and ecd==0 and Jab+Jcd==1 and Tzab==0 and Tzcd==0)
        {
          std::cout << "!!!!!!!!!!!!!!!!!! HERE !!!!!!!!!!!!!!!!!!!!" << std::endl;
          std::cout << "!! " << Lab_min << "<= Lab <= " << Lab_max  << std::endl;
          std::cout << "!! " << Lcd_min << "<= Lcd <= " << Lcd_max  << std::endl;
          std::cout << "!!!!!!!!!!!!!!!!!! HERE !!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
        for (int Lab=Lab_min; Lab<=Lab_max; ++Lab)
        {
          for (int Sab=std::max(0,std::abs(Lab-Jab)); Sab<=std::min(1,Lab+Jab); ++Sab)
          {
            double NormNineJab = sqrt((2*ja+1)*(2*jb+1)*(2*Lab+1)*(2*Sab+1)) * modelspace->GetNineJ(la,0.5,ja, lb,0.5,jb, Lab,Sab,Jab);
            for (int Lcd=Lcd_min; Lcd<=Lcd_max; ++Lcd)
            {
              for (int Scd=std::max(0,std::abs(Lcd-Jcd)); Scd<=std::min(1,Lcd+Jcd); ++Scd)
              {
                double NormNineJcd = sqrt((2*jc+1)*(2*jd+1)*(2*Lcd+1)*(2*Scd+1)) * modelspace->GetNineJ(lc,0.5,jc, ld,0.5,jd, Lcd,Scd,Jcd);

                // loop over rel/cm quantum numbers
                for (int n_ab=0; 2*n_ab<=eab; ++n_ab)
                {
                  for (int N_ab=0; 2*(N_ab+n_ab)<=eab; ++N_ab)
                  {
                    for (int lam_ab=0; 2*(N_ab+n_ab)+lam_ab<=eab; ++lam_ab)
                    {
                      int LAM_ab = eab-2*(N_ab+n_ab)-lam_ab;
                      if ((lam_ab + LAM_ab<Lab) or (std::abs(lam_ab-LAM_ab)>Lab)) continue;
                      double mosh_ab = modelspace->GetMoshinsky( N_ab, LAM_ab, n_ab, lam_ab, na, la, nb, lb, Lab);
                      int Tab = (lam_ab + Sab +1 )%2;
                      if (std::abs(Tzab)>Tab) continue;

                      for (int n_cd=0; 2*n_cd<=ecd; ++n_cd)
                      {
                        for (int N_cd=0; 2*(N_cd+n_cd)<=ecd; ++N_cd)
                        {
                          for (int lam_cd=0; 2*(N_cd+n_cd)+lam_cd<=ecd; ++lam_cd)
                          {
                            int LAM_cd = ecd-2*(N_cd+n_cd)-lam_cd;
                            if ((lam_cd + LAM_cd<Lcd) or (std::abs(lam_cd-LAM_cd)>Lcd)) continue;
                            double mosh_cd = modelspace->GetMoshinsky( N_cd, LAM_cd, n_cd, lam_cd, nc, lc, nd, ld, Lcd);
                            int Tcd = (lam_cd + Scd + 1)%2;
                            if (std::abs(Tzcd)>Tcd) continue;
                            double IsospinClebsch_ab = AngMom::CG(0.5,ta,0.5,tb, Tab,Tzab);
                            double IsospinClebsch_cd = AngMom::CG(0.5,tc,0.5,td, Tcd,Tzcd);
                            double coeff = NormNineJab*NormNineJcd*mosh_ab*mosh_cd*IsospinClebsch_ab*IsospinClebsch_cd;
                            size_t rel_index_bra = statemap[ javier_state_t(eab, n_ab, N_ab, Jab, Sab, Lab, lam_ab, LAM_ab, Tab, Tzab) ];
                            size_t rel_index_ket = statemap[ javier_state_t(ecd, n_cd, N_cd, Jcd, Scd, Lcd, lam_cd, LAM_cd, Tcd, Tzcd) ];
                            if (rel_index_bra<1) continue;
                            if (rel_index_ket<1) continue;
                            std::cout << "ab: " << javier_state_t(eab, n_ab, N_ab, Jab, Sab, Lab, lam_ab, LAM_ab, Tab, Tzab) << std::endl;
                            std::cout << "cd: " << javier_state_t(ecd, n_cd, N_cd, Jcd, Scd, Lcd, lam_cd, LAM_cd, Tcd, Tzcd) << std::endl;
                            std::cout << "ibra,iket = " << ibra << " , " << iket << "   rel_index_bra,rel_index_ket = " << rel_index_bra << " , " << rel_index_ket << std::endl;
                            if ((N_ab==N_cd) and (LAM_ab==LAM_cd))
                               itmat.second(ibra,iket) += coeff * matrel(rel_index_bra, rel_index_ket);
                            if ((n_ab==n_cd) and (lam_ab==lam_cd))
                               itmat.second(ibra,iket) += coeff * matcm( rel_index_bra, rel_index_ket);
                          }
                        }
                      }
                    }
                  }
                }

              }
            }
          }
        }
        
        if (eab==0 and ecd==0 and Jab+Jcd==1 and Tzab==0 and Tzcd==0) std::cout << "!!!!!!!!!!!!!!!!!! DONE !!!!!!!!!!!!!!!!!!!!" << std::endl;

      }
    }
  }
  std::cout << "Done!" << std::endl;

}


void ReadWrite::SetLECs(double c1, double c3, double c4, double cD, double cE)
{
   LECs = {c1,c3,c4,cD,cE};
}

void ReadWrite::SetLECs_preset(std::string key)
{
       if (key == "EM1.8_2.0")  LECs = {-0.81, -3.20, 5.40, 1.264, -0.120};
  else if (key == "EM2.0_2.0")  LECs = {-0.81, -3.20, 5.40, 1.271, -0.131};
  else if (key == "EM2.2_2.0")  LECs = {-0.81, -3.20, 5.40, 1.214, -0.137};
  else if (key == "EM2.8_2.0")  LECs = {-0.81, -3.20, 5.40, 1.278, -0.078};
  else if (key == "PWA2.0_2.0") LECs = {-0.76, -4.78, 3.96,-3.007, -0.686};
  else if (key == "N2LOSAT")    LECs = {-1.12152120, -3.92500586, 3.76568716, 0.861680589, -0.03957471}; // For testing purposes only. (This uses the wrong regulator).
}




