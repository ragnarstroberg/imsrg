#include "ReadWrite.hh"
#include "ModelSpace.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "string.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "H5Cpp.h"

#define LINESIZE 496
//#define HEADERSIZE 500
#define HEADERSIZE 255
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif
#ifndef HBARC
   #define HBARC 197.3269718 // hc in MeV * fm
#endif

using namespace std;
using namespace H5;

ReadWrite::~ReadWrite()
{
//  cout << "In ReadWrite destructor" << endl;
}

ReadWrite::ReadWrite()
: doCoM_corr(false), goodstate(true),LECs({-0.81,-3.20,5.40,1.271,-0.131}) // default to the EM2.0_2.0 LECs
{
}


void ReadWrite::ReadSettingsFile( string filename)
{
   char line[LINESIZE];
   string lstr;
   ifstream fin;
   fin.open(filename);
   if (! fin.good() )
   {
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "!!!!  Trouble reading Settings file " << filename << "!!!!!!" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      goodstate = false;
   }
   cout << "Reading settings file " << filename << endl;
   while (fin.getline(line,LINESIZE))
   {
      lstr = string(line);
      lstr = lstr.substr(0, lstr.find_first_of("#"));
      int colon = lstr.find_first_of(":");
      string param = lstr.substr(0, colon);
      if ( param.size() <1) continue;
      string value = lstr.substr(colon+1, lstr.length()-colon);
      value.erase(0,value.find_first_not_of(" \t\r\n"));
      value.erase(value.find_last_not_of(" \t\r\n")+1);
      if (value.size() < 1) continue;
      InputParameters[param] = value;
      cout << "parameter: [" << param << "] = [" << value << "]" << endl;
   }

}




/// Read two-body matrix elements from an Oslo-formatted file
void ReadWrite::ReadTBME_Oslo( string filename, Operator& Hbare)
{

  ifstream infile;
  char line[LINESIZE];
  int Tz,Par,J2,a,b,c,d;
  double fbuf[3];
  double tbme;
  int norbits = Hbare.GetModelSpace()->GetNumberOrbits();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
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

/// Read two-body matrix elements from an Oslo-formatted file
void ReadWrite::WriteTwoBody_Oslo( string filename, Operator& Op)
{

  ofstream outfile(filename);
  ModelSpace* modelspace = Op.GetModelSpace();
  int wint = 8;
  int wdouble = 12;
  int dprec = 6;
  outfile << fixed;
  outfile << setw(wint) << setprecision(dprec);

  if ( !outfile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble writing file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  outfile << "      Tz      Parity    J2        a        b        c        d     <ab|V|cd>    <ab|0|cd>    <ab|0|cd>    <ab|0|cd>" << endl;

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
        outfile << setw(wint) << Tz << " " << setw(wint) <<parity << " " << setw(wint) <<J2 << " " << setw(wint) << bra.p << " " << setw(wint) <<bra.q << " " << setw(wint) <<ket.p << " " << setw(wint) <<ket.q <<  " " << setw(wdouble) << setprecision(dprec) <<tbme << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << " " << setw(wdouble) << setprecision(dprec)<< 0.0 << endl;
      }
    }

  }

  return;
}


void ReadWrite::WriteOneBody_Oslo( string filename, Operator& Op)
{
  ofstream outfile(filename);
  ModelSpace* modelspace = Op.GetModelSpace();
  int wint = 3;
  int wdouble = 10;
  int dprec = 6;
  outfile << fixed;
  outfile << setw(wint) << setprecision(dprec);
  if ( !outfile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble writing file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  outfile << "  a     b      <a|V|b> " << endl;
  int norb = modelspace->GetNumberOrbits();
  for ( int i=0;i<norb;++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    for (int j : Op.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
    {
      outfile << setw(wint) << i << "   " << setw(wint) << j << "   " << setw(wdouble) << Op.OneBody(i,j) << endl;
    }
  }

}

void ReadWrite::ReadBareTBME_Jason( string filename, Operator& Hbare)
{

  ifstream infile;
  char line[LINESIZE];
  int Tz2,Par,J2,a,b,c,d;
  int na,nb,nc,nd;
  int la,lb,lc,ld;
  int ja,jb,jc,jd;
//  double fbuf[3];
  double tbme;
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norbits = modelspace->GetNumberOrbits();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }

//  infile.getline(line,LINESIZE);

//  while (!strstr(line,"<ab|V|cd>") && !infile.eof()) // Skip lines until we see the header
  //for (int& i : modelspace->holes ) // skip spe's at the top
  for (int i=0;i<6;++i ) // skip spe's at the top
  {
     infile.getline(line,LINESIZE);
     cout << "skipping line: " << line << endl;
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
void ReadWrite::ReadBareTBME_Navratil( string filename, Operator& Hbare)
{
  ifstream infile(filename);  
  if ( !infile.good() )
  {
    cerr << "************************************" << endl
         << "**    Trouble reading file  !!!   **" << filename << endl
         << "************************************" << endl;
    goodstate = false;
    return;
  }
  infile.close();

  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    ifstream infile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadBareTBME_Navratil_from_stream(zipstream, Hbare);
  }
  else
  {
    ifstream infile(filename);
    ReadBareTBME_Navratil_from_stream(infile, Hbare);
  }
}

void ReadWrite::ReadBareTBME_Navratil_from_stream( istream& infile, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
//  int emax = modelspace->GetNmax();
  int norb = modelspace->GetNumberOrbits();
//  int nljmax = norb/2;
//  int herm = Hbare.IsHermitian() ? 1 : -1 ;
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2;
     orbits_remap[nlj] = i;
  }

//  ifstream infile;
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

    if (abs(vpn)>1e-6)
    {
      Hbare.TwoBody.Set_pn_TBME_from_iso(J,T,0,a,b,c,d,vpn);
    }

  }
  
}



void ReadWrite::WriteTBME_Navratil( string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2;
//     orbits_remap[nlj] = i;
     orbits_remap[i]   = nlj+1;
//     orbits_remap[i+1] = nlj+1;
  }

  ofstream outfile(filename);

  double hw = modelspace->GetHbarOmega();
  hw = Hbare.GetModelSpace()->GetHbarOmega();
  double srg_lambda = 0;
  outfile << 0 << "    " << modelspace->GetNmax() << "    " << 2*modelspace->GetNmax() << "   " << hw << "     " << srg_lambda << endl;

  outfile << setiosflags(ios::fixed);

  double trel=0, h_ho_rel=0, vcoul=0;

  for (int a=0; a<norb; a+=2)
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b=0; b<=a; b+=2)
    {
      Orbit& ob = modelspace->GetOrbit(b);
      int jab_min = abs(oa.j2-ob.j2)/2;
      int jab_max = (oa.j2+ob.j2)/2;
      for (int c=0; c<=a; c+=2)
      {
        Orbit& oc = modelspace->GetOrbit(c);
        for (int d=0; d<=c; d+=2)
        {
          Orbit& od = modelspace->GetOrbit(d);
          if ( (oa.l + ob.l + oc.l + od.l)%2 > 0) continue;
          int jcd_min = abs(oc.j2-od.j2)/2;
          int jcd_max = (oc.j2+od.j2)/2;
          int jmin = max(jab_min,jcd_min);
          int jmax = min(jab_max,jcd_max);
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
            if (abs(vpp)>1e-7 or abs(vnn)>1e-7 or abs(v10)>1e-7)
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
                    << vnn << endl;
            }
            if (abs(v00)>1e-7)
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
                    << 0.0 << endl;
            }
          }
        }
      }
    }
  }
}








/// Decide if the file is gzipped or ascii, create a stream, then call ReadBareTBME_Darmstadt_from_stream().
void ReadWrite::ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int emax, int Emax, int lmax)
{

  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    ifstream infile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadBareTBME_Darmstadt_from_stream(zipstream, Hbare,  emax, Emax, lmax);
  }
  else if (filename.substr( filename.find_last_of(".")) == ".bin")
  {
    ifstream infile(filename, ios::binary);    
    if ( !infile.good() )
    {
      cerr << "problem opening " << filename << ". Exiting." << endl;
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
    
//    vector<double> v(n_elem);
//    infile.read((char*)&v[0], n_elem*sizeof(double));
    vector<float> v(n_elem);
    infile.read((char*)&v[0], n_elem*sizeof(float));
    infile.close();
    VectorStream vectorstream(v);
    cout << "n_elem = " << n_elem << endl;
    ReadBareTBME_Darmstadt_from_stream(vectorstream, Hbare,  emax, Emax, lmax);
  }
  else
  {
    ifstream infile(filename);
    ReadBareTBME_Darmstadt_from_stream(infile, Hbare,  emax, Emax, lmax);
  }
}


/// Decide the file format from the extension -- .me3j (Darmstadt group format, human-readable), .gz (gzipped me3j, less storage),
/// .bin (me3j converted to binary, faster to read), .h5 (HDF5 format). Default is to assume .me3j.
/// For the first three, the file is converted to a stream and sent to ReadDarmstadt_3body_from_stream().
/// For the HDF5 format, a separate function is called: Read3bodyHDF5().
void ReadWrite::Read_Darmstadt_3body( string filename, Operator& Hbare, int E1max, int E2max, int E3max)
{

  double start_time = omp_get_wtime();
  string extension = filename.substr( filename.find_last_of("."));

  if (extension == ".me3j")
  {
    ifstream infile(filename);
    Read_Darmstadt_3body_from_stream(infile, Hbare,  E1max, E2max, E3max);
  }
  else if ( extension == ".gz")
  {
    ifstream infile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    Read_Darmstadt_3body_from_stream(zipstream, Hbare,  E1max, E2max, E3max);
  }
  else if (extension == ".bin")
  {
    ifstream infile(filename, ios::binary);    
    if ( !infile.good() )
    {
      cerr << "problem opening " << filename << ". Exiting." << endl;
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
    
//    vector<float> v(n_elem);
//    infile.read((char*)&v[0], n_elem*sizeof(float));
    vector<float> v(n_elem);
    infile.read((char*)&v[0], n_elem*sizeof(float));
    infile.close();
    VectorStream vectorstream(v);
    cout << "n_elem = " << n_elem <<  endl;
    Read_Darmstadt_3body_from_stream(vectorstream, Hbare,  E1max, E2max, E3max);
  }
  else if (extension == ".h5")
  {
    Read3bodyHDF5_new( filename, Hbare );
//    Read3bodyHDF5( filename, Hbare );
  }
  else
  {
    cout << "assuming " << filename << " is of me3j format ... " << endl;
    ifstream infile(filename);
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
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  vector<int> orbits_remap;

  if (emax < 0)  emax = modelspace->Nmax;
  if (lmax < 0)  lmax = emax;

  for (int e=0; e<=min(emax,modelspace->Nmax); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = abs(2*l-1);
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
    if (e1 > modelspace->Nmax) break;

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
          int Jmin = max( abs(o1.j2 - o2.j2), abs(o3.j2 - o4.j2) )/2;
          int Jmax = min (o1.j2 + o2.j2, o3.j2+o4.j2)/2;
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

//             cout << a << " " << b << " " << c << " " << d << "   " << J << "   "
//                  << fixed << setprecision(7) << setw(11) << tbme_00 << " " << tbme_nn << " " << tbme_10 << " " << tbme_pp << endl;

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



/// Read me3j format three-body matrix elements. Pass in E1max, E2max, E3max for the file, so that it can be properly interpreted.
/// The modelspace truncation doesn't need to coincide with the file truncation. For example, you could have an emax=10 modelspace
/// and read from an emax=14 file, and the matrix elements with emax>10 would be ignored.
template <class T>
void ReadWrite::Read_Darmstadt_3body_from_stream( T& infile, Operator& Hbare, int E1max, int E2max, int E3max)
{
  if ( !infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  if (Hbare.particle_rank < 3)
  {
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << endl;
    cerr << " Oops. Looks like we're trying to read 3body matrix elements to a " << Hbare.particle_rank << "-body operator. For shame..." << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << endl;
    goodstate = false;
    return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int e1max = modelspace->GetNmax();
  int e2max = modelspace->GetN2max(); // not used yet
  int e3max = modelspace->GetN3max();
  cout << "Reading 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << endl;

  vector<int> orbits_remap(0);
  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.

  for (int e=0; e<=min(E1max,e1max); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size();



  // skip the first line
  char line[LINESIZE];
  infile.getline(line,LINESIZE);


  // begin giant nested loops
  int nread = 0;
  int nkept = 0;
  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
    int a =  orbits_remap[nlj1];
    Orbit & oa = modelspace->GetOrbit(a);
    int ea = 2*oa.n + oa.l;
    if (ea > E1max) break;
    if (ea > e1max) break;
//    cout << "a = " << a << endl;
    cout << setw(5) << setprecision(2) << nlj1*(nlj1+1.)/(nljmax*(nljmax+1))*100 << " % done" << '\r';
    cout.flush();

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
        int JabMin  = abs(oa.j2 - ob.j2)/2;

        int twoJCMindownbra;
        if (abs(oa.j2 - ob.j2) >oc.j2)
           twoJCMindownbra = abs(oa.j2 - ob.j2)-oc.j2;
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
              int JJabMin = abs(od.j2 - oe.j2)/2;

              int twoJCMindownket;
              if ( abs(od.j2 - oe.j2) > of.j2 )
                 twoJCMindownket = abs(od.j2 - oe.j2) - of.j2;
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
                int twoJCMin = max( abs(2*Jab - oc.j2), abs(2*JJab - of.j2));
                int twoJCMax = min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
       
                for(int twoJC = twoJCMin; twoJC <= twoJCMax; twoJC += 2)
                {
                 for(int tab = 0; tab <= 1; tab++) // the total isospin loop can be replaced by i+=5
                 {
                  for(int ttab = 0; ttab <= 1; ttab++)
                  {
                   //summation bounds
                   int twoTMin = 1; // twoTMin can just be used as 1
                   int twoTMax = min( 2*tab +1, 2*ttab +1);
       
                   for(int twoT = twoTMin; twoT <= twoTMax; twoT += 2)
                   {
//                    double V;
                    float V = 0;
                    infile >> V;
                    ++nread;
                    bool autozero = false;


                    if ( a==b and (tab+Jab)%2==0 ) autozero = true;
                    if ( d==e and (ttab+JJab)%2==0 ) autozero = true;
                    if ( a==b and a==c and twoT==3 and oa.j2<3 ) autozero = true;
                    if ( d==e and d==f and twoT==3 and od.j2<3 ) autozero = true;

                       if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
                          and (ea+eb+ec<=e3max) and (ed+ee+ef<=e3max) )
                       {
                         ++nkept;
                       }

                    if (not autozero and abs(V)>1e-5)
                    {
//                       double V0 = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
//                       V0 = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
                       if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
                          and (ea+eb+ec<=e3max) and (ed+ee+ef<=e3max) )
                       {
//                        cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << Jab << " " << JJab << " " << twoJC << " " << tab << " " << ttab << " " << twoT << " " << V << endl;
                        Hbare.ThreeBody.SetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f, V);
                       }
                    }

                    if (autozero)
                    {
//                       cout << " ( should be zero ) ";
                       if (abs(V) > 1e-6 and ea<=e1max and eb<=e1max and ec<=e1max)
                       {
                          cout << " <-------- AAAAHHHH!!!!!!!! Reading 3body file and this should be zero, but it's " << V << endl;
                          goodstate = false;
                          return;
                       }
                    }
 //                   cout << endl;
       
                   }//twoT
                  }//ttab
                 }//tab
//                 cout << " ------------------------" << endl;
                if (not infile.good() ) break;
                }//twoJ
               }//JJab
       
              }//Jab

            }
          }
        }
      }
    }
  }
  cout << "Read in " << nread << " floating point numbers (" << nread * sizeof(float)/1024./1024./1024. << " GB)" << endl;
  cout << "Stored " << nkept << " floating point numbers (" << nkept * sizeof(float)/1024./1024./1024. << " GB)" << endl;

}




/// Read three-body basis from HDF5 formatted file. This routine was ported to C++ from
/// a C routine by Heiko Hergert, with as little modification as possible.
void ReadWrite::GetHDF5Basis( ModelSpace* modelspace, string filename, vector<array<int,5>>& Basis)
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
     cerr << "Error: Failed to read dataset dimensions!" << endl;
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
      cerr << "Error. alpha != alphap " << endl;
      return;
    }

    Basis[alpha] = {o1,o2,o3,J12,twoJ};
    
  }

  delete[] dbuf[0];
  delete[] dbuf;

}


/// Read three-body matrix elements from HDF5 formatted file. This routine was ported to C++ from
/// a C routine by Heiko Hergert, with as little modification as possible.
void ReadWrite::Read3bodyHDF5( string filename,Operator& op )
{

  const int SLABSIZE = 10000000;

  ModelSpace* modelspace = op.GetModelSpace();
  vector<array<int,5>> Basis;
  GetHDF5Basis(modelspace, filename, Basis);


  H5File file(filename, H5F_ACC_RDONLY);
  DataSet label = file.openDataSet("vtnf_labels");
  DataSpace label_dspace = label.getSpace();
  DataSet value = file.openDataSet("vtnf");
  DataSpace value_dspace = value.getSpace();

  int label_nDim = label_dspace.getSimpleExtentNdims();
  if (label_nDim != 2)
  {
    cerr << "Error. Expected label_nDim==2, but got << label_nDim." << endl;
    return;
  }
  hsize_t label_maxDim[2];
  int label_status = label_dspace.getSimpleExtentDims(label_maxDim,NULL);
  if (label_status != label_nDim)
  {
    cerr << "Error. failed to read dataset dimensions for label." << endl;
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
    cerr << "Error. Expected value_nDim==2, but got << value_nDim." << endl;
    return;
  }
  hsize_t value_maxDim[2];
  int value_status = value_dspace.getSimpleExtentDims(value_maxDim,NULL);
  if (value_status != value_nDim)
  {
    cerr << "Error. failed to read dataset dimensions for value." << endl;
    return;
  }
  
  hsize_t value_curDim[2];
  value_curDim[0] = min(SLABSIZE, int(value_maxDim[0]));
  value_curDim[1] = 1 ;

  DataSpace value_buf_dspace(1,value_curDim);
  
  // Generate a 1d buffer in contiguous memory, also known as an array...
  double *value_buf = new double[value_curDim[0]];

  // break the file into slabs for reading
  int nSlabs = (int)((double)value_maxDim[0]/((double)SLABSIZE)) + 1;

  hsize_t stride[2] = {1,1};
  hsize_t count[2] = {1,1};

  // loop through the slabs
  for ( int n=0; n<nSlabs; ++n)
  {
    hsize_t start[2] = { n*value_curDim[0], 0};
    hsize_t label_block[2];
    hsize_t value_block[2];
    if (n==nSlabs-1)
    {
      label_block[0] = label_maxDim[0]-(nSlabs-1)*SLABSIZE;
      label_block[1] = label_maxDim[1];
      value_block[0] = value_maxDim[0]-(nSlabs-1)*SLABSIZE;
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
         cerr << "Error. Mismatching parity !  < "  << parity_abc << " " << parity_def << " " << Pi << "    " << endl;
       }
       if (J2 != twoJ or J2p != twoJ)
       {
         cerr << "Error. Mismatching total J! " << J2 << " " << J2p << " " << twoJ << "   alphas = " << alpha << ", " << alphap << endl;
       }

       me *= 0.5; // According to Heiko, this shouldn't be here. But comparing matrix elements with Petr's interaction suggests otherwise.
//       if (a!=d or b!=e or c!=f) me *=0.5;
//       if (alpha<50 and alphap<50)
//       if (a<5 and b<5 and c<5 and d<5 and e<5 and f<5)
       me *= modelspace->phase(oa.n+ob.n+oc.n+od.n+oe.n+of.n); // shamelessly copying Heiko. I don't understand this.
//       if (((alpha+1)/2==1 and (alphap/2==127)) or ( (alphap+1)/2==1 and alpha/2==127) )
//       if (((alpha+1)/2==1 or alpha/2==127) and (alphap/2==127 or (alphap+1)/2==1) )
//        {
//         double previous = op.ThreeBody.GetME(J12,JJ12,twoJ,T12,TT12,twoT,a,b,c,d,e,f);
//         double newme = op.ThreeBody.AddToME(J12,JJ12,twoJ,T12,TT12,twoT,a,b,c,d,e,f,me);
//         cout << "*** " << alpha << "  " << alphap << "  "
//              << a << "-" << b << "-" << c << " " << d << "-" << e << "-" << f << "  "
//              << T12 << "  " << TT12 << "  " << twoT << "  "
//              << J12 << "  " << JJ12 << "  " << twoJ << "  " << me << "     "
//              << previous << "  =>  "
//              <<  newme
//             << endl;
//       }
//       else
       op.ThreeBody.SetME(J12,JJ12,twoJ,T12,TT12,twoT,a,b,c,d,e,f, me);
       if (a==d and b==e and c==f and ( J12!=JJ12 ) )
//       if (a==d and b==e and c==f and ( J12!=JJ12 or T12 != TT12) )
          op.ThreeBody.SetME(JJ12,J12,twoJ,TT12,T12,twoT,a,b,c,d,e,f, me);

    } //loop over matrix elements
  } // loop over slabs
  delete[] label_buf[0];
  delete[] label_buf;
  delete[] value_buf;
  cout << "Writing me3j file..." << endl;
  Write_me3j(filename + "_to_me3j", op, 2, 24, 12);
  cout << "done" << endl;
}






// THIS ONE SEEMS TO WORK, SO FAR

void ReadWrite::Read3bodyHDF5_new( string filename,Operator& op )
{

  ModelSpace* modelspace = op.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();

  int t12p_list[5] = {0,0,1,1,1};
  int t12_list[5]  = {0,1,0,1,1};
  int twoT_list[5] = {1,1,1,1,3};

  H5File file(filename, H5F_ACC_RDONLY);

  DataSet basis = file.openDataSet("alphas");
  DataSpace basis_dspace = basis.getSpace();
  hsize_t iDim_basis[6];
  basis_dspace.getSimpleExtentDims(iDim_basis,NULL);


  // Generate a 2d buffer in contiguous memory
  //TODO: Do this with vectors so we don't have to worry about new and delete[]
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

  int i=-5;
  for (int alphaspp=0;alphaspp<alpha_max;++alphaspp)
  {
    int lap = dbuf[alphaspp][2];
    int lbp = dbuf[alphaspp][5];
    int lcp = dbuf[alphaspp][8];

    int ap = modelspace->GetOrbitIndex(dbuf[alphaspp][1],dbuf[alphaspp][2],dbuf[alphaspp][3],-1);
    int bp = modelspace->GetOrbitIndex(dbuf[alphaspp][4],dbuf[alphaspp][5],dbuf[alphaspp][6],-1);
    int cp = modelspace->GetOrbitIndex(dbuf[alphaspp][7],dbuf[alphaspp][8],dbuf[alphaspp][9],-1);
    int j12p = dbuf[alphaspp][10];
    int jtotp = dbuf[alphaspp][11];

    for (int alphasp=alphaspp; alphasp<alpha_max;++alphasp)
    {
      int la = dbuf[alphasp][2];
      int lb = dbuf[alphasp][5];
      int lc = dbuf[alphasp][8];
      int a = modelspace->GetOrbitIndex(dbuf[alphasp][1],dbuf[alphasp][2],dbuf[alphasp][3],-1);
      int b = modelspace->GetOrbitIndex(dbuf[alphasp][4],dbuf[alphasp][5],dbuf[alphasp][6],-1);
      int c = modelspace->GetOrbitIndex(dbuf[alphasp][7],dbuf[alphasp][8],dbuf[alphasp][9],-1);
      int j12 = dbuf[alphasp][10];
      int jtot = dbuf[alphasp][11];
      if (jtot != jtotp or (lap+lbp+lcp+la+lb+lc)%2>0) continue;
      i+=5;
      
      for (hsize_t k_iso=0;k_iso<5;++k_iso)
      {
       int T12  = t12p_list[k_iso];
       int TT12 = t12_list[k_iso];
       int twoT = twoT_list[k_iso];
       float *me = value_buf[i+k_iso];
       float summed_me = 0;
       for (int ii=0;ii<5;++ii) summed_me += LECs[ii] * me[ii] ;
       summed_me *= HBARC;
       summed_me *= modelspace->phase(dbuf[alphasp][1]+dbuf[alphasp][4]+dbuf[alphasp][7]+dbuf[alphaspp][1]+dbuf[alphaspp][4]+dbuf[alphaspp][7]);

       if ( ap==bp and (j12p+T12)%2 !=1 )
       {
         if ( abs(summed_me)>1.0e-6  )
         {
           cout << "AAHH!!  should be zero!" << endl;
         }
       }
       else if ( a==b  and (j12+TT12)%2 !=1 )
       {
         if ( abs(summed_me)>1.0e-6  )
         {
           cout << "AAHH!!  should be zero!" << endl;
         }
       }
       else if (a<norb and b<norb and c<norb and ap<norb and bp<norb and cp<norb)
       {
        op.ThreeBody.SetME(j12p,j12,jtot,T12,TT12,twoT,ap,bp,cp,a,b,c, summed_me);
        if (a==ap and b==bp and c==cp and j12 != j12p)
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




void ReadWrite::Write_me2j( string outfilename, Operator& Hbare, int emax, int Emax, int lmax)
{
  ofstream outfile(outfilename);
  if ( !outfile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble opening file  !!!   **" << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  vector<int> orbits_remap;

  if (emax < 0)  emax = modelspace->Nmax;
  if (lmax < 0)  lmax = emax;

  for (int e=0; e<=min(emax,modelspace->Nmax); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = abs(2*l-1);
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
//  char line[LINESIZE];
//  infile.getline(line,LINESIZE);
  outfile << "    blah blah blah header " << endl;
  int icount = 0;

  outfile << setiosflags(ios::fixed);
//  outfile << setprecision(7);
//  outfile << setw(12);

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
          int Jmin = max( abs(o1.j2 - o2.j2), abs(o3.j2 - o4.j2) )/2;
          int Jmax = min (o1.j2 + o2.j2, o3.j2+o4.j2)/2;
          if (Jmin > Jmax) continue;
          for (int J=Jmin; J<=Jmax; ++J)
          {
             // me2j format is unnormalized
             double norm_factor = 1;
             if (a==b)  norm_factor *= SQRT2;
             if (c==d)  norm_factor *= SQRT2;

             // Matrix elements are written in the file with (T,Tz) = (0,0) (1,1) (1,0) (1,-1)
             tbme_pp = Hbare.TwoBody.GetTBME(J,parity,-1,a,b,c,d);
             tbme_nn = Hbare.TwoBody.GetTBME(J,parity,1,a+1,b+1,c+1,d+1);
             tbme_10 = Hbare.TwoBody.Get_iso_TBME_from_pn(J,1,0,a,b,c,d);
             tbme_00 = Hbare.TwoBody.Get_iso_TBME_from_pn(J,0,0,a,b,c,d);

             
             outfile << setprecision(7) << setw(12) << tbme_00*norm_factor<< " "  ;
             if ((icount++)%10==9)
             {
               outfile << endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_nn<< " "  ;
             if ((icount++)%10==9)
             {
               outfile << endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_10*norm_factor << " "  ;
             if ((icount++)%10==9)
             {
               outfile << endl;
             }
             outfile << setprecision(7) << setw(12) << tbme_pp << " "  ;
             if ((icount++)%10==9)
             {
               outfile << endl;
             }


          }
        }
      }
    }
  }
  if (icount%10 !=9) outfile << endl;

}






void ReadWrite::Write_me3j( string ofilename, Operator& Hbare, int E1max, int E2max, int E3max)
{
  ofstream outfile(ofilename);
//  if ( !infile.good() )
//  {
//     cerr << "************************************" << endl
//          << "**    Trouble reading file  !!!   **" << endl
//          << "************************************" << endl;
//     goodstate = false;
//     return;
//  }
  if (Hbare.particle_rank < 3)
  {
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << endl;
    cerr << " Oops. Looks like we're trying to write 3body matrix elements from a " << Hbare.particle_rank << "-body operator. For shame..." << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! << " << endl;
    goodstate = false;
    return;
  }
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int e1max = modelspace->GetNmax();
  int e2max = modelspace->GetN2max(); // not used yet
  int e3max = modelspace->GetN3max();
  cout << "Writing 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << endl;

  vector<int> orbits_remap(0);
  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.

  for (int e=0; e<=min(E1max,e1max); ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size();



  // skip the first line
//  char line[LINESIZE];
//  infile.getline(line,LINESIZE);
//  int useless_counter=0;
  outfile << "         header nonsense... " << endl;
  outfile << setiosflags(ios::fixed);
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
        int JabMin  = abs(oa.j2 - ob.j2)/2;

        int twoJCMindownbra;
        if (abs(oa.j2 - ob.j2) >oc.j2)
           twoJCMindownbra = abs(oa.j2 - ob.j2)-oc.j2;
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
              int JJabMin = abs(od.j2 - oe.j2)/2;

              int twoJCMindownket;
              if ( abs(od.j2 - oe.j2) > of.j2 )
                 twoJCMindownket = abs(od.j2 - oe.j2) - of.j2;
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
                int twoJCMin = max( abs(2*Jab - oc.j2), abs(2*JJab - of.j2));
                int twoJCMax = min( 2*Jab + oc.j2 , 2*JJab + of.j2 );
       
                for(int twoJC = twoJCMin; twoJC <= twoJCMax; twoJC += 2)
                {
//                 ++useless_counter;
//                 if (a<5 and b<5 and c<5 and d<5 and e<5 and f<5)
//                 {
//                 cout << "#" << useless_counter << "  " << a << "-" << b << "-" << c << "-" << "  " << d << "-" << e << "-" << f << "  "
//                      << Jab << " " << JJab << " " << twoJC << endl;
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
                    float V = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
                    outfile << setprecision(7) << setw(12) << V << " "  ;
                    if ((icount++)%10==9)
                      outfile << endl;
                   }//twoT
                  }//ttab
                 }//tab
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
  if (icount%10 !=9) outfile << endl;

}







void ReadWrite::WriteNuShellX_op(Operator& op, string filename)
{
  WriteNuShellX_intfile(op,filename,"op");
}
void ReadWrite::WriteNuShellX_int(Operator& op, string filename)
{
  WriteNuShellX_intfile(op,filename,"int");
}

/// Write the valence space part of the interaction to a NuShellX *.int file.
/// Note that NuShellX assumes identical orbits for protons and neutrons,
/// so that the pnpn interaction should be equal to the pnnp interaction.
/// This is only approximately true for interactions generated with IMSRG,
/// so some averaging is required.
//void ReadWrite::WriteNuShellX_intfile(Operator& op, string filename, string mode)
//{
//   ofstream intfile;
//   intfile.open(filename, ofstream::out);
//   ModelSpace * modelspace = op.GetModelSpace();
//   int nvalence_proton_orbits = 0;
//   int proton_core_orbits = 0;
//   int neutron_core_orbits = 0;
//   int Acore = 0;
//   int wint = 4; // width for printing integers
//   int wdouble = 12; // width for printing doubles
//   int pdouble = 6; // precision for printing doubles
//   for (auto& i : modelspace->hole_qspace)
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      Acore += oi.j2 +1;
//      if (oi.tz2 < 0)
//      {
//         proton_core_orbits += 1;
//      }
//   }
//   neutron_core_orbits = modelspace->hole_qspace.size() - proton_core_orbits;
//
//   intfile << "! shell model effective interaction generated by IMSRG" << endl;
//   intfile << "! Zero body term: " << op.ZeroBody << endl;
//   intfile << "! Index   n l j tz" << endl;
//   // first do proton orbits
//   for (auto& i : modelspace->valence)
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      if (oi.tz2 > 0 ) continue;
//      int nushell_indx = i/2+1 -proton_core_orbits;
//      intfile << "!  " << nushell_indx << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << endl;
//      ++nvalence_proton_orbits;
//   }
//   // then do neutron orbits
//   for (auto& i : modelspace->valence)
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      if (oi.tz2 < 0 ) continue;
//      int nushell_indx = i/2+1 + nvalence_proton_orbits -neutron_core_orbits;
//      intfile << "!  " << nushell_indx << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << endl;
//   }
//
//   intfile << "!" << endl;
//   intfile << "-999  ";
//   for (auto& i : modelspace->valence)
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      if (oi.tz2 > 0 ) continue;
//      intfile << op.OneBody(i,i) << "  ";
//   }
//   for (auto& i : modelspace->valence)
//   {
//      Orbit& oi = modelspace->GetOrbit(i);
//      if (oi.tz2 < 0 ) continue;
//      intfile << op.OneBody(i,i) << "  ";
//   }
//   intfile << "  " << Acore << "  " << Acore+2 << "  0.00000 " << endl; // No mass dependence for now...
//
//
//   int nchan = modelspace->GetNumberTwoBodyChannels();
//   for (int ch=0;ch<nchan;++ch)
//   {
//      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
//      for (auto& ibra: tbc.GetKetIndex_vv() )
//      {
//         Ket &bra = tbc.GetKet(ibra);
//         int a = bra.p;
//         int b = bra.q;
//         Orbit& oa = modelspace->GetOrbit(a);
//         Orbit& ob = modelspace->GetOrbit(b);
//         for (auto& iket: tbc.GetKetIndex_vv())
//         {
//            if (iket<ibra) continue;
//            Ket &ket = tbc.GetKet(iket);
//            int c = ket.p;
//            int d = ket.q;
//            Orbit& oc = modelspace->GetOrbit(c);
//            Orbit& od = modelspace->GetOrbit(d);
//
//            // don't pull a_ind and b_ind out of this loop, on account of the swap below.
//            int a_ind = a/2+1 + ( oa.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int b_ind = b/2+1 + ( ob.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int c_ind = c/2+1 + ( oc.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int d_ind = d/2+1 + ( od.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int T = abs(tbc.Tz);
//
//            double tbme = op.TwoBody.GetTBME_norm(ch,a,b,c,d);
//            // NuShellX quirk: even though it uses pn formalism, it requires that Vpnpn = Vnpnp,
//            //  i.e. the spatial wf for protons and neutrons are assumed to be the same.
//            // This seems to only be an issue for expectation values of operators,
//            // so the mode "op" averages Vpnpn and Vnpnp to make them equal.
//            if (mode=="op" and oa.tz2!=ob.tz2)
//            {
//               int aa = a - oa.tz2;
//               int bb = b - ob.tz2;
//               int cc = c - oc.tz2;
//               int dd = d - od.tz2;
//               tbme += op.TwoBody.GetTBME_norm(ch,aa,bb,cc,dd);
//               tbme /= 2;
//            }
//
////            if ( abs(tbme) < 1e-6) tbme = 0;
//            if ( abs(tbme) < 1e-6) continue;
//            if (T==0)
//            {
//               if ( oa.j2 != ob.j2 or oa.l != ob.l or oa.n != ob.n ) tbme *= SQRT2; // pn TBMEs are unnormalized
//               if ( oc.j2 != od.j2 or oc.l != od.l or oc.n != od.n ) tbme *= SQRT2; // pn TBMEs are unnormalized
//               T = (tbc.J+1)%2;
//            }
//            // in NuShellX, the proton orbits must come first. This can be achieved by
//            // ensuring that the bra and ket indices are in ascending order.
//            if (a_ind > b_ind)
//            {
//               swap(a_ind,b_ind);
//               tbme *= bra.Phase(tbc.J);
//            }
//            if (c_ind > d_ind)
//            {
//               swap(c_ind,d_ind);
//               tbme *= ket.Phase(tbc.J);
//            }
//            if ((a_ind > c_ind) or (c_ind==a_ind and b_ind>d_ind) )
//            {
//              swap(a_ind,c_ind);
//              swap(b_ind,d_ind);
//            }
//
//            intfile  << setw(wint) << a_ind << setw(wint) << b_ind << setw(wint) << c_ind << setw(wint) << d_ind << "    "
//                     << setw(wint) << tbc.J << " " << setw(wint) << T     << "       "
//                     << setw(wdouble) << setiosflags(ios::fixed) << setprecision(pdouble) << tbme
//                     << endl;
//         }
//      }
//   }
//   intfile.close();
//}




void ReadWrite::WriteNuShellX_intfile(Operator& op, string filename, string mode)
{
   ofstream intfile;
   intfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();

   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   int pdouble = 6; // precision for printing doubles

   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   vector<int> valence_protons(modelspace->valence.size());
   vector<int> valence_neutrons(modelspace->valence.size());
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
   intfile << "! shell model effective interaction generated by IMSRG" << endl;
   intfile << "! Zero body term: " << op.ZeroBody << endl;
   intfile << "! Index   n l j tz" << endl;

   for ( auto& it : nushell2orb )
   {
      Orbit& oi = modelspace->GetOrbit(it.second);
      intfile << "!  " << it.first << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << endl;
   }
   intfile << "!" << endl;
   intfile << "-999  ";
   for ( auto& it : nushell2orb )
   {
      intfile << op.OneBody(it.second,it.second) << "  ";
   }

   intfile << "  " << Acore << "  " << Acore+2 << "  0.00000 " << endl; // No mass dependence for now...

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
//            int a_ind = a/2+1 + ( oa.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int b_ind = b/2+1 + ( ob.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int c_ind = c/2+1 + ( oc.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int d_ind = d/2+1 + ( od.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
            int a_ind = orb2nushell[a];
            int b_ind = orb2nushell[b];
            int c_ind = orb2nushell[c];
            int d_ind = orb2nushell[d];
            int T = abs(tbc.Tz);

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

            if ( abs(tbme) < 1e-6) continue;
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
                     << setw(wdouble) << setiosflags(ios::fixed) << setprecision(pdouble) << tbme
                     << endl;
         }
      }
   }
   intfile.close();
}





/// Write the valence model space to a NuShellX *.sp file.
void ReadWrite::WriteNuShellX_sps(Operator& op, string filename)
{
   ofstream spfile;
   spfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();

   // valence protons are the intersection of valence orbits and protons orbits. Likewise for neutrons.
   vector<int> valence_protons(modelspace->valence.size());
   vector<int> valence_neutrons(modelspace->valence.size());
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
   spfile << "! modelspace for IMSRG interaction" << endl;
   spfile << "pn" << endl; // proton-neutron formalism
   spfile << Acore << " " << Zcore << endl;
   spfile << modelspace->valence.size() << endl;
   spfile << "2 " << valence_protons.size() << " " << valence_neutrons.size() << endl;

   for (auto& it : nushell2orb )
   {
     Orbit& oi = modelspace->GetOrbit(it.second);
     spfile << it.first << " " << oi.n+1 << " " << oi.l << " " << oi.j2 << endl;
   }

   spfile << endl;
   spfile.close();

}




/// Not yet fully implemented.
void ReadWrite::WriteAntoine_int(Operator& op, string filename)
{
   ofstream intfile;
   intfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
//   int ncore_orbits = modelspace->holes.size();
   int nvalence_orbits = modelspace->valence.size();
//   int nvalence_proton_orbits = 0;
   int Acore = 0;
//   int wint = 4; // width for printing integers
//   int wdouble = 12; // width for printing doubles
   for (auto& i : modelspace->holes)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 +1;
   }
   intfile << "IMSRG INTERACTION" << endl;
   // 2 indicates pn treated separately
   intfile << "2 " << nvalence_orbits << " ";
   // write NLJ of the orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      intfile << oi.n*1000 + oi.l*100 + oi.j2 << " ";
   }
   intfile << endl;
   // Write proton SPE's
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      intfile << op.OneBody(i,i) << " ";
   }
   // Write neutron SPE's
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      intfile << op.OneBody(i,i) << " ";
   }
   intfile << endl;
   
   
}



/// Write an operator to a plain-text file
void ReadWrite::WriteOperatorHuman(Operator& op, string filename)
{
   ofstream opfile;
   opfile.open(filename, ofstream::out);
   if (not opfile.good() )
   {
     cout << "Trouble opening " << filename << ". Aborting WriteOperator." << endl;
     return;
   }
   ModelSpace * modelspace = op.GetModelSpace();

   if (op.IsHermitian() )
   {
      opfile << "Hermitian" << endl;
   }
   else if (op.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << endl;
   }
   else
   {
      opfile << "Non-Hermitian" << endl;
   }
   opfile << op.GetJRank() << "  " << op.GetTRank() << "  " << op.GetParity() << endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op.ZeroBody << endl;

   opfile << "$OneBody:\t" << endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (abs(op.OneBody(i,j)) > 0)
            opfile << i << "\t" << j << "\t" << setprecision(10) << op.OneBody(i,j) << endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << endl;

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
           if ( abs(tbme) > 1e-7 )
           {
             opfile << setw(4) << tbc_bra.J << " " << tbc_bra.parity << " " << tbc_bra.Tz  << "    "
                    << setw(4) << tbc_ket.J << " " << tbc_ket.parity << " " << tbc_ket.Tz  << "    "
                  << setw(4) << bra.p << " " << bra.q  << " " << ket.p << " " << ket.q  << "   "
                  << setw(10) << setprecision(6) << tbme << endl;
           }
        }
      }
   }

   opfile.close();

}





/// Write an operator to a plain-text file
void ReadWrite::WriteOperator(Operator& op, string filename)
{
   ofstream opfile;
   opfile.open(filename, ofstream::out);
   if (not opfile.good() )
   {
     cout << "Trouble opening " << filename << ". Aborting WriteOperator." << endl;
     return;
   }
   ModelSpace * modelspace = op.GetModelSpace();

   if (op.IsHermitian() )
   {
      opfile << "Hermitian" << endl;
   }
   else if (op.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << endl;
   }
   else
   {
      opfile << "Non-Hermitian" << endl;
   }
   opfile << op.GetJRank() << "  " << op.GetTRank() << "  " << op.GetParity() << endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op.ZeroBody << endl;

   opfile << "$OneBody:\t" << endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (abs(op.OneBody(i,j)) > 0)
            opfile << i << "\t" << j << "\t" << setprecision(10) << op.OneBody(i,j) << endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << endl;

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
           if ( abs(tbme) > 1e-7 )
           {
             opfile << setw(4) << chbra << " " << setw(4) << chket << "   "
                  << setw(4) << ibra  << " " << setw(4) << iket  << "   "
                  << setw(10) << setprecision(6) << tbme << endl;
           }
        }
      }
   }

   opfile.close();

}

/// Read an operator from a plain-text file
void ReadWrite::ReadOperator(Operator &op, string filename)
{
   ifstream opfile;
   opfile.open(filename);
   if (not opfile.good() )
   {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
   }

   string tmpstring;
   int i,j,chbra,chket;
   double v;

   opfile >> tmpstring;
   if (tmpstring == "Hermitian")
   {
      op.SetHermitian();
   }
   else if (tmpstring == "Anti-Hermitian")
   {
      op.SetAntiHermitian();
   }
   else
   {
      op.SetNonHermitian();
   }
   int jrank,trank,parity;
   opfile >> jrank >> trank >> parity;

   opfile >> tmpstring >> v;
   op.ZeroBody = v;

   getline(opfile, tmpstring);
   getline(opfile, tmpstring);
   getline(opfile, tmpstring);
   while (tmpstring[0] != '$')
   {
      stringstream ss(tmpstring);
      ss >> i >> j >> v;
      op.OneBody(i,j) = v;
      if ( op.IsHermitian() )
         op.OneBody(j,i) = v;
      else if ( op.IsAntiHermitian() )
         op.OneBody(j,i) = -v;
      getline(opfile, tmpstring);
   }

  while(opfile >> chbra >> chket >> i >> j >> v)
  {
    op.TwoBody.SetTBME(chbra,chket,i,j,v);
//    cout << "Just set mat el ij = " << op.TwoBody.GetTBME_norm(chbra,chket,i,j) << "  ji = " << op.TwoBody.GetTBME_norm(chket,chbra,j,i) << "  v = " << v << endl;
  }
//  if (op.IsHermitian())
//  {
//    op.Symmetrize();
//  }
//  else if (op.IsAntiHermitian())
//  {
//    op.AntiSymmetrize();
//  }

   opfile.close();
   
}


/// Write an operator to a plain-text file
void ReadWrite::CompareOperators(Operator& op1, Operator& op2, string filename)
{
   ofstream opfile;
   opfile.open(filename, ofstream::out);
   if (not opfile.good() )
   {
     cout << "Trouble opening " << filename << ". Aborting WriteOperator." << endl;
     return;
   }
   ModelSpace * modelspace = op1.GetModelSpace();

   if (op1.IsHermitian() )
   {
      opfile << "Hermitian" << endl;
   }
   else if (op1.IsAntiHermitian())
   {
      opfile << "Anti-Hermitian" << endl;
   }
   else
   {
      opfile << "Non-Hermitian" << endl;
   }
   opfile << op1.GetJRank() << "  " << op1.GetTRank() << "  " << op1.GetParity() << endl;

   opfile << "$ZeroBody:\t" << setprecision(10) << op1.ZeroBody << "   " << op2.ZeroBody << endl;

   opfile << "$OneBody:\t" << endl;

   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      int jmin = op1.IsNonHermitian() ? 0 : i;
      for (int j=jmin;j<norb;++j)
      {
         if (abs(op1.OneBody(i,j)) > 0 or abs(op2.OneBody(i,j))>0 )
            opfile << i << "\t" << j << "\t" << setprecision(10) << op1.OneBody(i,j) << "   " << op2.OneBody(i,j) << endl;
      }
   }

   opfile <<  "$TwoBody:\t"  << endl;

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
           if ( abs(tbme1) > 1e-7 or abs(tbme2>1e-7) )
           {
             opfile << setw(4) << tbc_bra.J << " " << tbc_bra.parity << " " << tbc_bra.Tz  << "    "
                  << setw(4) << bra.p << " " << bra.q  << " " << ket.p << " " << ket.q  << "   "
                  << setw(10) << setprecision(6) << tbme1 << "  "
                  << setw(10) << setprecision(6) << tbme2 << endl;
           }
        }
      }
   }

   opfile.close();

}



void ReadWrite::ReadOneBody_Takayuki(string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2 + 1;
     orbits_remap[nlj] = i;
  }


  ifstream infile(filename);
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

void ReadWrite::ReadTwoBody_Takayuki(string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2 + 1;
     orbits_remap[nlj] = i;
  }


  ifstream infile(filename);
  int a,b,c,d,tza,tzb,tzc,tzd,J;
  double me;
  while( infile >> tza >> a >> tzb >> b >> tzc >> c >> tzd >> d >> J >> me )
  {
    int aa = orbits_remap.at(a) + (tza+1)/2;
    int bb = orbits_remap.at(b) + (tzb+1)/2;
    int cc = orbits_remap.at(c) + (tzc+1)/2;
    int dd = orbits_remap.at(d) + (tzd+1)/2;
    if ( (aa==bb or cc==dd) and (J%2)>0 ) continue;
    if (abs(me)<1e-6) continue;
    Hbare.TwoBody.SetTBME_J(J,aa,bb,cc,dd,me);
  }
}


void ReadWrite::WriteOneBody_Takayuki(string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2 + 1;
//     orbits_remap[nlj] = i;
     orbits_remap[i]   = nlj;
     orbits_remap[i+1] = nlj;
  }

  ofstream outfile(filename);
  outfile << setiosflags(ios::fixed);

  for (int a=0; a<norb; ++a)
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
    {
      Orbit& ob = modelspace->GetOrbit(b);
      double me = Hbare.OneBody(a,b);
      if (abs(me) > 1e-7)
      {
      outfile << setw(3) << oa.tz2 << " " << setw(3) << orbits_remap.at(a) << " "
              << setw(3) << ob.tz2 << " " << setw(3) << orbits_remap.at(b) << " " 
              << setw(12) << setprecision(8) <<  me << endl;
      }
    }
  }

}



void ReadWrite::WriteTwoBody_Takayuki(string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  unordered_map<int,int> orbits_remap;
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2 + 1;
//     orbits_remap[nlj] = i;
     orbits_remap[i]   = nlj;
     orbits_remap[i+1] = nlj;
  }

  ofstream outfile(filename);
  outfile << setiosflags(ios::fixed);

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
        if (abs(tbme)<1e-8) continue;
        outfile << setw(3) << bra.op->tz2 << " " << setw(3) << orbits_remap.at(bra.p) << " "
                << setw(3) << bra.oq->tz2 << " " << setw(3) << orbits_remap.at(bra.q) << " "
                << setw(3) << ket.op->tz2 << " " << setw(3) << orbits_remap.at(ket.p) << " "
                << setw(3) << ket.oq->tz2 << " " << setw(3) << orbits_remap.at(ket.q) << " "
                << setw(3) << J << setw(12) << setprecision(8) << tbme << endl;
      }
    }
  }


}


void ReadWrite::WriteTensorOneBody(string filename, Operator& Op, string opname)
{
   ofstream outfile(filename);
   ModelSpace * modelspace = Op.GetModelSpace();
   int nvalence_proton_orbits = 0;
   int proton_core_orbits = 0;
   int neutron_core_orbits = 0;
   int Acore = 0;
   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   int pdouble = 6; // precision for printing doubles
   for (auto& i : modelspace->core)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 +1;
      if (oi.tz2 < 0)
      {
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->core.size() - proton_core_orbits;

   outfile << fixed << setprecision(pdouble);
   outfile << "!  One-body matrix elements for tensor operator: " << opname << "   generated with IM-SRG" << endl;
   outfile << "!  Rank_J :  " << Op.GetJRank() << endl;
   outfile << "!  Rank_T :  " << Op.GetTRank() << endl;
   outfile << "!  Parity :  " << showpos << 1-2*Op.GetParity() << noshowpos << endl;
   outfile << "!  Zero body term:  " << Op.ZeroBody << endl;
   outfile << "!  index   n   l   2j   2tz " << endl;
   // first do proton orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      int nushell_indx = i/2+1 -proton_core_orbits;
      outfile << "!     " << nushell_indx << "    " << oi.n << "   " << oi.l << "   " << oi.j2  << "   " << setw(wint) << oi.tz2 <<  endl;
      ++nvalence_proton_orbits;
   }
   // then do neutron orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      int nushell_indx = i/2+1 + nvalence_proton_orbits -neutron_core_orbits;
      outfile << "!     " << nushell_indx << "    " << oi.n << "   " << oi.l << "   " << oi.j2 << "   " << setw(wint) << oi.tz2 << endl;
   }

   outfile << "!" << endl << "!  a   b   < a || Op || b > " << endl;

   for ( auto a : modelspace->valence )
   {
     Orbit& oa = modelspace->GetOrbit(a);
     int a_ind = a/2+1 + ( oa.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//     for ( auto b : Op.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      for ( auto b : modelspace->valence )
     {
        Orbit& ob = modelspace->GetOrbit(b);
        double me = Op.OneBody(a,b);
        if ( abs(me) < 1e-7 ) continue;
        int b_ind = b/2+1 + ( ob.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
        outfile << setw(wint) << a_ind << " " << setw(wint) << b_ind << " " << setw(wdouble) << setprecision(pdouble) <<  me << endl;
     }
   }
}

void ReadWrite::WriteTensorTwoBody(string filename, Operator& Op, string opname)
{
   ofstream outfile(filename);
   ModelSpace * modelspace = Op.GetModelSpace();
   int nvalence_proton_orbits = 0;
   int proton_core_orbits = 0;
   int neutron_core_orbits = 0;
   int Acore = 0;
   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   int pdouble = 6; // precision for printing doubles
   for (auto& i : modelspace->core)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 +1;
      if (oi.tz2 < 0)
      {
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->core.size() - proton_core_orbits;

   outfile << fixed << setprecision(pdouble);
   outfile << "!  Two-body matrix elements for tensor operator: " << opname << "   generated with IM-SRG" << endl;
   outfile << "!  Rank_J :  " << Op.GetJRank() << endl;
   outfile << "!  Rank_T :  " << Op.GetTRank() << endl;
   outfile << "!  Parity :  " << showpos << 1-2*Op.GetParity() << noshowpos << endl;
   outfile << "!  Zero body term:  " << Op.ZeroBody << endl;
   outfile << "!  index   n   l   2j   2tz " << endl;
   // first do proton orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      int nushell_indx = i/2+1 -proton_core_orbits;
      outfile << "!     " << nushell_indx << "    " << oi.n << "   " << oi.l << "   " << oi.j2 << "   " << setw(wint) << oi.tz2 << endl;
      ++nvalence_proton_orbits;
   }
   // then do neutron orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      int nushell_indx = i/2+1 + nvalence_proton_orbits -neutron_core_orbits;
      outfile << "!     " << nushell_indx << "    " << oi.n << "   " << oi.l << "   " << oi.j2 << "   " << setw(wint) << oi.tz2 << endl;
   }

   outfile << "!  a    b    c    d     Jab  Jcd   <ab Jab || Op || cd Jcd>" << endl;
   for ( auto& itmat : Op.TwoBody.MatEl )
   {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(itmat.first[0]);
     TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(itmat.first[1]);
     auto& matrix = itmat.second;
//     int nbras = tbc_bra.GetNumberKets();
//     int nkets = tbc_ket.GetNumberKets();
//     for ( int ibra=0; ibra<nbras; ++ibra)
     for (auto& ibra: tbc_bra.GetKetIndex_vv() )
     {
       Ket& bra = tbc_bra.GetKet(ibra);
       Orbit& oa = modelspace->GetOrbit(bra.p);
       Orbit& ob = modelspace->GetOrbit(bra.q);
       int a_ind = oa.index/2+1 + ( oa.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
       int b_ind = ob.index/2+1 + ( ob.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//       for ( int iket=0; iket<nkets; ++iket)
       for (auto& iket: tbc_ket.GetKetIndex_vv() )
       {
         double me = matrix(ibra,iket);
         if (abs(me) < 1e-7) continue;
         Ket& ket = tbc_ket.GetKet(iket);
         Orbit& oc = modelspace->GetOrbit(ket.p);
         Orbit& od = modelspace->GetOrbit(ket.q);
         int c_ind = oc.index/2+1 + ( oc.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int d_ind = od.index/2+1 + ( od.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         outfile << setw(wint) << a_ind << " " << setw(wint) << b_ind << " " << setw(wint) << c_ind << " " << setw(wint) << d_ind << "   "
                 << setw(wint) << tbc_bra.J << " " << setw(wint) << tbc_ket.J << "   " << setw(wdouble) << setprecision(pdouble) << me << endl;
         
       }
     }

   }

}


void ReadWrite::ReadTwoBodyEngel(string filename, Operator& Op)
{
  if ( filename.substr( filename.find_last_of(".")) == ".gz")
  {
    ifstream infile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    ReadTwoBodyEngel_from_stream(zipstream, Op);
  }
  else
  {
    ifstream infile(filename);
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



void ReadWrite::SetLECs(double c1, double c3, double c4, double cD, double cE)
{
   LECs = {c1,c3,c4,cD,cE};
}

void ReadWrite::SetLECs_preset(string key)
{
       if (key == "EM1.8_2.0")  LECs = {-0.81, -3.20, 5.40, 1.264, -0.120};
  else if (key == "EM2.0_2.0")  LECs = {-0.81, -3.20, 5.40, 1.271, -0.131};
  else if (key == "EM2.2_2.0")  LECs = {-0.81, -3.20, 5.40, 1.214, -0.137};
  else if (key == "EM2.8_2.0")  LECs = {-0.81, -3.20, 5.40, 1.278, -0.078};
  else if (key == "PWA2.0_2.0") LECs = {-0.76, -4.78, 3.96,-3.007, -0.686};
  else if (key == "N2LOSAT")    LECs = {-1.12152120, -3.92500586, 3.76568716, 0.861680589, -0.03957471};
}




