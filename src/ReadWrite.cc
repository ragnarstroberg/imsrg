#include "ReadWrite.hh"
#include "ModelSpace.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "string.h"

#define LINESIZE 400
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif

using namespace std;

ReadWrite::ReadWrite()
{
   doCoM_corr = false;
   goodstate = true;
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




ModelSpace ReadWrite::ReadModelSpace( string filename)
{

   ModelSpace modelspace;
   ifstream infile;
   stringstream ss;
   char line[LINESIZE];
   char cbuf[10][20];
   int ibuf[2];
   int n,l,j2,tz2;
   double fbuf[2];
   double spe,hw;
   int A;
   
   infile.open(filename);
   if ( not infile.good() )
   {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
      cerr << "Trouble reading file " << filename << endl;
      goodstate = false;
      return modelspace;
   }
   
   infile.getline(line,LINESIZE);

   while (!strstr(line,"Mass number") && !infile.eof())
   {
      infile.getline(line,LINESIZE);
   }
   ss << line;
   for (int i=0;i<10;i++) ss >> cbuf[i];
   ss >> A;
   ss.str(string()); // clear up the stringstream
   ss.clear();
   
   while (!strstr(line,"Oscillator length") && !infile.eof())
   {
      infile.getline(line,LINESIZE);
   }
   ss << line;
   ss >> cbuf[0] >> cbuf[1] >> cbuf[2] >> cbuf[3] >> fbuf[0] >> hw;
   modelspace.SetHbarOmega(hw);
   modelspace.SetTargetMass(A);
//   cout << "Set hbar_omega to " << hw << endl;

   while (!strstr(line,"Legend:") && !infile.eof())
   {
      infile.getline(line,LINESIZE);
   }

   while (  infile >> cbuf[0] >> ibuf[0] >> n >> l >> j2 >> tz2
             >> ibuf[1] >> spe >> fbuf[0]  >> cbuf[1] >> cbuf[2])
   {
      int ph = 0; // 0=particle, 1=hole 
      int io = 0; // 0=inside, 1=outside projected model space
      if (strstr(cbuf[1],"hole")) ph++;
      if (strstr(cbuf[2],"outside")) io++;
      modelspace.AddOrbit( Orbit(n,l,j2,tz2,ph,io,spe) );

   }
   
   infile.close();
   modelspace.SetupKets();

   return modelspace;
}


void ReadWrite::ReadBareTBME( string filename, Operator& Hbare)
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


void ReadWrite::ReadBareTBME_Jason( string filename, Operator& Hbare)
{

  ifstream infile;
  char line[LINESIZE];
  int Tz2,Par,J2,a,b,c,d;
  int na,nb,nc,nd;
  int la,lb,lc,ld;
  int ja,jb,jc,jd;
  double fbuf[3];
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


// Read TBME's in Petr Navratil's format
// Setting Emax=-1 just uses the single-particle emax determined by the model space
//void ReadWrite::ReadBareTBME_Navratil( string filename, Operator3& Hbare)
//template<class OPERATOR>
//void ReadWrite::ReadBareTBME_Navratil( string filename, OPERATOR Hbare)
void ReadWrite::ReadBareTBME_Navratil( string filename, Operator& Hbare)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int emax = 0;
  int norb = modelspace->GetNumberOrbits();
  int nljmax = norb/2;
  int herm = Hbare.IsHermitian() ? 1 : -1 ;
  vector<int> orbits_remap(nljmax,-1);
  for (int i=0;i<norb;++i)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     if (oi.tz2 > 0 ) continue;
     int N = 2*oi.n + oi.l;
     int nlj = N*(N+1)/2 + max(oi.l-1,0) + (oi.j2 - abs(2*oi.l-1))/2;
     orbits_remap[nlj] = i;
     emax = max(emax,N);
  }
  emax *= 2;

  ifstream infile;
  char line[LINESIZE];
  infile.open(filename);
  double tbme_pp,tbme_nn,tbme_10,tbme_00;
  if ( !infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }

  // first line contains some information
  int total_number_elements;
  int N1max, N12max;
  double hw, srg_lambda;

  infile >> total_number_elements >> N1max >> N12max >> hw >> srg_lambda;

  // from Petr: The Trel and the H_HO_rel must be scaled by (2/A) hbar Omega
  int a,b,c,d,J,T;
  int nlj1,nlj2,nlj3,nlj4;
  double trel, h_ho_rel, vcoul, vpn, vpp, vnn;
  while( infile >> nlj1 >> nlj2 >> nlj3 >> nlj4 >> J >> T >> trel >> h_ho_rel >> vcoul >> vpn >> vpp >> vnn )
  {
    --nlj1;--nlj2;--nlj3;--nlj4;  // Fortran -> C indexing
    a = orbits_remap[nlj1];
    b = orbits_remap[nlj2];
    c = orbits_remap[nlj3];
    d = orbits_remap[nlj4];

    if (doCoM_corr)
    {
        vpn += trel * 2 * hw / Hbare.GetModelSpace()->GetTargetMass();
        vpp += trel * 2 * hw / Hbare.GetModelSpace()->GetTargetMass();
        vnn += trel * 2 * hw / Hbare.GetModelSpace()->GetTargetMass();
    }
    
    //if (a>=norb or b>=norb or c>=norb or d>=norb) continue;
    if (2*nlj1>=norb or 2*nlj2>=norb or 2*nlj3>=norb or 2*nlj4>=norb) continue;
    Orbit oa = modelspace->GetOrbit(a);
    Orbit ob = modelspace->GetOrbit(b);
    int parity = (oa.l + ob.l) % 2;


    // do pp and nn
    if (T==1)
    {
//      Hbare.SetTBME(J,parity,-1,a,b,c,d,vpp);
//      Hbare.SetTBME(J,parity,1,a+1,b+1,c+1,d+1,vnn);
      Hbare.TwoBody.SetTBME(J,parity,-1,a,b,c,d,vpp);
      Hbare.TwoBody.SetTBME(J,parity,1,a+1,b+1,c+1,d+1,vnn);
    }

    if (abs(vpn)<1e-6) continue;


    // Normalization
    if (a!=b)
    {
       vpn /= SQRT2;
    }
    if (c!=d)
    {
       vpn /= SQRT2;
    }
    if ( a==c and b==d )
    {
       vpn /= 2.0;
    }


   double vpnpn = vpn;
   double vpnnp = vpn;
   if (T==0) vpnnp *= -1;

  // These normalizations seem to work. Not sure why. There must be a better way.
  if ( a==d and b==c )
  {
     vpnnp /=2.0;
  }

    // now do pnpn, npnp, pnnp, nppn
  Hbare.TwoBody.AddToTBME(J,parity,0,a,b+1,c,d+1,vpnpn);
  if (a!=b and c!=d)
     Hbare.TwoBody.AddToTBME(J,parity,0,a+1,b,c+1,d,vpnpn);
  if ( c!=d)
    Hbare.TwoBody.AddToTBME(J,parity,0,a,b+1,c+1,d,vpnnp);
  if ( a!=b)
    Hbare.TwoBody.AddToTBME(J,parity,0,a+1,b,c,d+1,vpnnp);

  }
  
  return;
}





// Read TBME's in Darmstadt format
// Admittedly this is a mess. Emax is the 2body max energy
// Setting Emax=-1 just uses the single-particle emax determined by the model space
//void ReadWrite::ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax /*default=-1*/)
//void ReadWrite::ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int emax /*default=-1*/, int lmax)
void ReadWrite::ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int emax, int Emax, int lmax)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norb = modelspace->GetNumberOrbits();
  int herm = Hbare.IsHermitian() ? 1 : -1 ;
  vector<int> orbits_remap;

  if (emax < 0)  emax = modelspace->Nmax;
  if (lmax < 0)  lmax = emax;

  for (int e=0; e<emax; ++e)
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


  ifstream infile;
  char line[LINESIZE];
  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }

  double tbme_pp,tbme_nn,tbme_10,tbme_00;
  // skip the first line
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

             // convert isospin to pn formalism
             double tbme_pnpn = (tbme_10 + tbme_00)/2.0;
             double tbme_pnnp = (tbme_10 - tbme_00)/2.0;

             // Normalization
             if (a==b)
             {
                tbme_pp /= SQRT2;
                tbme_nn /= SQRT2;
             }
             if (c==d)
             {
                tbme_pp /= SQRT2;
                tbme_nn /= SQRT2;
             }

             // do pp and nn
             if (tbme_pp !=0)
                Hbare.TwoBody.SetTBME(J,parity,-1,a,b,c,d,tbme_pp);
             if (tbme_nn !=0)
                Hbare.TwoBody.SetTBME(J,parity,1,a+1,b+1,c+1,d+1,tbme_nn);

             // now do pnpn, npnp, pnnp, nppn
             if (tbme_pnpn !=0)
             {
               Hbare.TwoBody.SetTBME(J,parity,0,a,b+1,c,d+1,tbme_pnpn);
               Hbare.TwoBody.SetTBME(J,parity,0,a+1,b,c+1,d,tbme_pnpn);
             }
             if (tbme_pnnp !=0)
             {
               Hbare.TwoBody.SetTBME(J,parity,0,a,b+1,c+1,d,tbme_pnnp);
               Hbare.TwoBody.SetTBME(J,parity,0,a+1,b,c,d+1,tbme_pnnp);
             }
          }
        }
      }
    }
  }

}



void ReadWrite::Read_Darmstadt_3body( string filename, Operator& Hbare, int E1max, int E2max, int E3max)
{
  ModelSpace * modelspace = Hbare.GetModelSpace();
  vector<string> empty = {};
//  ModelSpace big_modelspace = ModelSpace(E1max,{},{});
  ModelSpace big_modelspace = ModelSpace(E1max,"skeleton");
  cout << "Done building big model space" << endl;
  int e1max = modelspace->GetNmax();
  int e2max = modelspace->GetN2max(); // not used yet
  int e3max = modelspace->GetN3max();
  cout << "Reading 3body file. emax limits for file: " << E1max << " " << E2max << " " << E3max << "  for modelspace: " << e1max << " " << e2max << " " << e3max << endl;

  vector<int> orbits_remap;
  int lmax = E1max; // haven't yet implemented the lmax truncation for 3body. Should be easy.

  for (int e=0; e<E1max; ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=min(e,lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
//         orbits_remap.push_back( modelspace->GetOrbitIndex(n,l,twoj,-1) );
         orbits_remap.push_back( big_modelspace.GetOrbitIndex(n,l,twoj,-1) );
      }
    }
  }
  int nljmax = orbits_remap.size()-1;
/*
  for (int i=0;i<nljmax;++i)
  {
    int a = orbits_remap[i];
    Orbit& obig = big_modelspace.GetOrbit(a);
    Orbit& osm = modelspace->GetOrbit(a);
    cout << i << "  -->  " << a << "  "
         << obig.n << " " << obig.l << " " << obig.j2 << " " << obig.tz2 << "    "
         << osm.n << " " << osm.l << " " << osm.j2 << " " << osm.tz2 << endl;
  }
*/



  ifstream infile;
  char line[LINESIZE];
  infile.open(filename);
  if ( not infile.good() )
  {
     cerr << "************************************" << endl
          << "**    Trouble reading file  !!!   **" << filename << endl
          << "************************************" << endl;
     goodstate = false;
     return;
  }
  infile.getline(line,LINESIZE);


  // begin giant nested loops
  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
    int a =  orbits_remap[nlj1];
//    Orbit & oa = modelspace->GetOrbit(a);
    Orbit & oa = big_modelspace.GetOrbit(a);
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
        {
           twoJCMindownbra = abs(oa.j2 - ob.j2)-oc.j2;
        }
        else if (oc.j2 < (oa.j2+ob.j2) )
        {
           twoJCMindownbra = 1;
        }
        else
        {
           twoJCMindownbra = oc.j2 - oa.j2 - ob.j2;
        }
        int twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2;


        // now loop over possible ket orbits
        for(int nnlj1=0; nnlj1<=nlj1; ++nnlj1)
        {
          int d =  orbits_remap[nnlj1];
          Orbit & od = modelspace->GetOrbit(d);
          int ed = 2*od.n + od.l;
          if (ed > E1max) break;

          for(int nnlj2=0; nnlj2 <= ((nlj1 == nnlj1) ? nlj2 : nnlj1); ++nnlj2)
          {
            int e =  orbits_remap[nnlj2];
            Orbit & oe = modelspace->GetOrbit(e);
            int ee = 2*oe.n + oe.l;
            if ( (ed+ee) > E2max) break;

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
              {
                 twoJCMindownket = abs(od.j2 - oe.j2) - of.j2;
              }
              else if ( of.j2 < (od.j2+oe.j2) )
              {
                 twoJCMindownket = 1;
              }
              else
              {
                 twoJCMindownket = of.j2 - od.j2 - oe.j2;
              }
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
                    double V;
                    infile >> V;
                    bool autozero = false;


                    if ( a==b and (tab+Jab)%2==0 ) autozero = true;
                    if ( d==e and (ttab+JJab)%2==0 ) autozero = true;
                    if ( a==b and a==c and twoT==3 and oa.j2<3 ) autozero = true;
                    if ( d==e and d==f and twoT==3 and od.j2<3 ) autozero = true;


                    if (not autozero and abs(V)>1e-5)
                    {
//                       double V0 = Hbare.ThreeBody.GetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
                       if(ea<=e1max and eb<=e1max and ec<=e1max and ed<=e1max and ee<=e1max and ef<=e1max
                          and ea+eb+ec<=e3max and ed+ee+ef<=e3max)
                       {
//                        cout << "abc " << a << " " << b << " " << c << "  def " << d << " " << e << " " << f << endl;
                        Hbare.ThreeBody.SetME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f, V);
/*
                        double Vout = Hbare.GetThreeBodyME(Jab,JJab,twoJC,tab,ttab,twoT,a,b,c,d,e,f);
                        if (abs(V-Vout)>1e-4)
                        {
                          cout << "STARTED WITH  " << V0 << "  PUT IN " << V << "   GOT OUT " << Vout << "     Jab,JJab,twoJC,tab,ttab,twoT = "
                               << Jab << " " << JJab << " " << twoJC << " " << tab << " " << ttab << " " << twoT << "   "
                               << " abcdef = "
                               <<  a << " " << b << " " << c << " " << d << " "<< e << " " << f << "  e3 = " << ea+eb+ec << " " << ed+ee+ef << endl;
                        }
*/
                       }
                    }

                    if (autozero)
                    {
//                       cout << " ( should be zero ) ";
                       if (abs(V) > 1e-6 and ea<=e1max and eb<=e1max and ec<=e1max)
                       {
                          cout << " <-------- AAAAHHHH!!!!!!!! ";
                       }
                    }
 //                   cout << endl;
       
                   }
                  }//ttab
                 }//tab
//                 cout << " ------------------------" << endl;
                if (not infile.good() ) break;
                }//twoJ
               }//JJab
       
              }



            }
          }
        }
      }
    }
  }

}






void ReadWrite::WriteOneBody(Operator& op, string filename)
{
   ofstream obfile;
   obfile.open(filename, ofstream::out);
   int norbits = op.GetModelSpace()->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      for (int j=0;j<norbits;j++)
      {
         if (abs(op.GetOneBody(i,j)) > 1e-6)
         {
            obfile << i << "    " << j << "       " << op.GetOneBody(i,j) << endl;

         }
      }
   }
   obfile.close();
}

void ReadWrite::WriteValenceOneBody(Operator& op, string filename)
{
   ofstream obfile;
   obfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int norbits = modelspace->GetNumberOrbits();
   obfile << " Zero body part: " << op.ZeroBody << endl;
   for (auto& i : modelspace->valence )
   {
      for (auto& j : modelspace->valence )
      {
         double obme = op.OneBody(i,j);
         if (abs(obme) > 1e-6)
         {
            obfile << i << "    " << j << "       " << obme << endl;
         }
      }
   }
   obfile.close();
}


void ReadWrite::WriteNuShellX_int(Operator& op, string filename)
{
   ofstream intfile;
   intfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int ncore_orbits = modelspace->holes.size();
   int nvalence_orbits = modelspace->valence.size();
   int nvalence_proton_orbits = 0;
   int proton_core_orbits = 0;
   int neutron_core_orbits = 0;
   int Acore = 0;
   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   int pdouble = 6; // precision for printing doubles
   for (auto& i : modelspace->hole_qspace)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 +1;
      if (oi.tz2 < 0)
      {
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->hole_qspace.size() - proton_core_orbits;

   intfile << "! shell model effective interaction generated by IMSRG" << endl;
   intfile << "! Zero body term: " << op.ZeroBody << endl;
   intfile << "! Index   nljtz" << endl;
   // first do proton orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      int nushell_indx = i/2+1 -proton_core_orbits;
      intfile << "!  " << nushell_indx << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << endl;
      ++nvalence_proton_orbits;
   }
   // then do neutron orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      int nushell_indx = i/2+1 + nvalence_proton_orbits -neutron_core_orbits;
      intfile << "!  " << nushell_indx << "   " << oi.n << " " << oi.l << " " << oi.j2 << "/2" << " " << oi.tz2 << "/2" << endl;
   }

   intfile << "!" << endl;
   intfile << "-999  ";
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      intfile << op.OneBody(i,i) << "  ";
   }
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      intfile << op.OneBody(i,i) << "  ";
   }
   intfile << "  " << Acore << " " << Acore+2 << "  0.00000 " << endl; // No mass dependence for now...

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
//      for (auto& ibra: tbc.KetIndex_vv)
      for (auto& ibra: tbc.GetKetIndex_vv() )
      {
         Ket &bra = tbc.GetKet(ibra);
         int a = bra.p;
         int b = bra.q;
         Orbit& oa = modelspace->GetOrbit(a);
         Orbit& ob = modelspace->GetOrbit(b);
         int a_ind = a/2+1 + ( oa.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int b_ind = b/2+1 + ( ob.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//         for (auto& iket: tbc.KetIndex_vv)
         for (auto& iket: tbc.GetKetIndex_vv())
         {
            if (iket<ibra) continue;
            Ket &ket = tbc.GetKet(iket);
            int c = ket.p;
            int d = ket.q;
            Orbit& oc = modelspace->GetOrbit(c);
            Orbit& od = modelspace->GetOrbit(d);
         int c_ind = c/2+1 + ( oc.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int d_ind = d/2+1 + ( od.tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
            int T = abs(tbc.Tz);
            //double tbme = op.TwoBody[ch](ibra,iket);
            double tbme = op.TwoBody.GetTBME(ch,bra,ket);
            if ( abs(tbme) < 1e-6) tbme = 0;
            if (T==0)
            {
               if (oa.j2 == ob.j2 and oa.l == ob.l and oa.n == ob.n) T = (tbc.J+1)%2;
               else tbme *= SQRT2; // pn TBMEs are unnormalized
               if (oc.j2 == od.j2 and oc.l == od.l and oc.n == od.n) T = (tbc.J+1)%2;
               else tbme *= SQRT2; // pn TBMEs are unnormalized
            }
            // in NuShellX, the proton orbits must come first.
            if (a_ind > b_ind)
            {
               intfile << setw(wint) << b_ind << setw(wint) << a_ind;
               tbme *= bra.Phase(tbc.J);
            }
            else
            {
               intfile << setw(wint) << a_ind << setw(wint) << b_ind;
            }
            if (c_ind > d_ind)
            {
               intfile << setw(wint) << d_ind << setw(wint) << c_ind;
               tbme *= ket.Phase(tbc.J);
            }
            else
            {
               intfile << setw(wint) << c_ind << setw(wint) << d_ind;
            }
            intfile
              << setw(wint) << tbc.J
              << "   "
              << setw(wint) << T
              << "       "
              << setw(wdouble) << setiosflags(ios::fixed) << setprecision(pdouble) << tbme
              << endl;
         }
      }
   }
   intfile.close();
}

void ReadWrite::WriteNuShellX_sps(Operator& op, string filename)
{
   ofstream spfile;
   spfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int ncore_orbits = modelspace->holes.size();
   int proton_core_orbits = 0;
   int neutron_core_orbits = 0;
   int nvalence_orbits = modelspace->valence.size();
   int nvalence_proton_orbits = 0;
   int Acore = 0;
   int Zcore = 0;
   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
   //for (int& i : modelspace->holes)
   for (auto& i : modelspace->hole_qspace)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      Acore += oi.j2 +1;
      if (oi.tz2 < 0)
      {
         Zcore += oi.j2+1;
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->hole_qspace.size() - proton_core_orbits;
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0)
      {
         ++nvalence_proton_orbits ;
      }
   }
   
   spfile << "! modelspace for IMSRG interaction" << endl;
   spfile << "pn" << endl; // proton-neutron formalism
   spfile << Acore << " " << Zcore << endl;
   spfile << nvalence_orbits << endl;
   spfile << "2 " << nvalence_proton_orbits << " " << nvalence_orbits-nvalence_proton_orbits << endl;
   // first do proton orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 > 0 ) continue;
      int nushell_indx = i/2-proton_core_orbits + 1;
      spfile << nushell_indx << " " << oi.n+1 << " " << oi.l << " " << oi.j2  << endl;
   }
   // then do neutron orbits
   for (auto& i : modelspace->valence)
   {
      Orbit& oi = modelspace->GetOrbit(i);
      if (oi.tz2 < 0 ) continue;
      int nushell_indx = i/2-neutron_core_orbits + nvalence_proton_orbits + 1;
      spfile << nushell_indx << " " << oi.n+1 << " " << oi.l << " " << oi.j2 << endl;
   }
   spfile << endl;
   spfile.close();

}




void ReadWrite::WriteAntoine_int(Operator& op, string filename)
{
   ofstream intfile;
   intfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int ncore_orbits = modelspace->holes.size();
   int nvalence_orbits = modelspace->valence.size();
   int nvalence_proton_orbits = 0;
   int Acore = 0;
   int wint = 4; // width for printing integers
   int wdouble = 12; // width for printing doubles
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


void ReadWrite::WriteTwoBody(Operator& op, string filename)
{
   ofstream tbfile;
   tbfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      for (int i=0;i<npq;++i)
      {
         Ket &bra = tbc.GetKet(i);
         Orbit &oa = modelspace->GetOrbit(bra.p);
         Orbit &ob = modelspace->GetOrbit(bra.q);
         for (int j=i;j<npq;++j)
         {
            Ket &ket = tbc.GetKet(j);
            double tbme = op.TwoBody.GetTBME(ch,bra,ket) / sqrt( (1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
            if ( abs(tbme)<1e-6 ) tbme = 0;
            Orbit &oc = modelspace->GetOrbit(ket.p);
            Orbit &od = modelspace->GetOrbit(ket.q);
            int wint = 4;
            int wdouble = 12;

            tbfile 
                   << setw(wint) << bra.p
                   << setw(wint) << bra.q
                   << setw(wint) << ket.p
                   << setw(wint) << ket.q
                   << setw(wint+3) << tbc.J << setw(wdouble) << std::fixed << tbme
                   << endl;
         }
      }
   }
   tbfile.close();
}


void ReadWrite::WriteValenceTwoBody(Operator& op, string filename)
{
   ofstream tbfile;
   tbfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;
//      for (auto& i : tbc.KetIndex_vv )
      for (auto& i : tbc.GetKetIndex_vv() )
      {
         Ket &bra = tbc.GetKet(i);
         Orbit &oa = modelspace->GetOrbit(bra.p);
         Orbit &ob = modelspace->GetOrbit(bra.q);
//         for (auto& j : tbc.KetIndex_vv )
         for (auto& j : tbc.GetKetIndex_vv() )
         {
            if (j<i) continue;
            Ket &ket = tbc.GetKet(j);
            double tbme = op.TwoBody.GetTBME(ch,bra,ket) / sqrt( (1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
            if ( abs(tbme)<1e-4 ) continue;
            Orbit &oc = modelspace->GetOrbit(ket.p);
            Orbit &od = modelspace->GetOrbit(ket.q);
            int wint = 4;
            int wdouble = 12;

            tbfile 
                   << setw(wint) << bra.p
                   << setw(wint) << bra.q
                   << setw(wint) << ket.p
                   << setw(wint) << ket.q
                   << setw(wint+3) << tbc.J << setw(wdouble) << std::fixed << tbme// << endl;
                   << "    < " << bra.p << " " << bra.q << " | V | " << ket.p << " " << ket.q << " >" << endl;
         }
      }
   }
   tbfile.close();
}


void ReadWrite::WriteOperator(Operator& op, string filename)
{
   ofstream opfile;
   opfile.open(filename, ofstream::out);
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

   opfile << "$ZeroBody:\t" << setprecision(10) << op.ZeroBody << endl;
   //opfile << "$OneBody:\t" << op.ZeroBody << endl;
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
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0; ch<nchan; ++ch)
   {
//      opfile << "ch: " << ch << endl;
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nkets; ++ibra)
      {
         Ket& bra = tbc.GetKet(ibra);
         int iket_min = op.IsNonHermitian() ? 0 : ibra;
         for (int iket=iket_min; iket<nkets; ++iket)
         {
            Ket& ket = tbc.GetKet(iket);
            double tbme = op.TwoBody.GetTBME(ch,ibra,iket);
            if (abs(tbme) > 1e-7)
              opfile << ch << "\t" << ibra << "\t" << iket << "\t" << setprecision(10) << tbme << endl;
         }
      }
   }
   opfile.close();

}


void ReadWrite::ReadOperator(Operator &op, string filename)
{
   ifstream opfile;
   opfile.open(filename);
   ModelSpace * modelspace = op.GetModelSpace();
   // Should put in some check for if the file exists

   string tmpstring;
   int i,j,ch;
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
   while (opfile >> ch >> i >> j >> v)
   {
      op.TwoBody.SetTBME(ch,i,j,v);
  }

   opfile.close();
   
}


