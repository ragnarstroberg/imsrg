#include "ReadWrite.hh"
#include "ModelSpace.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "string.h"

#define LINESIZE 400

using namespace std;

ReadWrite::ReadWrite()
{
   doCoM_corr = true;
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
      cout << "Trouble reading Settings file" << endl;
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
   float fbuf[2];
   float spe,hw;
   int A;
   
   infile.open(filename);
   if ( not infile.good() )
   {
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
   cout << "Set hbar_omega to " << hw << endl;

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
  float fbuf[3];
  float tbme;
  int norbits = Hbare.GetModelSpace()->GetNumberOrbits();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "Trouble reading file " << filename << endl;
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

     double com_corr = fbuf[2] * Hbare.GetModelSpace()->GetHbarOmega() / Hbare.GetModelSpace()->GetTargetMass();  // Some sort of COM correction. Check this

// NORMALIZATION: Read in normalized, antisymmetrized TBME's

     if (doCoM_corr)
     {
        //Hbare.SetTBME(J2/2,Par,Tz,a,b,c,d, tbme-com_corr );
        //Hbare.SetTBME(J2/2,Par,Tz,a,b,c,d, com_corr );
        Hbare.SetTBME(J2/2,Par,Tz,a,b,c,d, fbuf[2] );
     }
     else
     {
        Hbare.SetTBME(J2/2,Par,Tz,a,b,c,d, tbme ); // Don't do COM correction, for comparison with Darmstadt interaction.
     }

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
  float fbuf[3];
  float tbme;
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norbits = modelspace->GetNumberOrbits();

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "Trouble reading file " << filename << endl;
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
     Hbare.SetTBME(J2/2,Par,Tz2/2,a,b,c,d, tbme ); // Don't do COM correction, for comparison with Darmstadt interaction.

  }

  return;
}


// Read TBME's in Darmstadt format
void ReadWrite::ReadBareTBME_Darmstadt( string filename, Operator& Hbare, int Emax /*default=-1*/)
{

  ModelSpace * modelspace = Hbare.GetModelSpace();
  int emax = 0;
  int norb = modelspace->GetNumberOrbits();
  int nljmax = norb/2;
  int herm = Hbare.IsHermitian() ? 1 : -1 ;
  vector<int> orbits_darmstadt(nljmax,-1);
  for (int i=0;i<norb;++i)
  {
     Orbit* oi = modelspace->GetOrbit(i);
     if (oi->tz2 > 0 ) continue;
     int N = 2*oi->n + oi->l;
     int nlj = N*(N+1)/2 + max(oi->l-1,0) + (oi->j2 - abs(2*oi->l-1))/2;
     orbits_darmstadt[nlj] = i;
     if (N > emax) ++emax;
  }
  emax *= 2;
  if (Emax >=0)
     emax = Emax;

  ifstream infile;
  char line[LINESIZE];
  infile.open(filename);
  double tbme;
  double tbme_pp,tbme_nn,tbme_10,tbme_00;
  if ( !infile.good() )
  {
     cerr << "Trouble reading file " << filename << endl;
     goodstate = false;
     return;
  }

  // skip the first line
  infile.getline(line,LINESIZE);

  for(int nlj1=0; nlj1<nljmax; ++nlj1)
  {
    int a =  orbits_darmstadt[nlj1];
    Orbit * o1 = modelspace->GetOrbit(a);
    int e1 = 2*o1->n + o1->l;

    for(int nlj2=0; nlj2<=nlj1; ++nlj2)
    {
      int b =  orbits_darmstadt[nlj2];
      Orbit * o2 = modelspace->GetOrbit(b);
      int e2 = 2*o2->n + o2->l;
      if (e1+e2 > emax) break;
      int parity = (o1->l + o2->l) % 2;

      for(int nlj3=0; nlj3<=nlj1; ++nlj3)
      {
        int c =  orbits_darmstadt[nlj3];
        Orbit * o3 = modelspace->GetOrbit(c);
        int e3 = 2*o3->n + o3->l;

        for(int nlj4=0; nlj4<=(nlj3==nlj1 ? nlj2 : nlj3); ++nlj4)
        {
          int d =  orbits_darmstadt[nlj4];
          Orbit * o4 = modelspace->GetOrbit(d);
          int e4 = 2*o4->n + o4->l;
          if (e3+e4 > emax) break;
          if ( (o1->l + o2->l + o3->l + o4->l)%2 != 0) continue;
          int Jmin = max( abs(o1->j2 - o2->j2), abs(o3->j2 - o4->j2) )/2;
          int Jmax = min (o1->j2 + o2->j2, o3->j2+o4->j2)/2;
          if (Jmin > Jmax) continue;
          cout << endl << nlj1 << " " << nlj2 << " " << nlj3 << " " << nlj4 << "   -->  "
               << a << " "  
               << b << " "  
               << c << " "  
               << d << " "  
               << endl;
          for (int J=Jmin; J<=Jmax; ++J)
          {
             infile >> tbme_00 >> tbme_nn >> tbme_10 >> tbme_pp;

             // convert isospin to pn formalism
             double tbme_pnpn = (tbme_10 + tbme_00)/2;
             double tbme_pnnp = (tbme_10 - tbme_00)/2;

             cout << " J = " << J << " tbme: " << tbme_00 << " " << tbme_nn  << " " << tbme_10 << " " << tbme_pp << endl;
             if (a>=norb or b>=norb or c>=norb or d>=norb) continue;

             int ch_nn = modelspace->GetTwoBodyChannelIndex(J,parity,1);
             int ch_pn = modelspace->GetTwoBodyChannelIndex(J,parity,0);
             int ch_pp = modelspace->GetTwoBodyChannelIndex(J,parity,-1);
             TwoBodyChannel& tbc_nn = modelspace->GetTwoBodyChannel(ch_nn);
             TwoBodyChannel& tbc_pn = modelspace->GetTwoBodyChannel(ch_pn);
             TwoBodyChannel& tbc_pp = modelspace->GetTwoBodyChannel(ch_pp);

             // do pp first

             if (tbme_pp !=0)
             {
                   // normalize
                   if (a==b) tbme_pp /= sqrt(2);
                   if (c==d) tbme_pp /= sqrt(2);
                   int ibra = tbc_pp.GetLocalIndex(min(a,b),max(a,b));
                   int iket = tbc_pp.GetLocalIndex(min(c,d),max(c,d));
                   cout << "    PP: ibra,iket = " << ibra << "," << iket << endl;
                   Ket * bra = tbc_pp.GetKet(ibra);
                   Ket * ket = tbc_pp.GetKet(iket);
                   if (a>b) tbme_pp *= bra->Phase(J+0);
                   if (c>d) tbme_pp *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_pp](ibra,iket) = tbme_pp;
                   Hbare.TwoBody[ch_pp](iket,ibra) = tbme_pp;
             }

                   // then do nn
             if (tbme_nn !=0)
             {
                  // normalize
                   if (a==b) tbme_nn /= sqrt(2);
                   if (c==d) tbme_nn /= sqrt(2);
                   int ibra = tbc_nn.GetLocalIndex(min(a,b)+1,max(a,b)+1);
                   int iket = tbc_nn.GetLocalIndex(min(c,d)+1,max(c,d)+1);
                   cout << "    NN: ibra,iket = " << ibra << "," << iket << endl;
                   Ket * bra = tbc_nn.GetKet(ibra);
                   Ket * ket = tbc_nn.GetKet(iket);
                   if (a>b) tbme_nn *= bra->Phase(J+0);
                   if (c>d) tbme_nn *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_nn](ibra,iket) = tbme_nn;
                   Hbare.TwoBody[ch_nn](iket,ibra) = tbme_nn;
             }

                   // now do pnpn

             if (tbme_pnpn !=0)
             {
                   int ibra = tbc_pn.GetLocalIndex(min(a,b+1),max(a,b+1));
                   int iket = tbc_pn.GetLocalIndex(min(c,d+1),max(c,d+1));
                   cout << "    PNPN: ibra,iket = " << ibra << "," << iket << endl;
                   Ket * bra = tbc_pn.GetKet(ibra);
                   Ket * ket = tbc_pn.GetKet(iket);
                   int phase = 1;
                   if (a > b+1) phase *= bra->Phase(J+0);
                   if (c > d+1) phase *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_pn](ibra,iket) = phase * tbme_pnpn;
                   Hbare.TwoBody[ch_pn](iket,ibra) = phase * tbme_pnpn;


                   // now do npnp

                   ibra = tbc_pn.GetLocalIndex(min(a+1,b),max(a+1,b));
                   iket = tbc_pn.GetLocalIndex(min(c+1,d),max(c+1,d));
                   cout << "    NPNP: ibra,iket = " << ibra << "," << iket << endl;
                   bra = tbc_pn.GetKet(ibra);
                   ket = tbc_pn.GetKet(iket);
                   phase = 1;
                   if (a+1 > b) phase *= bra->Phase(J+0);
                   if (c+1 > d) phase *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_pn](ibra,iket) = phase * tbme_pnpn;
                   Hbare.TwoBody[ch_pn](iket,ibra) = phase * tbme_pnpn;
             }

                   // now do pnnp

             if (tbme_pnnp !=0)
             {
                   int ibra = tbc_pn.GetLocalIndex(min(a,b+1),max(a,b+1));
                   int iket = tbc_pn.GetLocalIndex(min(c+1,d),max(c+1,d));
                   cout << "    PNNP: ibra,iket = " << ibra << "," << iket << endl;
                   Ket * bra = tbc_pn.GetKet(ibra);
                   Ket * ket = tbc_pn.GetKet(iket);
                   int phase = 1;
                   if (a > b+1) phase *= bra->Phase(J+0);
                   if (c+1 > d) phase *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_pn](ibra,iket) = phase * tbme_pnnp;
                   Hbare.TwoBody[ch_pn](iket,ibra) = phase * tbme_pnnp;

                   // now do nppn

                   ibra = tbc_pn.GetLocalIndex(min(a+1,b),max(a+1,b));
                   iket = tbc_pn.GetLocalIndex(min(c,d+1),max(c,d+1));
                   cout << "    NPPN: ibra,iket = " << ibra << "," << iket << endl;
                   bra = tbc_pn.GetKet(ibra);
                   ket = tbc_pn.GetKet(iket);
                   phase = 1;
                   if (a+1 > b) phase *= bra->Phase(J+0);
                   if (c > d+1) phase *= ket->Phase(J+0);
                   Hbare.TwoBody[ch_pn](ibra,iket) = phase * tbme_pnnp;
                   Hbare.TwoBody[ch_pn](iket,ibra) = phase * tbme_pnnp;

             }


          }
        }
      }
    }
  }

  return;
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
   //for (int i=0;i<norbits;i++)
   obfile << " Zero body part: " << op.ZeroBody << endl;
   for (int& i : modelspace->valence )
   {
      for (int& j : modelspace->valence )
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
   int wfloat = 12; // width for printing floats
   //for (int& i : modelspace->holes)
   for (int& i : modelspace->hole_qspace)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      Acore += oi->j2 +1;
      if (oi->tz2 < 0)
      {
//         Zcore += oi->j2+1;
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->hole_qspace.size() - proton_core_orbits;

   intfile << "! shell model effective interaction generated by IMSRG" << endl;
   intfile << "! Zero body term: " << op.ZeroBody << endl;
   intfile << "! Index   nljtz" << endl;
   // first do proton orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      //int nushell_indx = (i-ncore_orbits)/2 + 1;
      int nushell_indx = i/2+1 -proton_core_orbits;
      intfile << "!  " << nushell_indx << "   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
//      intfile << "!  " << nushell_indx << "(" << i << ")   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
      ++nvalence_proton_orbits;
   }
   // then do neutron orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
      int nushell_indx = i/2+1 + nvalence_proton_orbits -neutron_core_orbits;
//      int nushell_indx = (i-ncore_orbits)/2 + nvalence_proton_orbits + 1;
      intfile << "!  " << nushell_indx << "   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
      //intfile << "!  " << nushell_indx << "(" << i << ")   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
   }

   intfile << "!" << endl;
   intfile << "-999  ";
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      intfile << op.OneBody(i,i) << "  ";
   }
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
      intfile << op.OneBody(i,i) << "  ";
   }
   intfile << "  " << Acore << " " << Acore+2 << "  0.00000 " << endl; // No mass dependence for now...

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      for (int& ibra: tbc.KetIndex_vv)
      {
         Ket *bra = tbc.GetKet(ibra);
         int a = bra->p;
         int b = bra->q;
         Orbit* oa = modelspace->GetOrbit(a);
         Orbit* ob = modelspace->GetOrbit(b);
//      int nushell_indx = i/2-proton_core_orbits + 1;
//      int nushell_indx = i/2-neutron_core_orbits + nvalence_proton_orbits + 1;
         int a_ind = a/2+1 + ( oa->tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int b_ind = b/2+1 + ( ob->tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//         int a_ind = (a-ncore_orbits)/2 + (oa->tz2+1)/2*(nvalence_proton_orbits) +1;
//         int b_ind = (b-ncore_orbits)/2 + (ob->tz2+1)/2*(nvalence_proton_orbits) +1;
         for (int& iket: tbc.KetIndex_vv)
         {
            if (iket<ibra) continue;
            Ket *ket = tbc.GetKet(iket);
            int c = ket->p;
            int d = ket->q;
            Orbit* oc = modelspace->GetOrbit(c);
            Orbit* od = modelspace->GetOrbit(d);
//         int c_ind = c/2+1 + ( tbc.Tz <1 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//         int d_ind = d/2+1 + ( tbc.Tz <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int c_ind = c/2+1 + ( oc->tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
         int d_ind = d/2+1 + ( od->tz2 <0 ? -proton_core_orbits : nvalence_proton_orbits - neutron_core_orbits);
//            int c_ind = (c-ncore_orbits)/2 + (oc->tz2+1)/2*(nvalence_proton_orbits) +1;
//            int d_ind = (d-ncore_orbits)/2 + (od->tz2+1)/2*(nvalence_proton_orbits) +1;
            int T = abs(tbc.Tz);
            double tbme = op.TwoBody[ch](ibra,iket);
            if (T==0)
            {
               if (oa->j2 == ob->j2 and oa->l == ob->l and oa->n == ob->n) T = (tbc.J+1)%2;
               else tbme *= sqrt(2); // pn TBMEs are unnormalized
               if (oc->j2 == od->j2 and oc->l == od->l and oc->n == od->n) T = (tbc.J+1)%2;
               else tbme *= sqrt(2); // pn TBMEs are unnormalized
            }
            // in NuShellX, the proton orbits must come first.
            if (a_ind > b_ind)
            {
               intfile << setw(wint) << b_ind << setw(wint) << a_ind;
               tbme *= bra->Phase(tbc.J);
            }
            else
            {
               intfile << setw(wint) << a_ind << setw(wint) << b_ind;
            }
            if (c_ind > d_ind)
            {
               intfile << setw(wint) << d_ind << setw(wint) << c_ind;
               tbme *= ket->Phase(tbc.J);
            }
            else
            {
               intfile << setw(wint) << c_ind << setw(wint) << d_ind;
            }
            intfile
              << setw(wint) << tbc.J
              << setw(wint) << T
              << setw(wfloat) << tbme
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
   int wfloat = 12; // width for printing floats
   //for (int& i : modelspace->holes)
   for (int& i : modelspace->hole_qspace)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      Acore += oi->j2 +1;
      if (oi->tz2 < 0)
      {
         Zcore += oi->j2+1;
         proton_core_orbits += 1;
      }
   }
   neutron_core_orbits = modelspace->hole_qspace.size() - proton_core_orbits;
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0)
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
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      //int nushell_indx = (i-ncore_orbits)/2 + 1;
      int nushell_indx = i/2-proton_core_orbits + 1;
      spfile << nushell_indx << " " << oi->n+1 << " " << oi->l << " " << oi->j2  << endl;
   }
   // then do neutron orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
      //int nushell_indx = (i-neutron_core_orbits) + nvalence_proton_orbits + 1;
      int nushell_indx = i/2-neutron_core_orbits + nvalence_proton_orbits + 1;
      spfile << nushell_indx << " " << oi->n+1 << " " << oi->l << " " << oi->j2 << endl;
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
   int wfloat = 12; // width for printing floats
   for (int& i : modelspace->holes)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      Acore += oi->j2 +1;
   }
   intfile << "IMSRG INTERACTION" << endl;
   // 2 indicates pn treated separately
   intfile << "2 " << nvalence_orbits << " ";
   // write NLJ of the orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      intfile << oi->n*1000 + oi->l*100 + oi->j2 << " ";
   }
   intfile << endl;
   // Write proton SPE's
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      intfile << op.OneBody(i,i) << " ";
   }
   // Write neutron SPE's
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
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
      //TwoBodyChannel * tbc = op.GetTwoBodyChannel(ch);
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      for (int i=0;i<npq;++i)
      {
         Ket *bra = tbc.GetKet(i);
         Orbit *oa = modelspace->GetOrbit(bra->p);
         Orbit *ob = modelspace->GetOrbit(bra->q);
         for (int j=i;j<npq;++j)
         {
            Ket *ket = tbc.GetKet(j);
            double tbme = op.GetTBME(ch,bra,ket) / sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()));
            if ( abs(tbme)<1e-4 ) continue;
            Orbit *oc = modelspace->GetOrbit(ket->p);
            Orbit *od = modelspace->GetOrbit(ket->q);
            int wint = 4;
            int wfloat = 12;

            tbfile 
                   << setw(wint) << bra->p
                   << setw(wint) << bra->q
                   << setw(wint) << ket->p
                   << setw(wint) << ket->q
                   << setw(wint+3) << tbc.J << setw(wfloat) << std::fixed << tbme
//                   << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >"
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
   //int nchan = op.GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      //TwoBodyChannel * tbc = op.GetTwoBodyChannel(ch);
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;
      //for (int i=0;i<npq;++i)
      for (int& i : tbc.KetIndex_vv )
      {
         Ket *bra = tbc.GetKet(i);
         Orbit *oa = modelspace->GetOrbit(bra->p);
         Orbit *ob = modelspace->GetOrbit(bra->q);
         //for (int j=i;j<npq;++j)
         for (int& j : tbc.KetIndex_vv )
         {
            if (j<i) continue;
            Ket *ket = tbc.GetKet(j);
            //double tbme = tbc.GetTBME(bra,ket);
            double tbme = op.GetTBME(ch,bra,ket) / sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()));
            if ( abs(tbme)<1e-4 ) continue;
            Orbit *oc = modelspace->GetOrbit(ket->p);
            Orbit *od = modelspace->GetOrbit(ket->q);
            int wint = 4;
            int wfloat = 12;

            tbfile 
                   << setw(wint) << bra->p
                   << setw(wint) << bra->q
                   << setw(wint) << ket->p
                   << setw(wint) << ket->q
                   << setw(wint+3) << tbc.J << setw(wfloat) << std::fixed << tbme// << endl;
                   << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
         }
      }
   }
   tbfile.close();
}




