#include "ReadWrite.hh"
#include "ModelSpace.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "string.h"

#define LINESIZE 400

using namespace std;


void ReadWrite::ReadSettingsFile( char* filename)
{
   char line[LINESIZE];
   string lstr;
   ifstream fin;
   fin.open(filename);
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




ModelSpace ReadWrite::ReadModelSpace( const char* filename)
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
   cout << "done reading header" << endl;

   while (  infile >> cbuf[0] >> ibuf[0] >> n >> l >> j2 >> tz2
             >> ibuf[1] >> spe >> fbuf[0]  >> cbuf[1] >> cbuf[2])
   {
      int ph = 0; // 0=particle, 1=hole 
      int io = 0; // 0=inside, 1=outside projected model space
      if (strstr(cbuf[1],"hole")) ph++;
      if (strstr(cbuf[2],"outside")) io++;
      modelspace.AddOrbit( Orbit(n,l,j2,tz2,ph,io,spe) );

   }
   cout << "done reading interaction" << endl;
   
   infile.close();
   modelspace.SetupKets();

   return modelspace;
}


void ReadWrite::ReadBareTBME( const char* filename, Operator& Hbare)
{

  ifstream infile;
  char line[LINESIZE];
  int Tz,Par,J2,a,b,c,d;
  float fbuf[3];
  float tbme;

  infile.open(filename);
  if ( !infile.good() )
  {
     cerr << "Trouble reading file " << filename << endl;
     return;
  }

  infile.getline(line,LINESIZE);

  while (!strstr(line,"<ab|V|cd>") && !infile.eof()) // Skip lines until we see the header
  {
     infile.getline(line,LINESIZE);
  }

  while ( infile >> Tz >> Par >> J2  // Read tbme
                 >> a >> b >> c >> d
                 >> tbme >> fbuf[0] >> fbuf[1] >> fbuf[2] )
  {
     a--; b--; c--; d--; // Fortran -> C  ==> 1 -> 0

     tbme -= fbuf[2] * Hbare.GetModelSpace()->GetHbarOmega() / Hbare.GetModelSpace()->GetTargetMass();  // Some sort of COM correction. Check this

// NORMALIZATION: Read in normalized, antisymmetrized TBME's

     Hbare.SetTBME(J2/2,Par,Tz,a,b,c,d,tbme);
  }

  return;
}


void ReadWrite::CalculateKineticEnergy(Operator *Hbare)
{
   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   int A = Hbare->GetModelSpace()->GetTargetMass();
   float hw = Hbare->GetModelSpace()->GetHbarOmega();
   for (int a=0;a<norbits;a++)
   {
      Orbit * orba = Hbare->GetModelSpace()->GetOrbit(a);
      Hbare->OneBody(a,a) = 0.5 * hw * (2*orba->n + orba->l +3./2); 
      for (int b=0;b<norbits;b++)  // make this better once OneBodyChannel is implemented
      {
         Orbit * orbb = Hbare->GetModelSpace()->GetOrbit(b);
         if (orba->l == orbb->l and orba->j2 == orbb->j2 and orba->tz2 == orbb->tz2)
         {
            if (orba->n == orbb->n)
               Hbare->OneBody(a,a) = 0.5 * hw * (2*orba->n + orba->l +3./2); 
            else if (orba->n == orbb->n+1)
               Hbare->OneBody(a,b) = 0.5 * hw * sqrt( (orba->n)*(orba->n + orba->l +1./2));
            else if (orba->n == orbb->n-1)
               Hbare->OneBody(a,b) = 0.5 * hw * sqrt( (orbb->n)*(orbb->n + orbb->l +1./2));
         }
      }
   }
   Hbare->OneBody *= (1-1./A);

}


void ReadWrite::WriteOneBody(Operator& op, const char* filename)
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

void ReadWrite::WriteValenceOneBody(Operator& op, const char* filename)
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


void ReadWrite::WriteNuShellX_int(Operator& op, const char* filename)
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

   intfile << "! shell model effective interaction generated by IMSRG" << endl;
   intfile << "! Zero body term: " << op.ZeroBody << endl;
   intfile << "! Index   nljtz" << endl;
   // first do proton orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 > 0 ) continue;
      int nushell_indx = (i-ncore_orbits)/2 + 1;
      intfile << "!  " << nushell_indx << "   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
      ++nvalence_proton_orbits;
   }
   // then do neutron orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
      int nushell_indx = (i-ncore_orbits)/2 + nvalence_proton_orbits + 1;
      intfile << "!  " << nushell_indx << "   " << oi->n << " " << oi->l << " " << oi->j2 << "/2" << " " << oi->tz2 << "/2" << endl;
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
         int a_ind = (a-ncore_orbits)/2 + (oa->tz2+1)/2*(nvalence_proton_orbits) +1;
         int b_ind = (b-ncore_orbits)/2 + (ob->tz2+1)/2*(nvalence_proton_orbits) +1;
         for (int& iket: tbc.KetIndex_vv)
         {
            if (iket<ibra) continue;
            Ket *ket = tbc.GetKet(iket);
            int c = ket->p;
            int d = ket->q;
            Orbit* oc = modelspace->GetOrbit(c);
            Orbit* od = modelspace->GetOrbit(d);
            int c_ind = (c-ncore_orbits)/2 + (oc->tz2+1)/2*(nvalence_proton_orbits) +1;
            int d_ind = (d-ncore_orbits)/2 + (od->tz2+1)/2*(nvalence_proton_orbits) +1;
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

void ReadWrite::WriteNuShellX_sps(Operator& op, const char* filename)
{
   ofstream spfile;
   spfile.open(filename, ofstream::out);
   ModelSpace * modelspace = op.GetModelSpace();
   int ncore_orbits = modelspace->holes.size();
   int nvalence_orbits = modelspace->valence.size();
   int nvalence_proton_orbits = 0;
   int Acore = 0;
   int Zcore = 0;
   int wint = 4; // width for printing integers
   int wfloat = 12; // width for printing floats
   for (int& i : modelspace->holes)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      Acore += oi->j2 +1;
      if (oi->tz2 < 0)
      {
         Zcore += oi->j2+1;
      }
   }
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
      int nushell_indx = (i-ncore_orbits)/2 + 1;
      spfile << nushell_indx << " " << oi->n+1 << " " << oi->l << " " << oi->j2  << endl;
   }
   // then do neutron orbits
   for (int& i : modelspace->valence)
   {
      Orbit* oi = modelspace->GetOrbit(i);
      if (oi->tz2 < 0 ) continue;
      int nushell_indx = (i-ncore_orbits)/2 + nvalence_proton_orbits + 1;
      spfile << nushell_indx << " " << oi->n+1 << " " << oi->l << " " << oi->j2 << endl;
   }
   spfile << endl;
   spfile.close();

}




void ReadWrite::WriteAntoine_int(Operator& op, const char* filename)
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


void ReadWrite::WriteTwoBody(Operator& op, const char* filename)
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
            //double tbme = tbc.GetTBME(bra,ket);
            double tbme = op.GetTBME(ch,bra,ket) / sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()));
            if ( abs(tbme)<1e-4 ) continue;
            Orbit *oc = modelspace->GetOrbit(ket->p);
            Orbit *od = modelspace->GetOrbit(ket->q);
            int wint = 4;
            int wfloat = 12;
/*                  tbfile 
                   << setw(wint)   << oa->n  << setw(wint) << oa->l  << setw(wint)<< oa->j2 << setw(wint) << oa->tz2 
                   << setw(wint+2) << ob->n  << setw(wint) << ob->l  << setw(wint)<< ob->j2 << setw(wint) << ob->tz2 
                   << setw(wint+2) << oc->n  << setw(wint) << oc->l  << setw(wint)<< oc->j2 << setw(wint) << oc->tz2 
                   << setw(wint+2) << od->n  << setw(wint) << od->l  << setw(wint)<< od->j2 << setw(wint) << od->tz2 
                   << setw(wint+3) << J << setw(wfloat) << std::fixed << tbme// << endl;
                   << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
                   //<< setw(wint+2) << p   << setw(wint) << Tz << setw(wint) << J << setw(wfloat) << std::fixed << tbme << endl;
*/
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


void ReadWrite::WriteValenceTwoBody(Operator& op, const char* filename)
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
/*                  tbfile 
                   << setw(wint)   << oa->n  << setw(wint) << oa->l  << setw(wint)<< oa->j2 << setw(wint) << oa->tz2 
                   << setw(wint+2) << ob->n  << setw(wint) << ob->l  << setw(wint)<< ob->j2 << setw(wint) << ob->tz2 
                   << setw(wint+2) << oc->n  << setw(wint) << oc->l  << setw(wint)<< oc->j2 << setw(wint) << oc->tz2 
                   << setw(wint+2) << od->n  << setw(wint) << od->l  << setw(wint)<< od->j2 << setw(wint) << od->tz2 
                   << setw(wint+3) << J << setw(wfloat) << std::fixed << tbme// << endl;
                   << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
                   //<< setw(wint+2) << p   << setw(wint) << Tz << setw(wint) << J << setw(wfloat) << std::fixed << tbme << endl;
*/
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




