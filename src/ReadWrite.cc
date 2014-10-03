#include "ReadWrite.hh"

#define LINESIZE 400

using namespace std;

//void ReadWrite::ReadModelSpace( char* filename, ModelSpace * modelspace)
ModelSpace ReadWrite::ReadModelSpace( char* filename)
{

   ModelSpace modelspace;
   ifstream infile;
   stringstream ss;
   char line[LINESIZE];
   char cbuf[10][20];
   int ibuf[2];
   int n,l,j2,tz2,hvq;
   float fbuf[2];
   float spe,hw;
   int A;
   
   infile.open(filename);
   if ( not infile.good() )
   {
      cerr << "Trouble reading file " << filename << endl;
      return 0;
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
//   cout << "I think that A is " << A << endl;
   
   while (!strstr(line,"Oscillator length") && !infile.eof())
   {
      infile.getline(line,LINESIZE);
   }
   ss << line;
   ss >> cbuf[0] >> cbuf[1] >> cbuf[2] >> cbuf[3] >> fbuf[0] >> hw;
   //modelspace->hbar_omega = hw;
   //modelspace->target_mass = A;
   modelspace.hbar_omega = hw;
   modelspace.target_mass = A;
   cout << "Set hbar_omega to " << hw << endl;

   while (!strstr(line,"Legend:") && !infile.eof())
   {
      infile.getline(line,LINESIZE);
   }
   cout << "done reading header" << endl;
   while (!infile.eof())
   {
      infile >> cbuf[0] >> ibuf[0] >> n >> l >> j2 >> tz2
             >> ibuf[1] >> spe >> fbuf[0]  >> cbuf[1] >> cbuf[2];
      hvq = 0; // 0=hole 1=valence 2=particle outside the valence space
      if (strstr(cbuf[1],"particle")) hvq++;
      if (strstr(cbuf[2],"outside")) hvq++;
      //modelspace->AddOrbit( Orbit(n,l,j2,tz2,hvq,spe) );
      modelspace.AddOrbit( Orbit(n,l,j2,tz2,hvq,spe) );
   
   }
   cout << "done reading interaction" << endl;
   
   infile.close();

   return modelspace;
}


//void ReadWrite::ReadBareInteraction( char* filename, Operator *Hbare)
//{
  // Not yet implemented.

//}



//void ReadWrite::ReadBareTBME( char* filename, Operator *Hbare)
void ReadWrite::ReadBareTBME( char* filename, Operator& Hbare)
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
     if (a==2 and b==6 and c==2 and d==6)
     {
        cout << "read TBME " << a << b << c << d << " J=" << J2/2 << " ("<<tbme<<")" << endl;
     }
     a--; b--; c--; d--; // Fortran -> C  => 1->0
     //TwoBodyChannel * h2body = &(Hbare->TwoBody[J2/2][Par][Tz+1]);
     //Ket * bra = Hbare->GetModelSpace()->GetKet(min(a,b),max(a,b));
     //Ket * ket = Hbare->GetModelSpace()->GetKet(min(c,d),max(c,d));
     TwoBodyChannel * h2body = Hbare.GetTwoBodyChannel(J2/2,Par,Tz);
     Ket * bra = Hbare.GetModelSpace()->GetKet(min(a,b),max(a,b));
     Ket * ket = Hbare.GetModelSpace()->GetKet(min(c,d),max(c,d));
     if (bra==NULL or ket==NULL) continue;
     //tbme -= fbuf[2] * Hbare->GetModelSpace()->hbar_omega / Hbare->GetModelSpace()->target_mass;  // Some sort of COM correction. Check this
     tbme -= fbuf[2] * Hbare.GetModelSpace()->hbar_omega / Hbare.GetModelSpace()->target_mass;  // Some sort of COM correction. Check this
     float phase = 1.0;
     if (a==b) phase *= sqrt(2.);   // Symmetry factor. Confirm that this should be here
     if (c==d) phase *= sqrt(2.);   // Symmetry factor. Confirm that this should be here
     if (a>b) phase *= bra->Phase(J2/2);
     if (c>d) phase *= ket->Phase(J2/2);
     //h2body->SetTBME(bra,ket, phase*tbme); // can later change this to [bra,ket] to speed up.
     h2body->SetTBME(bra,ket, phase*tbme);
     //if (a==0 and b==11 and c==1 and d==10)
     if (a==0 and b==11 and c==10 and d==1)
     {
        cout << "Setting TBME " << a << b << c << d << " J=" << J2/2 << " to " << phase*tbme << " ("<<tbme<<")" << endl;
     }
  }

  return;
}


void ReadWrite::CalculateKineticEnergy(Operator *Hbare)
{
   int norbits = Hbare->GetModelSpace()->GetNumberOrbits();
   int A = Hbare->GetModelSpace()->target_mass;
   float hw = Hbare->GetModelSpace()->hbar_omega;
   for (int a=0;a<norbits;a++)
   {
      Orbit * orba = Hbare->GetModelSpace()->GetOrbit(a);
      Hbare->OneBody(a,a) = 0.5 * hw * (2*orba->n + orba->l +3./2); 
      if (orba->n > 0)
      {
         int b = Hbare->GetModelSpace()->Index1(orba->n-1, orba->l, orba->j2, orba->tz2);
         Hbare->OneBody(a,b) = 0.5 * hw * sqrt( (orba->n)*(orba->n + orba->l +1./2));
         Hbare->OneBody(b,a) = Hbare->OneBody(a,b);
      }
   }
   Hbare->OneBody *= (1-1./A);

}



