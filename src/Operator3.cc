
#include "Operator3.hh"
#include <algorithm> // used for sort


/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator3::Operator3()
{
   modelspace = NULL;
   nChannels = 0;
   hermitian = true;
   antihermitian = false;
}

Operator3::Operator3(ModelSpace& ms) // Create a zero-valued operator in a given model space
{
  cout << "Constructing Operator3" << endl;
  modelspace = &ms;
  hermitian = true;
  antihermitian = false;
  ZeroBody = 0;
  cout << "ZeroBody went ok" << endl;
  int nOneBody = modelspace->GetNumberOrbits();
//  int nKets = modelspace->GetNumberKets();
  OneBody = arma::mat(nOneBody,nOneBody,arma::fill::zeros);
  cout << "OneBody went ok" << endl;
  nChannels = modelspace->GetNumberTwoBodyChannels();
  cout << "nchannels = " << nChannels << endl;
  TwoBody = vector<arma::mat>(modelspace->GetNumberTwoBodyChannels(), arma::mat());

  for (int ch=0;ch<nChannels;++ch)
  {
      cout << "ch = " << ch << endl;
      int npq = modelspace->GetTwoBodyChannel(ch).GetNumberKets();
      cout << "npq = " << npq << endl;
      TwoBody[ch] = arma::mat(npq,npq,arma::fill::zeros);
      cout << "success!" << endl;
  }
  cout << "TwoBody went ok" << endl;
  AllocateThreeBody();

}

/////////// COPY METHOD //////////////////////////
void Operator3::Copy(const Operator3& op)
{
   modelspace    = op.modelspace;
   nChannels     = op.nChannels;
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   TwoBody       = op.TwoBody;
   ThreeBody     = op.ThreeBody;
}


// Confusing nomenclature: J2 means 2 times the total J of the three body system
void Operator3::AllocateThreeBody()
{
   cout << "Begin AllocateThreeBody" << endl;
  ThreeBody.resize(5*2*modelspace->GetThreeBodyJmax());
  int norbits = modelspace->GetNumberOrbits();
  for (int a=0; a<norbits; ++a)
  {
   Orbit oa = modelspace->GetOrbit(a);
   for (int b=0; b<=a; ++b)
   {
    Orbit ob = modelspace->GetOrbit(b);
    int Jab_min = abs(oa.j2-ob.j2)/2;
    int Jab_max = (oa.j2+ob.j2)/2;
    int tab_min = abs(oa.tz2+ob.tz2)/2;
    int tab_max = 1;
     for (int c=0; c<=b; ++c)
     {
      Orbit oc = modelspace->GetOrbit(c);
      int J2_min = min(abs(2*Jab_min-oc.j2),abs(2*Jab_max-oc.j2));
      int J2_max = 2*Jab_max+oc.j2;
      int parity_abc = (oa.l+ob.l+oc.l)%2;
 

       // Begin loop over ket states
       for( int d=0; d<=a; ++d)
       {
        Orbit od = modelspace->GetOrbit(d);
        for (int e=0; e<=d; ++e)
        {
         Orbit oe = modelspace->GetOrbit(e);
         int Jde_min = abs(od.j2-oe.j2)/2;
         int Jde_max = (od.j2+oe.j2)/2;
         int tde_min = abs(od.tz2+oe.tz2)/2;
         int tde_max = 1;
         for (int f=0; f<=e; ++f)
         {
           Orbit of = modelspace->GetOrbit(f);

           // conserve parity
           int parity_def = (od.l+oe.l+of.l)%2;
           if (parity_def != parity_abc) continue;

           // Get hashing index for orbits
           long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);
  
           for (int J2=J2_min; J2<=J2_max; J2+=2)
           {
             // check triangle condition
             if ( (J2+of.j2)<Jde_min or abs(J2-of.j2)>Jde_max) continue;

             // Find actually allowed ranges of Jab and Jde for this J2
             int Jab_min_min = max(Jab_min,abs(J2-oc.j2)/2);
             int Jab_max_max = min(Jab_max,(J2+oc.j2)/2);
             int Jde_min_min = max(Jde_min,abs(J2-of.j2)/2);
             int Jde_max_max = min(Jde_max,(J2+of.j2)/2);
             

             int nJ_twobody = (Jab_max_max-Jab_min_min)*(Jde_max_max-Jde_min_min);
             if (nJ_twobody <1) continue;
             vector<double> zerovector(nJ_twobody,0.0);

             for (int tab=tab_min; tab<=tab_max; ++tab)
             {
              for (int tde=tde_min; tde<=tde_max; ++tde)
              {
                for (int T2=1;T2<=1+2*(tab*tde); T2+=2) // T can be 3/2 if tab=tde=1
                {
                  // A valid matrix element. Now set it to zero.
                  int ch = GetThreeBodyChannelIndex(J2,parity_abc,tab,tde,T2);
                  cout << "Allocating another zero vector" << endl;
                  ThreeBody[ch][orbit_index] = zerovector;
                } //T2
              } //tde 
             } //tab
        } //J2
       } //f
      } //e
     } //d
    } //c
   } //b
  } //a


}


double Operator3::GetThreeBodyME(int Jab, int Jde, int J2, int tab, int tde, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;
   int phase = 1;
   SortWithPhase(a,b,c,d,e,f,phase);

   Orbit oa = modelspace->GetOrbit(a);
   Orbit ob = modelspace->GetOrbit(b);
   Orbit oc = modelspace->GetOrbit(c);
   Orbit od = modelspace->GetOrbit(d);
   Orbit oe = modelspace->GetOrbit(e);
   Orbit of = modelspace->GetOrbit(f);

   // Flatten 12 indices to 3: (J,tab,tde,T2,parity), (a,b,c,d,e,f), (Jab,Jde)

   int parity = ( oa.l + ob.l + oc.l ) % 2;
   int ch = GetThreeBodyChannelIndex(J2,parity,tab,tde,T2);

   // Maybe it's better to have a map< vector<int>, vector<double> >
   long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);
   

   int J2index = (J2+oc.j2)*Jde + Jab;

   // remember to include phase
   return ThreeBody[ch][orbit_index][J2index] * phase;

//   vector<int> element = {Jab,Jde,J2,tab,tde,T2,a,b,c,d,e,f};
//   return ThreeBodyMatrixElements.at(element)*phase;

}

void Operator3::SetThreeBodyME(int Jab, int Jde, int J2, int tab, int tde, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V)
{
   // reorder so a>=b>=c and d>=e>=f
   cout << "SetThreeBodyME" << endl;
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;
   int phase = 1;
   SortWithPhase(a,b,c,d,e,f,phase);

   Orbit oa = modelspace->GetOrbit(a);
   Orbit ob = modelspace->GetOrbit(b);
   Orbit oc = modelspace->GetOrbit(c);
   Orbit od = modelspace->GetOrbit(d);
   Orbit oe = modelspace->GetOrbit(e);
   Orbit of = modelspace->GetOrbit(f);

   // Flatten 12 indices to 3: (J,tab,tde,T2,parity), (a,b,c,d,e,f), (Jab,Jde)

   int parity = ( oa.l + ob.l + oc.l ) % 2;
   int ch = GetThreeBodyChannelIndex(J2,parity,tab,tde,T2);

   // Maybe it's better to have a map< vector<int>, vector<double> >
   long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);
   

   int J2index = (J2+oc.j2)*Jde + Jab;

   // remember to include phase
   ThreeBody[ch][orbit_index][J2index] = V * phase;

//   vector<int> element = {Jab,Jde,J2,tab,tde,T2,a,b,c,d,e,f};
//   return ThreeBodyMatrixElements.at(element)*phase;

}







int Operator3::GetThreeBodyChannelIndex(int J2, int parity, int tab, int tde, int T2)
{
   int isospin_index = tde + 2*tab + (T2-1)/2; // ranges from 0-4
   return (2*isospin_index + parity) * modelspace->GetThreeBodyJmax() + J2;

}


long int Operator3::GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f)
{
   int orbit_index_bra =  a*(a+1)*(a+2)/6 + b*(b+1)/2 + c;
   int orbit_index_ket =  d*(d+1)*(d+2)/6 + e*(e+1)/2 + f;
   long int orb_indx = (orbit_index_bra * orbit_index_bra+1)/2 + orbit_index_ket;
   return orb_indx;
}

void Operator3::SortWithPhase(int& a, int& b, int& c, int& d, int& e, int& f, int& phase)
{

   /// NEED TO FIGURE OUT THE PHASES FOR EACH SWAP??
   if (a<b)
   {
     swap(a,b);
     phase *=1;
   } 
   if (b<c)
   {
     swap(b,c);
     phase *= 1;
   }
   if (a<b)
   {
     swap(a,b);
     phase *=1;
   } 
   if (d<e)
   {
     swap(d,e);
     phase *=1;
   } 
   if (e<f)
   {
     swap(e,f);
     phase *= 1;
   }
   if (d<e)
   {
     swap(d,e);
     phase *=1;
   }

   if (d>a)
   {
     swap(a,d);
     swap(b,e);
     swap(c,f);
   }
}

//   The normal ordered two body piece is 
//   Gamma(2)^J_ijkl = V(2)^J_ijkl + Sum_a n_a  Sum_K (2K+1)/(2J+1) V(3)^JJK_ijakla
//
//
Operator Operator3::DoNormalOrdering3()
{
   Operator3 opNO3 = Operator3(*modelspace);
//   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      for (int ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         for (int iket=ibra; iket<tbc.GetNumberKets(); ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            for (int& a : modelspace->holes)
            {
               Orbit & oa = modelspace->GetOrbit(a);
               int Tz2 = 2*tbc.Tz + oa.tz2;
               int kmin2 = abs(2*tbc.J-oa.j2);
               int kmax2 = 2*tbc.J+oa.j2;
               int parity = tbc.parity * (oa.l%2);
               for (int K2=kmin2; K2<=kmax2; ++K2)
               {
//                  #pragma omp critical
//                  opNO3.TwoBody[ch](ibra,iket) += (K2+1) * GetThreeBodyME(tbc.J,tbc.J,K2,Tz2,parity,i,j,a,k,l,a);
               }
            }
            opNO3.TwoBody[ch](ibra,iket) /= (2*tbc.J+1);
         }
      }
   }
   Operator opNO2 = opNO3.DoNormalOrdering();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);
   opNO2 += DoNormalOrdering();
   return opNO2;

}


