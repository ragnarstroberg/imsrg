
#include "Operator3.hh"
#include <algorithm> // used for sort


/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator3::Operator3()
{
   modelspace = NULL;
   nChannels = 0;
   hermitian = true;
   antihermitian = false;
   E3max = 6;
}

Operator3::Operator3(ModelSpace& ms) // Create a zero-valued operator in a given model space
{
  modelspace = &ms;
  E3max = 6;
  hermitian = true;
  antihermitian = false;
  ZeroBody = 0;
  int nOneBody = modelspace->GetNumberOrbits();
  OneBody = arma::mat(nOneBody,nOneBody,arma::fill::zeros);
  nChannels = modelspace->GetNumberTwoBodyChannels();
  TwoBody = vector<arma::mat>(modelspace->GetNumberTwoBodyChannels(), arma::mat());

  for (int ch=0;ch<nChannels;++ch)
  {
      int npq = modelspace->GetTwoBodyChannel(ch).GetNumberKets();
      TwoBody[ch] = arma::mat(npq,npq,arma::fill::zeros);
  }
  ThreeBody.resize(5*2*modelspace->GetThreeBodyJmax());
  AllocateThreeBody();

}

/////////// COPY METHOD //////////////////////////
void Operator3::Copy(const Operator3& rhs)
{
   modelspace    = rhs.modelspace;
   nChannels     = rhs.nChannels;
   hermitian     = rhs.hermitian;
   antihermitian = rhs.antihermitian;
   ZeroBody      = rhs.ZeroBody;
   OneBody       = rhs.OneBody;
   TwoBody       = rhs.TwoBody;
   ThreeBody     = rhs.ThreeBody;
}


// Confusing nomenclature: J2 means 2 times the total J of the three body system
void Operator3::AllocateThreeBody()
{
   int counter=0;
   cout << "Begin AllocateThreeBody" << endl;
  int norbits = modelspace->GetNumberOrbits();
  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   if ((2*oa.n + oa.l)>E3max) continue;
   for (int b=0; b<=a; b+=2)
   {
    Orbit& ob = modelspace->GetOrbit(b);
    if ((2*oa.n+oa.l + 2*ob.n+ob.l) > E3max) continue;
    int Jab_min = abs(oa.j2-ob.j2)/2;
    int Jab_max = (oa.j2+ob.j2)/2;
    int tab_min = 0;
//    int tab_min = abs(oa.tz2+ob.tz2)/2;
    int tab_max = 1;
     for (int c=0; c<=b; c+=2)
     {
      Orbit& oc = modelspace->GetOrbit(c);
      if ((2*oa.n+oa.l + 2*ob.n+ob.l + 2*oc.n+oc.l) > E3max) continue;
//      int J2_min = min(abs(2*Jab_min-oc.j2),abs(2*Jab_max-oc.j2));
//      int J2_max = 2*Jab_max+oc.j2;
      int J2_min = 1;
      int J2_max = oa.j2+ob.j2+oc.j2;
      int parity_abc = (oa.l+ob.l+oc.l)%2;
 

       // Begin loop over ket states
       for( int d=0; d<=a; d+=2)
       {
        Orbit& od = modelspace->GetOrbit(d);
        for (int e=0; e<=d; e+=2)
        {
         Orbit& oe = modelspace->GetOrbit(e);
         int Jde_min = abs(od.j2-oe.j2)/2;
         int Jde_max = (od.j2+oe.j2)/2;
         int tde_min = 0;
//         int tde_min = abs(od.tz2+oe.tz2)/2;
         int tde_max = 1;
         for (int f=0; f<=e; f+=2)
         {
           Orbit& of = modelspace->GetOrbit(f);

           // conserve parity
           int parity_def = (od.l+oe.l+of.l)%2;
           if (parity_def != parity_abc) continue;

           // Get hashing index for orbits
           long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);
  
           for (int J2=J2_min; J2<=J2_max; J2+=2)
           {
             // check triangle condition
//             if (orbit_index==560)
//             {
//                cout << "ORB_INDX = 560. J2 = " << J2 << "  Jde_min = " << Jde_min << " Jde_max  = " << Jde_max << endl;
//             }
             if ( (J2+of.j2)<(2*Jde_min) or abs(J2-of.j2)>(2*Jde_max)) continue;
//             if (orbit_index==49)
//             {
//                cout << "   survived the triangle condition" << endl;
//             }

             // Find actually allowed ranges of Jab and Jde for this J2
             int Jab_min_min = max(Jab_min,abs(J2-oc.j2)/2);
             int Jab_max_max = min(Jab_max,(J2+oc.j2)/2);
             int Jde_min_min = max(Jde_min,abs(J2-of.j2)/2);
             int Jde_max_max = min(Jde_max,(J2+of.j2)/2);
//             if (Jab_min_min > Jab_max_max or Jde_min_min > Jde_max_max) continue;           

             int nJ_twobody = (Jab_max_max-Jab_min_min+1)*(Jde_max_max-Jde_min_min+1);
//             if (orbit_index==560)
//             {
//                cout << "  nJ_twobody = " << nJ_twobody << endl;
//             }
             if (nJ_twobody <0) continue;
             vector<double> zerovector(nJ_twobody,0.0);

             for (int tab=tab_min; tab<=tab_max; ++tab)
             {
              for (int tde=tde_min; tde<=tde_max; ++tde)
              {
                for (int T2=1;T2<=1+2*(tab*tde); T2+=2) // T can be 3/2 if tab=tde=1
                {
                  // A valid matrix element. Now set it to zero.
                  int ch = GetThreeBodyChannelIndex(J2,parity_abc,tab,tde,T2);
//                   if (orbit_index==560)
//                   {
//                      cout << "ORBIT_INDEX 560,  ch=" << ch << " is ok.    "
//                           << J2_min << " " << J2_max << " "
//                           << tab_min << " " << tab_max << " "
//                           << tde_min << " " << tde_max << " "
//                           << "   " << J2 << " " << tab << " " << tde << " " << T2 << endl;
//                   }
//                  cout << "Allocating another zero vector" << endl;
//                  ThreeBody[ch][orbit_index] = vector<double>(nJ_twobody,0.0);
                  ThreeBody[ch][orbit_index] = zerovector;
//                  cout << "(" << J2 << "," << parity_abc << "," << tab << "," << tde << "," << T2 << ")"
//                       << "   ch = " << ch << " orbit_index = " << orbit_index << "  size = " << ThreeBody[ch].at(orbit_index).size() << endl;
                  counter += nJ_twobody;
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


double Operator3::GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;
   SortOrbits(a,b,c,d,e,f);

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);

   int Jab_min = max(abs(oa.j2-ob.j2),abs(J2-oc.j2))/2;
   int Jab_max = min(oa.j2+ob.j2,J2+oc.j2)/2;
   int Jde_min = max(abs(od.j2-oe.j2),abs(J2-of.j2))/2;
   int Jde_max = min(od.j2+oe.j2,J2+of.j2)/2;

   int tab_min = 0;
   int tab_max = 1;
   int tde_min = 0;
   int tde_max = 1;

   int parity = ( oa.l + ob.l + oc.l ) % 2;
   long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

   double V = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2);
        double Cj_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2);

        int J2index = (Jab_max-Jab_min)*(Jde-Jde_min) + (Jab-Jab_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2);
            double Ct_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2);
            int ch = GetThreeBodyChannelIndex(J2,parity,tab,tde,T2);

            V += Cj_abc * Cj_def * Ct_abc * Ct_def * ThreeBody[ch].at(orbit_index)[J2index];
         }
        }
      }
   }
   return V;

}

void Operator3::SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;
   SortOrbits(a,b,c,d,e,f);

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);

   int Jab_min = max(abs(oa.j2-ob.j2),abs(J2-oc.j2))/2;
   int Jab_max = min(oa.j2+ob.j2,J2+oc.j2)/2;
   int Jde_min = max(abs(od.j2-oe.j2),abs(J2-of.j2))/2;
   int Jde_max = min(od.j2+oe.j2,J2+of.j2)/2;

   int tab_min = 0;
   int tab_max = 1;
   int tde_min = 0;
   int tde_max = 1;

   int parity = ( oa.l + ob.l + oc.l ) % 2;
   long int orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

//   cout << Jab_min << "<= Jab <= " << Jab_max << "   " 
//        << Jde_min << "<= Jde <= " << Jde_max << endl;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2);
        double Cj_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2);

        int J2index = (Jab_max-Jab_min)*(Jde-Jde_min) + (Jab-Jab_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2);
            double Ct_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2);
            int ch = GetThreeBodyChannelIndex(J2,parity,tab,tde,T2);
      if (a==4 and b==0 and c==0 and d==4 and e==0 and f==0 and Cj_abc*Cj_def*Ct_abc*Ct_def !=0)
      {
           cout << "ch = " << ch << "  orbit_index = " << orbit_index << "  J2index = " << J2index << endl;
//         cout << "Cj_abc = " << Cj_abc << "   Cj_def = " << Cj_def << "  Ct_abc = " << Ct_abc << "  Ct_def = " << Ct_def << endl;
      }

            ThreeBody[ch].at(orbit_index)[J2index] += Cj_abc * Cj_def * Ct_abc * Ct_def * V;

        }
      }
    }
  }

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
   long int orb_indx = (orbit_index_bra * (orbit_index_bra+1))/2 + orbit_index_ket;
//   cout << "GETTING ORBIT INDEX " << a << b << c << d << e << f << "   "
//        << orbit_index_bra << " " << orbit_index_ket << "  " << orb_indx << endl;
   return orb_indx;
}

void Operator3::SortOrbits(int& a, int& b, int& c, int& d, int& e, int& f)
{

   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);

   if (d<e)  swap(d,e);
   if (e<f)  swap(e,f);
   if (d<e)  swap(d,e);

   if (d>a)
   {
     swap(a,d);
     swap(b,e);
     swap(c,f);
   }
}


double Operator3::RecouplingCoefficient(int a_in, int b_in, int c_in, int a, int b, int c, int Jab_in, int Jab, int J)
{
   double C;
   if (a==a_in)
   {
      if (b==b_in) // (abc)_in -> (abc)
      {
         C = Jab==Jab_in ? 1 : 0;
      }
      else  // (acb)_in -> (abc)
      {
         Orbit& oa = modelspace->GetOrbit(a_in);
         Orbit& ob = modelspace->GetOrbit(b_in);
         Orbit& oc = modelspace->GetOrbit(c_in);
         C = modelspace->phase( (ob.j2+oc.j2)/2+Jab-Jab_in) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ob.j2/2., oa.j2/2., Jab_in, oc.j2/2., J/2., Jab);
      }
   }
   else if (a==b_in)
   {
      Orbit& oa = modelspace->GetOrbit(a_in);
      Orbit& ob = modelspace->GetOrbit(b_in);
      Orbit& oc = modelspace->GetOrbit(c_in);
      if (b==a_in) // (bac)_in -> (abc)
      {
         C = Jab==Jab_in ? modelspace->phase((oa.j2+ob.j2)/2-Jab) : 0;
      }
      else // (bca)_in -> (abc)
      {
         C = modelspace->phase( (ob.j2+oc.j2)/2+Jab+1 ) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(oa.j2/2., ob.j2/2., Jab_in, oc.j2/2., J/2., Jab);
      }
   }
   else
   {
      Orbit& oa = modelspace->GetOrbit(a_in);
      Orbit& ob = modelspace->GetOrbit(b_in);
      Orbit& oc = modelspace->GetOrbit(c_in);
      if (b==a_in) // (cab)_in -> (abc)
      {
         C = modelspace->phase( (oa.j2+ob.j2)/2-Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ob.j2/2., oa.j2/2., Jab_in, oc.j2/2., J/2., Jab);
      }
      else // (cba)_in -> (abc)
      {
         C =  - sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(oa.j2/2., ob.j2/2., Jab_in, oc.j2/2., J/2., Jab);
      }
   }
   return C;

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


