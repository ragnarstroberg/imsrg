
#include "Operator3.hh"
#include <algorithm> // used for sort
#include "AngMom.hh" // used for CG


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
  AllocateThreeBody();

}

Operator3::Operator3(const Operator3& op)
{
   Copy(op);
}
Operator3::Operator3(const Operator& op)
{
   Copy(op);
}

/////////// COPY METHOD //////////////////////////
void Operator3::Copy(const Operator& rhs)
{
   modelspace    = rhs.modelspace;
   nChannels     = rhs.nChannels;
   hermitian     = rhs.hermitian;
   antihermitian = rhs.antihermitian;
   ZeroBody      = rhs.ZeroBody;
   OneBody       = rhs.OneBody;
   TwoBody       = rhs.TwoBody;
}
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

/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////

Operator3& Operator3::operator=(const Operator3& rhs)
{
   Copy(rhs);
   return *this;
}

// Set Operator3 to Operator2, and set 3body part to zero.
Operator3& Operator3::operator=(const Operator& rhs)
{
   modelspace    = rhs.modelspace;
   nChannels     = rhs.nChannels;
   hermitian     = rhs.hermitian;
   antihermitian = rhs.antihermitian;
   ZeroBody      = rhs.ZeroBody;
   OneBody       = rhs.OneBody;
   TwoBody       = rhs.TwoBody;
   AllocateThreeBody();
   return *this;
}

// multiply operator by a scalar
Operator3& Operator3::operator*=(const double rhs)
{
   ZeroBody *= rhs;
   OneBody *= rhs;
//   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
//   {
//      TwoBody[ch] *= rhs;
//   }
   for (arma::mat& two_body_matrix : TwoBody)
   {
      two_body_matrix *= rhs;
   }
   for ( auto it_orb : ThreeBody)
   {
      for (auto it_J : it_orb.second)
      {
         for (double& it_J2 : it_J)
         {
           it_J2 *= rhs;
         }
      }
   }
   return *this;
}

Operator3 Operator3::operator*(const double rhs) const
{
   Operator3 opout = Operator3(*this);
   opout *= rhs;
   return opout;
}

// Add non-member operator so we can multiply an operator
// by a scalar from the lhs, i.e. s*O = O*s
Operator3 operator*(const double lhs, const Operator3& rhs)
{
   return rhs * lhs;
}


// divide operator by a scalar
Operator3& Operator3::operator/=(const double rhs)
{
   return *this *=(1.0/rhs);
}

Operator3 Operator3::operator/(const double rhs) const
{
   Operator3 opout = Operator3(*this);
   opout /= rhs;
   return opout;
}

// Add Operator2 to Operator3
Operator3& Operator3::operator+=(const Operator& rhs)
{
   ZeroBody += rhs.ZeroBody;
   OneBody += rhs.OneBody;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBody[ch] += rhs.TwoBody[ch];
   }
/*   for ( auto it_orb : ThreeBody)
   {
      for (auto it_J : it_orb)
      {
         for (double& it_J2 : it_J)
         {
           it_J2 += rhs;
         }
      }
   }
*/
   return *this;
}

Operator3 Operator3::operator+(const Operator& rhs) const
{
   return ( Operator3(*this) += rhs );
}

// Subtract operator2 from operator3
Operator3& Operator3::operator-=(const Operator& rhs)
{
   ZeroBody -= rhs.ZeroBody;
   OneBody -= rhs.OneBody;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBody[ch] -= rhs.TwoBody[ch];
   }
/*
   for ( auto it_orb : ThreeBody)
   {
      for (auto it_J : it_orb)
      {
         for (double& it_J2 : it_J)
         {
           it_J2 -= rhs;
         }
      }
   }
*/
   return *this;
}

Operator3 Operator3::operator-(const Operator& rhs) const
{
   return ( Operator3(*this) -= rhs );
}

///////////////////// OTHER METHODS ///////////////////////////////

// Confusing nomenclature: J2 means 2 times the total J of the three body system
void Operator3::AllocateThreeBody()
{
  int norbits = modelspace->GetNumberOrbits();
  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   if ((2*oa.n + oa.l)>E3max) continue;
   for (int b=0; b<=a; b+=2)
   {
    Orbit& ob = modelspace->GetOrbit(b);
    if ((2*oa.n+oa.l + 2*ob.n+ob.l) > E3max) continue;
    for (int c=0; c<=b; c+=2)
    {
     Orbit& oc = modelspace->GetOrbit(c);
     if ((2*oa.n+oa.l + 2*ob.n+ob.l + 2*oc.n+oc.l) > E3max) continue;

     // Begin loop over ket states
     for( int d=0; d<=a; d+=2)
     {
      Orbit& od = modelspace->GetOrbit(d);
      for (int e=0; e<=d; e+=2)
      {
       Orbit& oe = modelspace->GetOrbit(e);
       for (int f=0; f<=e; f+=2)
       {
        Orbit& of = modelspace->GetOrbit(f);

        // conserve parity
        if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0) continue;

        int J2_min = max(max(1,oa.j2-ob.j2-oc.j2), od.j2-oe.j2-of.j2);
        int J2_max = min(oa.j2+ob.j2+oc.j2, od.j2+oe.j2+of.j2);
        if (J2_max < J2_min) continue;

        orbindx_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

        // Get hashing index for orbits
        ThreeBody[orbit_index].resize(J2_max-J2_min+1);

        for (int J2=J2_min; J2<=J2_max; J2+=2)
        {
          // Find actually allowed ranges of Jab and Jde for this J2
          int Jab_min = max(abs(oa.j2-ob.j2), abs(J2-oc.j2))/2;
          int Jab_max = min(oa.j2+ob.j2, J2+oc.j2)/2;
          int Jde_min = max(abs(od.j2-oe.j2), abs(J2-of.j2))/2;
          int Jde_max = min(od.j2+oe.j2, J2+of.j2)/2;
          int J_index = (J2-J2_min)/2;

          int nJ_twobody = (Jab_max-Jab_min+1)*(Jde_max-Jde_min+1);

          if (nJ_twobody <0) continue;
          vector<double> zerovector(nJ_twobody*5,0.0);
          ThreeBody[orbit_index][J_index] = zerovector;

        } //J2
       } //f
      } //e
     } //d
    } //c
   } //b
  } //a
//  for (auto it : ThreeBody)
//  {
//     cout << "it.first = " << it.first << "  it.second.size() = " << it.second.size() << endl;
//  }
}

// Get Three Body Matrix Element in proton-neutron formalism.
double Operator3::GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int Tz, int a, int b, int c, int d, int e, int f)
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);
   double tza = oa.tz2/2.0;
   double tzb = ob.tz2/2.0;
   double tzc = oc.tz2/2.0;
   double tzd = od.tz2/2.0;
   double tze = oe.tz2/2.0;
   double tzf = of.tz2/2.0;
   double Vpn=0;
   int Tmin = abs(Tz);
   for (int tab=abs(tza+tzb); tab<=1; ++tab)
   {
      // CG calculates the Clebsch-Gordan coefficient
      double CG1 = AngMom::CG(0.5,tza, 0.5,tzb, tab, tza+tzb);
      for (int tde=abs(tzd+tze); tde<=1; ++tde)
      {
         double CG2 = AngMom::CG(0.5,tzd, 0.5,tze, tde, tzd+tze);
         if (CG1*CG2==0) continue;
         for (int T=Tmin; T<=3; ++T)
         {
           double CG3 = AngMom::CG(tab,tza+tzb, 0.5,tzc, T/2., tza+tzb+tzc);
           double CG4 = AngMom::CG(tde,tzd+tze, 0.5,tzf, T/2., tzd+tze+tzf);
           if (CG3*CG4==0) continue;
           Vpn += CG1*CG2*CG3*CG4*GetThreeBodyME(Jab_in,Jde_in,J2,tab,tde,T,a-a%2,b-b%2,c-c%2,d-d%2,e-e%2,f-f%2);
         }
      }
   }
   return Vpn;
}

// Get Three Body Matrix Element in isospin formalism (which is how they're stored)
double Operator3::GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;

   SortOrbits(a,b,c);
   SortOrbits(d,e,f);

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);

   int Jindex = ( J2-max(1, max(oa.j2-ob.j2-oc.j2, od.j2-oe.j2-of.j2) ) )/2;

   int Jab_min = max(abs(oa.j2-ob.j2),abs(J2-oc.j2))/2;
   int Jab_max = min(oa.j2+ob.j2,J2+oc.j2)/2;
   int Jde_min = max(abs(od.j2-oe.j2),abs(J2-of.j2))/2;
   int Jde_max = min(od.j2+oe.j2,J2+of.j2)/2;

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;

   orbindx_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

   double V = 0;
   // Recouple J and T to get to the format in which it's stored.
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2,'j');
        double Cj_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2,'j');
        if (Cj_abc*Cj_def == 0) continue;

        int J2index = (a>=d) ? (Jab_max-Jab_min+1)*(Jde-Jde_min) + (Jab-Jab_min)
                             : (Jde_max-Jde_min+1)*(Jab-Jab_min) + (Jde-Jde_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2,'t');
            double Ct_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2,'t');

            int J2Tindex = (a>=d) ? J2index*5 + 2*tab + tde + (T2-1)/2
                                  : J2index*5 + 2*tde + tab + (T2-1)/2;

            V += Cj_abc * Cj_def * Ct_abc * Ct_def * ThreeBody.at(orbit_index).at(Jindex).at(J2Tindex);
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
   SortOrbits(a,b,c);
   SortOrbits(d,e,f);

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);

   int Jindex = ( J2-max(1, max(oa.j2-ob.j2-oc.j2, od.j2-oe.j2-of.j2) ) )/2;

   int Jab_min = max(abs(oa.j2-ob.j2),abs(J2-oc.j2))/2;
   int Jab_max = min(oa.j2+ob.j2,J2+oc.j2)/2;
   int Jde_min = max(abs(od.j2-oe.j2),abs(J2-of.j2))/2;
   int Jde_max = min(od.j2+oe.j2,J2+of.j2)/2;

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;

   orbindx_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2,'j');
        double Cj_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2,'j');

        int J2index = (a>=d) ? (Jab_max-Jab_min+1)*(Jde-Jde_min) + (Jab-Jab_min)
                             : (Jde_max-Jde_min+1)*(Jab-Jab_min) + (Jde-Jde_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = RecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2,'t');
            double Ct_def = RecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2,'t');

            int J2Tindex = (a>=d) ? J2index*5 + 2*tab + tde + (T2-1)/2
                                  : J2index*5 + 2*tde + tab + (T2-1)/2;
            
            ThreeBody.at(orbit_index).at(Jindex).at(J2Tindex) += Cj_abc * Cj_def * Ct_abc * Ct_def * V;

        }
      }
    }
  }

}


// Hashing function for compressing 6 orbit indices to one number
orbindx_t Operator3::GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f)
{
   unsigned char aa = a/2;
   unsigned char bb = b/2;
   unsigned char cc = c/2;
   unsigned char dd = d/2;
   unsigned char ee = e/2;
   unsigned char ff = f/2;

   orbindx_t orbit_index_bra = (aa << 16) + (bb << 8) + cc;
   orbindx_t orbit_index_ket = (dd << 16) + (ee << 8) + ff;

   // should change to uint64_t
   orbindx_t orb_indx = (aa>=dd) ? (orbit_index_bra << 32) + orbit_index_ket
                                : (orbit_index_ket << 32) + orbit_index_bra;

/*
   int aa = a-a%2;
   int bb = b-b%2;
   int cc = c-c%2;
   int dd = d-d%2;
   int ee = e-e%2;
   int ff = f-f%2;
   int orbit_index_bra =  aa*(aa+1)*(aa+2)/6 + bb*(bb+1)/2 + cc;
   int orbit_index_ket =  dd*(dd+1)*(dd+2)/6 + ee*(ee+1)/2 + ff;
   long int orb_indx = (aa>=dd) ? (orbit_index_bra * (orbit_index_bra+1))/2 + orbit_index_ket
                                : (orbit_index_ket * (orbit_index_ket+1))/2 + orbit_index_bra;
*/
   return orb_indx;
}

void Operator3::GetOrbitsFromIndex(orbindx_t orb_indx, int& a, int& b, int& c, int& d, int& e, int& f)
{
  a = ((orb_indx >> 48) & 0xF)*2;
  b = ((orb_indx >> 40) & 0xF)*2;
  c = ((orb_indx >> 32) & 0xF)*2;
  d = ((orb_indx >> 16) & 0xF)*2;
  e = ((orb_indx >> 8)  & 0xF)*2;
  f = ( orb_indx        & 0xF)*2;
  if (a<d)
  {
    swap(a,d);
    swap(b,e);
    swap(c,f);
  }
/*
  int orbit_index_bra = int(sqrt(1+8*orb_indx)-1)/2;
  int orbit_index_ket = orb_indx - orbit_index_bra;
  a = int (pow(orbit_index_bra,1./3.0));
  while (a*(a+1)*(a+2)<=6*orbit_index_bra) ++a;
  --a;
  orbit_index_bra -= a*(a+1)*(a+2)/6;
  b = int(sqrt(1+8*orbit_index_bra)-1)/2;
  c = orbit_index_bra-b*(b+1)/2;
  d = int (pow(orbit_index_ket,1./3.0));
  while (d*(d+1)*(d+2)<=6*orbit_index_ket) ++d;
  --d;
  orbit_index_ket -= d*(d+1)*(d+2)/6;
  e = int(sqrt(1+8*orbit_index_ket)-1)/2;
  f = orbit_index_bra-e*(e+1)/2;
*/
}


// Rearrange (abc) so that a>=b>=c
void Operator3::SortOrbits(int& a, int& b, int& c)
{
   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);
}

// Calculate the angular momentum or isospin recoupling coefficients needed to rearrange (abc)
double Operator3::RecouplingCoefficient(int a_in, int b_in, int c_in, int a, int b, int c, int Jab_in, int Jab, int J, char j_or_t)
{
   double C;
   double ja=0.5, jb=0.5, jc=0.5;
   if (j_or_t == 'j')
   {
      Orbit& oa = modelspace->GetOrbit(a_in);
      Orbit& ob = modelspace->GetOrbit(b_in);
      Orbit& oc = modelspace->GetOrbit(c_in);
      ja = oa.j2/2.0;
      jb = ob.j2/2.0;
      jc = oc.j2/2.0;
   }
   if (a==a_in)
   {
      if (b==b_in) // (abc)_in -> (abc)
      {
         C = Jab==Jab_in ? 1 : 0;
      }
      else  // (acb)_in -> (abc)
      {
         C = modelspace->phase( (jb+jc)+Jab-Jab_in) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab_in, jc, J/2., Jab);
      }
   }
   else if (a==b_in)
   {
      if (b==a_in) // (bac)_in -> (abc)
      {
         C = Jab==Jab_in ? modelspace->phase((ja+jb)/2-Jab) : 0;
      }
      else // (bca)_in -> (abc)
      {
         C = modelspace->phase( (jb+jc)+Jab+1 ) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab_in, jc, J/2., Jab);
      }
   }
   else
   {
      if (b==a_in) // (cab)_in -> (abc)
      {
         C = modelspace->phase( (ja+jb)-Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab_in, jc, J/2., Jab);
      }
      else // (cba)_in -> (abc)
      {
         C =  - sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab_in, jc, J/2., Jab);
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
   Operator opNO3 = Operator(*modelspace);
//   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      for (int ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         for (int iket=ibra; iket<tbc.GetNumberKets(); ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            for (int& a : modelspace->holes)
            {
               Orbit & oa = modelspace->GetOrbit(a);
               if ( (2*(oi.n+oj.n+oa.n)+oi.l+oj.l+oa.l)>E3max) continue;
               if ( (2*(ok.n+ol.n+oa.n)+ok.l+ol.l+oa.l)>E3max) continue;
               int Tz2 = 2*tbc.Tz + oa.tz2;
               int kmin2 = abs(2*tbc.J-oa.j2);
               int kmax2 = 2*tbc.J+oa.j2;
               for (int K2=kmin2; K2<=kmax2; ++K2)
               {
//                  #pragma omp critical
                  opNO3.TwoBody[ch](ibra,iket) += (K2+1) * GetThreeBodyME_pn(tbc.J,tbc.J,K2,Tz2,i,j,a,k,l,a);
               }
            }
            opNO3.TwoBody[ch](ibra,iket) /= (2*tbc.J+1);
         }
      }
   }
   Operator opNO2 = opNO3.DoNormalOrdering();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);

   // Also normal order the 1 and 2 body piece of the 3-body operator (if any)
   opNO2 += DoNormalOrdering();
   return opNO2;

}


