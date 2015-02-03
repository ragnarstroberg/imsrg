
#include "Operator.hh"
#include "AngMom.hh"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

double  Operator::bch_transform_threshold = 1e-6;
double  Operator::bch_product_threshold = 1e-4;

/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator() :
   modelspace(NULL), nChannels(0), hermitian(true), antihermitian(false),
    rank_J(0), rank_T(0), parity(0), particle_rank(2)
{}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    hermitian(true), antihermitian(false), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(ms.GetNumberTwoBodyChannels(), map<int,arma::mat>() ),
    TwoBodyTensorChannels(ms.GetNumberTwoBodyChannels(), vector<int>() ),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank)
{
  modelspace = &ms;
  AllocateTwoBody();
}



Operator::Operator(ModelSpace& ms) :
    hermitian(true), antihermitian(false), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(ms.GetNumberTwoBodyChannels(), map<int,arma::mat>() ),
    TwoBodyTensorChannels(ms.GetNumberTwoBodyChannels(), vector<int>() ),
    rank_J(0), rank_T(0), parity(0), particle_rank(2)
{
  modelspace = &ms;
  AllocateTwoBody();
}

Operator::Operator(const Operator& op)
{
   Copy(op);
}



void Operator::AllocateTwoBody()
{
  for (int ch_bra=0; ch_bra<nChannels;++ch_bra)
  {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
     for (int ch_ket=ch_bra; ch_ket<nChannels;++ch_ket)
     {
        TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
        if ( abs(tbc_bra.J-tbc_ket.J)>rank_J or (tbc_bra.J+tbc_ket.J)<rank_J
          or abs(tbc_bra.Tz-tbc_ket.Tz)>rank_T
           or (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;
        
        int n_bra = tbc_bra.GetNumberKets();
        TwoBodyTensorChannels[ch_bra].push_back(ch_ket);
        TwoBody[ch_bra][ch_ket] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
        
     }
  }
}


// Confusing nomenclature: J2 means 2 times the total J of the three body system
void Operator::AllocateThreeBody()
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

        orbindx3_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

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
}






/////////// COPY METHOD //////////////////////////
void Operator::Copy(const Operator& op)
{
   modelspace    = op.modelspace;
   nChannels     = op.nChannels;
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
   rank_J        = op.rank_J;
   rank_T        = op.rank_T;
   parity        = op.parity;
   particle_rank = op.particle_rank;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   TwoBody       = op.TwoBody;
   TwoBodyTensorChannels       = op.TwoBodyTensorChannels;
}

/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////
Operator& Operator::operator=(const Operator& rhs)
{
   Copy(rhs);
   return *this;
}

// multiply operator by a scalar
Operator& Operator::operator*=(const double rhs)
{
   ZeroBody *= rhs;
   OneBody *= rhs;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      for ( auto twobody : TwoBody[ch] )
      {
        arma::mat& matrix = twobody.second;
        matrix *= rhs;
      }
   }
   return *this;
}

Operator Operator::operator*(const double rhs) const
{
   Operator opout = Operator(*this);
   opout *= rhs;
   return opout;
}

// Add non-member operator so we can multiply an operator
// by a scalar from the lhs, i.e. s*O = O*s
Operator operator*(const double lhs, const Operator& rhs)
{
   return rhs * lhs;
}


// divide operator by a scalar
Operator& Operator::operator/=(const double rhs)
{
   return *this *=(1.0/rhs);
}

Operator Operator::operator/(const double rhs) const
{
   Operator opout = Operator(*this);
   opout *= (1.0/rhs);
   return opout;
}

// Add operators
Operator& Operator::operator+=(const Operator& rhs)
{
   ZeroBody += rhs.ZeroBody;
   OneBody += rhs.OneBody;
   for (int ch_bra=0; ch_bra<modelspace->GetNumberTwoBodyChannels();++ch_bra)
   {
      for (int& ch_ket : TwoBodyTensorChannels[ch_bra] )
      {
         TwoBody.at(ch_bra).at(ch_ket) += rhs.TwoBody.at(ch_bra).at(ch_ket);
      }
   }
   return *this;
}

Operator Operator::operator+(const Operator& rhs) const
{
   return ( Operator(*this) += rhs );
}

// Subtract operators
Operator& Operator::operator-=(const Operator& rhs)
{
   ZeroBody -= rhs.ZeroBody;
   OneBody -= rhs.OneBody;
   for (int ch_bra=0; ch_bra<modelspace->GetNumberTwoBodyChannels();++ch_bra)
   {
      for (int& ch_ket : TwoBodyTensorChannels[ch_bra] )
      {
         TwoBody.at(ch_bra).at(ch_ket) -= rhs.TwoBody.at(ch_bra).at(ch_ket);
      }
   }
   return *this;
}

Operator Operator::operator-(const Operator& rhs) const
{
   return ( Operator(*this) -= rhs );
}



///////// SETTER_GETTERS ///////////////////////////

//double Operator::GetTBME(int ch, int a, int b, int c, int d) const
double Operator::GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc_bra.GetNumberKets() or ket_ind > tbc_ket.GetNumberKets() )
     return 0;
   Ket & bra = tbc_bra.GetKet(bra_ind);
   Ket & ket = tbc_ket.GetKet(ket_ind);

   double phase = 1;
   if (a>b) phase *= bra.Phase(tbc_bra.J);
   if (c>d) phase *= ket.Phase(tbc_ket.J);
   if (a==b) phase *= sqrt(2.);
   if (c==d) phase *= sqrt(2.);
//   return phase * TwoBody[ch](bra_ind, ket_ind);
   return phase * TwoBody.at(ch_bra).at(ch_ket)(bra_ind, ket_ind);
}


void Operator::SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(min(c,d),max(c,d));
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   TwoBody.at(ch_bra).at(ch_ket)(bra_ind,ket_ind) = phase * tbme;
   if (hermitian) TwoBody.at(ch_bra).at(ch_ket)(ket_ind,bra_ind) = phase * tbme;
   if (antihermitian) TwoBody.at(ch_bra).at(ch_ket)(ket_ind,bra_ind) = - phase * tbme;
}


void Operator::AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(min(c,d),max(c,d));
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   TwoBody.at(ch_bra).at(ch_ket)(bra_ind,ket_ind) += phase * tbme;
   if (hermitian) TwoBody.at(ch_bra).at(ch_ket)(ket_ind,bra_ind) += phase * tbme;
   if (antihermitian) TwoBody.at(ch_bra).at(ch_ket)(ket_ind,bra_ind) -=  phase * tbme;
}

double Operator::GetTBME(int ch_bra, int ch_ket, Ket &bra, Ket &ket) const
{
   return GetTBME(ch_bra,ch_ket,bra.p,bra.q,ket.p,ket.q);
}
void Operator::SetTBME(int ch_bra, int ch_ket, Ket& bra, Ket& ket, double tbme)
{
   SetTBME(ch_bra, ch_ket, bra.p,bra.q,ket.p,ket.q,tbme);
}
void Operator::AddToTBME(int ch_bra, int ch_ket, Ket& bra, Ket& ket, double tbme)
{
   AddToTBME(ch_bra, ch_ket, bra.p,bra.q,ket.p,ket.q, tbme);
}

double Operator::GetTBME(int ch_bra, int ch_ket, int ibra, int iket) const
{
   // Use at() rather than [] because invalid access returns a useful error rather than a segfault
   return TwoBody.at(ch_bra).at(ch_ket)(ibra,iket);
}
void Operator::SetTBME(int ch_bra, int ch_ket, int iket, int ibra, double tbme)
{
   TwoBody.at(ch_bra).at(ch_ket)(ibra,iket) = tbme;
}
void Operator::AddToTBME(int ch_bra, int ch_ket, int iket, int ibra, double tbme)
{
   TwoBody.at(ch_bra).at(ch_ket)(ibra,iket) += tbme;
}

double Operator::GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d) const
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   return GetTBME(ch_bra,ch_ket,a,b,c,d);
}
void Operator::SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   SetTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}
void Operator::AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   AddToTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}
double Operator::GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket) const
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   return GetTBME(ch_bra,ch_ket,bra,ket);
}
void Operator::SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   SetTBME(ch_bra,ch_ket,bra,ket,tbme);
}
void Operator::AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   AddToTBME(ch_bra,ch_ket,bra,ket,tbme);
}



// for backwards compatibility...
double Operator::GetTBME(int ch, int a, int b, int c, int d) const
{
   return GetTBME(ch,ch,a,b,c,d);
}

void Operator::SetTBME(int ch, int a, int b, int c, int d, double tbme)
{
   return SetTBME(ch,ch,a,b,c,d,tbme);
}

double Operator::GetTBME(int ch, Ket &bra, Ket &ket) const
{
   return GetTBME(ch,ch,bra.p,bra.q,ket.p,ket.q);
}

void Operator::SetTBME(int ch, Ket& bra, Ket& ket, double tbme)
{
   SetTBME(ch,ch,bra.p,bra.q,ket.p,ket.q,tbme);
}

void Operator::AddToTBME(int ch, Ket& bra, Ket& ket, double tbme)
{
   AddToTBME(ch,ch,bra.p,bra.q,ket.p,ket.q,tbme);
}


double Operator::GetTBME(int ch, int ibra, int iket) const
{
   return GetTBME(ch,ch,ibra,iket);
}

void Operator::SetTBME(int ch, int ibra, int iket, double tbme)
{
   SetTBME(ch,ch,ibra,iket,tbme);
}

void Operator::AddToTBME(int ch, int ibra, int iket, double tbme)
{
   AddToTBME(ch,ch,ibra,iket,tbme);
}



double Operator::GetTBME(int j, int p, int t, int a, int b, int c, int d) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,ch,a,b,c,d);
}

void Operator::SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,ch,a,b,c,d,tbme);
}

void Operator::AddToTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   AddToTBME(ch,ch,a,b,c,d,tbme);
}

double Operator::GetTBME(int j, int p, int t, Ket& bra, Ket& ket) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,ch,bra,ket);
}

void Operator::SetTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,ch,bra,ket,tbme);
}

void Operator::AddToTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   AddToTBME(ch,ch,bra,ket,tbme);
}


double Operator::GetTBMEmonopole(int a, int b, int c, int d) const
{
   double mon = 0;
   Orbit &oa = modelspace->GetOrbit(a);
   Orbit &ob = modelspace->GetOrbit(b);
   Orbit &oc = modelspace->GetOrbit(c);
   Orbit &od = modelspace->GetOrbit(d);
   int Tzab = (oa.tz2 + ob.tz2)/2;
   int parityab = (oa.l + ob.l)%2;
   int Tzcd = (oc.tz2 + od.tz2)/2;
   int paritycd = (oc.l + od.l)%2;

   if (Tzab != Tzcd or parityab != paritycd) return 0;

   int jmin = abs(oa.j2 - ob.j2)/2;
   int jmax = (oa.j2 + ob.j2)/2;
   
   for (int J=jmin;J<=jmax;++J)
   {
      mon += (2*J+1) * GetTBME(J,parityab,Tzab,a,b,c,d);
   }
   mon /= (oa.j2 +1)*(ob.j2+1);
   return mon;
}


double Operator::GetTBMEmonopole(Ket & bra, Ket & ket) const
{
   return GetTBMEmonopole(bra.p,bra.q,ket.p,ket.q);
}





// Get Three Body Matrix Element in proton-neutron formalism.
double Operator::GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int Tz, int a, int b, int c, int d, int e, int f)
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
double Operator::GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;

   SortThreeBodyOrbits(a,b,c);
   SortThreeBodyOrbits(d,e,f);

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

   orbindx3_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

   double V = 0;
   // Recouple J and T to get to the format in which it's stored.
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = ThreeBodyRecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2,'j');
        double Cj_def = ThreeBodyRecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2,'j');
        if (Cj_abc*Cj_def == 0) continue;

        int J2index = (a>=d) ? (Jab_max-Jab_min+1)*(Jde-Jde_min) + (Jab-Jab_min)
                             : (Jde_max-Jde_min+1)*(Jab-Jab_min) + (Jde-Jde_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = ThreeBodyRecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2,'t');
            double Ct_def = ThreeBodyRecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2,'t');

            int J2Tindex = (a>=d) ? J2index*5 + 2*tab + tde + (T2-1)/2
                                  : J2index*5 + 2*tde + tab + (T2-1)/2;

            V += Cj_abc * Cj_def * Ct_abc * Ct_def * ThreeBody.at(orbit_index).at(Jindex).at(J2Tindex);
         }
        }
      }
   }
   return V;

}




void Operator::SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V)
{
   // reorder so a>=b>=c and d>=e>=f
   int a=a_in;
   int b=b_in;
   int c=c_in;
   int d=d_in;
   int e=e_in;
   int f=f_in;
   SortThreeBodyOrbits(a,b,c);
   SortThreeBodyOrbits(d,e,f);

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

   orbindx3_t orbit_index = GetThreeBodyOrbitIndex(a,b,c,d,e,f);

   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
        double Cj_abc = ThreeBodyRecouplingCoefficient(a_in,b_in,c_in,a,b,c,Jab_in,Jab,J2,'j');
        double Cj_def = ThreeBodyRecouplingCoefficient(d_in,e_in,f_in,d,e,f,Jde_in,Jde,J2,'j');

        int J2index = (a>=d) ? (Jab_max-Jab_min+1)*(Jde-Jde_min) + (Jab-Jab_min)
                             : (Jde_max-Jde_min+1)*(Jab-Jab_min) + (Jde-Jde_min);

        for (int tab=tab_min; tab<=tab_max; ++tab)
        {
          for (int tde=tde_min; tde<=tde_max; ++tde)
          {
            double Ct_abc = ThreeBodyRecouplingCoefficient(a_in,b_in,c_in,a,b,c,tab_in,tab,T2,'t');
            double Ct_def = ThreeBodyRecouplingCoefficient(d_in,e_in,f_in,d,e,f,tde_in,tde,T2,'t');

            int J2Tindex = (a>=d) ? J2index*5 + 2*tab + tde + (T2-1)/2
                                  : J2index*5 + 2*tde + tab + (T2-1)/2;
            
            ThreeBody.at(orbit_index).at(Jindex).at(J2Tindex) += Cj_abc * Cj_def * Ct_abc * Ct_def * V;

        }
      }
    }
  }

}



// Hashing function for compressing 6 orbit indices to one number
orbindx3_t Operator::GetThreeBodyOrbitIndex(int a, int b, int c, int d, int e, int f)
{
   unsigned char aa = a/2;
   unsigned char bb = b/2;
   unsigned char cc = c/2;
   unsigned char dd = d/2;
   unsigned char ee = e/2;
   unsigned char ff = f/2;

   orbindx3_t orbit_index_bra = (aa << 16) + (bb << 8) + cc;
   orbindx3_t orbit_index_ket = (dd << 16) + (ee << 8) + ff;

   // should change to uint64_t
   orbindx3_t orb_indx = (aa>=dd) ? (orbit_index_bra << 32) + orbit_index_ket
                                : (orbit_index_ket << 32) + orbit_index_bra;

   return orb_indx;
}

void Operator::GetOrbitsFromThreeBodyIndex(orbindx3_t orb_indx, int& a, int& b, int& c, int& d, int& e, int& f)
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
}

// Rearrange (abc) so that a>=b>=c
void Operator::SortThreeBodyOrbits(int& a, int& b, int& c)
{
   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);
}



// Calculate the angular momentum or isospin recoupling coefficients needed to rearrange (abc)
double Operator::ThreeBodyRecouplingCoefficient(int a_in, int b_in, int c_in, int a, int b, int c, int Jab_in, int Jab, int J, char j_or_t)
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





////////////////// MAIN INTERFACE METHODS //////////////////////////

Operator Operator::DoNormalOrdering()
{
   if (particle_rank==3)
   {
      return DoNormalOrdering3();
   }
   else
   {
      return DoNormalOrdering2();
   }
}

//
//  Normal ordering of a 2body operator
//  currently this only handles scalar operators
//
Operator Operator::DoNormalOrdering2()
{
   Operator opNO = *this;

   // Trivial parts
   opNO.ZeroBody = ZeroBody;
   opNO.OneBody = OneBody;
   opNO.TwoBody = TwoBody;



   for (int &k : modelspace->holes) // loop over hole orbits
   {
      opNO.ZeroBody += (modelspace->GetOrbit(k).j2+1) * OneBody(k,k);
   }



   int norbits = modelspace->GetNumberOrbits();

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;

      // Zero body part
      for (int& iket : tbc.KetIndex_hh) // loop over hole-hole kets in this channel
      {
//        opNO.ZeroBody += TwoBody[ch](iket,iket) * (2*J+1);  // <ab|V|ab>  (a,b in core)
        opNO.ZeroBody += GetTBME(ch,iket,iket) * (2*J+1);  // <ab|V|ab>  (a,b in core)
      }

      // One body part
      int ibra,iket;
      for (int a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         int bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for (int b=bstart;b<norbits;++b)
         {
            Orbit &ob = modelspace->GetOrbit(b);
            if (ob.j2 != oa.j2 or ob.tz2 != oa.tz2 or ob.l != oa.l) continue;
            for (int &h : modelspace->holes)  // C++11 syntax
            {
               if ( (ibra = tbc.GetLocalIndex(min(a,h),max(a,h)))<0) continue;
               if ( (iket = tbc.GetLocalIndex(min(b,h),max(b,h)))<0) continue;
               Ket & bra = tbc.GetKet(ibra);
               Ket & ket = tbc.GetKet(iket);
               double tbme = (2*J+1.0)/(oa.j2+1)  * GetTBME(ch,bra,ket);
               int phase = 1;
               if (a>h) phase *= bra.Phase(J);
               if (b>h) phase *= ket.Phase(J);

               opNO.OneBody(a,b) += phase * tbme;  // <ah|V|bh>
            }
            if (hermitian) opNO.OneBody(b,a) = opNO.OneBody(a,b);
            if (antihermitian) opNO.OneBody(b,a) = -opNO.OneBody(a,b);
         }
      }


   } // loop over channels
   return opNO;
}




//   Normal ordering of a 3body operator. Start by generating the normal ordered
//   2body piece, then use DoNormalOrdering2() to get the rest. (Note that there
//   are some numerical factors).
//   The normal ordered two body piece is 
//   Gamma(2)^J_ijkl = V(2)^J_ijkl + Sum_a n_a  Sum_K (2K+1)/(2J+1) V(3)^JJK_ijakla
//   Right now, this is only set up for scalar operators, but I don't anticipate
//   handling 3body tensor operators in the near future.
//
Operator Operator::DoNormalOrdering3()
{
   Operator opNO3 = Operator(*modelspace);
//   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& Gamma = (arma::mat&) opNO3.TwoBody[ch].at(ch);
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
//                  opNO3.TwoBody[ch](ibra,iket) += (K2+1) * GetThreeBodyME_pn(tbc.J,tbc.J,K2,Tz2,i,j,a,k,l,a);
                  Gamma(ibra,iket) += (K2+1) * GetThreeBodyME_pn(tbc.J,tbc.J,K2,Tz2,i,j,a,k,l,a);
               }
            }
//            opNO3.TwoBody[ch](ibra,iket) /= (2*tbc.J+1);
            Gamma(ibra,iket) /= (2*tbc.J+1);
         }
      }
   }
   Operator opNO2 = opNO3.DoNormalOrdering2();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);

   // Also normal order the 1 and 2 body piece of the 3-body operator (if any)
   opNO2 += DoNormalOrdering2();
   return opNO2;

}








void Operator::EraseTwoBody()
{
   for (int ch=0;ch<nChannels;++ch)
   {
      for (auto twobody : TwoBody[ch] )
      {
        arma::mat& matrix = twobody.second;
        matrix.zeros();
      }
   }
}

void Operator::ScaleZeroBody(double x)
{
   ZeroBody *= x;
}

void Operator::ScaleOneBody(double x)
{
   OneBody *= x;
}

void Operator::ScaleTwoBody(double x)
{
   for (int ch=0; ch<nChannels; ++ch)
   {
     for ( auto twobody : TwoBody[ch] )
     {
      arma::mat& matrix = twobody.second;
      matrix *= x;
     }
   }
}

void Operator::Eye()
{
   ZeroBody = 1;
   OneBody.eye();
   for (int ch=0; ch<nChannels; ++ch)
   {
     for ( auto twobody : TwoBody[ch] )
     {
        arma::mat& matrix = twobody.second;
        matrix.eye();
     }
   }
}


void Operator::CalculateKineticEnergy()
{
   OneBody.zeros();
   int norbits = modelspace->GetNumberOrbits();
   int A = modelspace->GetTargetMass();
   float hw = modelspace->GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace->GetOrbit(a);
      OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for (int b=a+1;b<norbits;++b)  // make this better once OneBodyChannel is implemented
      {
         Orbit & ob = modelspace->GetOrbit(b);
         if (oa.l == ob.l and oa.j2 == ob.j2 and oa.tz2 == ob.tz2)
         {
            if (oa.n == ob.n+1)
               OneBody(a,b) = 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
            else if (oa.n == ob.n-1)
               OneBody(a,b) = 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
            OneBody(b,a) = OneBody(a,b);
         }
      }
   }
}


//*****************************************************************************************
//  I think this explanation needs to be revisited...
//   Transform to a cross-coupled basis, for use in the 2body commutator relation
//    using a Pandya tranformation 
//                                               ____________________________________________________________________________
//                                              |                                                                            |
//   |a    |b  <ab|_J'    \a   /c  <ac|_J       |  <ac|V|bd>_J = sum_J' (2J'+1) (-1)^(ja+jd+J+J'+1) { ja jc J } <ab|V|cd>_J' |
//   |     |               \  /                 |                                                   { jd jb J'}              |
//   |_____|                \/_____             |____________________________________________________________________________|
//   |  V  |      ===>         V  /\ 
//   |     |                     /  \
//   |c    |d  |cd>_J'          /b   \d   |bd>_J
//
//   STANDARD                 CROSS  
//   COUPLING                 COUPLED  
//                                      
//void Operator::CalculateCrossCoupled()
void Operator::CalculateCrossCoupled(vector<arma::mat> &TwoBody_CC_left, vector<arma::mat> &TwoBody_CC_right) const
{
   // loop over cross-coupled channels
   #pragma omp parallel for
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.KetIndex_ph.size();
      int J_cc = tbc_cc.J;

//   These matrices don't actually need to be square, since we only care about
//   the particle-hole kets, which we sum over via matrix multiplication:
//   [ ]                  [     ]
//   |L|  x  [  R   ]  =  |  N  |
//   [ ]                  [     ]
//
      TwoBody_CC_left[ch_cc]  = arma::mat(2*nKets_cc, nph_kets,   arma::fill::zeros);
      TwoBody_CC_right[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);

      // loop over cross-coupled ph bras <ac| in this channel
      for (int i_ph=0; i_ph<nph_kets; ++i_ph)
      {
         Ket & bra_cc = tbc_cc.GetKet( tbc_cc.KetIndex_ph[i_ph] );
         int a = bra_cc.p;
         int c = bra_cc.q;
         Orbit & oa = modelspace->GetOrbit(a);
         Orbit & oc = modelspace->GetOrbit(c);
         double ja = oa.j2/2.0;
         double jc = oc.j2/2.0;

         // loop over cross-coupled kets |bd> in this channel
         // we go to 2*nKets to include |bd> and |db>
         for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
         {
            Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
            int b = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
            int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
            Orbit & ob = modelspace->GetOrbit(b);
            Orbit & od = modelspace->GetOrbit(d);
            double jb = ob.j2/2.0;
            double jd = od.j2/2.0;

            int phase_ad = modelspace->phase(ja+jd);

            // Get Tz,parity and range of J for <ab || cd > coupling
            int Tz_std = (oa.tz2 + ob.tz2)/2;
            int parity_std = (oa.l + ob.l)%2;
            int jmin = max(abs(ja-jb),abs(jc-jd));
            int jmax = min(ja+jb,jc+jd);
            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jc,J_cc,jd,jb,J_std);
               if (sixj == 0) continue;
               int phase = modelspace->phase(J_std);
               double tbme = GetTBME(J_std,parity_std,Tz_std,a,b,c,d);
               sm += (2*J_std+1) * phase * sixj * tbme ; 
            }
            TwoBody_CC_left[ch_cc](iket_cc,i_ph) = sm * phase_ad;


            // Get Tz,parity and range of J for <cb || ad > coupling
            Tz_std = (oa.tz2 + od.tz2)/2;
            parity_std = (oa.l + od.l)%2;
            jmin = max(abs(jc-jb),abs(ja-jd));
            jmax = min(jc+jb,ja+jd);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jc,J_cc,jb,jd,J_std);
               if (sixj == 0) continue;
               int phase = modelspace->phase(J_std);
               double tbme = GetTBME(J_std,parity_std,Tz_std,c,b,a,d);
               sm += (2*J_std+1) * phase * sixj * tbme ;
            }
            TwoBody_CC_right[ch_cc](i_ph,iket_cc) = - sm * phase_ad;

         }
      }
   }
}


//*****************************************************************************************
//  OpOut = exp(Omega) Op exp(-Omega)
Operator Operator::BCH_Transform( const Operator &Omega) const
{
//   double bch_transform_threshold = 1e-6;
   int max_iter = 20;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   Operator OpNested = *this;
   for (int i=1; i<max_iter; ++i)
   {
      OpNested = Omega.Commutator(OpNested) / i;
      
      OpOut += OpNested;

      if (OpNested.Norm() < (nx+ny)*bch_transform_threshold) return OpOut;
      if (i == warn_iter)
      {
         cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << endl;
      }

   }
   cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << endl;
   return OpOut;
}


//*****************************************************************************************
// Baker-Campbell-Hausdorff formula
//  returns Z, where
//  exp(Z) = exp(X) * exp(Y).
//  Z = X + Y + 1/2[X, Y]
//     + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]]
//     - 1/24 [Y,[X,[X,Y]]]
//     - 1/720 [Y,[Y,[Y,[Y,X]]]] - 1/720 [X,[X,[X,[X,Y]]]]
//     + ...
Operator Operator::BCH_Product( const Operator &Y) const
{
   const Operator& X = *this;
//   double bch_product_threshold = 1e-4;

   Operator Z = X + Y; 

   Operator XY = X.Commutator(Y);
   Z += XY*(1./2);    // [X,Y]
   double nx = X.Norm();
   double ny = Y.Norm();
   double nc1 = XY.Norm();

   if ( nc1/2 < (nx+ny)*bch_product_threshold ) return Z;

   Operator YYX = XY.Commutator(Y); // [[X,Y],Y] = [Y,[Y,X]]
   double nc2 = YYX.Norm();
   Z += YYX * (1/12.);      // [Y,[Y,X]]

   if ( nc2/12 < (nx+ny)*bch_product_threshold ) return Z;

   Operator XXY = X.Commutator(XY); // [X,[X,Y]]
   double nc3 = XXY.Norm();
   Z += XXY * (1/12.);      // [X,[X,Y]]

   if ( nc3/12 < (nx+ny)*bch_product_threshold ) return Z;

   cout << "Warning: BCH product expansion not converged after 3 nested commutators!" << endl;

   Operator YXXY = Y.Commutator(XXY); // [Y,[X,[X,Y]]]
   double nc4 = YXXY.Norm();
   Z += YXXY*(-1./24);

   Operator YYYX = Y.Commutator(YYX) ;
   Operator YYYYX = Y.Commutator(YYYX) ;
   double nc5 = YYYYX.Norm();
   Operator XXXY =  X.Commutator(XXY) ; // [X,[X,[X,[X,Y]]]]
   Operator XXXXY = X.Commutator(XXXY) ; // [X,[X,[X,[X,Y]]]]
   double nc6 = XXXXY.Norm();
   Z += (YYYYX + XXXXY)*(-1./720);

   if ( nc6/720. < (nx+ny)*bch_product_threshold ) return Z;
   cout << "Warning: BCH product expansion not converged after 5 nested commutators!" << endl;

   return Z;
}

// Frobenius norm of the operator
double Operator::Norm() const
{
   double nrm = 0;
   double n1 = OneBodyNorm();
   double n2 = TwoBodyNorm();
   return sqrt(n1*n1+n2*n2);
}

double Operator::OneBodyNorm() const
{
   return arma::norm(OneBody,"fro");
}

double Operator::TwoBodyNorm() const
{
   double nrm = 0;
   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      for (auto twobody : TwoBody[ch] )
      {
        arma::mat & matrix = twobody.second;
        double n2 = arma::norm(matrix,"fro");
        nrm += n2*n2;
      }
   }
   return sqrt(nrm);
}


Operator Operator::Commutator(const Operator& opright) const
{
   if (rank_J==0)
   {
      if (opright.rank_J==0)
      {
         return CommutatorScalarScalar(opright); // [S,S]
      }
      else
      {
         return CommutatorScalarTensor(opright); // [S,T]
      }
   }
   else if(opright.rank_J==0)
   {
      return (-1)*opright.CommutatorScalarTensor(*this); // [T,S]
   }
   else
   {
      cout << "In Tensor-Tensor because rank_J = " << rank_J << "  and opright.rank_J = " << opright.rank_J << endl;
      return CommutatorTensorTensor(opright); // [T,T]
   }
}


Operator Operator::CommutatorScalarScalar(const Operator& opright) const
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   if ( not out.IsAntiHermitian() )
   {
      comm110ss(opright, out);
      comm220ss(opright, out) ;
   }

   comm111ss(opright, out);
   comm121ss(opright, out);

   comm122ss(opright, out);
   comm222_pp_hh_221ss(opright, out);
   comm222_phss(opright, out);

   return out;
}


// Calculate [S,T]
Operator Operator::CommutatorScalarTensor(const Operator& opright) const
{
   Operator out = opright; // This ensures the commutator has the same tensor rank as opright
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   comm111st(opright, out);
   comm121st(opright, out);

   comm122st(opright, out);
   comm222_pp_hh_221st(opright, out);
   comm222_phst(opright, out);

   return out;
}

// This should really return a vector of operators, but I'm not going to write this until I need to.
Operator Operator::CommutatorTensorTensor(const Operator& opright) const
{
   Operator out = opright;
/*
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   comm111tt(opright, out);
   comm121tt(opright, out);

   comm122tt(opright, out);
   comm222_pp_hh_221tt(opright, out);
   comm222_phtt(opright, out);
*/
   cout << "Tensor-Tensor commutator not yet implemented" << endl;
   return out;
}



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
///////////////////////////////////////////////////////////////////////////////////////////

//*****************************************************************************************
//                ____Y    __         ____X
//          X ___(_)             Y___(_) 
//
//  [X1,Y1](0) = Sum_ab (2j_a+1) x_ab y_ba  (n_a-n_b) 
//             = Sum_a  (2j_a+1)  (xy-yx)_aa n_a
//
// -- AGREES WITH NATHAN'S RESULTS
void Operator::comm110ss(const Operator& opright, Operator& out) const
{
  if (IsHermitian() and opright.IsHermitian()) return ; // I think this is the case
  if (IsAntiHermitian() and opright.IsAntiHermitian()) return ; // I think this is the case

   int norbits = modelspace->GetNumberOrbits();
   arma::mat xyyx = OneBody*opright.OneBody - opright.OneBody*OneBody;
   for ( int& a : modelspace->holes) // C++11 range-for syntax: loop over all elements of the vector <particle>
   {
      out.ZeroBody += (modelspace->GetOrbit(a).j2+1) * xyyx(a,a);
   }
}


//*****************************************************************************************
//         __Y__       __X__
//        ()_ _()  -  ()_ _()
//           X           Y
//
//  [ X^(2), Y^(2) ]^(0) = 1/2 Sum_abcd  Sum_J (2J+1) x_abcd y_cdab (n_a n_b nbar_c nbar_d)
//                       = 1/2 Sum_J (2J+1) Sum_abcd x_abcd y_cdab (n_a n_b nbar_c nbar_d)  
//                       = 1/2 Sum_J (2J+1) Sum_ab  (X*P_pp*Y)_abab  P_hh
//
//  -- AGREES WITH NATHAN'S RESULTS (within < 1%)
void Operator::comm220ss( const Operator& opright, Operator& out) const
{
   if (IsHermitian() and opright.IsHermitian()) return; // I think this is the case
   if (IsAntiHermitian() and opright.IsAntiHermitian()) return; // I think this is the case
   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      #pragma omp critical
      out.ZeroBody += 2 * (2*tbc.J+1) * arma::trace( tbc.Proj_hh * TwoBody.at(ch).at(ch) * tbc.Proj_pp * opright.TwoBody.at(ch).at(ch) );

   }
}

//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// -- AGREES WITH NATHAN'S RESULTS
void Operator::comm111ss(const Operator & opright, Operator& out) const
{
   out.OneBody = OneBody*opright.OneBody - opright.OneBody*OneBody;
}

//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
void Operator::comm121ss(const Operator& opright, Operator& out) const
{
   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          Orbit &oj = modelspace->GetOrbit(j);
          if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue; // At some point, make a OneBodyChannel class...
          for (int &a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (int b=0;b<norbits;++b)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                out.OneBody(i,j) += (ob.j2+1) *  OneBody(a,b) * opright.GetTBMEmonopole(b,i,a,j) ;
                out.OneBody(i,j) -= (oa.j2+1) *  OneBody(b,a) * opright.GetTBMEmonopole(a,i,b,j) ;

                // comm211 part
                out.OneBody(i,j) -= (ob.j2+1) *  opright.OneBody(a,b) * GetTBMEmonopole(b,i,a,j) ;
                out.OneBody(i,j) += (oa.j2+1) *  opright.OneBody(b,a) * GetTBMEmonopole(a,i,b,j) ;
             }
          }
      }
   }
}



//*****************************************************************************************
//
//      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
//        |   (  )   _     |   (  )                                                                                      
//        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
//      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
//                                     
//
// -- AGREES WITH NATHAN'S RESULTS 
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm221ss(const Operator& opright, Operator& out) const
{

   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      arma::mat& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
      arma::mat& Matrixpp = (arma::mat&) Mpp.TwoBody.at(ch).at(ch);
      arma::mat& Matrixhh = (arma::mat&) Mpp.TwoBody.at(ch).at(ch);
      
      Matrixpp = ( LHS * tbc.Proj_pp * RHS - RHS * tbc.Proj_pp * LHS );
      Matrixhh = ( LHS * tbc.Proj_hh * RHS - RHS * tbc.Proj_hh * LHS );

      // If commutator is hermitian or antihermitian, we only
      // need to do half the sum. Add this.
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         for (int j=0;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,i,c,j,c);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (int &c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,i,c,j,c);
            }
            #pragma omp critical
            out.OneBody(i,j) += (2*tbc.J+1.0)/(oi.j2 +1.0) * cijJ;
         } // for j
      } // for i
   } //for ch
}





//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
// Right now, this is the slowest one...
void Operator::comm122ss(const Operator& opright, Operator& opout ) const
{
   int herm = opout.IsHermitian() ? 1 : -1;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      int npq = tbc.GetNumberKets();
      int norbits = modelspace->GetNumberOrbits();
      for (int ibra = 0;ibra<npq; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         int indx_ij = ibra;
         double pre_ij = i==j ? sqrt(2) : 1;
         for (int iket=ibra;iket<npq; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            int indx_kl = iket;
            double pre_kl = k==l ? sqrt(2) : 1;

            double cij = 0;
            double ckl = 0;
            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
               int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
               int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
               int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));

               if (indx_aj >= 0)
               {
                  double pre_aj = a>j ? tbc.GetKet(indx_aj).Phase(tbc.J) : 1;
                  if (a==j) pre_aj *= sqrt(2);
                  ckl += pre_aj  * OneBody(i,a) * RHS(indx_aj,indx_kl);
                  ckl -= pre_aj  * opright.OneBody(i,a) * LHS(indx_aj,indx_kl);
               }

               if (indx_ia >= 0)
               {
                  double pre_ia = i>a ? tbc.GetKet(indx_ia).Phase(tbc.J) : 1;
                  if (i==a) pre_ia *= sqrt(2);
                  ckl += pre_ia * OneBody(j,a) * RHS(indx_ia,indx_kl);
                  ckl -= pre_ia * opright.OneBody(j,a) * LHS(indx_ia,indx_kl);
               }

               if (indx_al >= 0)
               {
                  double pre_al = a>l ? tbc.GetKet(indx_al).Phase(tbc.J) : 1;
                  if (a==l) pre_al *= sqrt(2);
                  cij -= pre_al * OneBody(a,k) * RHS(indx_ij,indx_al);
                  cij += pre_al * opright.OneBody(a,k) * LHS(indx_ij,indx_al);
               }

               if (indx_ka >= 0)
               {
                  double pre_ka = k>a ? tbc.GetKet(indx_ka).Phase(tbc.J) : 1;
                  if (k==a) pre_ka *= sqrt(2);
                  cij -= pre_ka * OneBody(a,l) * RHS(indx_ij,indx_ka);
                  cij += pre_ka * opright.OneBody(a,l) * LHS(indx_ij,indx_ka);
               }

            }
            cijkl = (ckl*pre_kl + cij*pre_ij) / sqrt( (1.0+bra.delta_pq())*(1.0+ket.delta_pq()) );
            #pragma omp critical
            {
            OUT(ibra,iket) += cijkl;
            if (ibra != iket)
               OUT(iket,ibra) += herm * cijkl;
            }
         }
      }
   }
}





//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm222_pp_hhss(const Operator& opright, Operator& opout ) const
{
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      arma::mat Mpp = (LHS * tbc.Proj_pp * RHS - RHS * tbc.Proj_pp * LHS);
      arma::mat Mhh = (LHS * tbc.Proj_hh * RHS - RHS * tbc.Proj_hh * LHS);
      OUT += Mpp - Mhh;
   }
}







// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
void Operator::comm222_pp_hh_221ss(const Operator& opright, Operator& opout ) const 
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();

      arma::mat& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      arma::mat& Matrixpp = (arma::mat&) Mpp.TwoBody.at(ch).at(ch);
      arma::mat& Matrixhh = (arma::mat&) Mpp.TwoBody.at(ch).at(ch);
      
      Matrixpp =  LHS * tbc.Proj_pp * RHS;
      Matrixhh =  LHS * tbc.Proj_hh * RHS;

      
      if (opout.IsHermitian())
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else if (opout.IsAntiHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS * tbc.Proj_pp * LHS;
        Matrixhh -=  RHS * tbc.Proj_hh * LHS;
      }


      // The two body part
      OUT += Matrixpp - Matrixhh;
      //opout.TwoBody[ch] += Mpp.TwoBody[ch] - Mhh.TwoBody[ch];

      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         int jmin = opout.IsNonHermitian() ? 0 : i;
         for (int j=jmin;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,i,c,j,c);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (int &c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,i,c,j,c);
            }
            cijJ *= (2*tbc.J+1.0)/(oi.j2+1.0);
            opout.OneBody(i,j) += cijJ;
            if (! opout.IsNonHermitian() and i!=j)
            {
               opout.OneBody(j,i) += herm * cijJ;
            }
         } // for j
      } // for i
   } //for ch
}




//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.             [X2,Y2](2)_ph = -sum_ab (na-nb) sum_J' (-1)^(J'+ja+jb)
//
//                                                                             v----|    v----|                                    v----|    v----|
//   |          |      |          |              * [ (-1)^(ji+jl) { jk jl J } <bi|X|aj> <ak|Y|bl>   - (-1)^(jj+jl-J') { jk jl J } <bj|X|ai> <ak|Y|bl>
//   |     __Y__|      |     __X__|                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |                                            v----|    v----|                                    v----|    v----|
//   |____\/    |      |____\/    |            -  (-1)^(ji+jk-J') { jk jl J } <bi|X|aj> <al|Y|bk>   + (-1)^(jj+jk)    { jk jl J } <bj|X|ai> <al|Y|bk>  ]
//   |  X       |      |  Y       |                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//           
//            
// -- This appears to agree with Nathan's results
//
void Operator::comm222_phss(const Operator& opright, Operator& opout ) const
{

   int herm = opout.IsHermitian() ? 1 : -1;

   // Update Cross-coupled matrix elements
   vector<arma::mat> X_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> X_TwoBody_CC_right (nChannels, arma::mat() );

   vector<arma::mat> Y_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> Y_TwoBody_CC_right (nChannels, arma::mat() );

   CalculateCrossCoupled(X_TwoBody_CC_left, X_TwoBody_CC_right );
   //CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   opright.CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );

   // Construct the intermediate matrices N1 and N2
   vector<arma::mat> N1 (nChannels, arma::mat() );
   vector<arma::mat> N2 (nChannels, arma::mat() );
 // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
 // I think that N1 == N2.t() ??? Check this.
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] = X_TwoBody_CC_left[ch] *  Y_TwoBody_CC_right[ch];
      N2[ch] = Y_TwoBody_CC_left[ch] *  X_TwoBody_CC_right[ch];
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      int nKets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;

         for (int iket=ibra; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;


            double comm = 0;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi.l+ok.l)%2;
            int Tz_cc = abs(oi.tz2+ok.tz2)/2;
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+tbc.J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               double me2 = N2[ch_cc](indx_jl,indx_ik);
               double me3 = N2[ch_cc](indx_ik,indx_jl);
               double me4 = N1[ch_cc](indx_ik,indx_jl);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
            }

            parity_cc = (oi.l+ol.l)%2;
            Tz_cc = abs(oi.tz2+ol.tz2)/2;
            jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_il = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,l),max(i,l));
               int indx_jk = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,k),max(j,k));
               if (i>l)
                 indx_il += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>k)
                 indx_jk += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               double me2 = N2[ch_cc](indx_il,indx_jk);
               double me3 = N2[ch_cc](indx_jk,indx_il);
               double me4 = N1[ch_cc](indx_jk,indx_il);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
            }

            comm /= sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
            #pragma omp critical
            {
               OUT(ibra,iket) += comm;
//              opout.TwoBody[ch](ibra,iket) += comm;
              if (iket != ibra)
              {
                 OUT(iket,ibra) += herm*comm;
//                 opout.TwoBody[ch](iket,ibra) += herm*comm;
              }
            }
         }
      }
   }
}



//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-TENSOR COMMUTATORS      //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// -- AGREES WITH NATHAN'S RESULTS
// This is no different from the scalar-scalar version
void Operator::comm111st(const Operator & opright, Operator& out) const
{
   out.OneBody = OneBody*opright.OneBody - opright.OneBody*OneBody;
}

//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
void Operator::comm121st(const Operator& opright, Operator& out) const
{
   int norbits = modelspace->GetNumberOrbits();
   int nch = modelspace->GetNumberTwoBodyChannels();
   #pragma omp parallel for
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          Orbit &oj = modelspace->GetOrbit(j);
          if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue; // At some point, make a OneBodyChannel class...
          for (int &a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (int b=0;b<norbits;++b)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                for (int ch_bra=0;ch_bra<nch;++ch_bra)
                {
                   TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
                   int phase = modelspace->phase((oj.j2+ob.j2)/2+tbc_bra.J);
                   for (int  ch_ket : opright.TwoBodyTensorChannels[ch_bra])
                   {
                      // part with opleft one body and opright two body, i.e. [X(1),Y(2)](1)
                      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
                      double sixj_a = modelspace->GetSixJ(tbc_ket.J, tbc_bra.J, opright.rank_J, oi.j2/2.0, oj.j2/2.0, oa.j2/2.0);
                      double sixj_b = modelspace->GetSixJ(tbc_ket.J, tbc_bra.J, opright.rank_J, oi.j2/2.0, oj.j2/2.0, ob.j2/2.0);
                      double hatfactor = (2*tbc_bra.J+1)*sqrt((2*tbc_ket.J+1)/(oi.j2+1));

                      out.OneBody(i,j) += OneBody(a,b)*opright.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase*sixj_a*hatfactor;
                      out.OneBody(i,j) -= OneBody(b,a)*opright.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase*sixj_b*hatfactor;
                      
                      // part with opleft two body and opright one body, i.e. [X(2),Y(1)](1)
                      sixj_a = modelspace->GetSixJ(oi.j2/2.0, oa.j2/2.0, tbc_bra.J, ob.j2/2.0, oj.j2/2.0, opright.rank_J);
                      phase = modelspace->phase((oj.j2+ob.j2)/2);
                      hatfactor = (2*tbc_bra.J+1)*sqrt((oa.j2+1)/(oi.j2+1));
                      out.OneBody(i,j) -= opright.OneBody(a,b)*GetTBME(ch_bra,ch_ket,b,i,a,j) * phase*sixj_a*hatfactor;

                      sixj_b = modelspace->GetSixJ(oi.j2/2.0, ob.j2/2.0, tbc_bra.J, oa.j2/2.0, oj.j2/2.0, opright.rank_J);
                      phase = modelspace->phase((oj.j2+oa.j2)/2);
                      hatfactor = (2*tbc_bra.J+1)*sqrt((ob.j2+1)/(oi.j2+1));
                      out.OneBody(i,j) += opright.OneBody(b,a)*GetTBME(ch_bra,ch_ket,b,i,a,j) * phase*sixj_b*hatfactor;
                   }
                }

             }
          }
      }
   }
}




//*****************************************************************************************
//
//      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
//        |   (  )   _     |   (  )                                                                                      
//        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
//      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
//                                     
//
// -- AGREES WITH NATHAN'S RESULTS 
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm221st(const Operator& opright, Operator& out) const
{

   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
   {
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int npq = tbc_bra.GetNumberKets();
      for (int ch_ket : opright.TwoBodyTensorChannels[ch_bra] )
      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);

         arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& Matrixpp = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& Matrixhh = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);

         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);
         Matrixpp = (LHS * tbc_bra.Proj_pp * RHS  -  RHS * tbc_bra.Proj_pp * LHS);
         Matrixhh = (LHS * tbc_bra.Proj_hh * RHS  -  RHS * tbc_bra.Proj_hh * LHS);
      

      // If commutator is hermitian or antihermitian, we only
      // need to do half the sum. I should add this in later.
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         for (int j=0;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mpp.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (int &c : modelspace->particles)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mhh.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            }
            #pragma omp critical
            out.OneBody(i,j) += cijJ / (oi.j2 +1.0);
         } // for j
      } // for i
   } //for ch_ket
   } //for ch_bra
}





//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
// Right now, this is the slowest one...
// Haven't converted this one yet...
void Operator::comm122st(const Operator& opright, Operator& opout ) const
{
   int herm = opout.IsHermitian() ? 1 : -1;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      arma::mat& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      int norbits = modelspace->GetNumberOrbits();
      for (int ibra = 0;ibra<npq; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         int indx_ij = ibra;
         double pre_ij = i==j ? sqrt(2) : 1;
         for (int iket=ibra;iket<npq; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            int indx_kl = iket;
            double pre_kl = k==l ? sqrt(2) : 1;

            double cij = 0;
            double ckl = 0;
            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
               int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
               int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
               int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));

               if (indx_aj >= 0)
               {
                  double pre_aj = a>j ? tbc.GetKet(indx_aj).Phase(tbc.J) : 1;
                  if (a==j) pre_aj *= sqrt(2);
                  ckl += pre_aj  * OneBody(i,a) * RHS(indx_aj,indx_kl);
                  ckl -= pre_aj  * opright.OneBody(i,a) * LHS(indx_aj,indx_kl);
               }

               if (indx_ia >= 0)
               {
                  double pre_ia = i>a ? tbc.GetKet(indx_ia).Phase(tbc.J) : 1;
                  if (i==a) pre_ia *= sqrt(2);
                  ckl += pre_ia * OneBody(j,a) * RHS(indx_ia,indx_kl);
                  ckl -= pre_ia * opright.OneBody(j,a) * LHS(indx_ia,indx_kl);
               }

               if (indx_al >= 0)
               {
                  double pre_al = a>l ? tbc.GetKet(indx_al).Phase(tbc.J) : 1;
                  if (a==l) pre_al *= sqrt(2);
                  cij -= pre_al * OneBody(a,k) * RHS(indx_ij,indx_al);
                  cij += pre_al * opright.OneBody(a,k) * LHS(indx_ij,indx_al);
               }

               if (indx_ka >= 0)
               {
                  double pre_ka = k>a ? tbc.GetKet(indx_ka).Phase(tbc.J) : 1;
                  if (k==a) pre_ka *= sqrt(2);
                  cij -= pre_ka * OneBody(a,l) * RHS(indx_ij,indx_ka);
                  cij += pre_ka * opright.OneBody(a,l) * LHS(indx_ij,indx_ka);
               }

            }
            cijkl = (ckl*pre_kl + cij*pre_ij) / sqrt( (1.0+bra.delta_pq())*(1.0+ket.delta_pq()) );
            #pragma omp critical
            {
            OUT(ibra,iket) += cijkl;
            if (ibra != iket)
               OUT(iket,ibra) += herm * cijkl;
            }
         }
      }
   }
}





//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |              ( note that   1-n_a-n_b  =  nbar_a nbar_b - n_an_b )
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm222_pp_hhst(const Operator& opright, Operator& opout ) const
{
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
   {
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int npq = tbc_bra.GetNumberKets();
      for (int ch_ket : opright.TwoBodyTensorChannels[ch_bra] )
      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);


         arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_bra);
         arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch_bra).at(ch_ket);

         arma::mat Matrixpp = LHS * tbc_bra.Proj_pp * RHS;
         arma::mat Matrixhh = LHS * tbc_bra.Proj_hh * RHS;

         if (opout.IsHermitian())
         {
            Matrixpp -= Matrixpp.t();
            Matrixhh -= Matrixhh.t();
         }
         else if (opout.IsAntiHermitian())
         {
            Matrixpp += Matrixpp.t();
            Matrixhh += Matrixhh.t();
         }
         else
         {
            Matrixpp -=  RHS.t() * tbc_bra.Proj_pp * LHS.t();
            Matrixhh -=  RHS.t() * tbc_bra.Proj_hh * LHS.t();
         }

      
         OUT += Matrixpp - Matrixhh;
      }
   }
}







// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
void Operator::comm222_pp_hh_221st(const Operator& opright, Operator& opout ) const 
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;


   #pragma omp parallel for schedule(dynamic,5)
   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
   {
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);

      arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_bra);

      for (int ch_ket : opright.TwoBodyTensorChannels[ch_bra] )
      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);

         arma::mat& RHS  = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& OUT2 = (arma::mat&)   opout.TwoBody.at(ch_bra).at(ch_ket);

         arma::mat& Matrixpp = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& Matrixhh = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);
      
         Matrixpp =  LHS * tbc_bra.Proj_pp * RHS;
         Matrixhh =  LHS * tbc_bra.Proj_hh * RHS;

         if (opout.IsHermitian())
         {
            Matrixpp -= Matrixpp.t();
            Matrixhh -= Matrixhh.t();
         }
         else if (opout.IsAntiHermitian())
         {
            Matrixpp += Matrixpp.t();
            Matrixhh += Matrixhh.t();
         }
         else  // is this right?
         {
            Matrixpp -=  RHS.t() * tbc_bra.Proj_pp * LHS.t();
            Matrixhh -=  RHS.t() * tbc_bra.Proj_hh * LHS.t();
         }

         // The two body part
         OUT2 += Matrixpp - Matrixhh;


         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);

      // The one body part
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         for (int j=0;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mpp.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (int &c : modelspace->particles)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mhh.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            }
            #pragma omp critical
            opout.OneBody(i,j) += cijJ / (oi.j2 +1.0);
         } // for j
       } // for i
     } //for ch_ket
   } //for ch_bra
}




//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.             [X2,Y2](2)_ph = -sum_ab (na-nb) sum_J' (-1)^(J'+ja+jb)
//
//                                                                             v----|    v----|                                    v----|    v----|
//   |          |      |          |              * [ (-1)^(ji+jl) { jk jl J } <bi|X|aj> <ak|Y|bl>   - (-1)^(jj+jl-J') { jk jl J } <bj|X|ai> <ak|Y|bl>
//   |     __Y__|      |     __X__|                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |                                            v----|    v----|                                    v----|    v----|
//   |____\/    |      |____\/    |            -  (-1)^(ji+jk-J') { jk jl J } <bi|X|aj> <al|Y|bk>   + (-1)^(jj+jk)    { jk jl J } <bj|X|ai> <al|Y|bk>  ]
//   |  X       |      |  Y       |                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//           
//            
// -- This appears to agree with Nathan's results
//
// Haven't converted this one yet...
void Operator::comm222_phst(const Operator& opright, Operator& opout ) const
{
   int herm = opout.IsHermitian() ? 1 : -1;

   // Update Cross-coupled matrix elements
   vector<arma::mat> X_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> X_TwoBody_CC_right (nChannels, arma::mat() );

   vector<arma::mat> Y_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> Y_TwoBody_CC_right (nChannels, arma::mat() );

   CalculateCrossCoupled(X_TwoBody_CC_left, X_TwoBody_CC_right );
   //CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   opright.CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );

   // Construct the intermediate matrices N1 and N2
   vector<arma::mat> N1 (nChannels, arma::mat() );
   vector<arma::mat> N2 (nChannels, arma::mat() );
 // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
 // I think that N1 == N2.t() ??? Check this.
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] = X_TwoBody_CC_left[ch] *  Y_TwoBody_CC_right[ch];
      N2[ch] = Y_TwoBody_CC_left[ch] *  X_TwoBody_CC_right[ch];
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      int nKets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;

         for (int iket=ibra; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;


            double comm = 0;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi.l+ok.l)%2;
            int Tz_cc = abs(oi.tz2+ok.tz2)/2;
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+tbc.J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               double me2 = N2[ch_cc](indx_jl,indx_ik);
               double me3 = N2[ch_cc](indx_ik,indx_jl);
               double me4 = N1[ch_cc](indx_ik,indx_jl);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
            }

            parity_cc = (oi.l+ol.l)%2;
            Tz_cc = abs(oi.tz2+ol.tz2)/2;
            jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_il = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,l),max(i,l));
               int indx_jk = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,k),max(j,k));
               if (i>l)
                 indx_il += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>k)
                 indx_jk += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               double me2 = N2[ch_cc](indx_il,indx_jk);
               double me3 = N2[ch_cc](indx_jk,indx_il);
               double me4 = N1[ch_cc](indx_jk,indx_il);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
            }

            comm /= sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
            #pragma omp critical
            {
               OUT(ibra,iket) += comm;
//              opout.TwoBody[ch](ibra,iket) += comm;
              if (iket != ibra)
              {
                 OUT(iket,ibra) += herm*comm;
//                 opout.TwoBody[ch](iket,ibra) += herm*comm;
              }
            }
         }
      }
   }
}


