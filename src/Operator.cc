
#include "Operator.hh"
#include "AngMom.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif

using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

double  Operator::bch_transform_threshold = 1e-6;
double  Operator::bch_product_threshold = 1e-4;
map<string, double> Operator::timer;



/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator() :
   modelspace(NULL), nChannels(0), hermitian(true), antihermitian(false),
    rank_J(0), rank_T(0), parity(0), particle_rank(2)
{}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    hermitian(true), antihermitian(false), modelspace(&ms), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(ms.GetNumberTwoBodyChannels(), map<int,arma::mat>() ),
    TwoBodyTensorChannels(ms.GetNumberTwoBodyChannels(), vector<int>() ),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank)
{
//  modelspace = &ms;
  AllocateTwoBody();
  if (particle_rank >=3)
  {
    AllocateThreeBody();
  }
}



Operator::Operator(ModelSpace& ms) :
    hermitian(true), antihermitian(false), modelspace(&ms), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(ms.GetNumberTwoBodyChannels(), map<int,arma::mat>() ),
    TwoBodyTensorChannels(ms.GetNumberTwoBodyChannels(), vector<int>() ),
    rank_J(0), rank_T(0), parity(0), particle_rank(2)
{
//  modelspace = &ms;
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
        if ( abs(tbc_bra.J-tbc_ket.J)>rank_J ) continue;
        if ( (tbc_bra.J+tbc_ket.J)<rank_J ) continue;
        if ( abs(tbc_bra.Tz-tbc_ket.Tz)>rank_T ) continue;
        if ( (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;
        
        int n_bra = tbc_bra.GetNumberKets();
        TwoBodyTensorChannels[ch_bra].push_back(ch_ket);
        TwoBody[ch_bra][ch_ket] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
        
     }
  }
}


// Confusing nomenclature: J2 means 2 times the total J of the three body system
void Operator::AllocateThreeBody()
{
  E3max = modelspace->GetN3max();
  cout << "Begin AllocateThreeBody() with E3max = " << E3max << endl;
  vector<double> zerovector(5,0.0);
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
      for (int e=0; e<= (d==a ? b : d); e+=2)
      {
       Orbit& oe = modelspace->GetOrbit(e);
       for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
       {
        Orbit& of = modelspace->GetOrbit(f);

        // conserve parity
        if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0) continue;


        int Jab_min = abs(oa.j2-ob.j2)/2;
        int Jde_min = abs(od.j2-oe.j2)/2;
        int Jab_max = (oa.j2+ob.j2)/2;
        int Jde_max = (od.j2+oe.j2)/2;


        for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
        {
         for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
         {
           int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
           int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
           for (int J2=J2_min; J2<=J2_max; J2+=2)
           {
             ThreeBody[{a,b,c,d,e,f,J2,Jab,Jde}] = {0.,0.,0.,0.,0.};
           } //J2
         } //Jde
        } //Jab
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
   E2max         = op.E2max;
   E3max         = op.E3max;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   TwoBody       = op.TwoBody;
   ThreeBody     = op.ThreeBody;
   TwoBodyTensorChannels       = op.TwoBodyTensorChannels;
//   timer         = op.timer;
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
      for ( auto &twobody : TwoBody[ch] )
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
Operator operator*(const double lhs, const Operator&& rhs)
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

Operator Operator::operator-() const
{
   return (*this)*-1.0;
}



///////// SETTER_GETTERS ///////////////////////////

/// This returns the matrix element times a factor \f$ \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \f$
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
   if (a==b) phase *= SQRT2;
   if (c==d) phase *= SQRT2;
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
   if (IsHermitian())
      TwoBody.at(ch_bra).at(ch_ket)(iket,ibra) = tbme;
   else if(IsAntiHermitian())
      TwoBody.at(ch_bra).at(ch_ket)(iket,ibra) = -tbme;
}
void Operator::AddToTBME(int ch_bra, int ch_ket, int iket, int ibra, double tbme)
{
   TwoBody.at(ch_bra).at(ch_ket)(ibra,iket) += tbme;
   if (IsHermitian())
      TwoBody.at(ch_bra).at(ch_ket)(iket,ibra) += tbme;
   else if(IsAntiHermitian())
      TwoBody.at(ch_bra).at(ch_ket)(iket,ibra) -= tbme;
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

/// Returns an unnormalized monopole-like (angle-averaged) term
/// \f[ \bar{V}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \sum_{J}(2J+1) V_{ijkl}^J \f]
///
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





//*******************************************************************
/// Get three body matrix element in proton-neutron formalism.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
double Operator::GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int a, int b, int c, int d, int e, int f)
{

   double tza = modelspace->GetOrbit(a).tz2*0.5;
   double tzb = modelspace->GetOrbit(b).tz2*0.5;
   double tzc = modelspace->GetOrbit(c).tz2*0.5;
   double tzd = modelspace->GetOrbit(d).tz2*0.5;
   double tze = modelspace->GetOrbit(e).tz2*0.5;
   double tzf = modelspace->GetOrbit(f).tz2*0.5;

   double Vpn=0;
   int Tmin = min( abs(tza+tzb+tzc), abs(tzd+tze+tzf) );
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
           Vpn += CG1*CG2*CG3*CG4*GetThreeBodyME(Jab_in,Jde_in,J2,tab,tde,T,a,b,c,d,e,f);
         }
      }
   }
   return Vpn;
}


//*******************************************************************
/// Get three body matrix element in isospin formalism
/// \f$ V_{abcdef}^{J_{ab}J_{de}Jt_{ab}t_{de}T} \f$
///  (which is how they're stored).
/// The elements are stored with the following restrictions: \f$ a\geq b \geq c\f$,
/// \f$ d\geq e \geq f\f$, \f$ a\geq d\f$. If \f$ a=d\f$ then \f$ b \geq e\f$,
/// and if \f$ b=e \f$ then \f$ c \geq f \f$.
/// Other orderings are obtained by recoupling on the fly.
//*******************************************************************
double Operator::GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   return AddToThreeBodyME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in,0.0);
}

//*******************************************************************
/// Set a three body matrix element. Since only a subset of orbit
/// orderings are stored, we need to recouple if the input ordering
/// is different.
//*******************************************************************
void Operator::SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V)
{
   AddToThreeBodyME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in,V);
}

//*******************************************************************
/// Since setting and getting three body matrix elements requires
/// almost identical code, they are combined into one function
/// which adds \f$V_{in}\f$ to the matrix element and returns
/// its value \f$V_{out}\f$. To access a matrix element, we set
/// \f$V_{in}=0\f$. To set the matrix element, we simply
/// disregard $\V_{out}\f$.
//*******************************************************************
double Operator::AddToThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V_in)
{

   // Re-order so that a>=b>=c, d>=e>=f
   int a,b,c,d,e,f;
   int abc_recoupling_case = SortThreeBodyOrbits(a_in,b_in,c_in,a,b,c);
   int def_recoupling_case = SortThreeBodyOrbits(d_in,e_in,f_in,d,e,f);

   if (d>a or (d==a and e>b) or (d==a and e==b and f>c))
   {
      swap(a,d);
      swap(b,e);
      swap(c,f);
      swap(Jab_in,Jde_in);
      swap(tab_in,tde_in);
      swap(abc_recoupling_case, def_recoupling_case);
   }

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);
   if (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) return 0;
   if (2*(od.n+oe.n+of.n)+od.l+oe.l+of.l > E3max) return 0;


   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;
   double je = oe.j2*0.5;
   double jf = of.j2*0.5;


   int Jab_min = max( abs(ja-jb), abs(J2*0.5-jc) );
   int Jab_max = min( ja+jb, J2*0.5+jc );
   int Jde_min = max( abs(jd-je), abs(J2*0.5-jf) );
   int Jde_max = min( jd+je, J2*0.5+jf );

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;

   if (abc_recoupling_case==0 or abc_recoupling_case==2)
   {
     Jab_min = Jab_max = Jab_in;
     tab_min = tab_max = tab_in;
   }
   if (def_recoupling_case==0 or def_recoupling_case==2)
   {
      Jde_min = Jde_max = Jde_in;
      tde_min = tde_max = tde_in;
   }

   double V_out = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_abc = ThreeBodyRecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,J2);
       double Cj_def = ThreeBodyRecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,J2);

       array<double,5>& vj = ThreeBody.at({a,b,c,d,e,f,J2,Jab,Jde});
       for (int tab=tab_min; tab<=tab_max; ++tab)
       {
         for (int tde=tde_min; tde<=tde_max; ++tde)
         {
           double Ct_abc = ThreeBodyRecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,T2);
           double Ct_def = ThreeBodyRecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,T2);

           int Tindex = 2*tab + tde + (T2-1)/2;

           vj[Tindex] += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj[Tindex];

         }
       }
     }
   }
   return V_out;
}



//*******************************************************************
/// Coefficients for recoupling three body matrix elements
//*******************************************************************
double Operator::ThreeBodyRecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J)
{
   switch (recoupling_case)
   {
    case 0: return Jab==Jab_in ? 1 : 0;
    case 1: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 2: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case 3: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    case 4: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 5: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    default: return 0;
    }
}


//*******************************************************************
/// Rearrange orbits (abc) so that a>=b>=c
/// and return an int which reflects the required reshuffling
/// - 0: (abc)_in -> (abc)
/// - 1: (acb)_in -> (abc)
/// - 2: (bac)_in -> (abc)
/// - 3: (bca)_in -> (abc)
/// - 4: (cab)_in -> (abc)
/// - 5: (cba)_in -> (abc)
//*******************************************************************
int Operator::SortThreeBodyOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c)
{
   a=a_in;
   b=b_in;
   c=c_in;
   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);

   int recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? 0 : 1;
   else if (a_in==b)  recoupling_case = (b_in==a) ? 2 : 3;
   else               recoupling_case = (b_in==a) ? 4 : 5;

   a -= a%2;
   b -= b%2;
   c -= c%2;
   return recoupling_case;
}



void Operator::PrintTimes()
{
   cout << "==== TIMES ====" << endl;
   for ( auto it : timer )
   {
     cout << it.first << ":  " << it.second  << endl;
   }
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

//*************************************************************
///  Normal ordering of a 2body operator
///  currently this only handles scalar operators
//*************************************************************
Operator Operator::DoNormalOrdering2()
{
   Operator opNO = *this;

   // Trivial parts
   opNO.ZeroBody = ZeroBody;
   opNO.OneBody = OneBody;
   opNO.TwoBody = TwoBody;


   for (auto& k : modelspace->holes) // loop over hole orbits
   {
      opNO.ZeroBody += (modelspace->GetOrbit(k).j2+1) * OneBody(k,k);
   }


   int norbits = modelspace->GetNumberOrbits();

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;

      // Zero body part
//      for (auto& iket : tbc.KetIndex_hh) // loop over hole-hole kets in this channel
      for (auto& iket : tbc.GetKetIndex_hh() ) // loop over hole-hole kets in this channel
      {
        opNO.ZeroBody += GetTBME(ch,iket,iket) * (2*J+1);  // <ab|V|ab>  (a,b in core)
      }

      // One body part
      int ibra,iket;
      for (unsigned int a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         int bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for (unsigned int b=bstart;b<norbits;++b)
         {
            Orbit &ob = modelspace->GetOrbit(b);
            if (ob.j2 != oa.j2 or ob.tz2 != oa.tz2 or ob.l != oa.l) continue;
            for (auto& h : modelspace->holes)  // C++11 syntax
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



//*******************************************************************************
///   Normal ordering of a three body operator. Start by generating the normal ordered
///   two body piece, then use DoNormalOrdering2() to get the rest. (Note that there
///   are some numerical factors).
///   The normal ordered two body piece is 
///   \f[ \Gamma^J_{ijkl} = V^J_{ijkl} + \sum_a n_a  \sum_K \frac{2K+1}{2J+1} V^{(3)JJK}_{ijakla} \f]
///   Right now, this is only set up for scalar operators, but I don't anticipate
///   handling 3body tensor operators in the near future.
//*******************************************************************************
Operator Operator::DoNormalOrdering3()
{
   cout << "in DoNormalOrdering3() " << endl;
   Operator opNO3 = Operator(*modelspace);
   cout << "size of ThreeBody = " << ThreeBody.size() << endl;
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
            for (auto& a : modelspace->holes)
            {
               Orbit & oa = modelspace->GetOrbit(a);
               if ( (2*(oi.n+oj.n+oa.n)+oi.l+oj.l+oa.l)>E3max) continue;
               if ( (2*(ok.n+ol.n+oa.n)+ok.l+ol.l+oa.l)>E3max) continue;
               int kmin2 = abs(2*tbc.J-oa.j2);
               int kmax2 = 2*tbc.J+oa.j2;
               for (int K2=kmin2; K2<=kmax2; K2+=2)
               {
                  Gamma(ibra,iket) += (K2+1) * GetThreeBodyME_pn(tbc.J,tbc.J,K2,i,j,a,k,l,a);
               }
            }
            Gamma(ibra,iket) /= (2*tbc.J+1);
         }
      }
   }
   cout << "Done looping over channels" << endl;
   Operator opNO2 = opNO3.DoNormalOrdering2();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);

   // Also normal order the 1 and 2 body pieces
   opNO2 += DoNormalOrdering2();
   return opNO2;

}




void Operator::Erase()
{
  EraseZeroBody();
  EraseOneBody();
  EraseTwoBody();
  if (particle_rank >=3)
    EraseThreeBody();
}

void Operator::EraseOneBody()
{
   OneBody.zeros();
}

void Operator::EraseTwoBody()
{
   for (int ch=0;ch<nChannels;++ch)
   {
      for (auto &twobody : TwoBody[ch] )
      {
        arma::mat& matrix = twobody.second;
        matrix.zeros();
      }
   }
}

void Operator::EraseThreeBody()
{
 for (auto& it_Orb: ThreeBody)
 {
    it_Orb.second = {0.,0.,0.,0.,0.,};
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
     for ( auto &twobody : TwoBody[ch] )
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
     for ( auto &twobody : TwoBody[ch] )
     {
        arma::mat& matrix = twobody.second;
        matrix.eye();
     }
   }
}


//***********************************************
/// Calculates the kinetic energy operator in the 
/// harmonic oscillator basis.
/// \f[ t_{ab} = \frac{1}{2}\hbar\omega
/// \delta_{\ell_a \ell_b} \delta_{j_aj_b} \delta_{t_{za}t_{zb}}
/// \left\{
/// \begin{array}{ll}
/// 2n_a + \ell_a + \frac{3}{2} &: n_a=n_b\\
/// \sqrt{n_{a}(n_{a}+\ell_a + \frac{1}{2})} &: n_a=n_b+1\\
/// \end{array} \right. \f]
//***********************************************
void Operator::CalculateKineticEnergy()
{
   OneBody.zeros();
   int norbits = modelspace->GetNumberOrbits();
   int A = modelspace->GetTargetMass();
   double hw = modelspace->GetHbarOmega();
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

//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
//
//
//void Operator::DoPandyaTransformation( vector<arma::mat>& TBMECC)
// SCALAR VARIETY
/// \f[
///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  X^{J}_{ilkj}
/// \f]
void Operator::DoPandyaTransformation( Operator& opCC)
{

   for (int ch_bra=0; ch_bra<nChannels; ++ch_bra)
   {
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int nBras = tbc_bra.GetNumberKets();
      int J = tbc_bra.J;
      for (int& ch_ket : opCC.TwoBodyTensorChannels[ch_bra] )
      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
         int nKets = tbc_ket.GetNumberKets();
         arma::mat& MatrixCC = opCC.TwoBody.at(ch_bra).at(ch_ket);
         MatrixCC.zeros();
         
         for (int ibra=0; ibra<nKets; ++ibra)
         {
            Ket& bra = modelspace->GetKet(ibra);
            int i = bra.p;
            int j = bra.q;
            Orbit& oi = modelspace->GetOrbit(i);
            Orbit& oj = modelspace->GetOrbit(j);
            double ji = oi.j2/2.0;
            double jj = oj.j2/2.0;
            for (int iket=0; iket<nKets; ++iket) // Note: if X is hermitian, then X_CC is also hermitian.
            {
               Ket& ket = modelspace->GetKet(iket);
               int k = ket.p;
               int l = ket.q;
               Orbit& ok = modelspace->GetOrbit(k);
               Orbit& ol = modelspace->GetOrbit(l);
               if ( (oi.tz2+ol.tz2)!=(oj.tz2+ok.tz2)) continue;
               int Tz_std = (oi.tz2+ol.tz2)/2;
               int parity_std = (oi.l+ol.l)%2;
               double jk = ok.j2/2.0;
               double jl = ol.j2/2.0;
               int minJprime = max( abs(ji-jj), abs(jk-jl) );
               int maxJprime = min( ji+jj, jk+jl );
               for (int Jprime=minJprime; Jprime<=maxJprime; ++Jprime)
               {
                  double sixj = modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
                  MatrixCC(ibra,iket) -= (2*Jprime+1) * sixj * GetTBME(Jprime,parity_std,Tz_std,i,l,k,j);
               }
            }
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
void Operator::CalculateCrossCoupled(vector<arma::mat> &TwoBody_CC_left, vector<arma::mat> &TwoBody_CC_right)
{
   // loop over cross-coupled channels
   #pragma omp parallel for
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      arma::uvec kets_ph = tbc_cc.GetKetIndex_ph();
      int nph_kets = kets_ph.n_rows;
//      int nph_kets = tbc_cc.KetIndex_ph.size();
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
//         Ket & bra_cc = tbc_cc.GetKet( tbc_cc.KetIndex_ph[i_ph] );
         Ket & bra_cc = tbc_cc.GetKet( kets_ph[i_ph] );
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
Operator Operator::BCH_Transform(  Operator &Omega)
{
//   double bch_transform_threshold = 1e-6;
   int max_iter = 20;
//   int max_iter = 6;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   Operator OpNested = *this;
//   cout << "Begin BCH_Transform: ||H1|| = " << OpOut.OneBodyNorm()
//        << "  ||Nested1|| = " << OpNested.OneBodyNorm()
//        << "  ||Omega1|| = " << Omega.OneBodyNorm() << endl;
   for (int i=1; i<max_iter; ++i)
   {
      OpNested = Omega.Commutator(OpNested) / i;

//     cout << "Before adding: ||H1|| = " << OpOut.OneBodyNorm()
//        << "  ||Nested1|| = " << OpNested.OneBodyNorm()
//        << "  ||Omega1|| = " << Omega.OneBodyNorm() << endl;
      
      OpOut += OpNested;

//     cout << "In BCH_Transform: ||H1|| = " << OpOut.OneBodyNorm()
//        << "  ||Nested1|| = " << OpNested.OneBodyNorm()
//        << "  ||Omega1|| = " << Omega.OneBodyNorm() << endl;
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
Operator Operator::BCH_Product(  Operator &Y)
{
   Operator& X = *this;
//   double bch_product_threshold = 1e-4;

   Operator Z = X + Y; 
//   return Z; // Remove this!!
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
   for (int ch=0;ch<nChannels;++ch)
   {
      for (auto &twobody : TwoBody[ch] )
      {
        arma::mat& matrix = (arma::mat&) twobody.second;
        double n2 = arma::norm(matrix,"fro");
        nrm += n2*n2;
      }
   }
   return sqrt(nrm);
}

void Operator::Symmetrize()
{
   OneBody = arma::symmatu(OneBody);
   for (int ch=0;ch<nChannels;++ch)
   {
     for (int chket : TwoBodyTensorChannels[ch] )
     {
       arma::mat& X = TwoBody[ch].at(chket);
       X = arma::symmatu(X);
     }
   }
}

void Operator::AntiSymmetrize()
{
   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      for(int j=i+1;j<norb;++j)
      {
        OneBody(j,i) = -OneBody(i,j);
      }
   }
   for (int ch=0;ch<nChannels;++ch)
   {
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
     int nbras = tbc.GetNumberKets();
     for (int chket : TwoBodyTensorChannels[ch] )
     {
       TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(chket);
       int nkets = tbc_ket.GetNumberKets();
       for (int ibra=0; ibra<nbras; ++ibra)
       {
         for (int iket=ibra+1; iket<nkets; ++iket)
         {
            arma::mat& X = TwoBody[ch].at(chket);
            X(iket,ibra) = -X(ibra,iket);
         }
       }
     }
   }
}

Operator Operator::Commutator( Operator& opright)
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


Operator Operator::CommutatorScalarScalar( Operator& opright) 
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
      if (particle_rank>1 and opright.particle_rank>1)
        comm220ss(opright, out) ;
   }

    double t = omp_get_wtime();
   comm111ss(opright, out);
     t = omp_get_wtime() - t;
    timer["comm111ss"] += t;

    t = omp_get_wtime();
//   comm111st(opright, out);  // << equivalent in scalar case
   comm121ss(opright, out);
//   comm121st(opright, out);  // << equivalent in scalar case
     t = omp_get_wtime() - t;
    timer["comm121ss"] += t;

    t = omp_get_wtime();
   comm122ss(opright, out); //  This is the slow one for some reason.
//   comm122st(opright, out); // << equivalent in scalar case
     t = omp_get_wtime() - t;
    timer["comm122ss"] += t;

   if (particle_rank>1 and opright.particle_rank>1)
   {
    t = omp_get_wtime();
    comm222_pp_hh_221ss(opright, out);
     t = omp_get_wtime() - t;
    timer["comm222_pp_hh_221ss"] += t;
     
////   comm222_pp_hh_221st(opright, out); // << equivalent in scalar case

    t = omp_get_wtime();
    comm222_phss(opright, out);
//   comm222_phst_pandya(opright, out);
////   comm222_phst(opright, out);
     t = omp_get_wtime() - t;
    timer["comm222_phss"] += t;
   }


   if ( out.IsHermitian() )
   {
      out.Symmetrize();
   }
   else if (out.IsAntiHermitian() )
   {
      out.AntiSymmetrize();
   }


   return out;
}


// Calculate [S,T]
Operator Operator::CommutatorScalarTensor( Operator& opright) 
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
Operator Operator::CommutatorTensorTensor( Operator& opright) 
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
/// \f[
///  [X_{1)},Y_{(1)}]_{(0)} = \sum_{a} n_a (2j_a+1) \left(X_{(1)}Y_{(1)}-Y_{(1)}X_{(1)}\right)_{aa}
/// \f]
void Operator::comm110ss( Operator& opright, Operator& out) 
{
  if (IsHermitian() and opright.IsHermitian()) return ; // I think this is the case
  if (IsAntiHermitian() and opright.IsAntiHermitian()) return ; // I think this is the case

   int norbits = modelspace->GetNumberOrbits();
   arma::mat xyyx = OneBody*opright.OneBody - opright.OneBody*OneBody;
   for ( auto& a : modelspace->holes) // C++11 range-for syntax: loop over all elements of the vector <particle>
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
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) X_{abcd}^{J} Y_{cdab}^{J}
/// \f]
/// may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{ab} (\mathcal{P}_{hh} X_{(2)}^{J} \mathcal{P}_{pp} Y_{(2)}^{J})_{abab}
/// \f] where \f$ \mathcal{P}_{hh} \f$ is a projector onto hole-hole two body states.
void Operator::comm220ss(  Operator& opright, Operator& out) 
{
   if (IsHermitian() and opright.IsHermitian()) return; // I think this is the case
   if (IsAntiHermitian() and opright.IsAntiHermitian()) return; // I think this is the case
//   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//      #pragma omp critical  // Now this is just plain stupid.
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
/// \f[
/// [X_{(1)},Y_{(1)}]_{(1)} = X_{(1)}Y_{(1)} - Y_{(1)}X_{(1)}
/// \f]
void Operator::comm111ss( Operator & opright, Operator& out) 
{
   out.OneBody += OneBody*opright.OneBody - opright.OneBody*OneBody;
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
//                     (note: I think this should actually be)
//                                                = sum_ab (n_a nbar_b) sum_J (2J+1)/(2j_i+1)
//                                                      * y_ab xbiag - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
/// Returns \f$ [X_{(1)},Y_{(2)}] - [Y_{(1)},X_{(2)}] \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1}\sum_{ab} (n_a \bar{n}_b) \sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj})
/// \f]
void Operator::comm121ss( Operator& opright, Operator& out) 
{
   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for 
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
//      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      int jmin = out.IsNonHermitian() ? 0 : i;
      for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
//          Orbit &oj = modelspace->GetOrbit(j);
//          if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue; // At some point, make a OneBodyChannel class...
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (auto& b : modelspace->particles)
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
// \f[
// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{c}
// \left( \mathcal{P}_{pp} (X \mathcal{P}_{hh} Y^{J} 
// - Y^{J} \mathcal{P}_{hh} X^{J}) \mathcal{P}_{pp}
//  - \mathcal{P}_{hh} (X^{J} \mathcal{P}_{pp} Y^{J} 
//  -  Y^{J} \mathcal{P}_{pp} X^{J}) \mathcal{P}_{hh} \right)_{cicj}
// \f]
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{abc} (\bar{n}_a\bar{n}_bn_c + n_an_b\bar{n}_c)
///  (X^{J}_{ciab} Y^{J}_{abcj} - Y^{J}_{ciab}X^{J}_{abcj})
/// \f]
/// This may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1} \sum_{c} \sum_{J} (2J+1) \left( n_c \mathcal{M}^{J}_{pp,icjc} + \bar{n}_c\mathcal{M}^{J}_{hh,icjc} \right)
/// \f]
/// With the intermediate matrix \f[ \mathcal{M}^{J}_{pp} \equiv \frac{1}{2} (X^{J}\mathcal{P}_{pp} Y^{J} - Y^{J}\mathcal{P}_{pp}X^{J}) \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$
void Operator::comm221ss( Operator& opright, Operator& out) 
{

   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;

   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = (arma::mat&) TwoBody.at(ch).at(ch);
      auto& RHS = (arma::mat&) opright.TwoBody.at(ch).at(ch);
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
            for (auto& c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,i,c,j,c);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,i,c,j,c);
            }
            cijJ *= (2*tbc.J+1.0)/(oi.j2 +1.0);
            #pragma omp critical
            out.OneBody(i,j) +=  cijJ;
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
/// Returns \f$ [X_{(1)},Y_{(2)}]_{(2)} - [Y_{(1)},X_{(2)}]_{(2)} \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]^{J}_{ijkl} = \sum_{a} ( X_{ia}Y^{J}_{ajkl} + X_{ja}Y^{J}_{iakl} - X_{ak} Y^{J}_{ijal} - X_{al} Y^{J}_{ijka} )
/// \f]
void Operator::comm122ss( Operator& opright, Operator& opout ) 
{
   auto& L1 = OneBody;
   auto& R1 = opright.OneBody;

   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& L2 = TwoBody.at(ch).at(ch);
      auto& R2 = opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT = opout.TwoBody.at(ch).at(ch);


      int npq = tbc.GetNumberKets();
      int norbits = modelspace->GetNumberOrbits();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         double pre_ij = i==j ? SQRT2 : 1;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         arma::Row<double> L2_ij = L2.row(indx_ij); // trying this to better use the cache. not sure if it helps.
         arma::Row<double> R2_ij = R2.row(indx_ij);
         int klmin = opout.IsNonHermitian() ? 0 : indx_ij;
         for (int indx_kl=klmin;indx_kl<npq; ++indx_kl)
         {
            Ket & ket = tbc.GetKet(indx_kl);
            int k = ket.p;
            int l = ket.q;
            double pre_kl = k==l ? SQRT2 : 1;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            arma::vec L2_kl = L2.unsafe_col(indx_kl);
            arma::vec R2_kl = R2.unsafe_col(indx_kl);

            double cijkl = 0;


            for (int a : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
                 int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
                 if (indx_aj < 0) continue;
                 double pre_aj = a>j ? tbc.GetKet(indx_aj).Phase(tbc.J) : (a==j ? SQRT2 : 1);
//                 cijkl += pre_kl * pre_aj  * ( L1(i,a) * R2(indx_aj,indx_kl) - R1(i,a) * L2(indx_aj,indx_kl) );
                 cijkl += pre_kl * pre_aj  * ( L1(i,a) * R2_kl(indx_aj) - R1(i,a) * L2_kl(indx_aj) );
            }

            for (int a : modelspace->OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
            {
                 int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
                 if (indx_ia < 0) continue;
                 double pre_ia = i>a ? tbc.GetKet(indx_ia).Phase(tbc.J) : (i==a ? SQRT2 : 1);
//                 cijkl += pre_kl * pre_ia * ( L1(j,a) * R2(indx_ia,indx_kl) - R1(j,a) * L2(indx_ia,indx_kl) );
                 cijkl += pre_kl * pre_ia * ( L1(j,a) * R2_kl(indx_ia) - R1(j,a) * L2_kl(indx_ia) );
             }

            for (int a : modelspace->OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
                int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
                if (indx_al < 0) continue;
                double pre_al = a>l ? tbc.GetKet(indx_al).Phase(tbc.J) : (a==l ? SQRT2 : 1);
//                cijkl += pre_ij * pre_al * ( R1(a,k) * L2(indx_ij,indx_al) - L1(a,k) * R2(indx_ij,indx_al) );
                cijkl += pre_ij * pre_al * ( R1(a,k) * L2_ij(indx_al) - L1(a,k) * R2_ij(indx_al) );
            }

            for (int a : modelspace->OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
            {
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));
               if (indx_ka < 0) continue;
               double pre_ka = k>a ? tbc.GetKet(indx_ka).Phase(tbc.J) : (k==a ? SQRT2 : 1);
//               cijkl += pre_ij * pre_ka * ( R1(a,l) * L2(indx_ij,indx_ka) - L1(a,l) * R2(indx_ij,indx_ka) );
               cijkl += pre_ij * pre_ka * ( R1(a,l) * L2_ij(indx_ka) - L1(a,l) * R2_ij(indx_ka) );
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            OUT(indx_ij,indx_kl) += cijkl / norm;
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
/// Calculates the part of the commutator \f$ [X_{(2)},Y_{(2)}]_{(2)} \f$ which involves particle-particle
/// or hole-hole intermediate states.
/// \f[
/// [X_{(2)},Y_{(2)}]^{J}_{ijkl} = \frac{1}{2} \sum_{ab} (\bar{n}_a\bar{n}_b - n_an_b) (X^{J}_{ijab}Y^{J}_{ablk} - Y^{J}_{ijab}X^{J}_{abkl})
/// \f]
/// This may be written as
/// \f[
/// [X_{(2)},Y_{(2)}]^{J} = \mathcal{M}^{J}_{pp} - \mathcal{M}^{J}_{hh}
/// \f]
/// With the intermediate matrices
/// \f[
/// \mathcal{M}^{J}_{pp} \equiv \frac{1}{2}(X^{J} \mathcal{P}_{pp} Y^{J} - Y^{J} \mathcal{P}_{pp} X^{J})
/// \f]
/// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$.
void Operator::comm222_pp_hhss( Operator& opright, Operator& opout ) 
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







/// Since comm222_pp_hhss() and comm221ss() both require the ruction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
void Operator::comm222_pp_hh_221ss( Operator& opright, Operator& opout )  
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();

   Operator Mpp = opout;
   Operator Mhh = opout;

   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo, and the one-body part isn't thread-safe.
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = TwoBody.at(ch).at(ch);
      auto& RHS = opright.TwoBody.at(ch).at(ch);
      arma::mat& OUT =  opout.TwoBody.at(ch).at(ch);

      arma::mat & Matrixpp = Mpp.TwoBody.at(ch).at(ch);
      arma::mat & Matrixhh = Mhh.TwoBody.at(ch).at(ch);

      arma::uvec kets_pp = tbc.GetKetIndex_pp();
      arma::uvec kets_hh = tbc.GetKetIndex_hh();
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * RHS.rows(kets_hh);

      if (opout.IsHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (opout.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * LHS.rows(kets_hh);
      }

      // The two body part
      OUT += Matrixpp - Matrixhh;

      // The one body part
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         int jmin = opout.IsNonHermitian() ? 0 : i;
         double Jfactor = (2*tbc.J+1.0)/(oi.j2+1.0);
         for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
         {
            if (j<jmin) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,c,i,c,j);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,c,i,c,j);
            }
            opout.OneBody(i,j) += cijJ * Jfactor;
         } // for j
      } // for i
   } //for ch
/*
   if (not opout.IsNonHermitian())
   {
    if (opout.IsHermitian())
    {
//       opout.OneBody = arma::symmatu(opout.OneBody);
    }
    else
    {
     for (int i=0;i<norbits;++i)
     {
       for (int j=i;j<norbits;++j)
        {
           opout.OneBody(j,i) = herm* opout.OneBody(i,j);
        }
     }
    }
   }
*/
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
void Operator::comm222_phss( Operator& opright, Operator& opout ) 
{

   // Update Cross-coupled matrix elements
   vector<arma::mat> X_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> X_TwoBody_CC_right (nChannels, arma::mat() );

   vector<arma::mat> Y_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> Y_TwoBody_CC_right (nChannels, arma::mat() );

   double t = omp_get_wtime();
   CalculateCrossCoupled(X_TwoBody_CC_left, X_TwoBody_CC_right );
   opright.CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   t = omp_get_wtime() - t;
   timer["CalculateCrossCoupled"] += t;

   // Construct the intermediate matrices N1 and N2
   vector<arma::mat> N1 (nChannels, arma::mat() );
 // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
   for (int ch=0;ch<nChannels;++ch)
   {

      N1[ch] = X_TwoBody_CC_left[ch] * Y_TwoBody_CC_right[ch] - Y_TwoBody_CC_left[ch] * X_TwoBody_CC_right[ch];
      N1[ch] += N1[ch].t(); // maybe this works?  Yes !
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);
      int J = tbc.J;

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
         int ketmin = opout.IsNonHermitian() ? 0 : ibra;
         for (int iket=ketmin; iket<nKets; ++iket)
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
               double sixj = modelspace->GetSixJ(jk,jl,J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
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
               double sixj = modelspace->GetSixJ(jk,jl,J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            comm /= norm;
            OUT(ibra,iket) += comm;
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
// This is no different from the scalar-scalar version
void Operator::comm111st( Operator & opright, Operator& out) 
{
   out.OneBody += OneBody*opright.OneBody - opright.OneBody*OneBody;
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
// X is scalar, Y is tensor
// 
void Operator::comm121st( Operator& opright, Operator& out) 
{
   int norbits = modelspace->GetNumberOrbits();
   int nch = modelspace->GetNumberTwoBodyChannels();
   int Lambda = opright.rank_J;
   #pragma omp parallel for
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          Orbit &oj = modelspace->GetOrbit(j);
          double jj = oj.j2/2.0;
          if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue; // At some point, make a OneBodyChannel class...

          double zij = 0;
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             double ja = oa.j2/2.0;
             for (auto& b : modelspace->particles)  // C++11 syntax
             {
                Orbit &ob = modelspace->GetOrbit(b);
                double jb = ob.j2/2.0;
                for (int ch_bra=0;ch_bra<nch;++ch_bra)
                {
                   TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
                   double J = tbc_bra.J;

                   // part with opleft two body and opright one body, i.e. [X(2),Y(1)]  ==  -[Y(1),X(2)]
                   for (int ch_ket : TwoBodyTensorChannels[ch_bra])
                   {
                      double sixj_ab = modelspace->GetSixJ(jj, ja, J, jb, ji, Lambda);
                      double sixj_ba = modelspace->GetSixJ(jj, jb, J, ja, ji, Lambda);
                      int phase_ab = modelspace->phase(ji+jb+J);
                      int phase_ba = modelspace->phase(ji+ja+J);
                      double hatfactor_ab = (2*J+1)*sqrt((oa.j2+1));
                      double hatfactor_ba = (2*J+1)*sqrt((ob.j2+1));
   
                      zij -= opright.OneBody(a,b)*GetTBME(ch_bra,ch_ket,b,i,a,j) * phase_ab * sixj_ab * hatfactor_ab;
                      zij += opright.OneBody(b,a)*GetTBME(ch_bra,ch_ket,a,i,b,j) * phase_ba * sixj_ba * hatfactor_ba;
                   }

                   // part with opleft one body and opright two body, i.e. [X(1),Y(2)]
                   if (oa.j2 != ob.j2) continue;  // X(1) is scalar -> ja==jb
                   for (int  ch_ket : opright.TwoBodyTensorChannels[ch_bra])
                   {
                      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
                      double Jprime = tbc_ket.J;
                      double sixj = modelspace->GetSixJ(Jprime, J, Lambda, ji, jj, ja);
                      int phase = modelspace->phase(jj+ja+J);
                      double hatfactor = (2*J+1)*sqrt((2*Jprime+1));

                      zij += OneBody(a,b)*opright.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase * sixj * hatfactor;
                      zij -= OneBody(b,a)*opright.GetTBME(ch_bra,ch_ket,a,i,b,j) * phase * sixj * hatfactor; // Here hatfactor_ab == hatfactor_ba
                   }
                }
             }
          }
          out.OneBody(i,j) += zij/sqrt(oi.j2+1);
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
void Operator::comm221st( Operator& opright, Operator& out) 
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
            for (auto& c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mpp.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
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
// Agrees with previous code in the scalar-scalar limit
void Operator::comm122st( Operator& opright, Operator& opout ) 
{
   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();
   int Lambda = opout.rank_J;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch_bra=0; ch_bra<nChannels; ++ch_bra)
   {

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int J = tbc_bra.J;
      int npq = tbc_bra.GetNumberKets();



      for (int ch_ket : opright.TwoBodyTensorChannels[ch_bra] )
      {
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int Jprime = tbc_ket.J;
//         arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_bra);
//         arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch_bra).at(ch_ket);

      for (int ibra = 0;ibra<npq; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.0;
         double jj = oj.j2/2.0;
         for (int iket=ibra;iket<npq; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.0;
            double jl = ol.j2/2.0;

            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
              Orbit& oa = modelspace->GetOrbit(a);
              double ja = oa.j2/2.0;

               cijkl += OneBody(i,a) * opright.GetTBME(ch_bra,ch_ket,a,j,k,l);
               cijkl += OneBody(j,a) * opright.GetTBME(ch_bra,ch_ket,i,a,k,l);
               cijkl -= OneBody(a,k) * opright.GetTBME(ch_bra,ch_ket,i,j,a,l);
               cijkl -= OneBody(a,l) * opright.GetTBME(ch_bra,ch_ket,i,j,k,a);


               double sixj_ia = modelspace->GetSixJ(J,Jprime,Lambda,ja,ji,jj);
               int phase_ia = modelspace->phase(ji+jj+J);
               cijkl -= opright.OneBody(i,a) * GetTBME(ch_bra,ch_bra,a,j,k,l) * sixj_ia * phase_ia * sqrt((oi.j2+1)*(2*Jprime+1));

               double sixj_ja = modelspace->GetSixJ(J,Jprime,Lambda,ja,jj,ji);
               int phase_ja = modelspace->phase(ji+jj+J);
               cijkl -= opright.OneBody(j,a) * GetTBME(ch_bra,ch_bra,i,a,k,l) * sixj_ja * phase_ja * sqrt((oj.j2+1)*(2*Jprime+1));

               double sixj_ak = modelspace->GetSixJ(J,Jprime,Lambda,jk,ja,jl);
               int phase_ak = modelspace->phase(ja+jl+J+Lambda);
               cijkl += opright.OneBody(a,k) * GetTBME(ch_bra,ch_bra,i,j,a,l) * sixj_ak * phase_ak * sqrt((oa.j2+1)*(2*Jprime+1));

               double sixj_al = modelspace->GetSixJ(J,Jprime,Lambda,jl,ja,jk);
               int phase_al = modelspace->phase(ja+jk+J+Lambda);
               cijkl += opright.OneBody(a,l) * GetTBME(ch_bra,ch_bra,i,j,k,a) * sixj_al * phase_al * sqrt((oa.j2+1)*(2*Jprime+1));

            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            cijkl /= norm;
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
void Operator::comm222_pp_hhst( Operator& opright, Operator& opout ) 
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







// Since comm222_pp_hh and comm211 both require the ruction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
void Operator::comm222_pp_hh_221st( Operator& opright, Operator& opout )  
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();

   Operator Mpp = opout;
   Operator Mhh = opout;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
   {
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);

      auto& LHS = TwoBody.at(ch_bra).at(ch_bra);

      for (int ch_ket : opright.TwoBodyTensorChannels[ch_bra] )
      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
         if ( ch_ket != ch_bra )
         {
            cout << " AAH! ch_bra = " << ch_bra << "  ch_ket = " << ch_ket << endl;
         }

         auto& RHS  =  opright.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& OUT2 =    opout.TwoBody.at(ch_bra).at(ch_ket);

         arma::mat& Matrixpp =  Mpp.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& Matrixhh =  Mhh.TwoBody.at(ch_bra).at(ch_ket);
     
//         arma::uvec kets_pp = arma::uvec(tbc_bra.KetIndex_pp);
//         arma::uvec kets_hh = arma::uvec(tbc_bra.KetIndex_hh);
         arma::uvec kets_pp = tbc_bra.GetKetIndex_pp();
         arma::uvec kets_hh = tbc_bra.GetKetIndex_hh();
      
         Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
         Matrixhh =  LHS.cols(kets_hh) * RHS.rows(kets_hh);
 
         if (opout.IsHermitian())
         {
            Matrixpp += Matrixpp.t();
            Matrixhh += Matrixhh.t();
         }
         else if (opout.IsAntiHermitian()) // i.e. LHS and RHS are hermitian
         {
            Matrixpp -= Matrixpp.t();
            Matrixhh -= Matrixhh.t();
         }
         else 
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
         int jmin = opout.IsNonHermitian() ? 0 : i;
         for (int j=jmin;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ -=   Mpp.GetTBME(ch_bra,ch_ket,c,i,j,c)*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ -=   Mhh.GetTBME(ch_bra,ch_ket,c,i,j,c)*sixj;
            }
            #pragma omp critical
            opout.OneBody(i,j) += cijJ *hatfactor / sqrt(oi.j2 +1.0);
            if (not opout.IsNonHermitian())
            {
               opout.OneBody(j,i) = herm * opout.OneBody(i,j);
            }
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
void Operator::comm222_phst( Operator& opright, Operator& opout ) 
{
   int herm = opout.IsHermitian() ? 1 : -1;


   Operator Xcc = Operator(*modelspace, 0, 1, 0, 2);
   Operator Ycc = Operator(Xcc);

   DoPandyaTransformation(Xcc);
   opright.DoPandyaTransformation(Ycc);


   // Update Cross-coupled matrix elements
/*
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
 // I think that N1 == N2.t() ??? Check this. No.
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] = X_TwoBody_CC_left[ch] * Y_TwoBody_CC_right[ch] - Y_TwoBody_CC_left[ch] * X_TwoBody_CC_right[ch];
      N1[ch] += N1[ch].t(); // maybe this works?  Yes !
//      N2[ch] = Y_TwoBody_CC_left[ch] *  X_TwoBody_CC_right[ch];
//      N1[ch] -= N2[ch]; // maybe this works? Yes!
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
               comm -= (2*Jprime+1) * phase * sixj * (me1);
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
               comm -= (2*Jprime+1) * phase * sixj * (me1);
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            comm /= norm;
            #pragma omp critical
            {
               OUT(ibra,iket) += comm;
              if (iket != ibra)
              {
                 OUT(iket,ibra) += herm*comm;
              }
            }
            if (abs(comm) > 1e-7)
            {
               cout << modelspace->phase(ji+jj+jk+jl) << ":  "  << comm << endl;
            }

         }
      }
   }
*/
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
///
/// \f[
///  [X_{(2)}Y_{(2)}]^{J}_{ijkl} = \sum_{abJ'} (n_a - n_b) (2J'+1)
///  \left( \left\{ \begin{array}{lll}
///  j_i & j_ j & J \\
///  j_k & j_ l & J' \\
///  \end{array} \right\}
///  ( \bar{X}^{J'}_{i\bar{l}a\bar{b}} \bar{Y}^{J'}_{a\bar{b}k\bar{j}}
///  -\bar{X}^{J'}_{i\bar{l}a\bar{b}} \bar{X}^{J'}_{a\bar{b}k\bar{j}})
///  - (-1)^{J-j_1-j_2} [ i \leftrightarrow j ] \right)
/// \f]
/// This may be rewritten as
/// \f[
///   [X_{(2)}Y_{(2)}]^{J}_{ijkl} = \sum_{J'} (2J'+1)
///  \left\{ \begin{array}{lll}
///  j_i & j_ j & J \\
///  j_k & j_ l & J' \\
///  \end{array} \right\}
///  ( \bar{Z}^{J'}_{i\bar{l}k\bar{j}} -
/// (-1)^{J-j_i-j_j} 
///  \left\{ \begin{array}{lll}
///  j_j & j_ i & J \\
///  j_k & j_ l & J' \\
///  \end{array} \right\}
/// \bar{Z}^{J'}_{j\bar{l}k\bar{i}} )
/// \f]
///  with the Pandya transformed commutator defined as
/// \f[
///  \bar{Z}^{J'}_{i\bar{j}k\bar{l}} \equiv \sum_{ab} (n_a\bar{n}_b - \bar{n}_a n_b )
/// ( \bar{X}^{J'}_{i\bar{j}a\bar{b}} \bar{Y}^{J'}_{a\bar{b}k\bar{j}}
///  - \bar{Y}^{J'}_{i\bar{j}a\bar{b}} \bar{X}^{J'}_{a\bar{b}k\bar{j}} )
/// \f]
void Operator::comm222_phst_pandya( Operator& opright, Operator& opout ) 
{
   int herm = opout.IsHermitian() ? 1 : -1;
   if (opout.IsNonHermitian() ) herm = 0;
   Operator& X = *this;
   Operator& Y = opright;
   

   Operator XCC = Operator(*modelspace, 0,1,0,2);
   Operator YCC = Operator(*modelspace, 0,1,0,2);
   Operator WCC = Operator(*modelspace, 0,1,0,2);
   Operator W   = Operator(*modelspace, 0,0,0,2);

   X.DoPandyaTransformation(XCC);
   Y.DoPandyaTransformation(YCC);

   for (int chbra=0; chbra<nChannels; ++chbra)
   {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(chbra);
//     arma::uvec phbras = arma::uvec(tbc_bra.KetIndex_ph);
     arma::uvec phbras = tbc_bra.GetKetIndex_ph();
     for ( int chket : XCC.TwoBodyTensorChannels[chbra] )
     {
       TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(chket);
       arma::mat& wcc = WCC.TwoBody.at(chbra).at(chket);
       arma::mat& xcc = XCC.TwoBody.at(chbra).at(chket);
       arma::mat& ycc = YCC.TwoBody.at(chbra).at(chket);
//       arma::uvec phkets = arma::uvec(tbc_ket.KetIndex_ph);
       arma::uvec phkets = tbc_ket.GetKetIndex_ph();
       wcc = xcc.rows(phbras) * ycc.cols(phkets);
       wcc -= wcc.t();
     }
   }   
   WCC.DoPandyaTransformation(W);
/*   
   ZCC = 

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);
      int J = tbc.J;

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
               double sixj = modelspace->GetSixJ(jk,jl,J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
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
               double sixj = modelspace->GetSixJ(jk,jl,J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            comm /= norm;
            OUT(ibra,iket) += comm;
            if (iket != ibra and herm !=0)
            {
               OUT(iket,ibra) = herm * OUT(ibra,iket);
            }
         }
      }
   }
*/

}

