
#include "TwoBodyME.hh"
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif

// destructor defined for debugging purposes
TwoBodyME::~TwoBodyME()
{}

TwoBodyME::TwoBodyME()
: modelspace(NULL), nChannels(0), hermitian(true),antihermitian(false),
  rank_J(0), rank_T(0), parity(0)
{
//  cout << "Default TwoBodyME constructor" << endl;
}


TwoBodyME::TwoBodyME(ModelSpace* ms)
: modelspace(ms), nChannels(ms->GetNumberTwoBodyChannels()),
  hermitian(true), antihermitian(false), rank_J(0), rank_T(0), parity(0)
{
  Allocate();
}


TwoBodyME::TwoBodyME(ModelSpace* ms, int rJ, int rT, int p)
: modelspace(ms), nChannels(ms->GetNumberTwoBodyChannels()),
  hermitian(true), antihermitian(false), rank_J(rJ), rank_T(rT), parity(p)
{
  Allocate();
}


 TwoBodyME& TwoBodyME::operator*=(const double rhs)
 {
//   for ( auto& itmat : MatEl )
   for ( auto& itmat : MtxEl )
   {
      itmat *= rhs;
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator+=(const TwoBodyME& rhs)
 {
//   for ( auto& itmat : MatEl )
   for ( auto& itindex : rhs.MtxIndex )
   {
//      int ch_bra = itindex.first[0];
//      int ch_ket = itindex.first[1];
//      MtxEl[itindex.second] += rhs.GetMatrix(ch_bra,ch_ket);
      GetMatrix(itindex.first) += rhs.GetMatrix(itindex.first);
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator+=(const double rhs)
 {
   for ( auto& itmat : MtxEl )
   {
      itmat += rhs;
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator-=(const TwoBodyME& rhs)
 {
//   for ( auto& itmat : rhs.MatEl )
   for ( auto& itindex : rhs.MtxIndex )
   {
//      int ch_bra = itmat.first[0];
//      int ch_ket = itmat.first[1];
//      GetMatrix(ch_bra,ch_ket) -= itmat.second;
      GetMatrix(itindex.first) -= rhs.GetMatrix(itindex.first);
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator/=(const TwoBodyME& rhs)
 {
   for ( auto& itindex : MtxIndex )
   {
      GetMatrix(itindex.first) /= rhs.GetMatrix(itindex.first);
   }
   return *this;
 }

 TwoBodyME abs(const TwoBodyME& in)
 {
   TwoBodyME out(in);
   for ( auto& itmat : out.MtxEl )   itmat = arma::abs(itmat);
   return out;
 }

 arma::mat& TwoBodyME::GetMatrix(int chbra, int chket)
 {
    return MtxEl[ MtxIndex.at({chbra,chket}) ];
 }

 const arma::mat& TwoBodyME::GetMatrix(int chbra, int chket) const
 {
    return MtxEl[ MtxIndex.at({chbra,chket}) ];
 }


void TwoBodyME::Allocate()
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
//        MatEl[{ch_bra,ch_ket}] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
        MtxIndex[{ch_bra,ch_ket}] =  MtxEl.size();
        MtxEl.emplace_back( tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros );
     }
  }
}

void TwoBodyME::SetHermitian()
{
  hermitian = true;
  antihermitian = false;
}

void TwoBodyME::SetAntiHermitian()
{
  hermitian = false;
  antihermitian = true;
}

void TwoBodyME::SetNonHermitian()
{
  hermitian = false;
  antihermitian = false;
}



/// This returns the matrix element times a factor \f$ \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \f$
double TwoBodyME::GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const
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
   return phase * GetMatrix(ch_bra,ch_ket)(bra_ind, ket_ind);
}

/// This returns the normalized matrix element 
double TwoBodyME::GetTBME_norm(int ch_bra, int ch_ket, int a, int b, int c, int d) const
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
   return phase * GetMatrix(ch_bra,ch_ket)(bra_ind, ket_ind);
}

void TwoBodyME::SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(min(c,d),max(c,d));
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   GetMatrix(ch_bra,ch_ket)(bra_ind,ket_ind) = phase * tbme;
   if (ch_ket != ch_bra) return;
   if (hermitian) GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind) = phase * tbme;
   if (antihermitian) GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind) = - phase * tbme;
}


void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(min(c,d),max(c,d));
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   GetMatrix(ch_bra,ch_ket)(bra_ind,ket_ind) += phase * tbme;
   if (ch_bra==ch_ket and ket_ind==bra_ind) return;
   if (hermitian) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) += phase * tbme;
   if (antihermitian) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) -=  phase * tbme;
}

double TwoBodyME::GetTBME(int ch_bra, int ch_ket, Ket &bra, Ket &ket) const
{
   return GetTBME(ch_bra,ch_ket,bra.p,bra.q,ket.p,ket.q);
}
void TwoBodyME::SetTBME(int ch_bra, int ch_ket, Ket& bra, Ket& ket, double tbme)
{
   SetTBME(ch_bra, ch_ket, bra.p,bra.q,ket.p,ket.q,tbme);
}
void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, Ket& bra, Ket& ket, double tbme)
{
   AddToTBME(ch_bra, ch_ket, bra.p,bra.q,ket.p,ket.q, tbme);
}

double TwoBodyME::GetTBME_norm(int ch_bra, int ch_ket, int ibra, int iket) const
{
   return GetMatrix(ch_bra,ch_ket)(ibra,iket);
}

void TwoBodyME::SetTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme)
{
   GetMatrix(ch_bra,ch_ket)(ibra,iket) = tbme;
   if (IsHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) = tbme;
   else if(IsAntiHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) = -tbme;
}
void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme)
{
   GetMatrix(ch_bra,ch_ket)(ibra,iket) += tbme;
   if (IsHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) += tbme;
   else if(IsAntiHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) -= tbme;
}

double TwoBodyME::GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d) const
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   return GetTBME(ch_bra,ch_ket,a,b,c,d);
}
void TwoBodyME::SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   SetTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}
void TwoBodyME::AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   AddToTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}
double TwoBodyME::GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket) const
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   return GetTBME(ch_bra,ch_ket,bra,ket);
}
void TwoBodyME::SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   SetTBME(ch_bra,ch_ket,bra,ket,tbme);
}
void TwoBodyME::AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme)
{
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,p_bra,t_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,p_ket,t_ket);
   AddToTBME(ch_bra,ch_ket,bra,ket,tbme);
}
double TwoBodyME::GetTBME_J(int j_bra, int j_ket, int a, int b, int c, int d) const
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   int parity_bra = (oa.l+ob.l)%2;
   int parity_ket = (oc.l+od.l)%2;
   int Tz_bra = (oa.tz2+ob.tz2)/2;
   int Tz_ket = (oc.tz2+od.tz2)/2;
   if ( (parity+parity_bra+parity_ket)%2 > 0) return 0;
   if ( abs(Tz_bra-Tz_ket)>rank_T) return 0;
   if ( abs(j_bra-j_ket) > rank_J) return 0;
   if ( j_bra + j_ket < rank_J) return 0;
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,parity_bra,Tz_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,parity_ket,Tz_ket);
   if (ch_bra <= ch_ket)
     return GetTBME(ch_bra,ch_ket,a,b,c,d);
   return modelspace->phase(j_bra - j_ket) * GetTBME(ch_ket,ch_bra,c,d,a,b);
}
void TwoBodyME::SetTBME_J(int j_bra, int j_ket, int a, int b, int c, int d, double tbme)
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,(oa.l+ob.l)%2,(oa.tz2+ob.tz2)/2);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,(oc.l+od.l)%2,(oc.tz2+od.tz2)/2);
   SetTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}
void TwoBodyME::AddToTBME_J(int j_bra, int j_ket, int a, int b, int c, int d, double tbme)
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,(oa.l+ob.l)%2,(oa.tz2+ob.tz2)/2);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,(oc.l+od.l)%2,(oc.tz2+od.tz2)/2);
   AddToTBME(ch_bra,ch_ket,a,b,c,d,tbme);
}

double TwoBodyME::GetTBME_J_norm(int j_bra, int j_ket, int a, int b, int c, int d) const
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   int parity_bra = (oa.l+ob.l)%2;
   int parity_ket = (oc.l+od.l)%2;
   int Tz_bra = (oa.tz2+ob.tz2)/2;
   int Tz_ket = (oc.tz2+od.tz2)/2;
   if ( (parity+parity_bra+parity_ket)%2 > 0) return 0;
   if ( abs(Tz_bra-Tz_ket)>rank_T) return 0;
   if ( abs(j_bra-j_ket) > rank_J) return 0;
   if ( j_bra + j_ket < rank_J) return 0;
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,parity_bra,Tz_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,parity_ket,Tz_ket);
   if (ch_bra <= ch_ket)
     return GetTBME_norm(ch_bra,ch_ket,a,b,c,d);
   return modelspace->phase(j_bra - j_ket) * GetTBME_norm(ch_ket,ch_bra,c,d,a,b);
}

// for backwards compatibility...
double TwoBodyME::GetTBME(int ch, int a, int b, int c, int d) const
{
   return GetTBME(ch,ch,a,b,c,d);
}
double TwoBodyME::GetTBME_norm(int ch, int a, int b, int c, int d) const
{
   return GetTBME_norm(ch,ch,a,b,c,d);
}

void TwoBodyME::SetTBME(int ch, int a, int b, int c, int d, double tbme)
{
   return SetTBME(ch,ch,a,b,c,d,tbme);
}

double TwoBodyME::GetTBME(int ch, Ket &bra, Ket &ket) const
{
   return GetTBME(ch,ch,bra.p,bra.q,ket.p,ket.q);
}

double TwoBodyME::GetTBME_norm(int ch, Ket &bra, Ket &ket) const
{
   return GetTBME_norm(ch,ch,bra.p,bra.q,ket.p,ket.q);
}

void TwoBodyME::SetTBME(int ch, Ket& bra, Ket& ket, double tbme)
{
   SetTBME(ch,ch,bra.p,bra.q,ket.p,ket.q,tbme);
}

void TwoBodyME::AddToTBME(int ch, Ket& bra, Ket& ket, double tbme)
{
   AddToTBME(ch,ch,bra.p,bra.q,ket.p,ket.q,tbme);
}


double TwoBodyME::GetTBME_norm(int ch, int ibra, int iket) const
{
   return GetTBME_norm(ch,ch,ibra,iket);
}

void TwoBodyME::SetTBME(int ch, int ibra, int iket, double tbme)
{
   SetTBME(ch,ch,ibra,iket,tbme);
}

void TwoBodyME::AddToTBME(int ch, int ibra, int iket, double tbme)
{
   AddToTBME(ch,ch,ibra,iket,tbme);
}



double TwoBodyME::GetTBME(int j, int p, int t, int a, int b, int c, int d) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,ch,a,b,c,d);
}

double TwoBodyME::GetTBME_norm(int j, int p, int t, int a, int b, int c, int d) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME_norm(ch,ch,a,b,c,d);
}

void TwoBodyME::SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,ch,a,b,c,d,tbme);
}

void TwoBodyME::AddToTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   AddToTBME(ch,ch,a,b,c,d,tbme);
}

double TwoBodyME::GetTBME(int j, int p, int t, Ket& bra, Ket& ket) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,ch,bra,ket);
}

void TwoBodyME::SetTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,ch,bra,ket,tbme);
}

void TwoBodyME::AddToTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   AddToTBME(ch,ch,bra,ket,tbme);
}

double TwoBodyME::GetTBME_J(int j, int a, int b, int c, int d) const
{
   return GetTBME_J(j,j,a,b,c,d);
}

void TwoBodyME::SetTBME_J(int j, int a, int b, int c, int d, double tbme)
{
   SetTBME_J(j,j,a,b,c,d,tbme);
}

void TwoBodyME::AddToTBME_J(int j, int a, int b, int c, int d, double tbme)
{
   AddToTBME_J(j,j,a,b,c,d,tbme);
}

double TwoBodyME::GetTBME_J_norm(int j, int a, int b, int c, int d) const
{
   return GetTBME_J_norm(j,j,a,b,c,d);
}
///
/// Useful for reading in files in isospin formalism
/// \f[
/// \left\langle ab | V | cd \right\rangle_{pnpn} = \frac{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}}{2}
///       \left( \left\langle ab | V | cd \right\rangle_{10} + \left\langle ab | V | cd \right\rangle_{00}   \right)
/// \f]
/// \f[
/// \left\langle ab | V | cd \right\rangle_{pnnp} = \frac{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}}{2}
///       \left( \left\langle ab | V | cd \right\rangle_{10} - \left\langle ab | V | cd \right\rangle_{00}   \right)
/// \f]
///
void TwoBodyME::Set_pn_TBME_from_iso(int j, int T, int tz, int a, int b, int c, int d, double tbme)
{
   // convert everyting to proton labels. Incrementing by 1 gets the neutron label
   a -= a%2;
   b -= b%2;
   c -= c%2;
   d -= d%2;
   int parity = (modelspace->GetOrbit(a).l + modelspace->GetOrbit(b).l)%2;
   int isospin_phase = 2*T-1;
   tbme *= 0.5;
   if (a==b) tbme *= SQRT2;
   if (c==d) tbme *= SQRT2;
   AddToTBME(j,parity,tz,a,b+1,c  ,d+1,tbme);
   if (c!=d)
     AddToTBME(j,parity,tz,a,b+1,c+1,d  ,tbme*isospin_phase);
   if (a!=b and c!=d)
     AddToTBME(j,parity,tz,a+1,b,c+1,d,tbme);
   if (a!=b and (a!=c or b!=d) )
     AddToTBME(j,parity,tz,a+1,b,c,d+1,tbme*isospin_phase);

}

double TwoBodyME::Get_iso_TBME_from_pn(int j, int T, int tz, int a, int b, int c, int d)
{
   // convert everyting to proton labels. Incrementing by 1 gets the neutron label
   a -= a%2;
   b -= b%2;
   c -= c%2;
   d -= d%2;
   int parity = (modelspace->GetOrbit(a).l + modelspace->GetOrbit(b).l)%2;
   int isospin_phase = 2*T-1;
   double tbme = GetTBME(j,parity,tz,a,b+1,c,d+1) + GetTBME(j,parity,tz,a+1,b,c+1,d) + isospin_phase * ( GetTBME(j,parity,tz,a,b+1,c+1,d) + GetTBME(j,parity,tz,a+1,b,c,d+1) );
   tbme *= 0.5;
   if (a==b) tbme /= SQRT2;
   if (c==d) tbme /= SQRT2;
   return tbme;
}


/// Returns an unnormalized monopole-like (angle-averaged) term
/// \f[ \bar{V}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \frac{\sum_{J}(2J+1) V_{ijkl}^J}{(2j_i+1)(2j_j+1)} \f]
double TwoBodyME::GetTBMEmonopole(int a, int b, int c, int d) const
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

double TwoBodyME::GetTBMEmonopole_norm(int a, int b, int c, int d) const
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

      mon += (2*J+1) * GetTBME_norm(J,parityab,Tzab,a,b,c,d);
   }
   mon /= (oa.j2 +1)*(ob.j2+1);
   return mon;
}


double TwoBodyME::GetTBMEmonopole(Ket & bra, Ket & ket) const
{
   return GetTBMEmonopole(bra.p,bra.q,ket.p,ket.q);
}

void TwoBodyME::Erase()
{
//  for ( auto& itmat : MatEl )
  for ( auto& itmat : MtxEl )
  {
//     arma::mat& matrix = itmat.second;
//     matrix.zeros();
    itmat.zeros();
  }
}


double TwoBodyME::Norm() const
{
   double nrm = 0;
//   for ( auto& itmat : MatEl )
   for ( const auto& itmat : MtxEl )
   {
//      const arma::mat& matrix = itmat.second;
//      double n2 = arma::norm(matrix,"fro");
      double n2 = arma::norm(itmat,"fro");
      nrm += n2*n2;
   }
   return sqrt(nrm);
}


void TwoBodyME::Symmetrize()
{
  if (rank_J>0 or rank_T>0 or parity>0) return;
//  for (auto& itmat : MatEl )
  for (auto& itmat : MtxEl )
  {
//      arma::mat& matrix = itmat.second;
//      matrix = arma::symmatu(matrix);
      itmat = arma::symmatu(itmat);
  }
}

void TwoBodyME::AntiSymmetrize()
{
  if (rank_J>0) return;
//  for (auto& itmat : MatEl )
  for (auto& itmat : MtxEl )
  {
//    arma::mat& matrix = itmat.second;
//    matrix = arma::trimatu(matrix) - arma::trimatu(matrix).t();
    itmat = arma::trimatu(itmat) - arma::trimatu(itmat).t();
  }


}


void TwoBodyME::Scale(double x)
{
//   for ( auto& itmat : MatEl )
   for ( auto& itmat : MtxEl )
   {
//      arma::mat& matrix = itmat.second;
//      matrix *= x;
      itmat *= x;
   }

}

void TwoBodyME::Eye()
{
//   for ( auto& itmat : MatEl )
   for ( auto& itmat : MtxEl )
   {
//      arma::mat& matrix = itmat.second;
//      matrix.eye();
     itmat.eye();
   }
}


int TwoBodyME::Dimension()
{
   int dim = 0;
//   for ( auto& itmat : MatEl )
   for ( auto& itmat : MtxEl )
   {
//      int N = itmat.second.n_cols;
      int N = itmat.n_cols;
      dim += N*(N+1)/2;
   }
   return dim;
}

int TwoBodyME::size()
{
  int size=0;
  for ( auto& itmat : MtxEl )
     size += itmat.size();
//  for ( auto& itmat : MatEl )
//     size += itmat.second.size();
  return size*sizeof(double);
}



