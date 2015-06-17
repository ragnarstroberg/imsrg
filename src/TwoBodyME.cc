
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

TwoBodyME::TwoBodyME(const TwoBodyME& tbme)
: modelspace(tbme.modelspace), MatEl(tbme.MatEl), nChannels(tbme.nChannels),
  hermitian(tbme.hermitian), antihermitian(tbme.antihermitian),
  rank_J(tbme.rank_J), rank_T(tbme.rank_T), parity(tbme.parity)
{}


 TwoBodyME& TwoBodyME::operator=(const TwoBodyME& rhs)
 {
   Copy(rhs);
   return *this;
 }

 TwoBodyME& TwoBodyME::operator*=(const double rhs)
 {
   for ( auto& itmat : MatEl )
   {
      itmat.second *= rhs;
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator+=(const TwoBodyME& rhs)
 {
   for ( auto& itmat : MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      itmat.second += rhs.GetMatrix(ch_bra,ch_ket);
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator-=(const TwoBodyME& rhs)
 {
   for ( auto& itmat : rhs.MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
//      itmat.second -= rhs.GetMatrix(ch_bra,ch_ket);
      GetMatrix(ch_bra,ch_ket) -= itmat.second;
   }
   return *this;
 }


void TwoBodyME::Copy(const TwoBodyME& rhs)
{
   modelspace = rhs.modelspace;
   MatEl      = rhs.MatEl;
   nChannels  = rhs.nChannels;
   hermitian  = rhs.hermitian;
   antihermitian  = rhs.antihermitian;
   rank_J     = rhs.rank_J;
   rank_T     = rhs.rank_T;
   parity     = rhs.parity;

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
        
        MatEl[{ch_bra,ch_ket}] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
     }
  }
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

/// This returns the matrix element times a factor \f$ \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \f$
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
   if (hermitian) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) = phase * tbme;
   if (antihermitian) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) = - phase * tbme;
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

void TwoBodyME::SetTBME(int ch_bra, int ch_ket, int iket, int ibra, double tbme)
{
   GetMatrix(ch_bra,ch_ket)(ibra,iket) = tbme;
   if (IsHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) = tbme;
   else if(IsAntiHermitian())
      GetMatrix(ch_bra,ch_ket)(iket,ibra) = -tbme;
}
void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, int iket, int ibra, double tbme)
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
/// \f[ \bar{V}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \sum_{J}(2J+1) V_{ijkl}^J \f]
///
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
  for ( auto& itmat : MatEl )
  {
     arma::mat& matrix = itmat.second;
     matrix.zeros();
  }
}


double TwoBodyME::Norm() const
{
   double nrm = 0;
   for ( auto& itmat : MatEl )
   {
      const arma::mat& matrix = itmat.second;
      double n2 = arma::norm(matrix,"fro");
      nrm += n2*n2;
   }
   return sqrt(nrm);
}


void TwoBodyME::Symmetrize()
{
  for (auto& itmat : MatEl )
  {
    arma::mat& matrix = itmat.second;
    matrix = arma::symmatu(matrix);
  }
}

void TwoBodyME::AntiSymmetrize()
{
  for (auto& itmat : MatEl )
  {
    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];
    int nbras = modelspace->GetTwoBodyChannel(ch_bra).GetNumberKets();
    int nkets = modelspace->GetTwoBodyChannel(ch_ket).GetNumberKets();
    arma::mat& matrix = itmat.second;
    for (int ibra=0; ibra<nbras; ++ibra)
    {
      for (int iket=ibra+1; iket<nkets; ++iket)
      {
         matrix(iket,ibra) = -matrix(ibra,iket);
      }
    }
  }


}


void TwoBodyME::Scale(double x)
{
   for ( auto& itmat : MatEl )
   {
      arma::mat& matrix = itmat.second;
      matrix *= x;
   }

}

void TwoBodyME::Eye()
{
   for ( auto& itmat : MatEl )
   {
      arma::mat& matrix = itmat.second;
      matrix.eye();
   }
}


int TwoBodyME::Dimension()
{
   int dim = 0;
   for ( auto& itmat : MatEl )
   {
      int N = itmat.second.n_cols;
      dim += N*(N+1)/2;
   }
   return dim;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TwoBodyME_ph::TwoBodyME_ph( TwoBodyME& TBME )
: TwoBodyME(TBME)
{
  Allocate();
  


}





void TwoBodyME_ph::DoPandyaTransformation(TwoBodyME& TwoBody_pp, string option)
{
//   int n_nonzero = modelspace->SortedTwoBodyChannels_CC.size();
//   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
/*
   for ( auto& itmat : MatEl )
//   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      arma::uvec& kets_ph = tbc_cc.GetKetIndex_ph();
      int nph_kets = kets_ph.n_rows;
      int J_cc = tbc_cc.J;

      TwoBody_CC_hp[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);
      TwoBody_CC_ph[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);

      // loop over cross-coupled ph bras <ac| in this channel
      for (int ibra=0; ibra<nph_kets; ++ibra)
      {
         Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
         int a = bra_cc.p;
         int b = bra_cc.q;
         Orbit & oa = modelspace->GetOrbit(a);
         Orbit & ob = modelspace->GetOrbit(b);
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;

         // loop over cross-coupled kets |cd> in this channel
         // we go to 2*nKets to include |cd> and |dc>
         for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
         {
            Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
            int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
            int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
            Orbit & oc = modelspace->GetOrbit(c);
            Orbit & od = modelspace->GetOrbit(d);
            double jc = oc.j2*0.5;
            double jd = od.j2*0.5;


            int jmin = max(abs(ja-jd),abs(jc-jb));
            int jmax = min(ja+jd,jc+jb);
            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
               if (abs(sixj) < 1e-8) continue;
               double tbme = TwoBody.GetTBME_J(J_std,a,d,c,b);
               sm -= (2*J_std+1) * sixj * tbme ;
            }
            TwoBody_CC_hp[ch_cc](ibra,iket_cc) = sm;


            // Exchange (a <-> b) to account for the (n_a - n_b) term
            // Get Tz,parity and range of J for <bd || ca > coupling
            jmin = max(abs(jb-jd),abs(jc-ja));
            jmax = min(jb+jd,jc+ja);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(jb,ja,J_cc,jc,jd,J_std);
               if (abs(sixj) < 1e-8) continue;
               double tbme = TwoBody.GetTBME_J(J_std,b,d,c,a);
               sm -= (2*J_std+1) * sixj * tbme ;
            }
            TwoBody_CC_ph[ch_cc](ibra,iket_cc) = sm;

         }
      }
   }
*/
}







void TwoBodyME_ph::Allocate()
{
  for (int ch_bra=0; ch_bra<nChannels;++ch_bra)
  {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
     for (int ch_ket=ch_bra; ch_ket<nChannels;++ch_ket)
     {
        TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
        if ( abs(tbc_bra.J-tbc_ket.J)>rank_J ) continue;
        if ( (tbc_bra.J+tbc_ket.J)<rank_J ) continue;
        if ( abs(tbc_bra.Tz)==abs(tbc_bra.Tz) ) continue; // this is the only difference
        if ( (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;
        
        MatEl[{ch_bra,ch_ket}] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
     }
  }
}








