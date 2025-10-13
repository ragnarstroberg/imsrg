
#include "TwoBodyME.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh" // for SQRT2
//#ifndef SQRT2
//  #define SQRT2 1.4142135623730950488
//#endif

// destructor defined for debugging purposes
TwoBodyME::~TwoBodyME()
{}

TwoBodyME::TwoBodyME()
: modelspace(NULL), nChannels(0), hermitian(true),antihermitian(false),allocated(false),
  rank_J(0), rank_T(0), parity(0)
{
//  cout << "Default TwoBodyME constructor" << endl;
}


TwoBodyME::TwoBodyME(ModelSpace* ms)
: modelspace(ms), nChannels(ms->GetNumberTwoBodyChannels()),
  hermitian(true), antihermitian(false),allocated(false), rank_J(0), rank_T(0), parity(0)
{
  Allocate();
}


TwoBodyME::TwoBodyME(ModelSpace* ms, int rJ, int rT, int p)
: modelspace(ms), nChannels(ms->GetNumberTwoBodyChannels()),
  hermitian(true), antihermitian(false), allocated(false), rank_J(rJ), rank_T(rT), parity(p)
{
  Allocate();
}


 TwoBodyME& TwoBodyME::operator*=(const double rhs)
 {
   for ( auto& itmat : MatEl )
   {
      itmat.second *= rhs;
   }
   return *this;
 }

 TwoBodyME TwoBodyME::operator*(const double rhs) const
 {
   TwoBodyME out(*this);
   out *=rhs;
   return out;
 }

 TwoBodyME& TwoBodyME::operator+=(const TwoBodyME& rhs)
 {
   if (rhs.IsAllocated() )
   {
     if (not this->IsAllocated() )
     {
        *this = rhs;
     }
     else
     {
       for ( auto& itmat : MatEl )
       {
          auto ch_bra = itmat.first[0];
          auto ch_ket = itmat.first[1];
          itmat.second += rhs.GetMatrix(ch_bra,ch_ket);
       }
     }
   }
   return *this;
 }

 TwoBodyME& TwoBodyME::operator-=(const TwoBodyME& rhs)
 {
   if (not this->IsAllocated() )
   {
      *this = rhs * (-1.0);
   }
   else
   {
      for ( auto& itmat : rhs.MatEl )
      {
         auto ch_bra = itmat.first[0];
         auto ch_ket = itmat.first[1];
         GetMatrix(ch_bra,ch_ket) -= itmat.second;
      }
   }
   return *this;
 }



void TwoBodyME::Allocate()
{
  nChannels = modelspace->GetNumberTwoBodyChannels();
  MatEl.clear();
  for (size_t ch_bra=0; ch_bra<nChannels;++ch_bra)
  {
     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
     for (size_t ch_ket=ch_bra; ch_ket<nChannels;++ch_ket)
     {
        TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
        if ( not AngMom::Triangle(tbc_bra.J,tbc_ket.J,rank_J) ) continue;
//        if ( std::abs(tbc_bra.J-tbc_ket.J)>rank_J ) continue;
//        if ( (tbc_bra.J+tbc_ket.J)<rank_J ) continue;
        if ( std::abs(tbc_bra.Tz-tbc_ket.Tz)!=rank_T ) continue; // we don't couple to T, so rank_T really means |delta Tz|
        if ( (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;
        MatEl[{ch_bra,ch_ket}] =  arma::mat(tbc_bra.GetNumberKets(), tbc_ket.GetNumberKets(), arma::fill::zeros);
     }
  }
  allocated = true;
}

void TwoBodyME::Deallocate()
{
  MatEl.clear();
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


bool TwoBodyME::IsAllocated()const
{
   return allocated;
}


/// This returns the matrix element times a factor \f$ \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \f$
double TwoBodyME::GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const
{
  double norm = 1;
   if (a==b) norm *= PhysConst::SQRT2;
   if (c==d) norm *= PhysConst::SQRT2;
   return norm * GetTBME_norm(ch_bra,ch_ket,a,b,c,d);
}

/// This returns the normalized matrix element 
double TwoBodyME::GetTBME_norm(int ch_bra, int ch_ket, int a, int b, int c, int d) const
{
   if (not allocated) return 0;
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   auto bra_ind = tbc_bra.GetLocalIndex(std::min(a,b),std::max(a,b));
   auto ket_ind = tbc_ket.GetLocalIndex(std::min(c,d),std::max(c,d));

   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc_bra.GetNumberKets() or ket_ind > tbc_ket.GetNumberKets() )
   {
      return 0;
   }
     
   Ket & bra = tbc_bra.GetKet(bra_ind);
   Ket & ket = tbc_ket.GetKet(ket_ind);

   double phase = 1;
   if (a>b) phase *= bra.Phase(tbc_bra.J);
   if (c>d) phase *= ket.Phase(tbc_ket.J);
   if (ch_bra > ch_ket)
   {
     return hermitian ?   phase * modelspace->phase(tbc_bra.J-tbc_ket.J) * GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind)
                      : - phase * modelspace->phase(tbc_bra.J-tbc_ket.J) * GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind);
   }
   return phase * GetMatrix(ch_bra,ch_ket)(bra_ind, ket_ind);
}

void TwoBodyME::SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   if (not allocated)  return;
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   int bra_ind = tbc_bra.GetLocalIndex(std::min(a,b),std::max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(std::min(c,d),std::max(c,d));
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);

// new lines suggested by Takayuki Miyagi
   if( ch_bra > ch_ket )                        
   {                               
     std::swap(ch_bra, ch_ket);     
     std::swap(bra_ind, ket_ind);     
     phase *= modelspace->phase(tbc_bra.J-tbc_ket.J);   
     if (antihermitian) phase *= -1; //additional fix from Takayuki
   }       
// end new lines

   GetMatrix(ch_bra,ch_ket)(bra_ind,ket_ind) = phase * tbme;
   if (ch_ket != ch_bra) return;
   if (hermitian and ch_bra==ch_ket) GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind) = phase * tbme;
   if (antihermitian and ch_bra==ch_ket) GetMatrix(ch_ket,ch_bra)(ket_ind,bra_ind) = - phase * tbme;
}


void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
//   cout << "Tzbra = " << tbc_bra.Tz << "   Tzket = " << tbc_ket.Tz << "  rank_T = " << rank_T << endl;
   int bra_ind = tbc_bra.GetLocalIndex(std::min(a,b),std::max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(std::min(c,d),std::max(c,d));
//   cout << "bra_ind = " << bra_ind << " = local index of " << min(a,b) << " " << max(a,b) << " from channel " << ch_bra << endl;
//   cout << "ket_ind = " << ket_ind << " = local index of " << min(c,d) << " " << max(c,d) << " from channel " << ch_ket << endl;
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   if (ch_bra > ch_ket)
   {
    std::swap(ch_bra,ch_ket);
//    swap(tbc_bra,tbc_ket);
    std::swap(bra_ind,ket_ind);
    phase *= modelspace->phase(tbc_bra.J-tbc_ket.J);
   }
//   cout << "Getting Matrix " << ch_bra << "," << ch_ket << "(" << bra_ind << "," << ket_ind
//        << "), dimension = " << GetMatrix(ch_bra,ch_ket).n_rows << "x" << GetMatrix(ch_bra,ch_ket).n_cols << endl;
   GetMatrix(ch_bra,ch_ket)(bra_ind,ket_ind) += phase * tbme;
   if (ch_bra!=ch_ket or ket_ind==bra_ind) return;
//   if (ch_ket != ch_bra) return;
   if (hermitian and ch_bra==ch_ket) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) += phase * tbme;
   if (antihermitian and ch_bra==ch_ket) GetMatrix(ch_bra,ch_ket)(ket_ind,bra_ind) -=  phase * tbme;
}

void TwoBodyME::AddToTBMENonHermNonNormalized(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
//   cout << "Tzbra = " << tbc_bra.Tz << "   Tzket = " << tbc_ket.Tz << "  rank_T = " << rank_T << endl;
   int bra_ind = tbc_bra.GetLocalIndex(std::min(a,b),std::max(a,b));
   int ket_ind = tbc_ket.GetLocalIndex(std::min(c,d),std::max(c,d));
//   cout << "bra_ind = " << bra_ind << " = local index of " << min(a,b) << " " << max(a,b) << " from channel " << ch_bra << endl;
//   cout << "ket_ind = " << ket_ind << " = local index of " << min(c,d) << " " << max(c,d) << " from channel " << ch_ket << endl;
   double phase = 1;
   if (a>b) phase *= tbc_bra.GetKet(bra_ind).Phase(tbc_bra.J);
   if (c>d) phase *= tbc_ket.GetKet(ket_ind).Phase(tbc_ket.J);
   if (a==b) phase *= 1.0 / PhysConst::SQRT2;
   if (c==d) phase *= 1.0 / PhysConst::SQRT2;
//   cout << "Getting Matrix " << ch_bra << "," << ch_ket << "(" << bra_ind << "," << ket_ind
//        << "), dimension = " << GetMatrix(ch_bra,ch_ket).n_rows << "x" << GetMatrix(ch_bra,ch_ket).n_cols << endl;
   GetMatrix(ch_bra,ch_ket)(bra_ind,ket_ind) += phase * tbme;
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
   if (IsHermitian() and ch_bra==ch_ket)
      GetMatrix(ch_ket,ch_bra)(iket,ibra) = tbme;
   else if(IsAntiHermitian() and ch_bra==ch_ket)
      GetMatrix(ch_ket,ch_bra)(iket,ibra) = -tbme;
}
void TwoBodyME::AddToTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme)
{
   if (ch_bra>ch_ket)
   {
     std::swap(ch_bra,ch_ket);
     std::swap(ibra,iket);
     tbme *= modelspace->phase( modelspace->GetTwoBodyChannel(ch_bra).J - modelspace->GetTwoBodyChannel(ch_ket).J);
   }
   GetMatrix(ch_bra,ch_ket)(ibra,iket) += tbme;

   if (ch_bra==ch_ket and ibra!=iket)
   {
     if (IsHermitian() and ch_bra==ch_ket)
       GetMatrix(ch_bra,ch_ket)(iket,ibra) += tbme;
     else if(IsAntiHermitian() and ch_bra==ch_ket)
       GetMatrix(ch_bra,ch_ket)(iket,ibra) -= tbme;
   }
}

void TwoBodyME::AddToTBMENonHerm(int ch_bra, int ch_ket, int ibra, int iket, double tbme)
{
   GetMatrix(ch_bra,ch_ket)(ibra,iket) += tbme;
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
   if ( std::abs(Tz_bra-Tz_ket)!=rank_T) return 0;
   if ( not AngMom::Triangle(j_bra,j_ket,rank_J) ) return 0;
   if ( not AngMom::Triangle(oa.j2,ob.j2,2*j_bra) ) return 0;
   if ( not AngMom::Triangle(oc.j2,od.j2,2*j_ket) ) return 0;
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,parity_bra,Tz_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,parity_ket,Tz_ket);
   return GetTBME(ch_bra,ch_ket,a,b,c,d);
}
void TwoBodyME::SetTBME_J(int j_bra, int j_ket, int a, int b, int c, int d, double tbme)
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   // std::cout<<"here"<<std::endl;
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
   if ( std::abs(Tz_bra-Tz_ket)!=rank_T) return 0;
   if ( not AngMom::Triangle(j_bra,j_ket,rank_J) ) return 0;
//   if ( std::abs(j_bra-j_ket) > rank_J) return 0;
//   if ( j_bra + j_ket < rank_J) return 0;
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,parity_bra,Tz_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,parity_ket,Tz_ket);
   return GetTBME_norm(ch_bra,ch_ket,a,b,c,d);
}

// This is intended for use in the Pandya transformation needed for the 222 ph commutator, but it might be useful in other cases.
// The idea is that we want to access the same element from two different operators, so we do all the phases and index lookup only once.
void TwoBodyME::GetTBME_J_norm_twoOps(const TwoBodyME& OtherTBME, int j_bra, int j_ket, int a, int b, int c, int d, double& tbme_this, double& tbme_other) const
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   int parity_bra = (oa.l+ob.l)%2;
   int parity_ket = (oc.l+od.l)%2;
   int Tz_bra = (oa.tz2+ob.tz2)/2;
   int Tz_ket = (oc.tz2+od.tz2)/2;
   tbme_this  =0;
   tbme_other =0;
   if ( (parity+parity_bra+parity_ket)%2 > 0) return;
   if ( std::abs(Tz_bra-Tz_ket)!=rank_T) return;
   if ( not AngMom::Triangle(j_bra,j_ket,rank_J) ) return;
//   if ( std::abs(j_bra-j_ket) > rank_J) return;
//   if ( j_bra + j_ket < rank_J) return;
   int ch_bra = modelspace->GetTwoBodyChannelIndex(j_bra,parity_bra,Tz_bra);
   int ch_ket = modelspace->GetTwoBodyChannelIndex(j_ket,parity_ket,Tz_ket);

   TwoBodyChannel& tbc_bra =  modelspace->GetTwoBodyChannel(ch_bra);
   TwoBodyChannel& tbc_ket =  modelspace->GetTwoBodyChannel(ch_ket);
   auto bra_ind = tbc_bra.GetLocalIndex(std::min(a,b),std::max(a,b));
   auto ket_ind = tbc_ket.GetLocalIndex(std::min(c,d),std::max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc_bra.GetNumberKets() or ket_ind > tbc_ket.GetNumberKets() )
   {
     return;
   }
   Ket & bra = tbc_bra.GetKet(bra_ind);
   Ket & ket = tbc_ket.GetKet(ket_ind);

   double phase = 1;
   if (a>b) phase *= bra.Phase(tbc_bra.J);
   if (c>d) phase *= ket.Phase(tbc_ket.J);
   if (ch_bra > ch_ket) // Phase for conjugating a reduced matrix element.
   {
     phase *=  modelspace->phase(tbc_bra.J-tbc_ket.J) ;
     std::swap(ch_bra,ch_ket);
     std::swap(bra_ind,ket_ind);
   }
   tbme_this =  phase * GetMatrix(ch_bra,ch_ket)(bra_ind, ket_ind);
   tbme_other =  phase * OtherTBME.GetMatrix(ch_bra,ch_ket)(bra_ind, ket_ind);
}


void TwoBodyME::GetTBME_J_twoOps(const TwoBodyME& OtherTBME, int j_bra, int j_ket, int a, int b, int c, int d, double& tbme_this, double& tbme_other) const
{
   GetTBME_J_norm_twoOps(OtherTBME, j_bra,j_ket, a,b,c,d, tbme_this, tbme_other);
   if (a==b)
   {
      tbme_this *= PhysConst::SQRT2;
      tbme_other *= PhysConst::SQRT2;
   }
   if (c==d)
   {
      tbme_this *= PhysConst::SQRT2;
      tbme_other *= PhysConst::SQRT2;
   }
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
   if (a==b) tbme *= PhysConst::SQRT2;
   if (c==d) tbme *= PhysConst::SQRT2;
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
   if (a==b) tbme /= PhysConst::SQRT2;
   if (c==d) tbme /= PhysConst::SQRT2;
   return tbme;
}

/*
/// Returns an unnormalized monopole-like (angle-averaged) term
/// \f[ \bar{V}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \frac{\sum_{J}(2J+1) V_{ijkl}^J}{(2j_i+1)(2j_j+1)} \f]
double TwoBodyME::GetTBMEmonopole(int a, int b, int c, int d) const
{
  double norm = 1;
   if (a==b) norm *= PhysConst::SQRT2;
   if (c==d) norm *= PhysConst::SQRT2;
   return norm * GetTBMEmonopole_norm(a,b,c,d);



//   double mon = 0;
//   Orbit &oa = modelspace->GetOrbit(a);
//   Orbit &ob = modelspace->GetOrbit(b);
//   Orbit &oc = modelspace->GetOrbit(c);
//   Orbit &od = modelspace->GetOrbit(d);
//   int Tzab = (oa.tz2 + ob.tz2)/2;
//   int parityab = (oa.l + ob.l)%2;
//   int Tzcd = (oc.tz2 + od.tz2)/2;
//   int paritycd = (oc.l + od.l)%2;
//
////   if (Tzab != Tzcd or parityab != paritycd) return 0;
//
//   int jmin = std::abs(oa.j2 - ob.j2)/2;
//   int jmax = (oa.j2 + ob.j2)/2;
//   
//   for (int J=jmin;J<=jmax;++J)
//   {
//      size_t chab = modelspace->GetTwoBodyChannelIndex(J,parityab,Tzab);
//      size_t chcd = modelspace->GetTwoBodyChannelIndex(J,paritycd,Tzcd);
//      mon += (2*J+1) * GetTBME(chab,chcd,a,b,c,d);
////      mon += (2*J+1) * GetTBME(J,parityab,Tzab,a,b,c,d);
//   }
//   mon /= (oa.j2 +1)*(ob.j2+1);
//   return mon;
}
*/

/// \f[ \bar{V}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \frac{\sum_{J}(2J+1) V_{ijkl}^J}{(2j_i+1)(2j_j+1)} \f]
double TwoBodyME::GetTBMEmonopole(int a, int b, int c, int d) const
{

  double norm = 1;
  if (a==b) norm *= PhysConst::SQRT2;
  if (c==d) norm *= PhysConst::SQRT2;
  return norm * GetTBMEmonopole_norm( a,b,c,d);
//   double mon = 0;
//   Orbit &oa = modelspace->GetOrbit(a);
//   Orbit &ob = modelspace->GetOrbit(b);
//   Orbit &oc = modelspace->GetOrbit(c);
//   Orbit &od = modelspace->GetOrbit(d);
//   int Tzab = (oa.tz2 + ob.tz2)/2;
//   int parityab = (oa.l + ob.l)%2;
//   int Tzcd = (oc.tz2 + od.tz2)/2;
//   int paritycd = (oc.l + od.l)%2;
//
////   if (Tzab != Tzcd or parityab != paritycd) return 0;
//
//
//   int jmin = std::abs(oa.j2 - ob.j2)/2;
//   int jmax = (oa.j2 + ob.j2)/2;
//   
//   for (int J=jmin;J<=jmax;++J)
//   {
//      size_t chab = modelspace->GetTwoBodyChannelIndex(J,parityab,Tzab);
//      size_t chcd = modelspace->GetTwoBodyChannelIndex(J,paritycd,Tzcd);
////      mon += (2*J+1) * GetTBME(J,parityab,Tzab,a,b,c,d);
//      mon += (2*J+1) * GetTBME(chab,chcd,a,b,c,d);
////      mon += (2*J+1) * GetTBME_norm(chab,chcd,a,b,c,d);
////      std::cout << " :::::  J= " << J << " " << GetTBME(chab,chcd,a,b,c,d) << "  =>    mon = " << mon << std::endl;
//   }
//   mon /= (oa.j2 +1)*(ob.j2+1);
////   std::cout << " ;;;;;  mon /= " << (oa.j2 +1)*(ob.j2+1) << "   => " << mon << std::endl;
//   return mon;
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

//   if (Tzab != Tzcd or parityab != paritycd) return 0;
   if (  (parityab + paritycd + parity)%2 !=0 ) return 0;
   if ( std::abs( Tzab - Tzcd) > rank_T ) return 0;

   int jmin = std::abs(oa.j2 - ob.j2)/2;
   int jmax = (oa.j2 + ob.j2)/2;
   
   for (int J=jmin;J<=jmax;++J)
   {
      size_t chab = modelspace->GetTwoBodyChannelIndex(J,parityab,Tzab);
      size_t chcd = modelspace->GetTwoBodyChannelIndex(J,paritycd,Tzcd);
      mon += (2*J+1) * GetTBME_norm(chab,chcd,a,b,c,d);
//      mon += (2*J+1) * GetTBME_norm(J,parityab,Tzab,a,b,c,d);
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
   if (not allocated) return 0;
   for ( auto& itmat : MatEl )
   {
      const arma::mat& matrix = itmat.second;
      int Jbra = modelspace->GetTwoBodyChannel( itmat.first[0] ).J;
      int Jket = modelspace->GetTwoBodyChannel( itmat.first[1] ).J;
      int degeneracy = (2*Jket+1) * (std::min(Jbra,Jket+rank_J) - std::max(-Jbra,Jket-rank_J)+1);
//      int degeneracy = 1;//(2*Jket+1) * (std::min(Jbra,Jket+rank_J) - std::max(-Jbra,Jket-rank_J)+1);
//      double n2 = arma::norm(matrix,"fro");
      double n2 = arma::norm(matrix,"fro") * degeneracy;
//      nrm += n2*n2;
      // If bra and ket are different channels, we only store
      // one ordering. For the norm, we then need a factor of 2.
      nrm += (itmat.first[0]==itmat.first[1]) ? n2*n2 : 2*n2*n2;
//      if (itmat.first[0] != itmat.first[1])
//         nrm += n2*n2;
   }
   return sqrt(nrm);
}


// arma::symmatu reflects the upper triangle to the lower triangle
// the upper triangle are the elements Mat(ibra,iket) with ibra<=iket
void TwoBodyME::Symmetrize()
{
  if (rank_J>0 or rank_T>0 or parity>0) return;
  for (auto& itmat : MatEl )
  {
      arma::mat& matrix = itmat.second;
      matrix = arma::symmatu(matrix);
  }
}

// arma::trimatu takes the upper triangle of a matrix
// we then set the lower triangle to the negtive of the upper triangle
void TwoBodyME::AntiSymmetrize()
{
  if (rank_J>0) return;
  for (auto& itmat : MatEl )
  {
    arma::mat& matrix = itmat.second;
    matrix = arma::trimatu(matrix) - arma::trimatu(matrix).t();
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

void TwoBodyME::PrintMatrix(size_t chbra,size_t chket) const
{
   std::cout.precision(12); 
   MatEl.at({chbra,chket}).raw_print();
//   MatEl.at({chbra,chket}).print();
}

void TwoBodyME::PrintAllMatrices() const
{
  for ( auto& itmat : MatEl )
  {
    
    arma::uvec subscript = itmat.second.is_empty()  ?  arma::uvec({0,0}) 
                           :   arma::ind2sub( arma::size(itmat.second),  arma::abs(itmat.second).index_max() ) ; // get row,column of maximum entry
    std::cout << "ch_bra, ch_ket : " << itmat.first[0] << " " << itmat.first[1] << "     norm = " << arma::norm( itmat.second, "fro")
              << "  max entry at ( " << subscript(0) << " , " << subscript(1) << " ) ";
    if ( not itmat.second.is_empty() )  std::cout << "     " << itmat.second(subscript(0),subscript(1));
    std::cout << std::endl  << itmat.second << std::endl << std::endl;

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

int TwoBodyME::size()
{
  int size=0;
  for ( auto& itmat : MatEl )
     size += itmat.second.size();
  return size*sizeof(double);
}



void TwoBodyME::WriteBinary( std::ofstream& of )
{
  of.write((char*)&nChannels,sizeof(nChannels));
  of.write((char*)&hermitian,sizeof(hermitian));
  of.write((char*)&antihermitian,sizeof(antihermitian));
  of.write((char*)&rank_J,sizeof(rank_J));
  of.write((char*)&rank_T,sizeof(rank_T));
  of.write((char*)&parity,sizeof(parity));
  for ( auto& itmat : MatEl )
    of.write((char*)itmat.second.memptr(),itmat.second.size()*sizeof(double));

}


void TwoBodyME::ReadBinary( std::ifstream& of )
{
  of.read((char*)&nChannels,sizeof(nChannels));
  of.read((char*)&hermitian,sizeof(hermitian));
  of.read((char*)&antihermitian,sizeof(antihermitian));
  of.read((char*)&rank_J,sizeof(rank_J));
  of.read((char*)&rank_T,sizeof(rank_T));
  of.read((char*)&parity,sizeof(parity));
  Allocate();
  for ( auto& itmat : MatEl )
    of.read((char*)itmat.second.memptr(),itmat.second.size()*sizeof(double));

}



// non-member operator overloads


TwoBodyME operator+(const TwoBodyME& lhs, const TwoBodyME& rhs)
{
  TwoBodyME TBout(lhs);
  TBout += rhs;
  return TBout;
}

TwoBodyME operator-(const TwoBodyME& lhs, const TwoBodyME& rhs)
{
  TwoBodyME TBout(lhs);
  TBout -= rhs;
  return TBout;
}


TwoBodyME operator*(const double lhs, const TwoBodyME& rhs)
{
  return rhs * lhs;
}
