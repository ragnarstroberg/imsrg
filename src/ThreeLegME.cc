
#include "ThreeLegME.hh"
//#include "TwoBodyChannel.hh"
#include "PhysicalConstants.hh" // for SQRT2
#include <armadillo>



void ThreeLegME::Allocate()
{
  MatEl.clear();
  size_t nChannels = modelspace->GetNumberTwoBodyChannels();
  for (size_t ch=0; ch<nChannels;++ch)
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    MatEl[{ch}] =  arma::mat(tbc.GetNumberKets(), modelspace->GetNumberOrbits(), arma::fill::zeros);
  }
}

void ThreeLegME::Deallocate()
{
  MatEl.clear();
}


/// ThreeLegME matrix element access functions
double ThreeLegME::GetME(size_t ch, size_t a, size_t b, size_t c) const
{
   double norm = 1;
   if (a==b) norm *= PhysConst::SQRT2;
   return norm * GetME_norm(ch,a,b,c);
}

double ThreeLegME::GetME_norm(size_t ch, size_t a, size_t b, size_t c) const
{
   TwoBodyChannel& tbc =  modelspace->GetTwoBodyChannel(ch);
   auto bra_ind = tbc.GetLocalIndex(std::min(a,b),std::max(a,b));
   if (bra_ind < 0 or bra_ind > tbc.GetNumberKets() )
     return 0;
   Ket & bra = tbc.GetKet(bra_ind);

   double phase = (a>b) ? bra.Phase(tbc.J) : 1;
//   double phase = 1;
//   if (a>b) phase *= bra.Phase(tbc.J);

   return phase * GetMatrix(ch)(bra_ind, c);
}

double ThreeLegME::GetME_J(int J, size_t a, size_t b, size_t c) const
{
   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   int parity_bra = (oa.l+ob.l)%2;
   int Tz_bra = (oa.tz2+ob.tz2)/2;
//   if ( (parity+parity_bra+parity_ket)%2 > 0) return 0;
//   if ( std::abs(Tz_bra-Tz_ket)!=rank_T) return 0;
//   if ( std::abs(j_bra-j_ket) > rank_J) return 0;
//   if ( j_bra + j_ket < rank_J) return 0;
   int ch = modelspace->GetTwoBodyChannelIndex(J,parity_bra,Tz_bra);
   return GetME(ch,a,b,c);
}

// Assume we're passing in a normalized matrix element
void ThreeLegME::SetME(size_t ch, size_t a, size_t b, size_t c, double me)
{
   TwoBodyChannel& tbc =  modelspace->GetTwoBodyChannel(ch);
   auto bra_ind = tbc.GetLocalIndex(std::min(a,b),std::max(a,b));
   if (bra_ind < 0 or bra_ind > tbc.GetNumberKets() )  return;

   Ket & bra = tbc.GetKet(bra_ind);
   double phase = 1;
   if (a>b) phase *= bra.Phase(tbc.J);

   GetMatrix(ch)(bra_ind, c) = me * phase;
}

// Assume we're passing in a normalized matrix element
void ThreeLegME::AddToME(size_t ch, size_t a, size_t b, size_t c, double me)
{
   TwoBodyChannel& tbc =  modelspace->GetTwoBodyChannel(ch);
   auto bra_ind = tbc.GetLocalIndex(std::min(a,b),std::max(a,b));
   if (bra_ind < 0 or bra_ind > tbc.GetNumberKets() )  return;

   Ket & bra = tbc.GetKet(bra_ind);
   double phase = 1;
   if (a>b) phase *= bra.Phase(tbc.J);

   GetMatrix(ch)(bra_ind, c) += me * phase;
}

void ThreeLegME::AddToME_J( int J, size_t a, size_t b, size_t c, double me)
{

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   int parity_bra = (oa.l+ob.l)%2;
   int Tz_bra = (oa.tz2+ob.tz2)/2;
   int ch = modelspace->GetTwoBodyChannelIndex(J,parity_bra,Tz_bra);

   AddToME( ch, a,b,c, me);
}



// Overloaded operators for ThreeLegME

ThreeLegME& ThreeLegME::operator*=(const double rhs)
{
  for ( auto& iter : MatEl ) iter.second *= rhs;
  return *this;
}

ThreeLegME  ThreeLegME::operator*(const double rhs) const
{
  ThreeLegME lhs(*this);
  return lhs *= rhs;
}

ThreeLegME& ThreeLegME::operator+=(const ThreeLegME& rhs)
{
  for (auto& iter : MatEl)  iter.second += rhs.GetMatrix(iter.first);
  return *this;
}

ThreeLegME  ThreeLegME::operator+(const ThreeLegME& rhs) const
{
  ThreeLegME lhs(*this);
  return lhs += rhs;
}

ThreeLegME& ThreeLegME::operator-=(const ThreeLegME& rhs)
{
  for (auto& iter : MatEl)  iter.second -= rhs.GetMatrix(iter.first);
  return *this;
}

ThreeLegME  ThreeLegME::operator-(const ThreeLegME& rhs) const
{
  ThreeLegME lhs(*this);
  return lhs -= rhs;
}


ThreeLegME ThreeLegME::operator-() const
{
  ThreeLegME lhs(*this);
  return lhs *=-1;
}



double ThreeLegME::Norm() const
{
  double norm = 0;
  for (auto& iter : MatEl )
  {
    double n = arma::norm( iter.second, "fro");
    norm += n*n;
  }
  return sqrt(norm);

}

void ThreeLegME::Erase()
{
  for (auto& iter : MatEl) iter.second.zeros();
}
