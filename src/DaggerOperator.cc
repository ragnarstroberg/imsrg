
#include "DaggerOperator.hh"
#include "PhysicalConstants.hh" // for SQRT2


DaggerOperator::DaggerOperator()
{}



DaggerOperator::DaggerOperator(ModelSpace& ms)
: Operator(ms) , ThreeLeg(ms)
{
  SetNumberLegs(3);
  SetNonHermitian();
}

DaggerOperator::DaggerOperator(ModelSpace& ms, size_t Qorb)
: Operator(ms),  ThreeLeg(ms)
{
  SetNumberLegs(3);
  SetNonHermitian();
  SetQSpaceOrbit(Qorb);
}





/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////

DaggerOperator& DaggerOperator::operator=(const DaggerOperator& rhs) = default;
//DaggerOperator& DaggerOperator::operator=(DaggerOperator&& rhs) = default;


// multiply operator by a scalar
DaggerOperator& DaggerOperator::operator*=(const double rhs)
{
   ZeroBody *= rhs;
   OneBody *= rhs;
   ThreeLeg *= rhs;
   return *this;
}

DaggerOperator DaggerOperator::operator*(const double rhs) const
{
   DaggerOperator opout(*this);
   opout *= rhs;
   return opout;
}

// Add non-member operator so we can multiply an operator
// by a scalar from the lhs, i.e. s*O = O*s
DaggerOperator operator*(const double lhs, const DaggerOperator& rhs)
{
   return rhs * lhs;
}
DaggerOperator operator*(const double lhs, const DaggerOperator&& rhs)
{
   return rhs * lhs;
}


// divide operator by a scalar
DaggerOperator& DaggerOperator::operator/=(const double rhs)
{
   return *this *=(1.0/rhs);
}

DaggerOperator DaggerOperator::operator/(const double rhs) const
{
   DaggerOperator opout = DaggerOperator(*this);
   opout *= (1.0/rhs);
   return opout;
}

// Add operators
DaggerOperator& DaggerOperator::operator+=(const DaggerOperator& rhs)
{
   ZeroBody += rhs.ZeroBody;
   OneBody  += rhs.OneBody;
   if (rhs.GetParticleRank() > 1)
     ThreeLeg  += rhs.ThreeLeg;
   return *this;
}

DaggerOperator DaggerOperator::operator+(const DaggerOperator& rhs) const
{
   if (GetParticleRank() >= rhs.GetParticleRank())
     return ( DaggerOperator(*this) += rhs );
   else
     return ( DaggerOperator(rhs) += *this );
}



// Subtract operators
DaggerOperator& DaggerOperator::operator-=(const DaggerOperator& rhs)
{
   OneBody -= rhs.OneBody;
   if (rhs.GetParticleRank() > 1)
     ThreeLeg -= rhs.ThreeLeg;
   return *this;
}

DaggerOperator DaggerOperator::operator-(const DaggerOperator& rhs) const
{
   return ( DaggerOperator(*this) -= rhs );
}


// Negation operator
DaggerOperator DaggerOperator::operator-() const
{
   return (*this)*-1.0;
}









/// ThreeLegME matrix element access functions



void ThreeLegME::Allocate()
{
  MatEl.clear();
  size_t nChannels = modelspace->GetNumberTwoBodyChannels();
  for (size_t ch=0; ch<nChannels;++ch)
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//    if ( std::abs(tbc_bra.J-tbc_ket.J)>rank_J ) continue;
//    if ( (tbc_bra.J+tbc_ket.J)<rank_J ) continue;
//    if ( std::abs(tbc_bra.Tz-tbc_ket.Tz)!=rank_T ) continue; // we don't couple to T, so rank_T really means |delta Tz|
//    if ( (tbc_bra.parity + tbc_ket.parity + parity)%2>0 ) continue;
    MatEl[{ch}] =  arma::mat(tbc.GetNumberKets(), modelspace->GetNumberOrbits(), arma::fill::zeros);
  }

}


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




double DaggerOperator::ThreeLegNorm()
{
  double norm = 0;
  for (auto& iter : ThreeLeg.MatEl )
  {
    double n = arma::norm( iter.second, "fro");
    norm += n*n;
  }
  return sqrt(norm);

}


void DaggerOperator::Erase()
{

  OneBody.zeros();
  for (auto& iter : ThreeLeg.MatEl) iter.second.zeros();

}



