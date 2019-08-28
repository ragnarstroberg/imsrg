
#include "ThreBodyMEpn.hh"







ThreeBodyMEpn::ThreeBodyMEpn(ModelSpace* ms)
 : modelspace(ms)
{}

ThreeBodyME::ThreeBodyMEpn(ModelSpace* ms, int e3max)
: modelspace(ms), E3max(e3max)
{
  Allocate();
}


ThreeBodyME::ThreeBodyMEpn(const ThreeBodyME& tbme)
 : modelspace(tbme.modelspace), MatEl(tbme.MatEl), emax(tbme.emax), E3max(tbme.E3max),
   herm(tbme.herm)
{
  Allocate();
}





// This will need to be more elaborate if we want to use tensor 3-body.
void ThreeBodyME::Allocate()
{
  total_dimension = 0;
  for (auto ch : modelspace->ThreeBodyChannels )
  {
    size_t channel_dim = modelspace->GetNumberKet3s( ch );
    MatEl[{ch,ch}] = SymmMatrix<ME_type>( channel_dim );
    total_dimension += channel_dim * (channel_dim+1)/2;
  }
}





/// These are the ones that eventually get called, but typically the other methods will
/// be more convenient to call.
ME_type ThreeBodyMEpn::GetME_pn(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket) const
{
  return MatEl[{ch_bra,ch_ket}](ibra,iket);
}

// We have this here in case we want to set a matrix element, but we store it in a different
// coupling order. In that case, we need to add to multiple matrix elements with the appropriate
// recoupling coefficients included (see below).
void ThreeBodyMEpn::AddToME_pn(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type matel)
{
  MatEl[{ch_bra,ch_ket}](ibra,iket) += matel;
}

void ThreeBodyMEpn::SetME_pn(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ME_type matel)
{
  MatEl[{ch_bra,ch_ket}](ibra,iket) = matel;
}







void ThreeBodyMEpn::AddToME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type me_add )
{

  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> ibra;
  std::vector<size_t> iket;
  size_t ch_bra, ch_ket;
  size_t ch_bra = GetKetIndex_withRecoupling( twoJ, Jab, a,b,c, ibra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( twoJ, Jde, d,e,f, iket, recouple_ket );
  if ( ch_bra != ch_ket) return 0;

  double me_out = 0;
  for ( int i=0; i<ibra.size(); i++)
  {
    for (int j=0; j<iket.size(); j++)
    {
       AddToME_pn( ch_bra, ch_ket, ibra[i], iket[i], recouple_bra[i] * recouple_ket[i] * me_add );
    }
  }
}


// This does the recoupling twice, so it could certainly be more efficient.
// But it's easier to code and understand this way, so I'll change it if need be.
void ThreeBodyMEpn::SetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ME_type me_set )
{
  double me_previous = GetMEpn( Jab_in, Jde_in, twoJ, a,b,c,d,e,f);
  AddToME_pn( Jab_in, Jde_in, twoJ, a,b,c,d,e,f,  me_set-me_previous );
}



//// Let's not implement this until we actually need it...
////
//ME_type ThreBodyMEpn::GetME_pn(size_t ch_bra, size_t ch_ket, Ket3& bra, Ket3& ket) const
//{
//  ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
//  ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
//  size_t ibra = Tbc_bra.GetLocalIndex( bra );
//  size_t iket = Tbc_bra.GetLocalIndex( ket );
//  return GetME_pn( chbra, chket, ibra, iket );
//}






ME_type ThreBodyMEpn::GetME_pn(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f) const
{

  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> ibra;
  std::vector<size_t> iket;
  size_t ch_bra, ch_ket;
  size_t ch_bra = GetKetIndex_withRecoupling( twoJ, Jab, a,b,c, ibra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( twoJ, Jde, d,e,f, iket, recouple_ket );
  if ( ch_bra != ch_ket) return 0;

  double me_out = 0;
  for ( int i=0; i<ibra.size(); i++)
  {
    for (int j=0; j<iket.size(); j++)
    {
      me_out += recouple_bra[i] * recouple_ket[i] * GetME_pn( ch_bra, ch_ket, ibra[i], iket[i] );
    }
  }

  return me_out;
}







size_t GetKetIndex_withRecoupling( int twoJ, int Jab_in, size_t a_in, size_t b_in, size_t c_in, std::vector<size_t>& iket , std::vector<double>& recouple)
{

  int recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);

  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);
  int parity = ( oa.l + ob.l + oc.l )%2;
  int twoTz = ( oa.tz2 + ob.tz2 + oc.tz2 );
  int ch_bra = modelspace->GetThreeBodyChannelIndex( twoJ, parity_bra, twoTz_bra );

  int Jab_min = std::max( std::abs(oa.j2-ob.j2), std::abs(oc.j2-twoJ) )/2;
  int Jab_max = std::min( oa.j2+ob.j2, oc.j2+twoJ)/2;

  double ja = oa.j2*0.5;
  double jb = ob.j2*0.5;
  double jc = oc.j2*0.5;

  // Loop over possible values of Jab with the new ordering of a,b,c and
  // fill a vector of the index of where each of those |a,b,c,Jab> states lives
  // as well as the recouplng coefficient.
  for (int Jab=Jab_min; Jab<=Jab; Jab++)
  {
    iket.push_back( modelspace->GetThreeBodyChannel(ch_bra).GetLocalIndex( a,b,c,Jab ) );
    recouple.push_back( RecouplingCoefficient( recoupling_case, ja,jb,jc, Jab_in, Jab, twoJ );
  }

  return ch_bra;

}






// Define some constants for the various permutations of three indices
// for use in RecouplingCoefficient and SortOrbits
const int ThreeBodyME::ABC = 0;
const int ThreeBodyME::BCA = 1;
const int ThreeBodyME::CAB = 2;
const int ThreeBodyME::ACB = 3;
const int ThreeBodyME::BAC = 4;
const int ThreeBodyME::CBA = 5;



//*******************************************************************
/// Coefficients for recoupling three body matrix elements.
/// Note that this does not include the -1 factor for an odd
/// permutation of fermionic operators. That is handled in ThreeBodyME::AccessME.
/// Here, we just only with the angular momentum / isospin recoupling factors
//*******************************************************************
double ThreeBodyME::RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const
{
   if ( std::abs(int(ja-jb))>Jab  or int(ja+jb)<Jab) return 0;
   if ( std::abs(int(jc-twoJ/2.))>Jab  or int(jc+twoJ/2.)<Jab) return 0;
   switch (recoupling_case)
   {
    case ABC: return Jab==Jab_in ? 1 : 0;
    case BCA: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    case CAB: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case ACB: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case BAC: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case CBA: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    default: return 0;
    }
}




//*******************************************************************
/// Rearrange orbits (abc) so that a>=b>=c
/// and return an int which reflects the required reshuffling
/// - 0: (abc)_in -> (abc)
/// - 1: (bca)_in -> (abc)
/// - 2: (cab)_in -> (abc)
/// - 3: (acb)_in -> (abc)  -- odd permutation
/// - 4: (bac)_in -> (abc)  -- odd permutation
/// - 5: (cba)_in -> (abc)  -- odd permutation
//*******************************************************************
int ThreeBodyMEpn::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
{
   a_in -= a_in%2;
   b_in -= b_in%2;
   c_in -= c_in%2;
   a=a_in;
   b=b_in;
   c=c_in;
   if (a<b)  std::swap(a,b);
   if (b<c)  std::swap(b,c);
   if (a<b)  std::swap(a,b);

   int recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? ABC : ACB;
   else if (a_in==b)  recoupling_case = (b_in==a) ? BAC : BCA;
   else               recoupling_case = (b_in==a) ? CAB : CBA;

   return recoupling_case;
}
  
