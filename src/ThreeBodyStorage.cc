

#include "ThreeBodyStorage.hh"



ThreeBodyStorage::ThreeBodyStorage()
{}

ThreeBodyStorage::ThreeBodyStorage(ModelSpace* ms)
 : modelspace(ms), emax(ms->GetEMax3Body()), E2max(ms->GetE3max()), E3max(ms->GetE3max()), lmax(ms->GetLmax())
{}

ThreeBodyStorage::ThreeBodyStorage(ModelSpace* ms, int e3max)
 : modelspace(ms),  emax(ms->GetEMax3Body()), E2max(ms->GetE2max()), E3max(ms->GetE3max()), lmax(ms->GetLmax())
{}

ThreeBodyStorage::ThreeBodyStorage( const ThreeBodyStorage& TBS_in )
: modelspace(TBS_in.modelspace), emax(TBS_in.emax), E2max(TBS_in.E2max), E3max(TBS_in.E3max), lmax(TBS_in.lmax), herm(TBS_in.herm), 
 rank_J(TBS_in.rank_J), rank_T(TBS_in.rank_T), parity(TBS_in.parity), ISOSPIN_BLOCK_DIMENSION(TBS_in.ISOSPIN_BLOCK_DIMENSION),
 is_allocated(TBS_in.is_allocated), ch_start(TBS_in.ch_start), ch_dim(TBS_in.ch_dim)
{}

ThreeBodyStorage::ThreeBodyStorage(ModelSpace* ms, int rJ, int rT, int p)
 : modelspace(ms),  emax(ms->GetEMax3Body()), E2max(ms->GetE2max()), E3max(ms->GetE3max()), lmax(ms->GetLmax()), rank_J(rJ), rank_T(rT), parity(p)
{}

ThreeBodyStorage::ThreeBodyStorage(ModelSpace* ms, int e3max , int rJ, int rT, int p)
 : modelspace(ms), emax(ms->GetEMax3Body()), E2max(ms->GetE2max()), E3max(e3max), lmax(ms->GetLmax()), rank_J(rJ), rank_T(rT), parity(p)
{}





// In general, not all storage types can/must implement all ways of accessing 3-body matrix elments.
// If we try to call a method that's not implemented, we make some noise and die.
// This will presumably be due to a mistake in the code, but if we really want to access the matrix
// element in that particular way using that particular storage mode, then we'd better go ahead and implement it.
//
//ThreeBodyStorage::ME_type ThreeBodyStorage::NotImplemented() const
ThreeBodyStorage::ME_type ThreeBodyStorage::NotImplemented(std::string funcname) const
{
  if ( not is_allocated ) // if we're not allocated, then just return zero. No harm, no foul.
  return 0;

  std::cout << "ERROR: " << funcname << " is not yet implemented for ThreeBodyStorage mode [" << GetStorageMode()  << "]. Dying now." << std::endl;
  std::exit(EXIT_FAILURE);
  return 0;
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
//int ThreeBodyME::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
ThreeBodyStorage::Permutation ThreeBodyStorage::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
{
   a=a_in;
   b=b_in;
   c=c_in;
   if (a>b)  std::swap(a,b); // the a <= b <= c ordering matches with the 2body storage, which have if (q<p) continue;
   if (b>c)  std::swap(b,c);
   if (a>b)  std::swap(a,b);
//   if (a<b)  std::swap(a,b);
//   if (b<c)  std::swap(b,c);
//   if (a<b)  std::swap(a,b);

//   int recoupling_case;
   Permutation recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? ABC : ACB;
   else if (a_in==b)  recoupling_case = (b_in==a) ? BAC : BCA;
   else               recoupling_case = (b_in==a) ? CAB : CBA;

   return recoupling_case;
}





//*******************************************************************
/// Coefficients for recoupling three body matrix elements.
/// Note that this does not include the -1 factor for an odd
/// permutation of fermionic operators. That is given by ThreeBodyStorage::PermutationPhase. 
/// Here, we just treat the angular momentum / isospin recoupling factors
//*******************************************************************
double ThreeBodyStorage::RecouplingCoefficient(ThreeBodyStorage::Permutation recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const
{
   if ( std::abs(int(ja-jb))>Jab  or int(ja+jb)<Jab) return 0;
   if ( std::abs(int(jc-twoJ/2.))>Jab  or int(jc+twoJ/2.)<Jab) return 0;
   switch (recoupling_case)
   {
    case ABC: return Jab==Jab_in ? 1 : 0;
    case BAC: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case BCA: return modelspace->phase( jb+jc+Jab_in+1)   * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    case CAB: return modelspace->phase( ja+jb-Jab+1)      * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case ACB: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case CBA: return    -1                                * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    default: return 0;// shouldn't ever get here?
    }
}

/// If we make an odd permutation of indices, we get a fermionic minus sign.
int ThreeBodyStorage::PermutationPhase( ThreeBodyStorage::Permutation recoupling_case ) const
{
  return  (recoupling_case==ABC or recoupling_case==BCA or recoupling_case==CAB) ?  +1 : -1;
}


/// Apply a permutation to a set of indices a,b,c
void ThreeBodyStorage::Permute( Permutation perm, size_t a_in, size_t b_in, size_t c_in, size_t& a_out, size_t& b_out, size_t& c_out )
{
  switch (perm)
  {
    case ABC: a_out=a_in; b_out=b_in; c_out=c_in; break;
    case BAC: a_out=b_in; b_out=a_in; c_out=c_in; break;
    case BCA: a_out=b_in; b_out=c_in; c_out=a_in; break;
    case CAB: a_out=c_in; b_out=a_in; c_out=b_in; break;
    case ACB: a_out=a_in; b_out=c_in; c_out=b_in; break;
    case CBA: a_out=c_in; b_out=b_in; c_out=a_in; break;
  }
}


std::vector<ThreeBodyStorage::Permutation> ThreeBodyStorage::UniquePermutations( size_t a, size_t b, size_t c ) const
{
  std::vector<ThreeBodyStorage::Permutation> unique_perms = {ABC};
  if ( a!=b ) unique_perms.push_back( BAC );
  if ( a!=c ) unique_perms.push_back( CBA );
  if ( b!=c ) unique_perms.push_back( ACB );
  if ( a!=b and a!=c and b!=c)
  {
    unique_perms.push_back( BCA );
    unique_perms.push_back( CAB );
  }
//  if ( a!=b and a!=c and b!=c) unique_perms.push_back( CAB );
  return unique_perms;
}


bool ThreeBodyStorage::IsKetValid( int Jab, int twoJ, size_t a, size_t b, size_t c) const 
{
  if (!IsKetInEMaxTruncations(a, b, c)) {
    return false;
  }

  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);

  int jj_a = oa.j2;
  int jj_b = ob.j2;
  int jj_c = oc.j2;

  int Jab_min = std::abs((jj_a - jj_b) / 2);
  int Jab_max = (jj_a + jj_b) / 2;

  int twoJ_min = std::abs(2 * Jab - jj_c);
  int twoJ_max = 2 * Jab + jj_c;

  // Check against j_a, j_b coupling range.
  if ((Jab > Jab_max) || (Jab < Jab_min)) {
    return false;
  }

  // Check against J_ab, j_c coupling range.
  if ((twoJ > twoJ_max) || (twoJ < twoJ_min)) {
    return false;
  }

  return true;
}

bool ThreeBodyStorage::IsKetInEMaxTruncations(size_t a, size_t b, size_t c) const 
{
  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);

  int e_a = oa.n * 2 + oa.l;
  int e_b = ob.n * 2 + ob.l;
  int e_c = oc.n * 2 + oc.l;

  int E3 = e_a + e_b + e_c;

  // Check against emax 3-body cut.
  if ((e_a > modelspace->GetEMax3Body()) 
      || (e_b > modelspace->GetEMax3Body())
      || (e_c > modelspace->GetEMax3Body())) {
        return false;
  }

  // Check against E3max cut.
  if (E3 > modelspace->GetE3max()) {
    return false;
  }

  return true;
}

bool ThreeBodyStorage::IsOrbitIn3BodyEMaxTruncation(size_t a) const 
{
  const Orbit& oa = modelspace->GetOrbit(a);

  return IsOrbitIn3BodyEMaxTruncation(oa);
}

bool ThreeBodyStorage::IsOrbitIn3BodyEMaxTruncation(const Orbit& oa) const 
{

  int e_a = oa.n * 2 + oa.l;

  // Check against emax 3-body cut.
  return e_a <= modelspace->GetEMax3Body();
}

size_t ThreeBodyStorage::GetKetIndex_withRecoupling( int Jab_in, int twoJ, size_t a_in, size_t b_in, size_t c_in, std::vector<size_t>& iket , std::vector<double>& recouple) const
{
  int a,b,c;
  if (!IsKetInEMaxTruncations(a_in, b_in, c_in)) {
    std::cout << "Warning: Accessing matrix element that is 0 by truncations.\n";
    std::cout << a_in << ", " << b_in << ", " << c_in << "\n";
  }
  Permutation recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);

  int permutation_phase = PermutationPhase( recoupling_case );

  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);
  int parity = ( oa.l + ob.l + oc.l )%2;
  int twoTz = ( oa.tz2 + ob.tz2 + oc.tz2 );
  int ch = modelspace->GetThreeBodyChannelIndex( twoJ, parity, twoTz );
  if ( (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l) > modelspace->E3max) return ch;

  if (ch < 0 ) return ch;
  auto Tbc = modelspace->GetThreeBodyChannel(ch);

  int Jab_min = std::max( std::abs(oa.j2-ob.j2), std::abs(oc.j2-twoJ) )/2;
  int Jab_max = std::min( oa.j2+ob.j2, oc.j2+twoJ)/2;

  if (  ( a_in==b_in and (Jab_in%2)>0 ) 
     or ( a_in==b_in and (Jab_in > (modelspace->GetOrbit(a_in).j2-1)) )
     or ( a_in==b_in and a_in==c_in and ( twoJ > (3*oa.j2-3)) )
     or ( twoJ==(oa.j2+ob.j2+oc.j2) and (a==b or a==c or b==c) ) )
  {
    Jab_max = Jab_min-1;
  }

  double ja = oa.j2*0.5;
  double jb = ob.j2*0.5;
  double jc = oc.j2*0.5;

  // Loop over possible values of Jab with the new ordering of a,b,c and
  // fill a vector of the index of where each of those |a,b,c,Jab> states lives
  // as well as the recouplng coefficient.
  for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
  {
    if (a==b and (Jab%2)>0) continue;
    size_t index = Tbc.GetLocalIndex( a,b,c,Jab );
    if ( index == size_t(-1) ) continue;

    double coefficient =  permutation_phase * RecouplingCoefficient( recoupling_case, ja,jb,jc, Jab_in, Jab, twoJ );
    if (std::abs(coefficient)<1e-10) continue;

    iket.push_back( index );
    recouple.push_back( coefficient );
  }

  return ch;

}




std::vector<ThreeBodyStorage::ME_type> ThreeBodyStorage::GetME_pn_TwoOps(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f, const ThreeBodyStorage& X, const ThreeBodyStorage& Y) const
{
  std::vector<double> me_out( 2, 0.0 );
  if (!IsKetValid(Jab, twoJ, a, b, c) || !IsKetValid(Jde, twoJ, d, e, f)) return me_out;
  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> ibra;
  std::vector<size_t> iket;
  size_t ch_bra = GetKetIndex_withRecoupling( Jab, twoJ, a,b,c, ibra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( Jde, twoJ, d,e,f, iket, recouple_ket );
  if ( ch_bra != ch_ket) return me_out;
  //TODO: Should we also throw an exception if twoJ is even?

  for ( size_t i=0; i<ibra.size(); i++)
  {
    for (size_t j=0; j<iket.size(); j++)
    {
        me_out[0] += recouple_bra[i] * recouple_ket[j] * X.GetME_pn_ch( ch_bra, ch_ket, ibra[i], iket[j] );
        me_out[1] += recouple_bra[i] * recouple_ket[j] * Y.GetME_pn_ch( ch_bra, ch_ket, ibra[i], iket[j] );

    }
  }
  return me_out;
}




ThreeBodyStorage::ME_type ThreeBodyStorage::GetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2b, int twoT) const 
{
   Orbit& oc = modelspace->GetOrbit(c);
   int twoJmin = std::abs( 2*J2b - oc.j2 );
   int twoJmax = 2*J2b + oc.j2; 
   ME_type me_no2b = 0;
   if ( rank_J==0 and rank_T==0 and parity==0)
   {
     for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
     {
       me_no2b += GetME_iso(J2b,J2b,twoJ, Tab,Tde,twoT,twoT, a,b,c,d,e,f) * (twoJ+1.);
     }
   }
   else
   {
     for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
     {
       me_no2b += GetME_iso(J2b,J2b,twoJ, Tab,Tde,twoT,twoT, a,b,c,d,e,f) * sqrt(twoJ+1.);
     }
   }
   return me_no2b;
}


// Get ME with Jab=Jde=J2b, summed over total J with weight 2J+1
// For most storage modes which implement GetME_pn, this is all that's needed.
// For the NO2B storage mode, this will be overloaded.
ThreeBodyStorage::ME_type ThreeBodyStorage::GetME_pn_no2b(int a, int b, int c, int d, int e, int f,  int J2b) const
{
   Orbit& oc = modelspace->GetOrbit(c);
   int twoJmin = std::abs( 2*J2b - oc.j2 );
   int twoJmax = 2*J2b + oc.j2; 
   ME_type me_no2b = 0;

   if ( rank_J==0 and rank_T==0 and parity==0)
   {
     for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
     {
       me_no2b += GetME_pn(J2b,J2b,twoJ, a,b,c,d,e,f) * (twoJ+1.);
     }
   }
   else
   {
     for (int twoJ=twoJmin; twoJ<=twoJmax; twoJ+=2)
     {
       me_no2b += GetME_pn(J2b,J2b,twoJ, a,b,c,d,e,f) * sqrt(twoJ+1.);
     }
   }
   return me_no2b;
}

ThreeBodyStorage::ME_type ThreeBodyStorage::GetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int T3) const 
{
    int j2a = modelspace->GetOrbit(a).j2;
    int j2b = modelspace->GetOrbit(b).j2;
//    int j2c = modelspace->GetOrbit(c).j2;
    double v = 0;

    int j2min =  std::abs(j2a-j2b) /2;
    int j2max = (j2a+j2b)/2;
    for (int j2=j2min; j2<=j2max; ++j2)
    {
      v += GetME_iso_no2b(a,b,c,Tab,d,e,f,Tde,j2,T3);
   }
//   v /= j2c+1.0;
   return v;
}

ThreeBodyStorage::ME_type ThreeBodyStorage::GetME_pn_mono(int a, int b, int c, int d, int e, int f) const
{
    int j2a = modelspace->GetOrbit(a).j2;
    int j2b = modelspace->GetOrbit(b).j2;
//    int j2c = modelspace->GetOrbit(c).j2;
    double v = 0;

    int j2min =  std::abs(j2a-j2b) /2;
    int j2max = (j2a+j2b)/2;
    for (int j2=j2min; j2<=j2max; ++j2)
    {
      v += GetME_pn_no2b(a,b,c,d,e,f,j2);
   }
//   v /= j2c+1.0;
   return v;
}



