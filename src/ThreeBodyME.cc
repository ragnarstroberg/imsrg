#include "ThreeBodyME.hh"
//#include "ThreeBodyStorage_pn.hh"
#include "ModelSpace.hh"
#include "ThreeBodyStorage.hh"
#include "ThreeBodyStorage_iso.hh"
#include "IMSRGProfiler.hh"
#include "AngMom.hh"

#include <omp.h>
#include <unordered_map>


// Keeps track of whether any ThreeBodyME has been allocated
// so that we can print out details just on the first time.
//bool ThreeBodyME::none_allocated = true;

//ThreeBodyME::~ThreeBodyME()
//{}

ThreeBodyME::ThreeBodyME()
: threebody_storage(new ThreeBodyStorage_iso()),   modelspace(NULL),E3max(0),herm(1)
{
}

ThreeBodyME::ThreeBodyME(ModelSpace* ms)
: threebody_storage(new ThreeBodyStorage_iso()), modelspace(ms), E3max(ms->E3max), emax(ms->GetEMax3Body()), herm(1)
{}

ThreeBodyME::ThreeBodyME(ModelSpace* ms, int rJ, int rT, int p)
 : threebody_storage(new ThreeBodyStorage_iso(ms,rJ,rT,p)),
 modelspace(ms), E3max(ms->E3max), emax(ms->GetEMax3Body()), herm(1),  rank_J(rJ), rank_T(rT), parity(p)
{
  ISOSPIN_BLOCK_DIMENSION = 5;
  if (rank_T==1) ISOSPIN_BLOCK_DIMENSION = 9; 
  else if (rank_T==3) ISOSPIN_BLOCK_DIMENSION = 1;
}

ThreeBodyME::ThreeBodyME(ModelSpace* ms, int e3max)
: threebody_storage(new ThreeBodyStorage_iso(ms,e3max)), modelspace(ms),E3max(e3max), emax(ms->GetEMax3Body()), herm(1)
{}

ThreeBodyME::ThreeBodyME(ModelSpace* ms, int e3max, int rJ, int rT, int p)
: threebody_storage(new ThreeBodyStorage_iso(ms,e3max,rJ,rT,p)), modelspace(ms),E3max(e3max), emax(ms->GetEMax3Body()),
    herm(1), rank_J(rJ), rank_T(rT), parity(p)
{
  ISOSPIN_BLOCK_DIMENSION = 5;
  if (rank_T==1) ISOSPIN_BLOCK_DIMENSION = 9; 
  else if (rank_T==3) ISOSPIN_BLOCK_DIMENSION = 1;
}

ThreeBodyME::ThreeBodyME(const ThreeBodyME& Tbme)
: threebody_storage( Tbme.threebody_storage->Clone()), modelspace(Tbme.modelspace), 
//:  modelspace(Tbme.modelspace), 
    E3max(Tbme.E3max), emax(Tbme.emax), herm(Tbme.herm),
   rank_J(Tbme.rank_J), rank_T(Tbme.rank_T), parity(Tbme.parity)
{
//   *threebody_storage = *(Tbme.threebody_storage);
}




 ThreeBodyME& ThreeBodyME::operator=(const ThreeBodyME& rhs)
 {
   modelspace = rhs.modelspace;
   E3max = rhs.E3max;
   emax = rhs.emax;
   herm = rhs.herm;
   rank_J = rhs.rank_J;
   rank_T = rhs.rank_T;
   parity = rhs.parity;
//   *threebody_storage = *(rhs.threebody_storage);
   threebody_storage = std::unique_ptr<ThreeBodyStorage>( rhs.threebody_storage->Clone());
   return *this;
 }

 ThreeBodyME& ThreeBodyME::operator*=(const double rhs)
 {
   threebody_storage->Multiply(rhs);
   return *this;
 }

 ThreeBodyME& ThreeBodyME::operator+=(const ThreeBodyME& rhs)
 {
   if ( IsAllocated() and rhs.IsAllocated() )
   {
      threebody_storage->Add(*rhs.threebody_storage);
   }
   else if ( rhs.IsAllocated() )
   {
     threebody_storage = rhs.threebody_storage->Clone();
   }
   return *this;
 }

 ThreeBodyME& ThreeBodyME::operator-=(const ThreeBodyME& rhs)
 {
   if ( IsAllocated()  and  rhs.IsAllocated())
   {
      threebody_storage->Subtract(*rhs.threebody_storage);
   }
   else if ( rhs.IsAllocated() )
   {
     threebody_storage = rhs.threebody_storage->Clone();
     threebody_storage->Multiply(-1);
   }
   return *this;
 }


void ThreeBodyME::Allocate()
{
  if (not threebody_storage->IsAllocated())
  {
    threebody_storage->Allocate();
  }
//  std::cout << "Done calling ThreeBodyME::Allocate() " << std::endl;
}

// Interface setter-getters


//*******************************************************************
/// Get three body matrix element in proton-neutron formalism, regardless of how it is stored.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
ThreeBME_type ThreeBodyME::GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const
{
  return threebody_storage->GetME_pn( Jab_in, Jde_in, twoJ, a, b, c, d, e, f);
}

ThreeBME_type ThreeBodyME::GetME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket) const
{
   return threebody_storage->GetME_pn_ch(ch_bra, ch_ket, ibra, iket);
}

// Get a 3-body matrix element in isospin formalism, regardless of how it is stored.
ThreeBME_type ThreeBodyME::GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const
{
  return threebody_storage->GetME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoTabc, twoTdef, a, b, c, d, e, f);
}

// If we only provide one value of total isospin, assume it's the same T for bra and ket.
ThreeBME_type ThreeBodyME::GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f) const
{
    if (rank_T > 0)
    {
      std::cout << "TROUBLE  !!! " << __FILE__ << "  line " << __LINE__ << std::endl;
    }
    return threebody_storage->GetME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoT, twoT, a, b, c, d, e, f);
}

ThreeBME_type ThreeBodyME::GetME_pn(  int Jab_in, int j0, int Jde_in, int j1, int a, int b, int c, int d, int e, int f) const
{
  return threebody_storage->GetME_pn( Jab_in, j0, Jde_in, j1, a, b, c, d, e, f);
}


void ThreeBodyME::SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    threebody_storage->SetME_pn(Jab_in, Jde_in, twoJ, a, b, c, d, e, f,  V);
}


void ThreeBodyME::SetME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBME_type V)
{
   threebody_storage->SetME_pn_ch(ch_bra, ch_ket, ibra, iket,  V);
}

void ThreeBodyME::SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    threebody_storage->SetME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoTabc, twoTdef, a, b, c, d, e, f,  V);
}

// If we only provide one value of total isospin, assume it's the same T for bra and ket.
void ThreeBodyME::SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    if (rank_T > 0)
    {
      std::cout << "TROUBLE  !!! " << __FILE__ << "  line " << __LINE__ << std::endl;
    }
    threebody_storage->SetME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoT, twoT, a, b, c, d, e, f,  V);
}






void ThreeBodyME::AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    threebody_storage->AddToME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoTabc, twoTdef, a, b, c, d, e, f,  V);
}

void ThreeBodyME::AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBME_type V)
{
   threebody_storage->AddToME_pn_ch(ch_bra, ch_ket, ibra, iket,  V);
}

void ThreeBodyME::AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    if (rank_T > 0)
    {
      std::cout << "TROUBLE  !!! " << __FILE__ << "  line " << __LINE__ << std::endl;
    }
    threebody_storage->AddToME_iso(Jab_in, Jde_in, twoJ, tab_in, tde_in, twoT, twoT, a, b, c, d, e, f,  V);
}

void ThreeBodyME::AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBME_type V)
{
    threebody_storage->AddToME_pn(Jab_in, Jde_in, twoJ, a, b, c, d, e, f,  V);
}







std::vector<double> ThreeBodyME::GetME_pn_TwoOps(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f, const ThreeBodyME& X, const ThreeBodyME& Y) const
{
  return threebody_storage->GetME_pn_TwoOps(Jab, Jde, twoJ, a, b, c, d, e, f, *(X.threebody_storage), *(Y.threebody_storage) );
}



ThreeBME_type ThreeBodyME::GetME_pn_no2b(int a, int b, int c, int d, int e, int f,  int J2b) const
{
  return threebody_storage->GetME_pn_no2b(a,b,c,d,e,f,J2b);
}

ThreeBME_type ThreeBodyME::GetME_pn_mono(int a, int b, int c, int d, int e, int f) const
{
  return threebody_storage->GetME_pn_mono(a,b,c,d,e,f);
}


////////////////////////////////////////////////////////////////////////////////////////////////
/// Begin under-the-hood implementation of the matrix-element setter-getter-adders
//////////////////////////////////////////////////////////////////////////////////////////////









//*****************************************************************************
// We were storing the matrix elements in isospin format.
// Switch to proton-neutron format. This will take a little while, but accessing
// matrix elements will be faster once it's done. So this will potentially be
// useful when doing IMSRG(3) calculations.
// In the future, we can maybe implement some cuts here so that we only
// transform matrix elements that will be used later, or impose some E3max-type cuts.
//*****************************************************************************
void ThreeBodyME::TransformToPN()
{
//  ThreeBodyStorage_pn TBS_pn( modelspace, E3max, rank_J, rank_T, parity  );
//  std::shared_ptr<ThreeBodyStorage> TBS_pn( new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );
  std::unique_ptr<ThreeBodyStorage> TBS_pn( new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );
  TBS_pn->SetHerm( this->herm );

//  TBS_pn.Allocate();
  TBS_pn->Allocate();
  

//  size_t nch = modelspace->GetNumberThreeBodyChannels();
  std::vector<size_t> chbra_list,chket_list;
  for (auto iter : TBS_pn->ch_start)
  {
    chbra_list.push_back(iter.first.ch_bra);
    chket_list.push_back(iter.first.ch_ket);
  }
  size_t nch = chbra_list.size();
//  std::cout << "BEGIN LOOP" << std::endl;

  #pragma omp parallel for schedule(dynamic,1)
  for (size_t ich=0; ich<nch; ich++)
  {
    size_t ch_bra = chbra_list[ich];
    size_t ch_ket = chket_list[ich];
    ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
    ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
    int twoJ = Tbc_bra.twoJ;
    size_t nbras = Tbc_bra.GetNumber3bKets();
    size_t nkets = Tbc_ket.GetNumber3bKets();
    for (size_t ibra=0; ibra<nbras; ibra++)
    {
//      std::cout << "   ibra = " << ibra << std::endl;
      Ket3& bra = Tbc_bra.GetKet(ibra);
//      size_t iketmax = (ch_bra==ch_ket) ? ibra : nkets-1;
      size_t iketmin = (ch_bra==ch_ket) ? ibra : 0;
//      std::cout << "   ibra = " << ibra << "  iketmax = " << iketmax << std::endl;
//      for (size_t iket=0; iket<=ibra; iket++)
      for (size_t iket=iketmin; iket<nkets; iket++)
      {
//        std::cout << "   iket = " << iket << std::endl;
        if ( (ch_bra==ch_ket) and (ibra==iket) and (herm==-1) ) continue;
        Ket3& ket = Tbc_ket.GetKet(iket);
//        std::cout << "Call GetME_pn with isospin mode  JJJ "
//                  << bra.Jpq << " " << ket.Jpq << " " << twoJ << "    abc def "
//                  << bra.p << " " << bra.q << " " << bra.r << "   "
//                  << ket.p << " " << ket.q << " " << ket.r 
//                  << std::endl;
        double me_pn = threebody_storage->GetME_pn( bra.Jpq, ket.Jpq, twoJ,  bra.p, bra.q, bra.r, ket.p, ket.q, ket.r );
//        std::cout << "Call SetME_pn_ch with pn mode " << std::endl;
        TBS_pn->SetME_pn_ch( ch_bra, ch_ket, ibra, iket, me_pn);
//        std::cout << "Set " << TBS_pn->GetME_pn_ch( ch_bra, ch_ket, ibra, iket ) << "  and it should be " << me_pn << std::endl;
      }
    }
  }
//  std::cout << "DONE WITH LOOP" << std::endl;

  // Assign the new ThreeBodyStorage to our shared pointer.
  // The old data will be cleaned up automatically
//  threebody_storage = std::shared_ptr<ThreeBodyStorage_pn>( TBS_pn );
//  threebody_storage = TBS_pn;
  *threebody_storage = *TBS_pn;
//  std::cout << "DONE ASSIGNING" << std::endl;

//  storage_mode = pn;
}






//*****************************************************************************
// This is the one to call if we haven't read in any matrix elements and
// we just want to switch, but don't need to transform from isospin.
//*****************************************************************************
void ThreeBodyME::SwitchToPN_and_discard()
{
  double t_start = omp_get_wtime();
//  threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );
  threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );

  threebody_storage->SetHerm( this->herm );
  threebody_storage->Allocate();
//  threebody_storage = TBS_pn ;
//  threebody_storage = std::shared_ptr<ThreeBodyStorage_pn>( TBS_pn );
//  storage_mode = pn;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}

void ThreeBodyME::SetMode(std::string mode)
{
  double t_start = omp_get_wtime();
  if (mode == "isospin" )
  {
//      threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_iso( modelspace, E3max, rank_J, rank_T, parity)  );
      threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_iso( modelspace, E3max, rank_J, rank_T, parity)  );
  }
  else if ( mode == "pn" )
  {
//      threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );
      threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_pn( modelspace, E3max, rank_J, rank_T, parity)  );
  }
  else if (mode == "no2b" )
  {
//    threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_no2b<ME_single_type>( modelspace, E3max, rank_J, rank_T, parity)  );
    threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_no2b<ME_single_type>( modelspace, E3max, rank_J, rank_T, parity)  );
  }
  else if (mode == "no2bhalf" )
  {
//    threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_no2b<ME_half_type>( modelspace, E3max, rank_J, rank_T, parity)  );
    threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_no2b<ME_half_type>( modelspace, E3max, rank_J, rank_T, parity)  );
  }
  else if (mode == "mono" )
  {
//    threebody_storage = std::shared_ptr<ThreeBodyStorage>(new ThreeBodyStorage_mono<ME_single_type>( modelspace, E3max, rank_J, rank_T, parity)  );
    threebody_storage = std::unique_ptr<ThreeBodyStorage>(new ThreeBodyStorage_mono<ME_single_type>( modelspace, E3max, rank_J, rank_T, parity)  );
  }
  else
  {
    std::cout << " WARNING: I don't know 3N mode " << mode << " .  Continuing with what we had before which is " << threebody_storage->GetStorageMode() << std::endl;
  }
  threebody_storage->SetHerm( this->herm );
  threebody_storage->Allocate();
//  storage_mode = pn;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}

bool ThreeBodyME::IsKetValid( int Jab_in, int twoJ, size_t a_in, size_t b_in, size_t c_in) const
{
   return threebody_storage->IsKetValid(  Jab_in, twoJ, a_in, b_in, c_in);
}

bool ThreeBodyME::IsKetInEMaxTruncations(size_t a_in, size_t b_in, size_t c_in) const
{
   return threebody_storage->IsKetInEMaxTruncations(a_in, b_in, c_in);
}

bool ThreeBodyME::IsOrbitIn3BodyEMaxTruncation(size_t a) const
{
   return threebody_storage->IsOrbitIn3BodyEMaxTruncation(a);
}

bool ThreeBodyME::IsOrbitIn3BodyEMaxTruncation(const Orbit& oa) const
{
   return threebody_storage->IsOrbitIn3BodyEMaxTruncation(oa);
}

size_t ThreeBodyME::GetKetIndex_withRecoupling( int Jab_in, int twoJ, size_t a_in, size_t b_in, size_t c_in, std::vector<size_t>& iket , std::vector<double>& recouple) const
{
   return threebody_storage->GetKetIndex_withRecoupling(  Jab_in, twoJ, a_in, b_in, c_in, iket ,recouple);
}




//*******************************************************************
/// Coefficients for recoupling three body matrix elements.
/// Note that this does not include the -1 factor for an odd
/// permutation of fermionic operators. That is handled in ThreeBodyME::AccessME.
/// Here, we just only with the angular momentum / isospin recoupling factors
//*******************************************************************
double ThreeBodyME::RecouplingCoefficient(ThreeBodyStorage::Permutation recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const
{
  return threebody_storage->RecouplingCoefficient(recoupling_case, ja, jb, jc, Jab_in, Jab, twoJ);
}

int ThreeBodyME::PermutationPhase( ThreeBodyStorage::Permutation recoupling_case ) const
{
  return threebody_storage->PermutationPhase( recoupling_case );
}


void ThreeBodyME::Permute(ThreeBodyStorage::Permutation perm, size_t a_in, size_t b_in, size_t c_in, size_t& a_out, size_t& b_out, size_t& c_out )
{
  threebody_storage->Permute(perm, a_in,b_in,c_in, a_out,b_out,c_out);
}

std::vector<ThreeBodyStorage::Permutation> ThreeBodyME::UniquePermutations( size_t a, size_t b, size_t c ) const
{
  return threebody_storage->UniquePermutations(a,b,c);
}

std::unordered_map<ThreeBodyStorageChannel,size_t,ThreeBodyStorageChannelHash>& ThreeBodyME::Get_ch_start() const
{
  return threebody_storage->ch_start;
}

std::vector<size_t>& ThreeBodyME::Get_ch_dim()
{
  return threebody_storage->ch_dim;
}


void ThreeBodyME::SetHermitian()
{
  herm = +1;
  threebody_storage->SetHerm(herm);
}

void ThreeBodyME::SetAntiHermitian()
{
  herm = -1;
  threebody_storage->SetHerm(herm);
}

void ThreeBodyME::SetE3max(int e)
{
  E3max = e;
  threebody_storage->SetE3max(E3max);
}
void ThreeBodyME::Setemax(int e)
{
  emax=  e;
  threebody_storage->SetEmax(emax);
}


double ThreeBodyME::Norm() const
{
   return threebody_storage->Norm();
}


void ThreeBodyME::Erase()
{
  threebody_storage->Erase();
}

void ThreeBodyME::Deallocate()
{
  threebody_storage->Deallocate();
}


size_t ThreeBodyME::size()
{
  return threebody_storage->size();
}

void ThreeBodyME::WriteBinary(std::ofstream& f)
{
  threebody_storage->WriteBinary(f);
}

void ThreeBodyME::ReadBinary(std::ifstream& f)
{
  threebody_storage->ReadBinary(f);
}


void ThreeBodyME::WriteFile(std::vector<std::string> StringInputs, std::vector<int> IntInputs ) const
{
  threebody_storage->WriteFile( StringInputs, IntInputs );
}

void ThreeBodyME::ReadFile( std::vector<std::string> StringInputs, std::vector<int> IntInputs )
{
  threebody_storage->ReadFile( StringInputs, IntInputs );
}


std::string ThreeBodyME::GetStorageMode() const
{
  return threebody_storage->GetStorageMode();
}



void ThreeBodyME::Print() const
{
   threebody_storage->Print();
}
