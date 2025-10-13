
#include "AngMom.hh"
#include "ThreeBodyStorage.hh"
#include "ThreeBodyStorage_pn.hh"



#include <omp.h>

bool ThreeBodyStorage_pn::none_allocated = true;
int ThreeBodyStorage_pn::number_allocated = 0;


ThreeBodyStorage_pn::ThreeBodyStorage_pn( const ThreeBodyStorage_pn& TBS_in )
//: ThreeBodyStorage( TBS_in ), MatEl(TBS_in.MatEl), total_dimension(0)
: ThreeBodyStorage( TBS_in ), MatEl(TBS_in.MatEl), total_dimension(TBS_in.total_dimension)
{
  if (is_allocated)   number_allocated++;
}


ThreeBodyStorage_pn::~ThreeBodyStorage_pn( )
{
  if (is_allocated)  number_allocated--;
}

//std::shared_ptr<ThreeBodyStorage> ThreeBodyStorage_pn::Clone() const { return std::shared_ptr<ThreeBodyStorage>( new ThreeBodyStorage_pn( *this)); };
std::unique_ptr<ThreeBodyStorage> ThreeBodyStorage_pn::Clone() const { return std::unique_ptr<ThreeBodyStorage>( new ThreeBodyStorage_pn( *this)); };


  void ThreeBodyStorage_pn::Multiply(const double rhs) 
  {
    for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] *= rhs;
  }

  void ThreeBodyStorage_pn::Add(const ThreeBodyStorage& rhs) 
  {
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_pn*)(&rhs))->MatEl;
       for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] += rhsMatEl[i];
    }
    else
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << "  allocated?  " << this->IsAllocated() << "  " << rhs.IsAllocated()  << "  norms " << this->Norm() << "  " << rhs.Norm() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void ThreeBodyStorage_pn::Subtract(const ThreeBodyStorage& rhs)
  {
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
      auto& rhsMatEl = ((ThreeBodyStorage_pn*)(&rhs))->MatEl;
      for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] -= rhsMatEl[i];
    }
    else
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

//  ThreeBodyStorage_pn& ThreeBodyStorage_pn::operator+=(const ThreeBodyStorage_pn& rhs)
//  {
//    for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] += rhs.MatEl[i];
//  }
//
//  ThreeBodyStorage_pn& ThreeBodyStorage_pn::operator-=(const ThreeBodyStorage_pn& rhs)
//  {
//    for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] -= rhs.MatEl[i];
//  }


void ThreeBodyStorage_pn::Allocate()
{
//  std::cout << "BEGIN " << __func__ << "  in  " << __FILE__ << std::endl;
  double tstart = omp_get_wtime();
  MatEl.clear();
  total_dimension = 0;
  size_t nch = modelspace->GetNumberThreeBodyChannels();
  for (size_t ch_bra=0; ch_bra<nch; ch_bra++)
  {
    ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel( ch_bra );
    size_t nkets_bra = Tbc_bra.GetNumber3bKets(); // Number of kets in this 3body J,p,Tz channel
    ch_dim.push_back( nkets_bra );
    for (size_t ch_ket=ch_bra; ch_ket<nch; ch_ket++)
    {
      ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel( ch_ket );
      if (  ( std::abs(Tbc_bra.twoJ-Tbc_ket.twoJ)<=2*rank_J ) and ( (Tbc_bra.twoJ+Tbc_ket.twoJ)>=2*rank_J )
          and ( (Tbc_bra.parity+Tbc_ket.parity)%2==parity ) and ( std::abs(Tbc_bra.twoTz-Tbc_ket.twoTz)==2*rank_T )  )
      {
         size_t nkets_ket = Tbc_ket.GetNumber3bKets(); // Number of kets in this 3body J,p,Tz channel
         if (nkets_ket > 0) // If there aren't any kets, don't add it to the list.
            ch_start[{ch_bra,ch_ket}] = total_dimension;
         if (ch_bra==ch_ket)
         {
           total_dimension += nkets_bra * (nkets_bra+1)/2; // only need to store half the matrix
         }
         else
         {
           total_dimension += nkets_bra * nkets_ket;
         }
      }
    }
  }
//  std::cout << "In ThreeBodyStorage_pn::Allocate()  total dimension = " << total_dimension 
//            << "  -> " << total_dimension*sizeof(pnME_type) / (1024.*1024.*1024.) << " GB"  << std::endl;
//  std::cout << "  number of ThreeBodyStorage_pn currently allocated: " << CountAllocations() << "  ->  " << CountAllocations() * total_dimension * sizeof(pnME_type) / (1024.*1024.*1024.) << " GB " <<  std::endl;
  MatEl.resize(total_dimension,0.0);
  if (none_allocated)
  {
     std::cout << "DONE ALLOCATING PN 3-body, size of MatEl is " << MatEl.size()
               << "  ->  " <<MatEl.size()*sizeof(pnME_type) / (1024.*1024.*1024.) << " GB" << std::endl;
  }
//  std::cout << " ... allocation successful." << std::endl;
  is_allocated = true;
  none_allocated = false;
  number_allocated ++;
//  PN_mode = true;
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - tstart;
}






//////////////////////////////////////////////////////////////////////////////////////
//////   Interface methods
///////////////////////////////////////////////////////////////////////////////////

ThreeBodyStorage::ME_type ThreeBodyStorage_pn::GetME_pn(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const
{
  if (!IsKetValid(Jab_in, twoJ, a, b, c) || !IsKetValid(Jde_in, twoJ, d, e, f)) return 0.0;
  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> index_bra;
  std::vector<size_t> index_ket;
  size_t ch_bra = GetKetIndex_withRecoupling( Jab_in, twoJ, a,b,c, index_bra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( Jde_in, twoJ, d,e,f, index_ket, recouple_ket );
  if ( rank_J==0 and rank_T==0 and parity==0 and  (ch_bra != ch_ket) ) return 0;
  if (ch_bra==-1 or ch_ket==-1) return 0;
  ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
  ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
  if ( (Tbc_bra.parity + Tbc_ket.parity)%2 != parity) return 0;
  if ( std::abs( Tbc_bra.twoTz - Tbc_ket.twoTz) != 2*rank_T ) return 0;
  //TODO: Should we also throw an exception if twoJ is even?



  double me_out = 0;
  for ( size_t I=0; I<index_bra.size(); I++)
  {
    for (size_t J=0; J<index_ket.size(); J++)
    {
      me_out += recouple_bra[I] * recouple_ket[J] * GetME_pn_ch( ch_bra, ch_ket, index_bra[I], index_ket[J] );
    }
  }

  return me_out;
}

// tensor getter
ThreeBodyStorage::ME_type ThreeBodyStorage_pn::GetME_pn(  int Jab_in, int j0, int Jde_in, int j1, int a, int b, int c, int d, int e, int f) const
{
  if (!IsKetValid(Jab_in, j0, a, b, c) || !IsKetValid(Jde_in, j1, d, e, f)) return 0.0;
  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> index_bra;
  std::vector<size_t> index_ket;
  size_t ch_bra = GetKetIndex_withRecoupling( Jab_in, j0, a,b,c, index_bra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( Jde_in, j1, d,e,f, index_ket, recouple_ket );
  if ( rank_J==0 and rank_T==0 and parity==0 and  (ch_bra != ch_ket) ) return 0;
  if (ch_bra==-1 or ch_ket==-1) return 0;
  if ( std::abs(j0 - j1)  > rank_J * 2 or  (j0 + j1) < rank_J * 2  ) return 0;
  ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
  ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
  if ( (Tbc_bra.parity + Tbc_ket.parity)%2 != parity) return 0;
  if ( std::abs( Tbc_bra.twoTz - Tbc_ket.twoTz) != 2*rank_T ) return 0;
  //TODO: Should we also throw an exception if j0 and j1 are even?

  double me_out = 0;
  for ( size_t I=0; I<index_bra.size(); I++)
  {
    for (size_t J=0; J<index_ket.size(); J++)
    {
      me_out += recouple_bra[I] * recouple_ket[J] * GetME_pn_ch( ch_bra, ch_ket, index_bra[I], index_ket[J] );
    }
  }

  return me_out;
}






ThreeBodyStorage::ME_type ThreeBodyStorage_pn::GetME_iso( int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const
{
  double me_iso = 0;
//  int twoTz =  twoT; // it should be independent of Tz, so we just pick one
  int twoTz =  1; // it should be independent of Tz, so we just pick one
  if (a==b and (Jab_in+tab_in)%2==0) return 0;
  if (d==e and (Jde_in+tde_in)%2==0) return 0;
  for (int tz2a : {-1,1} )
  {
    for (int tz2b : {-1,1} )
    {
      double clebsch_ab = AngMom::CG(0.5,0.5*tz2a, 0.5,0.5*tz2b, tab_in, 0.5*(tz2a+tz2b) );
      if ( std::abs(clebsch_ab)<1e-7) continue;
      for (int tz2c : {-1,1} )
      {
        if ( (tz2a + tz2b + tz2c) != twoTz) continue;
        double clebsch_abc = AngMom::CG(tab_in, 0.5*(tz2a+tz2b), 0.5, 0.5*tz2c, 0.5*twoTabc, 0.5*twoTz );
        if ( std::abs(clebsch_abc)<1e-7) continue;
        for (int tz2d : {-1,1} )
        {
          for (int tz2e : {-1,1} )
          {
            int tz2f = twoTz - tz2d - tz2e;
            if (std::abs(tz2f)!=1) continue;
            double clebsch_de = AngMom::CG(0.5,0.5*tz2d, 0.5,0.5*tz2e, tde_in, 0.5*(tz2d+tz2e) );
            if ( std::abs(clebsch_de)<1e-7) continue;
            double clebsch_def = AngMom::CG(tde_in, 0.5*(tz2d+tz2e), 0.5, 0.5*tz2f, 0.5*twoTdef, 0.5*twoTz );
            if ( std::abs(clebsch_def)<1e-7) continue;
            size_t a_pn = 2*(a/2) + (tz2a+1)/2;
            size_t b_pn = 2*(b/2) + (tz2b+1)/2;
            size_t c_pn = 2*(c/2) + (tz2c+1)/2;
            size_t d_pn = 2*(d/2) + (tz2d+1)/2;
            size_t e_pn = 2*(e/2) + (tz2e+1)/2;
            size_t f_pn = 2*(f/2) + (tz2f+1)/2;
            double me_pn = GetME_pn( Jab_in, Jde_in, twoJ, a_pn,b_pn,c_pn,d_pn,e_pn,f_pn);
            me_iso += me_pn  * clebsch_ab * clebsch_abc * clebsch_de * clebsch_def ;
          }
        }
      }
    }
  }
  return me_iso;
}


/*
void ThreeBodyStorage_pn::SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}


void ThreeBodyStorage_pn::AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}



void ThreeBodyStorage_pn::SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBodyStorage_pn::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}

void ThreeBodyStorage_pn::AddToME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBodyStorage_pn::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}
*/




ThreeBodyStorage::ME_type ThreeBodyStorage_pn::GetME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket) const
{
  if ( not is_allocated ) return 0;
  if (ch_bra==ch_ket and ibra==iket and herm==-1) return 0;
  size_t index;
  int herm_flip;
  AccessME(ch_bra,ch_ket,ibra,iket,index,herm_flip);
  //  std::cout << "   got index = " << index << "   MatEl size is " << MatEl.size() << std::hex << "  address = " << &(MatEl[index]) << std::dec << "   ->  " << MatEl[index]  << "   herm_flip = " << herm_flip << "   because herm is " << herm << std::endl;
  return MatEl.at(index)*herm_flip;
}



void ThreeBodyStorage_pn::SetME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyStorage::ME_type V)
{
//  std::cout << " IN " << __FILE__ << " " << __func__ << "  is_allocated = " << is_allocated << std::endl;
  size_t index;
  int herm_flip;
  AccessME(ch_bra,ch_ket,ibra,iket,index,herm_flip);
  //  std::cout << "   got index = " << index << "   MatEl size is " << MatEl.size() << std::endl;
  MatEl.at(index) = herm_flip * V;
  //  std::cout << " all is well. " << std::endl;
}


void ThreeBodyStorage_pn::AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyStorage::ME_type V)
{
  size_t index;
  int herm_flip;
  AccessME(ch_bra,ch_ket,ibra,iket,index,herm_flip);
//  std::cout << " IN " << __FILE__ << "  " << __func__ << "  index = " << index << "  size of MatEl = " << MatEl.size()
//            << "   adding " << herm_flip << " * " << V << "  to " << MatEl.at(index) << std::endl;
  MatEl.at(index) += herm_flip * V;
}


double ThreeBodyStorage_pn::Norm() const
{
  double norm = 0;
  for (auto me : MatEl)  norm += me*me;
  return norm;
}


// Set all elements to zero
void ThreeBodyStorage_pn::Erase()
{
   MatEl.assign( MatEl.size(), 0. );
}

void ThreeBodyStorage_pn::Deallocate()
{
    std::vector<ThreeBodyStorage_pn::pnME_type>().swap( MatEl);
    is_allocated = false;
}

size_t ThreeBodyStorage_pn::size() const
{
  return total_dimension * sizeof(pnME_type);
}


void ThreeBodyStorage_pn::WriteBinary(std::ofstream& f)
{
  f.write((char*)&E3max,sizeof(E3max));
  f.write((char*)&total_dimension,sizeof(total_dimension));
  f.write((char*)&MatEl[0],total_dimension);
}

void ThreeBodyStorage_pn::ReadBinary(std::ifstream& f)
{
  f.read((char*)&E3max,sizeof(E3max));
  f.read((char*)&total_dimension,sizeof(total_dimension));
  Allocate();
  f.read((char*)&MatEl[0],total_dimension*sizeof(ThreeBodyStorage_pn::pnME_type));
}


//////////////////////////////////////////////////////////////////////////////////////
//////   Internal implementation methods
///////////////////////////////////////////////////////////////////////////////////





//  Here's how we fold a symmetric/antisymmetric matrix into a 1D array. We could either use column-major or row-major ordering.
// 
// ( 0,0   0,1  0,2  0,3  )        0    1    2    3     4    5    6     7    8    9        dimension D=4
// ( 1,0   1,1  1,2  1,3  )  -> [ 0,0  1,0  2,0  3,0   1,1  2,1  3,1   2,2  3,2  3,3 ]     index(a,b) = a + (2*D-b-1)*b/2   (column-major)
// ( 2,0   2,1  2,2  2,3  )  -> [ 0,0  0,1  0,2  0,3   1,1  1,2  1,3   2,2  2,3  3,3 ]     index(a,b) = b + (2*D-a-1)*a/2   (row-major)
// ( 3,0   3,1  3,2  3,3  )        0   1+0  2+0  3+0   1+3  2+3  3+3   2+5  3+5  3+6  
//
//  For the storage of 3N matrix elements, we go with row-major ordering. This is because it seems more natural to structure loops as
//  for bra, for ket,   which means the ket (second index) is incremented in the inner loop, and so that data should be adjacent,
//  which will happen with row-major ordering.
void ThreeBodyStorage_pn::AccessME(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, size_t& index, int& herm_flip) const
{
  herm_flip = (  (ch_bra < ch_ket ) or ((ch_ket==ch_bra) and (ibra>=iket))) ? 1 : herm; // we store bra < ket
  size_t ch_1 = std::min(ch_bra,ch_ket);
  size_t ch_2 = std::max(ch_bra,ch_ket);
  size_t iket_1 = (ch_bra==ch_ket) ? std::min(ibra,iket) : (  (ch_bra<ch_ket) ? ibra : iket   );
  size_t iket_2 = (ch_bra==ch_ket) ? std::max(ibra,iket) : (  (ch_bra<ch_ket) ? iket : ibra   );
  // so now ch_1,ch_2 and iket_1,iket_2 are ordered the way we store them

  auto iter_ch_start = ch_start.find({ch_1,ch_2});
  if ( ( iter_ch_start == ch_start.end() )
      or    iket_1>ch_dim[ch_1] or iket_2>ch_dim[ch_2])
  {
    std::ostringstream oss;
    oss << __FILE__ << " " << __func__ << " ch_bra,ch_ket " << ch_bra << " " << ch_ket << "  ibra,iket " << ibra << " " << iket;
    throw std::domain_error( oss.str() );
  }
  // ch_start points to where the matrix for this channel starts, and the rest
  // folds two indices into one, assuming we only store the half-triangular matrix.
  if (ch_1==ch_2)
  {
    index = (iter_ch_start->second) +   (2*ch_dim[ch_2] - iket_1 - 1)*iket_1/2 + iket_2 ;
  }
  else
  {
    index = (iter_ch_start->second) + ch_dim[ch_2] * iket_1 + iket_2;
    // include phase factor for 3b tensor operator
    // < i || T || j > = herm * (-)^(i-j) < j || T || i >*
    if ( ch_bra > ch_ket )  
    {
      ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
      herm_flip *= modelspace->phase((Tbc_bra.twoJ - Tbc_ket.twoJ)/2);
    }
  }
  // check for trouble
  if (index>=MatEl.size())
  {
    std::ostringstream oss;
    oss << __func__ << " ch_bra,ch_ket " << ch_bra << " " << ch_ket << "  ibra,iket " << ibra << " " << iket << "  index= " << index << " > MatEl.size() = " << MatEl.size() ;
    throw std::domain_error( oss.str() );
  }
}





void ThreeBodyStorage_pn::Print()
{
    for (auto &it : ch_start )
    {
      size_t ch_bra = it.first.ch_bra;
      size_t ch_ket = it.first.ch_ket;
      ThreeBodyChannel &Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
      ThreeBodyChannel &Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
      size_t nbras = Tbc_bra.GetNumberKets();
      size_t nkets = Tbc_ket.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras; ibra++)
      {
        Ket3& bra = Tbc_bra.GetKet( ibra);
        for ( size_t iket = 0; iket <nkets; iket++)
        {
           Ket3& ket = Tbc_ket.GetKet( iket);
           double ME = GetME_pn_ch( ch_bra,  ch_ket, ibra,  iket);
           std::cout << "  " << bra.p << " " << bra.q << " " << bra.r << " " << bra.Jpq << "   " << ket.p << " " << ket.q << " " << ket.r << " " << ket.Jpq << " , " << Tbc_bra.twoJ << " " << Tbc_ket.twoJ
                     << "( " << bra.op->cvq <<  bra.oq->cvq << bra.oR->cvq << ket.op->cvq << ket.oq->cvq << ket.oR->cvq  << " )" << "  => " << ME << std::endl;
        }
      }
    }


}



