
#include "AngMom.hh"
#include <sstream>
#include "ThreeBodyStorage_iso.hh"



bool ThreeBodyStorage_iso::none_allocated = true;

ThreeBodyStorage_iso::ThreeBodyStorage_iso( const ThreeBodyStorage_iso& TBS_in )
:  ThreeBodyStorage( TBS_in ), MatEl(TBS_in.MatEl), OrbitIndexHash(TBS_in.OrbitIndexHash), total_dimension(0)
{}


//ThreeBodyStorage_iso::~ThreeBodyStorage_iso( )
//{
//}


//std::shared_ptr<ThreeBodyStorage> ThreeBodyStorage_iso::Clone() const { return std::shared_ptr<ThreeBodyStorage>( new ThreeBodyStorage_iso( *this)); };
std::unique_ptr<ThreeBodyStorage> ThreeBodyStorage_iso::Clone() const { return std::unique_ptr<ThreeBodyStorage>( new ThreeBodyStorage_iso( *this)); };



  void ThreeBodyStorage_iso::Multiply(const double rhs) 
  {
    for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] *= rhs;
  }

  void ThreeBodyStorage_iso::Add(const ThreeBodyStorage& rhs) 
  {
    if ( not rhs.IsAllocated() ) return;
//    if ( not this->IsAllocated() ) 
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_iso*)(&rhs))->MatEl;
       for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] += rhsMatEl[i];
    }
    else 
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << "  allocated?  " << this->IsAllocated() << "  " << rhs.IsAllocated()  << "  norms " << this->Norm() << "  " << rhs.Norm() << "   dimensions " << this->size() << "   " << rhs.size() << " dying." << std::endl;
//      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void ThreeBodyStorage_iso::Subtract(const ThreeBodyStorage& rhs)
  {
    if ( not rhs.IsAllocated() )return;
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_iso*)(&rhs))->MatEl;
       for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] += rhsMatEl[i];
    }
    else
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }



//ThreeBodyStorage_iso& ThreeBodyStorage_iso::operator*=(const double rhs)
//{
//  for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] *= rhs;
//}
//
//ThreeBodyStorage_iso& ThreeBodyStorage_iso::operator+=(const ThreeBodyStorage_iso& rhs)
//{
//  for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] += rhs.MatEl[i];
//}
//
//ThreeBodyStorage_iso& ThreeBodyStorage_iso::operator-=(const ThreeBodyStorage_iso& rhs)
//{
//  for ( size_t i=0; i<MatEl.size();i++ )   MatEl[i] -= rhs.MatEl[i];
//}


void ThreeBodyStorage_iso::Allocate()
{
  MatEl.clear();
  OrbitIndexHash.clear();
  E3max = modelspace->GetE3max();
  int norbits = modelspace->GetNumberOrbits();
  std::cout << "Begin AllocateThreeBody() in " << __FILE__ << " with E3max = " << E3max << " norbits = " << norbits << std::endl;
  int lmax = 50000; // maybe do something with this later...
  total_dimension = 0;

  ISOSPIN_BLOCK_DIMENSION = 5; // If rank_T is 0 or 2, then there are 5 possible isospin combinations
  if (rank_T==1) ISOSPIN_BLOCK_DIMENSION = 9; 
  else if (rank_T==3) ISOSPIN_BLOCK_DIMENSION = 1;

  // loop through even orbits, since odd ones are neutrons and we're doing isospin here.
  // we store matrix elements for a>=b>=c and d>=e>=f
  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   int ea = 2*oa.n+oa.l;
   if (ea>E3max) break;
   if (ea>emax) break;
   for (int b=0; b<=a; b+=2)
   {
     if (oa.l > lmax) break;
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;

     int Jab_min = std::abs(oa.j2-ob.j2)/2;
     int Jab_max = (oa.j2+ob.j2)/2;
     for (int c=0; c<=b; c+=2)
     {
       if (ob.l > lmax) break;
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       for (int d=0; d<=a; d+=2)
       {
         if (oc.l > lmax) break;
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           if (od.l > lmax) break;
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             if (oe.l > lmax) break;
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             if (((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2 !=parity) or (of.l > lmax)) continue;

             OrbitIndexHash[ KeyHash(a,b,c,d,e,f) ] = total_dimension;
//             std::cout << "  hashing " << a << " " << b << " " << c << " " << d << " " << e << " " << f << "  : " << KeyHash(a,b,c,d,e,f) << "  -> " << total_dimension << "  check:  " << OrbitIndexHash.at( KeyHash(a,b,c,d,e,f) ) << std::endl;
             int Jde_min = std::abs(od.j2-oe.j2)/2;
             int Jde_max = (od.j2+oe.j2)/2;

             for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
             {
              for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
              {
                int twoJ_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
                int twoJ_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
                for (int twoJ=twoJ_min; twoJ<=twoJ_max; twoJ+=2)
                {
                  total_dimension += ISOSPIN_BLOCK_DIMENSION; // no need to loop over the isospin combinations if we're just counting. 
                } //J2
              } //Jde
             } //Jab
           } //f
         } //e
       } //d
     } //c
   } //b
  } //a

  MatEl.resize(total_dimension,0.0);
  std::cout << "Allocated " << total_dimension << " three body matrix elements (" <<  total_dimension * sizeof(isoME_type)/1024./1024./1024. << " GB), "
       << std::endl << "  number of buckets in hash table: " << OrbitIndexHash.bucket_count() << "  and load factor = " << OrbitIndexHash.load_factor()
       << "  estimated storage ~ " << ((OrbitIndexHash.bucket_count()+OrbitIndexHash.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB"
       << std::endl;
  is_allocated = true;
}








//////////////////////////////////////////////////////////////////////////////////////
//////   Interface methods
///////////////////////////////////////////////////////////////////////////////////



//*******************************************************************
/// Get three body matrix element in isospin formalism, when we store in isospin formalism
/// \f$ V_{abcdef}^{J_{ab}J_{de}Jt_{ab}t_{de}T} \f$
///  (which is how they're stored).
/// The elements are stored with the following restrictions: \f$ a\geq b \geq c\f$,
/// \f$ d\geq e \geq f\f$, \f$ a\geq d\f$. If \f$ a=d\f$ then \f$ b \geq e\f$,
/// and if \f$ b=e \f$ then \f$ c \geq f \f$.
/// Other orderings are obtained by recoupling on the fly.
//*******************************************************************
ThreeBodyStorage::ME_type ThreeBodyStorage_iso::GetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const
{
   if ((a_in/2==b_in/2) and (Jab_in+tab_in)%2==0) return 0; // Make sure this is ok
   if ((d_in/2==e_in/2) and (Jde_in+tde_in)%2==0) return 0; // Make sure this is ok
   auto elements =  AccessME(Jab_in,Jde_in,twoJ,tab_in,tde_in,twoTabc,twoTdef,a_in,b_in,c_in,d_in,e_in,f_in);
   double me = 0;
   for (auto elem : elements) me += MatEl.at(elem.first) * elem.second;
   return me;
}


//*******************************************************************
/// Get three body matrix element in proton-neutron formalism.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
ThreeBodyStorage::ME_type ThreeBodyStorage_iso::GetME_pn(int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f) const
{

   if (a==b and a==c and modelspace->GetOrbit(a).j2<3) return 0;
   if (d==e and d==f and modelspace->GetOrbit(d).j2<3) return 0;
   if (a==b and Jab_in%2>0) return 0;
   if (d==e and Jde_in%2>0) return 0;
   double tza = modelspace->GetOrbit(a).tz2*0.5;
   double tzb = modelspace->GetOrbit(b).tz2*0.5;
   double tzc = modelspace->GetOrbit(c).tz2*0.5;
   double tzd = modelspace->GetOrbit(d).tz2*0.5;
   double tze = modelspace->GetOrbit(e).tz2*0.5;
   double tzf = modelspace->GetOrbit(f).tz2*0.5;
   double dTz =  (tza+tzb+tzc) - (tzd+tze+tzf);
   if ( std::abs(dTz) > rank_T ) return 0;


   double Vpn=0;
   int twoTabc_min = std::abs(tza+tzb+tzc)*2;
   for (int tab=std::abs(tza+tzb); tab<=1; ++tab)
   {
      // CG calculates the Clebsch-Gordan coefficient  TODO: There are only a few CG cases, and we can probably use a specific formula rather than the general one.
      double CG1 = AngMom::CG(0.5,tza, 0.5,tzb, tab, tza+tzb);
      for (int tde=std::abs(tzd+tze); tde<=1; ++tde)
      {
         double CG2 = AngMom::CG(0.5,tzd, 0.5,tze, tde, tzd+tze);
         if (CG1*CG2==0) continue;
         for (int twoTabc=twoTabc_min; twoTabc<=3; twoTabc+=2)
         {
           double CG3 = AngMom::CG(tab,tza+tzb, 0.5,tzc, twoTabc/2., tza+tzb+tzc);
           int twoTdef_min = std::max( (int)std::abs(tzd+tze+tzf)*2, std::abs(twoTabc-2*rank_T));
           int twoTdef_max = twoTabc + 2*rank_T;
           for (int twoTdef=twoTdef_min; twoTdef<=twoTdef_max; twoTdef+=2)
           {
             double CG4 = AngMom::CG(tde,tzd+tze, 0.5,tzf, twoTdef/2., tzd+tze+tzf);
             if (CG3*CG4==0) continue;
             double Viso = GetME_iso(Jab_in,Jde_in,twoJ,tab,tde,twoTabc,twoTdef,a,b,c,d,e,f);
             if (rank_T > 0 ) // the matrix elements are reduced in T, so apply Wigner-Eckart
             {
               Viso *= AngMom::CG(0.5*twoTdef,tzd+tze+tzf,  rank_T,dTz,  0.5*twoTabc,tza+tzb+tzc  ) / sqrt( twoTabc+1);
             }
             Vpn += CG1*CG2*CG3*CG4* Viso;
//             std::cout << __FILE__ << " " << __LINE__ << "  " << Jab_in << " " << Jde_in << " " << twoJ << "  " << tab << " " << tde << " " << twoTabc << " " << twoTdef << "   " << a << " " << b << " " << c << " " << d << " "<< e << " " << f << "     " << CG1 << " " << CG2 << " " << CG3 << " " << CG4 << " " << Viso << "  -> " << Vpn << std::endl;
           }
         }
      }
   }
   return Vpn;
}


//*******************************************************************
/// Set a three body matrix element. Since only a subset of orbit
/// orderings are stored, we need to recouple if the input ordering
/// is different.
//*******************************************************************
void ThreeBodyStorage_iso::SetME_iso(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, ThreeBodyStorage::ME_type V)
{
   auto elements = AccessME(Jab_in,Jde_in,twoJ,tab_in,tde_in,twoTabc,twoTdef,a_in,b_in,c_in,d_in,e_in,f_in);
   double me = 0;
   for (auto elem : elements)  me += MatEl.at(elem.first) * elem.second;
   for (auto elem : elements)  MatEl.at(elem.first) += (V-me)*elem.second;
}



/*
void ThreeBodyStorage_iso::SetME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}
*/



void ThreeBodyStorage_iso::AddToME_iso(int Jab, int Jde, int twoJ, int tab, int tde, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ThreeBodyStorage::ME_type V)
{
   auto elements = AccessME(Jab,Jde,twoJ,tab,tde,twoTabc,twoTdef,a,b,c,d,e,f);
   for (auto elem : elements)  MatEl.at(elem.first) += V * elem.second;
}


/*
void ThreeBodyStorage_iso::AddToME_pn( int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}


// Most likely this will not be called. Lets not bother implementing it unless we actually need it.
ThreeBodyStorage::ME_type ThreeBodyStorage_iso::GetME_pn_ch(size_t chbra, size_t chket, size_t ibra, size_t iket) const
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}


void ThreeBodyStorage_iso::AddToME_pn_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}

void ThreeBodyStorage_iso::SetME_pn_ch(  size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyStorage::ME_type V)
{
  std::cout << __func__ << "   on " << __FILE__ << "  line " << __LINE__ << "  is not yet implemented." << std::endl;
  std::exit(EXIT_FAILURE);
}
*/

double ThreeBodyStorage_iso::Norm() const
{
  double norm = 0;
  for (auto me : MatEl)  norm += me*me;
  return norm;
}

// Set all elements to zero
void ThreeBodyStorage_iso::Erase()
{
   MatEl.assign( MatEl.size(), 0. );
}

void ThreeBodyStorage_iso::Deallocate()
{
   std::vector<ThreeBodyStorage_iso::isoME_type>().swap(MatEl);
   std::unordered_map<size_t, size_t>().swap( OrbitIndexHash  );
}

size_t ThreeBodyStorage_iso::size() const
{
  return total_dimension * sizeof(isoME_type);
}


void ThreeBodyStorage_iso::WriteBinary(std::ofstream& f)
{
  f.write((char*)&E3max,sizeof(E3max));
  f.write((char*)&total_dimension,sizeof(total_dimension));
//  f.write((char*)&MatEl[0],total_dimension);
  f.write((char*)&MatEl[0],total_dimension);
}

void ThreeBodyStorage_iso::ReadBinary(std::ifstream& f)
{
  f.read((char*)&E3max,sizeof(E3max));
  f.read((char*)&total_dimension,sizeof(total_dimension));
  Allocate();
  f.read((char*)&MatEl[0],total_dimension*sizeof(ThreeBodyStorage_iso::isoME_type));
}




//////////////////////////////////////////////////////////////////////////////////////
//////   Internal implementation methods
///////////////////////////////////////////////////////////////////////////////////



/// Hash function to map six indices to a single long unsigned int.
/// Each index gets 10 bits, for a maximum of 1024 indices.
/// If we have good isospin, this means we're ok up to emax=43.
// In the future, if this isn't sufficient, we could change
// the key type to a bitset<128>, which would allow up to emax=1446.
// I'm guessing I chose the weird ordering because the matrix elements are stored
// such that a>=b>=c  and d>=e>=f, but I don't immediately see why there should be a benefit from this.
size_t ThreeBodyStorage_iso::KeyHash(size_t a,size_t b,size_t c,size_t d,size_t e,size_t f) const
{
 return (   (a/2)
         + ((d/2) << 10)
         + ((b/2) << 20)
         + ((e/2) << 30)
         + ((c/2) << 40)
         + ((f/2) << 50) );

}

void ThreeBodyStorage_iso::KeyUnhash(size_t& key, size_t& a, size_t& b, size_t& c, size_t& d, size_t& e, size_t& f) const
{
  size_t Lowest_ten_bits = 0x3FF;
  a = 2* ((key >>  0) & Lowest_ten_bits);
  d = 2* ((key >> 10) & Lowest_ten_bits);
  b = 2* ((key >> 20) & Lowest_ten_bits);
  e = 2* ((key >> 30) & Lowest_ten_bits);
  c = 2* ((key >> 40) & Lowest_ten_bits);
  f = 2* ((key >> 50) & Lowest_ten_bits);

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
ThreeBodyStorage::Permutation ThreeBodyStorage_iso::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
{
   a_in -= a_in%2;  // We make the indices even, because in isospin formalism we only use proton (even) orbits and so less sorting is needed
   b_in -= b_in%2;
   c_in -= c_in%2;
   a=a_in;
   b=b_in;
   c=c_in;
//   if (a>b)  std::swap(a,b); // the a <= b <= c ordering matches with the 2body storage, which have if (q<p) continue;
//   if (b>c)  std::swap(b,c);
//   if (a>b)  std::swap(a,b);
  // we store matrix elements for a>=b>=c and d>=e>=f  so this is the correct sorting for the isospin matrix elements
   if (a<b)  std::swap(a,b);
   if (b<c)  std::swap(b,c);
   if (a<b)  std::swap(a,b);


//   int recoupling_case;
   Permutation recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? ABC : ACB;
   else if (a_in==b)  recoupling_case = (b_in==a) ? BAC : BCA;
   else               recoupling_case = (b_in==a) ? CAB : CBA;

   return recoupling_case;
}









//*******************************************************************
/// Since the code for setting or getting matrix elements is almost
/// identical, do all the work here to pull out a list of indices
/// and coefficients which are needed for setting or getting.
//*******************************************************************
//std::vector<std::pair<size_t,double>> ThreeBodyME::AccessME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const
std::vector<std::pair<size_t,double>> ThreeBodyStorage_iso::AccessME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const
{

   std::vector<std::pair<size_t,double>> elements;
   // Re-order so that a>=b>=c, d>=e>=f
   int a,b,c,d,e,f;
   Permutation abc_recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);
   Permutation def_recoupling_case = SortOrbits(d_in,e_in,f_in,d,e,f);
   int herm_flip = +1;

   if (d>a or (d==a and e>b) or (d==a and e==b and f>c))
   {
	std::swap(a,d);
      	std::swap(b,e);
      	std::swap(c,f);
      	std::swap(Jab_in,Jde_in);
      	std::swap(tab_in,tde_in);
        std::swap(twoTabc,twoTdef);
      	std::swap(abc_recoupling_case, def_recoupling_case);
        herm_flip *= herm * AngMom::phase((twoTabc-twoTdef)/2);
//        herm_flip *= herm;
   }

   auto it_hash = OrbitIndexHash.find(KeyHash(a,b,c,d,e,f));
//   if ( it_hash == end(OrbitIndexHash) )   return elements;

   if ( it_hash == end(OrbitIndexHash) )
   {
      std::cout << " IN " << __FILE__ << " " << __LINE__ << "  orbits not in hash for " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
      return elements;
   }


   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);
   if (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) return elements;
   if (2*(od.n+oe.n+of.n)+od.l+oe.l+of.l > E3max) return elements;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;
   double je = oe.j2*0.5;
   double jf = of.j2*0.5;

   int Jab_min = std::abs(ja-jb);
   int Jde_min = std::abs(jd-je);
   int Jab_max = (ja+jb);
   int Jde_max = (jd+je);

   int tab_min = twoTabc==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = twoTdef==3 ? 1 : 0;
   int tde_max = 1;


   auto indx = it_hash->second;
   if (indx > MatEl.size())
   {
     std::cout << "!!!! Troubling index in " << __func__ << ". abcdef = " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
     std::ostringstream oss;
     oss << "ThreeBodyStorage_iso::AccessME() --  AAAAHHH indx = " << indx << "  but MatEl.size() = " << MatEl.size() << std::endl;
     throw std::domain_error( oss.str() );
   }


//    std::cout << __FILE__ << " " << __func__ << "  ---------------------------------------"   << std::endl;

   int J_index = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {

     double Cj_abc = RecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,twoJ);
     // Pick up a -1 for odd permutations
//     if ( abc_recoupling_case/3 != def_recoupling_case/3 ) Cj_abc *= -1;
      Cj_abc *= PermutationPhase( abc_recoupling_case );

     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_def = RecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,twoJ);
       Cj_def *= PermutationPhase( def_recoupling_case );

       int twoJ_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
       int twoJ_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
       if (twoJ_min>twoJ_max) continue;
       J_index += (twoJ-twoJ_min)/2 * ISOSPIN_BLOCK_DIMENSION;
//       J_index += (twoJ-twoJ_min)/2*5;

       if (twoJ>=twoJ_min and twoJ<=twoJ_max and std::abs(Cj_abc*Cj_def)>1e-8)
       {
         for (int tab=tab_min; tab<=tab_max; ++tab)
         {
           double Ct_abc = RecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,twoTabc);
           for (int tde=tde_min; tde<=tde_max; ++tde)
           {
             double Ct_def = RecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,twoTdef);
             if (std::abs(Ct_abc*Ct_def)<1e-8) continue;
             if (herm==-1 and a==d and b==e and c==f and Jab==Jde and tab==tde) continue; // TODO: check this is ok

             // Index mapping for the various isospin cases
             // rank_T=0, 2*tab+tde+twoTabc/2:                 (0,1;0,1)=>0  (0,1;1,1)=>1  (1,1;0,1)=>2 (1,1;1,1)=>3  (1,3;1,3)=>4
             //
             // rank_T=1, tab+twoTabc/2 + 3(tde+twoTdef/2):    (0,1;0,1)=>0 (1,1;0,1)=>1 (1,3;0,1)=>2
             //                                                (0,1;1,1)=>3 (1,1;1,1)=>4 (1,3;1,1)=>5
             //                                                (0,1;1,3)=>6 (1,1;1,3)=>7 (1,3;1,3)=>8
             //
             // rank_T=2, (tab-tde) + (twoTabc-twoTdef)/2 + 2: (0,1;1,3)=>0  (1,1;1,3)=>1  (1,3;1,3)=>2  (1,3;1,1)=>3  (1,3,0,1)=>4
             //
             int Tindex = 2*tab + tde + (twoTabc)/2;    // (0,1) (1,1) (1,3) => 0,1,2   twoTabc/2 +2*tab or   tab + (twoTabc)/2
             if ( rank_T==1 ) Tindex = ( tab + twoTabc/2) + 3*( tde + twoTdef/2);
             else if ( rank_T==2 ) Tindex = (tab-tde) + (twoTabc-twoTdef)/2 + 2;
             else if ( rank_T==3 ) Tindex = 0;

//              std::cout << __FILE__ << " " << __func__ << "   " << indx << " + " << J_index << " + " << Tindex << " * "
//                        << "  " << Cj_abc << " " << Cj_def << " " << Ct_abc << " " << Ct_def << "  * " << herm_flip << std::endl;

              if (indx+J_index+Tindex > MatEl.size())
              {
                 std::cout << "ThreeBodyStorage_iso::AccessME() " << __LINE__ <<  " --  AAAAHHH indx = "
                           << indx  << " + " << J_index << " + " << Tindex << " = " << indx+J_index+Tindex
                           << "  but MatEl.size() = " << MatEl.size() << "    t values: "
                           << tab << " " << tde << " " << twoTabc << " " << twoTdef << std::endl;
                 std::ostringstream oss;
                 oss << "ThreeBodyStorage_iso::AccessME() line " << __LINE__ << std::endl;
                 throw std::domain_error( oss.str() );
              }

             elements.emplace_back( std::make_pair(indx + J_index + Tindex, Cj_abc * Cj_def * Ct_abc * Ct_def * herm_flip )) ;
           }
         }
       }
       J_index += (twoJ_max-twoJ+2)/2 * ISOSPIN_BLOCK_DIMENSION;
     }
   }
   return elements;
}



