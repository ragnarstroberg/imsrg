
#include "ThreeBodyStorage_no2b.hh"
#include "AngMom.hh"

#include <omp.h>




//OrbitIsospin::OrbitIsospin(int idx, int n, int l, int j)
//  : idx(idx), n(n), l(l), j(j), e(2*n+l)
//{}


ThreeBodyChannelNO2B::ThreeBodyChannelNO2B()
  : J2(-1), P2(-1), J1(-1), P1(-1), T3(-1), Ndim(-1)
{}

ThreeBodyChannelNO2B::ThreeBodyChannelNO2B(int J2, int P2, int J1, int P1, int T3,  ThreeBodySpaceNO2B& thr)
  : J2(J2), P2(P2), J1(J1), P1(P1), T3(T3)
{
  Ndim = 0;
  int emax  = thr.Emax;
  int e2max = thr.E2max;
  int e3max = thr.E3max;
  int Norbs = thr.GetNumberIsospinOrbits();
  for (int ia=0; ia < Norbs; ia++)
  {
    OrbitIsospin & oa = thr.GetIsospinOrbit(ia);
    if(oa.e > emax) continue;
    for (int ib=ia; ib < Norbs; ib++)
    {
      OrbitIsospin & ob = thr.GetIsospinOrbit(ib);
      if(ob.e > emax) continue;
      if(oa.e + ob.e > e2max) continue;
      if( not AngMom::Triangle( oa.j, ob.j, 2*J2) ) continue;
//      if(std::abs(oa.j-ob.j) > 2*J2) continue;
//      if(        (oa.j+ob.j) < 2*J2) continue;
      if((oa.l+ob.l)%2 != P2) continue;
      for (int ic=0; ic < Norbs; ic++)
      {
        OrbitIsospin & oc = thr.GetIsospinOrbit(ic);
        if(oc.j != J1) continue;
        if(oc.l%2 != P1) continue;
        if(oc.e > emax) continue;
        if(oa.e + oc.e > e2max) continue;
        if(ob.e + oc.e > e2max) continue;
        if(oa.e + ob.e + oc.e > e3max) continue;
        for (int T2: {0,1})
        {
//          if( std::abs(2*T2-1) > T3 ) continue;
//          if(         (2*T2+1) < T3 ) continue;
          if( not AngMom::Triangle( 2*T2, 1, T3) ) continue;
          if(((J2+T2)%2 == 0) and (ia == ib) ) continue;
          int ph = 1 - 2*(((oa.j + ob.j)/2 + J2 + T2)%2);
          int key_ab = Hash_abct(ia, ib, ic, T2);
          int key_ba = Hash_abct(ib, ia, ic, T2);
          abct2n[key_ab] = Ndim;
          abct2n[key_ba] = Ndim;
          iphase[key_ab] = 1;
          iphase[key_ba] = ph;
          Ndim += 1;
        }// for T2
      }// for ic
    }// for ib
  }// for ia
}

int ThreeBodyChannelNO2B::Hash_abct(int a, int b, int c, int Tab) const
{
  return ( a + (b<<10) + (c<<20) + (Tab<<30) );
}

void ThreeBodyChannelNO2B::UnHash_abct(int key, int & a, int & b, int & c, int & Tab) const
{
  int Lowest_ten_bits = 0x3FF;
  a = ( (key>> 0) & Lowest_ten_bits);
  b = ( (key>>10) & Lowest_ten_bits);
  c = ( (key>>20) & Lowest_ten_bits);
  Tab = ( (key>>30) & Lowest_ten_bits);
}



//ThreeBodySpaceNO2B::~ThreeBodySpaceNO2B()
//{}

ThreeBodySpaceNO2B::ThreeBodySpaceNO2B()
  :  ThreeBodyChannels() , NChannels(-1)
{}


ThreeBodySpaceNO2B::ThreeBodySpaceNO2B( int emax, int e2max, int e3max, int lmax)
 : Emax(emax), E2max(e2max), E3max(e3max), Lmax(lmax)
{
  
  double t_start = omp_get_wtime();

// Allocate all the 1b orbits
  int idx = 0;
  for (int e=0; e<=Emax; ++e)
  {
    int lmin = e%2;
    for (int l=lmin; l<=std::min(e,Lmax); l+=2)
    {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2)
      {
        OrbitIsospin orb(idx,n,l,twoj);
        iOrbits.push_back(orb);
        nlj2idx[{n,l,twoj}]=idx;
        idx += 1;
      }
    }
  }


// Allocate the ThreeBody channels
  int J2max = std::min(2*lmax+1, e2max+1);
  int count = 0;
  for (int J2=0; J2<=J2max; ++J2)
  {
    for (int P2: {0,1})
    {
      for (int J1=1; J1<=2*lmax+1; J1+=2)
      {
        for (int P1: {0,1})
        {
          for (int T3 : {1,3})
          {
            ThreeBodyChannelNO2B channel(J2, P2, J1, P1, T3, *this);
            if(channel.Ndim < 1) continue;

            ThreeBodyChannels.push_back(channel);
            int key = Hash_Channel(J2,P2,J1,P1,T3);
            idcs2ch[key] = count;
            count += 1;
          }
        }
      }
    }
  }
  NChannels = ThreeBodyChannels.size();
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}



int ThreeBodySpaceNO2B::Hash_Channel(int J2, int P2, int J1, int P1, int T3) const
{
  return ( J2 + (J1<<8) + (P2<<16) + (P1<<17) + (T3<<18) );
}

void ThreeBodySpaceNO2B::UnHash_Channel(int key, int& J2, int& P2, int& J1, int& P1, int& T3) const
{
  J2 = ( (key>> 0) & 255 );
  J1 = ( (key>> 8) & 255 );
  P2 = ( (key>>16) & 1 );
  P1 = ( (key>>17) & 1 );
  T3 = ( (key>>18) & 255 );
}

// Take Jab,Pab,Jc,Pc,T and return the channel number for subsequent lookup in the vector ThreeBodyChannels
int ThreeBodySpaceNO2B::GetChannelIndex(int Jab, int Pab, int Jc, int Pc, int T)const
{
  auto key = Hash_Channel(Jab, Pab, Jc, Pc, T);
  auto iter = idcs2ch.find(key);
  int ch =-1;
  if ( iter!=idcs2ch.end() ) ch = iter->second;
  return ch;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////  Begin implementation of ThreeBodyStorage_no2b methods. This is the main event.
////////////////////////////////////////////////////////////////////////////////////////

template<class StoreType>
ThreeBodyStorage_no2b<StoreType>::ThreeBodyStorage_no2b( const ThreeBodyStorage_no2b<StoreType>& TBS_in )
//: MatEl(TBS_in.MatEl), threebodyspace(TBS_in.threebodyspace),    ThreeBodyStorage( TBS_in )
:  ThreeBodyStorage( TBS_in ),  MatEl(TBS_in.MatEl), threebodyspace(TBS_in.threebodyspace)
{}


template<class StoreType>
std::unique_ptr<ThreeBodyStorage> ThreeBodyStorage_no2b<StoreType>::Clone() const { return std::unique_ptr<ThreeBodyStorage>( new ThreeBodyStorage_no2b<StoreType>( *this)); };
//std::shared_ptr<ThreeBodyStorage> ThreeBodyStorage_no2b<StoreType>::Clone() const { return std::shared_ptr<ThreeBodyStorage>( new ThreeBodyStorage_no2b<StoreType>( *this)); };


template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Allocate()
{
  // Create the three body space
  threebodyspace = ThreeBodySpaceNO2B(emax,E2max,E3max, std::min(emax,lmax) );

  // Allocate the vectors for the matrix element storage
  for (int ch=0; ch<threebodyspace.NChannels; ch++)
  {
    ThreeBodyChannelNO2B& ch_no2b = threebodyspace.ThreeBodyChannels.at(ch);
    size_t n = ch_no2b.Ndim;

    MatEl[ch] = std::vector<StoreType>( n*(n+1)/2, StoreType(0.0));
  }
  is_allocated = true;
}





template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Multiply(const double rhs) 
{
  for ( auto& itmat : MatEl )
  {
    for ( auto& it : itmat.second )
    {
      it *= rhs;
    }
  }
}


template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Add(const ThreeBodyStorage& rhs) 
{
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_no2b*)(&rhs))->MatEl;

       for ( auto& itmat : rhsMatEl )
       {
         auto ch = itmat.first;
         for ( size_t i=0; i<itmat.second.size(); i++)
         {
           MatEl[ch][i] += itmat.second[i];
         }
       }
    }
    else
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
}


template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Subtract(const ThreeBodyStorage& rhs) 
{
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_no2b*)(&rhs))->MatEl;

       for ( auto& itmat : rhsMatEl )
       {
         auto ch = itmat.first;
         for ( size_t i=0; i<itmat.second.size(); i++)
         {
           MatEl[ch][i] -= itmat.second[i];
         }
       }
    }
    else
    {
      std::cout << "OOPS!!! Tried to " << __func__ << "  with incompatible storage modes  " << rhs.GetStorageMode() << " and " << this->GetStorageMode() << " dying." << std::endl;
      std::exit(EXIT_FAILURE);
    }
}



template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::SetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int twoT, ThreeBodyStorage::ME_type V)
{
   OrbitIsospin & oa = threebodyspace.iOrbits[a];
   OrbitIsospin & ob = threebodyspace.iOrbits[b];
   OrbitIsospin & oc = threebodyspace.iOrbits[c];
   OrbitIsospin & od = threebodyspace.iOrbits[d];
   OrbitIsospin & oe = threebodyspace.iOrbits[e];
   OrbitIsospin & of = threebodyspace.iOrbits[f];
 
   int P1 = oc.l%2;
   if(P1 != of.l%2) return;
 
   int J1 = oc.j;
   if(J1 != of.j) return;
 
   int P2 = (oa.l + ob.l)%2;
   if(P2 != (od.l + oe.l)%2) return;
 
   int ch = threebodyspace.GetChannelIndex(J2,P2,J1,P1,twoT);
   if (ch==-1) return;
   ThreeBodyChannelNO2B & ch_no2b = threebodyspace.ThreeBodyChannels[ch];
   
   int key_bra = ch_no2b.GetIndex(a,b,c,Tab);
   int key_ket = ch_no2b.GetIndex(d,e,f,Tde);
 
   int ph = ch_no2b.GetPhase(key_bra) * ch_no2b.GetPhase(key_ket);
   int ibra = ch_no2b.abct2n.at(key_bra);
   int iket = ch_no2b.abct2n.at(key_ket);
   auto index = threebodyspace.idx1d(ibra,iket);
 
   MatEl[ch][index] = StoreType(V*ph);
}


// Return a three-body matrix element where a,b are coupled to J2 and Tab, and d,e are coupled to J2 and Tde
// Tab,c are coupled to T3, and Tde,c are coupled to T3
// There is no total J quantum number, because it has been summed over with a weight 2J+1
template<class StoreType>
ThreeBodyStorage::ME_type ThreeBodyStorage_no2b<StoreType>::GetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int twoT) const
{

  ThreeBodyStorage::ME_type vout = 0;
//  if ( (a==b) and (Tab+J2)%2==0 ) return vout;
//  if ( (d==e) and (Tde+J2)%2==0 ) return vout;

  const OrbitIsospin & oa = threebodyspace.iOrbits[a];
  const OrbitIsospin & ob = threebodyspace.iOrbits[b];
  const OrbitIsospin & oc = threebodyspace.iOrbits[c];
  const OrbitIsospin & od = threebodyspace.iOrbits[d];
  const OrbitIsospin & oe = threebodyspace.iOrbits[e];
  const OrbitIsospin & of = threebodyspace.iOrbits[f];


  int P1 = oc.l%2;
  if(P1 != of.l%2) return vout;  // return 0.

  int J1 = oc.j;
  if(J1 != of.j) return vout;  // return 0.

  int P2 = (oa.l + ob.l)%2;
  if(P2 != (od.l + oe.l)%2) return vout;  // return 0.

  int ch = threebodyspace.GetChannelIndex(J2,P2,J1,P1,twoT);
  if (ch==-1) return vout;
  const ThreeBodyChannelNO2B & ch_no2b = threebodyspace.ThreeBodyChannels[ch];
  int key_bra = ch_no2b.GetIndex(a,b,c,Tab);
  int key_ket = ch_no2b.GetIndex(d,e,f,Tde);
  if(ch_no2b.iphase.find(key_bra) == ch_no2b.iphase.end()) return vout;  // return 0.
  if(ch_no2b.iphase.find(key_ket) == ch_no2b.iphase.end()) return vout; // return 0.
  int ph = ch_no2b.GetPhase(key_bra) * ch_no2b.GetPhase(key_ket);
  int ibra = ch_no2b.abct2n.at(key_bra);
  int iket = ch_no2b.abct2n.at(key_ket);

  auto index = threebodyspace.idx1d(ibra,iket);
  vout = ph * MatEl.at(ch)[ index ];

  return vout;
}

//  Get a 3b matrix element in pn mode, summed over total J with weight 2J+1
template<class StoreType>
ThreeBodyStorage::ME_type ThreeBodyStorage_no2b<StoreType>::GetME_pn_no2b(int a, int b, int c, int d, int e, int f, int J2) const
{
  Orbit & oa = modelspace->GetOrbit(a);
  Orbit & ob = modelspace->GetOrbit(b);
  Orbit & oc = modelspace->GetOrbit(c);
  Orbit & od = modelspace->GetOrbit(d);
  Orbit & oe = modelspace->GetOrbit(e);
  Orbit & of = modelspace->GetOrbit(f);


  ThreeBodyStorage::ME_type vout = 0;

  int P1 = oc.l%2;
  if(P1 != of.l%2) return vout; // return 0.

  int J1 = oc.j2;
  if(J1 != of.j2) return  vout; // return 0.

  int P2 = (oa.l + ob.l)%2;
  if(P2 != (od.l + oe.l)%2) return   vout; // return 0.

  int Z3 = oa.tz2 + ob.tz2 + oc.tz2; // total isospin projection, times 2
  if(Z3 != od.tz2 + oe.tz2 + of.tz2) return   vout; // return 0.

  // convert to the indices used interally
  int i1 = threebodyspace.nlj2idx.at({oa.n, oa.l, oa.j2});
  int i2 = threebodyspace.nlj2idx.at({ob.n, ob.l, ob.j2});
  int i3 = threebodyspace.nlj2idx.at({oc.n, oc.l, oc.j2});
  int i4 = threebodyspace.nlj2idx.at({od.n, od.l, od.j2});
  int i5 = threebodyspace.nlj2idx.at({oe.n, oe.l, oe.j2});
  int i6 = threebodyspace.nlj2idx.at({of.n, of.l, of.j2});

  for (int T12 : {0,1})
  {
    if(std::abs(oa.tz2 + ob.tz2) > 2*T12) continue;
    for (int T45 : {0,1})
    {
      if(std::abs(od.tz2 + oe.tz2) > 2*T45) continue;
//      for (int twoT=std::max( std::abs(2*T12-1), std::abs(2*T45-1) );
      for (int twoT=1;  twoT<= std::min( 2*T12+1, 2*T45+1 ); twoT+=2)
      {
        vout += GetME_iso_no2b(i1, i2, i3, T12, i4, i5, i6, T45, J2, twoT) *
             AngMom::CG(0.5, oa.tz2*0.5, 0.5, ob.tz2*0.5, T12, (oa.tz2+ob.tz2)*0.5) *
             AngMom::CG(0.5, od.tz2*0.5, 0.5, oe.tz2*0.5, T45, (od.tz2+oe.tz2)*0.5) *
             AngMom::CG(T12, (oa.tz2+ob.tz2)*0.5, 0.5, oc.tz2*0.5, twoT*0.5, Z3*0.5) *
             AngMom::CG(T45, (od.tz2+oe.tz2)*0.5, 0.5, of.tz2*0.5, twoT*0.5, Z3*0.5);
      }
    }
  }
  return vout;
}







template<class StoreType>
double ThreeBodyStorage_no2b<StoreType>::Norm() const 
{
  double norm = 0;
  for ( auto& itmat : MatEl )
  {
    for ( auto& it : itmat.second )
    {
      norm += it*it;
    }
  }
  return sqrt(norm);
}

template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Erase()  // Set all elements to zero
{
  for ( auto& itmat : MatEl )
  {
    for ( auto& it : itmat.second )
    {
      it = 0.;
    }
  }
}

template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::Deallocate() 
{
    std::map<int, std::vector<StoreType>>().swap( MatEl);
}

template<class StoreType>
size_t ThreeBodyStorage_no2b<StoreType>::size() const
{
   size_t thesize =0;
   for ( auto& itmat : MatEl )
   {
     thesize += itmat.second.size();
   }
   return thesize;
}



//////////////////////////////////////////////////////////////////////
////// File IO stuff
/////////////////////////////////////////////////////////////////


template<class StoreType>
size_t ThreeBodyStorage_no2b<StoreType>::CountME(int Emax_file, int E2max_file, int E3max_file, int Lmax_file, std::vector<OrbitIsospin>& file_Orbits) const
{
  double t_start = omp_get_wtime();
  size_t counter=0;
  int Norbs = file_Orbits.size();
  #pragma omp parallel for schedule(dynamic,1) reduction (+:counter)
  for (int i1=0; i1 < Norbs; i1++) {
    OrbitIsospin & o1 = file_Orbits[i1];
    int j1 = o1.j;
    int l1 = o1.l;
    int e1 = o1.e;
    if(e1 > emax) continue;
    if(e1 > Emax_file) continue;
    for (int i2=0; i2 <= i1; i2++) {
      OrbitIsospin & o2 = file_Orbits[i2];
      int j2 = o2.j;
      int l2 = o2.l;
      int e2 = o2.e;
      if(e2 > Emax_file) continue;
      if(e1 + e2 > E2max_file) continue;
      for (int i3=0; i3 < Norbs; i3++) {
        OrbitIsospin & o3 = file_Orbits[i3];
        int j3 = o3.j;
        int l3 = o3.l;
        int e3 = o3.e;
        if(e3 > Emax_file) continue;
        if(e2 + e3 > E2max_file) continue;
        if(e1 + e3 > E2max_file) continue;
        if(e1 + e2 + e3 > E3max_file) continue;

        for (int i4=0; i4 <= i1; i4++) {
          OrbitIsospin & o4 = file_Orbits[i4];
          int j4 = o4.j;
          int l4 = o4.l;
          int e4 = o4.e;
          if(e4 > Emax_file) continue;
          for (int i5=0; i5 <= i4; i5++) {
            OrbitIsospin & o5 = file_Orbits[i5];
            int j5 = o5.j;
            int l5 = o5.l;
            int e5 = o5.e;
            if(e5 > Emax_file) continue;
            if(e4 + e5 > E2max_file) continue;
            for (int i6=0; i6 < Norbs; i6++) {
              OrbitIsospin & o6 = file_Orbits[i6];
              int j6 = o6.j;
              int l6 = o6.l;
              int e6 = o6.e;
              if(j6 != j3) continue;
              if(l6 != l3) continue;
              if(e6 > Emax_file) continue;
              if(e4 + e6 > E2max_file) continue;
              if(e5 + e6 > E2max_file) continue;
              if(e4 + e5 + e6 > E3max_file) continue;
              if( (l1+l2+l3)%2 != (l4+l5+l6)%2 ) continue;


              int Jmin = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
              int Jmax =std::min( j1+j2, j4+j5 )/2;
              if (Jmin>Jmax) continue;
              size_t JT_block_size = (Jmax+1-Jmin) * ISOSPIN_BLOCK_DIMENSION; // ISOSPIN_BLOCK_DIMENSION inherited form base class

              counter += JT_block_size;

            }
          }
        }
      }
    }
  }
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return counter;
}





//  The inputs should be  (  {filename} ,  {Emax_file, E2max_file, E3max_file, Lmax_file} )   // where the last int argument is optional
template<class StoreType>
void ThreeBodyStorage_no2b<StoreType>::ReadFile( std::vector<std::string>& StringInputs, std::vector<int>& IntInputs )
{
   double t_start = omp_get_wtime();
   const size_t MAX_READ = 8* 1024*1024*1024L; // Read in 8 GB chunks
 
   if ( ( StringInputs.size()<1 ) or ( IntInputs.size()<2) )
   {
     std::cout << "ERROR! " << __FILE__ << "  " << __LINE__ << " " << __func__ << std::endl;
     std::cout << "Expected 1 string input and >=2 int input, but found " << StringInputs.size() << " and " << IntInputs.size() << std::endl;
     std::cout << "Dying." << std::endl;
     std::exit(EXIT_FAILURE);
   }
 
   std::string FileName = StringInputs[0];
   int Emax_file = IntInputs[0] ;
   int E2max_file = IntInputs[1] ;
   int E3max_file = IntInputs[2] ;
   int Lmax_file = Emax_file;
   if ( IntInputs.size() > 3 ) Lmax_file = IntInputs[3];
 
   std::cout << __func__ << ". from the input, I extracted " << FileName << "  " << Emax_file << " " << E2max_file << " " << E3max_file << " " << Lmax_file << std::endl;
 
   std::ifstream infile;
   boost::iostreams::filtering_istream zipstream;
   std::string filemode = "me3j";
   if(FileName.find("stream.bin") != std::string::npos)  filemode = "bin";
   else if (FileName.find(".gz") != std::string::npos) filemode = "gz";
 
   size_t nwords = sizeof(StoreType);
   std::cout << __func__ << "  reading/storing with " << 8*nwords << "  bit floats. filemode is " << filemode << std::endl;
 
 
   // Generate a vector of all the orbits present in the file.
   // This may be different from what's present in our model space
   std::vector<OrbitIsospin> file_Orbits;
   for (int e=0; e<=Emax_file; ++e)
   {
     int lmin = e%2;
     for (int l=lmin; l<=std::min(e,Lmax_file); l+=2) 
     {
       int n = (e-l)/2;
       int twojMin = std::abs(2*l-1);
       int twojMax = 2*l+1;
       for (int twoj=twojMin; twoj<=twojMax; twoj+=2) 
       {
         int idx = file_Orbits.size();
         file_Orbits.push_back( OrbitIsospin(idx,n,l,twoj));
       }
     }
   }
 
 
   size_t n_elem_to_read = CountME(Emax_file, E2max_file, E3max_file, Lmax_file, file_Orbits); // number of elements we want
 
   if ( filemode == "bin" ) // binary. check how big the file is.
   {
     infile = std::ifstream(FileName, std::ios::binary);
     infile.seekg(0,infile.end);
     size_t n_elem_in_file = infile.tellg();
     infile.seekg(0, infile.beg);
     n_elem_in_file -= infile.tellg();
     n_elem_in_file /= nwords;
     n_elem_to_read = std::min(n_elem_to_read, n_elem_in_file);
   }
   else if ( filemode == "gz")  // a gzipped me3j file
   {
     infile = std::ifstream(FileName, std::ios_base::in | std::ios_base::binary);
     zipstream.push(boost::iostreams::gzip_decompressor());
     zipstream.push(infile);
   }
   else if (filemode == "me3j" ) // a plain-text me3j file
   {
     infile = std::ifstream(FileName);
   }
 
 
   size_t buffer_size = std::min( MAX_READ,  n_elem_to_read);
 
   std::vector<StoreType> vbuf(buffer_size);
   std::cout << "Allocated a vector of size " << buffer_size << std::endl;
 
 
   int Norbs_file = file_Orbits.size();
 
   // Note that this is NOT a parallel for, but just a parallel. All threads do all the iterations
   #pragma omp parallel
   {
     size_t counter = 0;
     size_t total_counter = 0;
     size_t orbit_counter = 0;
     size_t num_threads = omp_get_num_threads();
     size_t this_thread = omp_get_thread_num();

     for (int iF1=0; iF1 < Norbs_file; iF1++)
     {
       OrbitIsospin & o1 = file_Orbits[iF1];
       int j1 = o1.j;
       int l1 = o1.l;
       int e1 = o1.e;
       if(e1 > emax) continue;
       if(e1 > Emax_file) continue;
       for (int iF2=0; iF2 <= iF1; iF2++) 
       {
         OrbitIsospin & o2 = file_Orbits[iF2];
         int j2 = o2.j;
         int l2 = o2.l;
         int e2 = o2.e;
         if(e2 > Emax_file) continue;
         if(e1 + e2 > E2max_file) continue;
         for (int iF3=0; iF3 < Norbs_file; iF3++)
         {
           OrbitIsospin & o3 = file_Orbits[iF3];
           int j3 = o3.j;
           int l3 = o3.l;
           int e3 = o3.e;
           if(e3 > Emax_file) continue;
           if(e2 + e3 > E2max_file) continue;
           if(e1 + e3 > E2max_file) continue;
           if(e1 + e2 + e3 > E3max_file) continue;
   
           for (int iF4=0; iF4 <= iF1; iF4++)
           {
             OrbitIsospin & o4 = file_Orbits[iF4];
             int j4 = o4.j;
             int l4 = o4.l;
             int e4 = o4.e;
             if(e4 > Emax_file) continue;
             for (int iF5=0; iF5 <= iF4; iF5++)
             {
               OrbitIsospin & o5 = file_Orbits[iF5];
               int j5 = o5.j;
               int l5 = o5.l;
               int e5 = o5.e;
               if(e5 > Emax_file) continue;
               if(e4 + e5 > E2max_file) continue;
               for (int iF6=0; iF6 < Norbs_file; iF6++)
               {
                 OrbitIsospin & o6 = file_Orbits[iF6];
                 int j6 = o6.j;
                 int l6 = o6.l;
                 int e6 = o6.e;
                 if(j6 != j3) continue;
                 if(l6 != l3) continue;
                 if(e6 > Emax_file) continue;
                 if(e4 + e6 > E2max_file) continue;
                 if(e5 + e6 > E2max_file) continue;
                 if(e4 + e5 + e6 > E3max_file) continue;
                 if( (l1+l2+l3)%2 != (l4+l5+l6)%2 ) continue;
   
                 orbit_counter +=1;

                 int Jmin = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
                 int Jmax =std::min( j1+j2, j4+j5 )/2;
//                 int JT_block_size = (Jmax+1-Jmin) * 5; // 5 comes from the 5 possible isospin combinations 001 011 101 111 113
//                 int JT_block_size = (Jmax+1-Jmin) * ISOSPIN_BLOCK_DIMENSION ; 
//                 if (Jmin>Jmax) JT_block_size = 0;
//   
//                 // If this thread isn't going to do the storing, and if we're not going to reload the buffer, just skip ahead
//                 if ( (orbit_counter % num_threads != this_thread) and ( (counter + JT_block_size) < buffer_size ) )
//                 {
//                     counter += JT_block_size;
//                     continue;
//                 }
//   
   
//                 for (int J = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
//                     J <= std::min( j1+j2, j4+j5 )/2; J++)
                 for (int J = Jmin; J <= Jmax; J++)
                 {
                   for (int T12: {0,1})
                   {
                     for (int T45: {0,1})
                     {
                       for (int T3 = 1; T3 <= std::min( 2*T12+1, 2*T45+1 ); T3+=2)
                       {
                         if(counter == buffer_size) counter = 0;
                         if(counter == 0)
                         {
   
                           // All threads stop when we've run out of buffer, and thread 0 reads new data into the buffer
                           #pragma omp barrier 
                             if ( this_thread==0 )
                             {
                             // fill the buffer
                             size_t n_this_pass = std::min( buffer_size, n_elem_to_read-total_counter );
                             if (filemode == "bin")        infile.read((char*)&vbuf[0], n_this_pass*nwords);
                             else if (filemode == "gz")    for (size_t iread=0;iread<n_this_pass;iread++) zipstream >> vbuf[iread];
                             else if (filemode == "me3j")  for (size_t iread=0;iread<n_this_pass;iread++) infile >> vbuf[iread];
                             }
                           #pragma omp barrier
                           // back to work, little threads...
   
                         }
                         counter += 1;
                         total_counter += 1;
   
                         if (orbit_counter % num_threads != this_thread) continue; // This is how we divide the jobs among the threads
   
                         if( (e1 > emax) or (e2 > emax) or (e3 > emax) or (e4 > emax) or (e5 > emax) or (e6 > emax)) continue;
                         if( std::max({e1,e2,e3,e4,e5,e6}) > emax) continue;
                         if( ((e1+e2) > E2max) or ((e1+e3) > E2max) or ((e2+e3) > E2max) or ((e4+e5) > E2max) or ((e4+e6) > E2max) or ((e5+e6) > E2max)) continue;
                         if( (e1+e2+e3 > E3max) or ( e4+e5+e6 > E3max) ) continue;
                         if( total_counter > n_elem_to_read) break;
   
                         int i1 = threebodyspace.nlj2idx.at({o1.n,o1.l,o1.j});
                         int i2 = threebodyspace.nlj2idx.at({o2.n,o2.l,o2.j});
                         int i3 = threebodyspace.nlj2idx.at({o3.n,o3.l,o3.j});
                         int i4 = threebodyspace.nlj2idx.at({o4.n,o4.l,o4.j});
                         int i5 = threebodyspace.nlj2idx.at({o5.n,o5.l,o5.j});
                         int i6 = threebodyspace.nlj2idx.at({o6.n,o6.l,o6.j});
   
                         StoreType vset = vbuf[counter-1];
   
   
                         if( (i1==i2 and (J+T12)%2 ==0 ) or ( i4==i5 and (J+T45)%2 ==0 ) ) {
                           if( abs(vset) > 1.e-6 ){
                             std::cout << "Warning: something wrong, this three-body matrix element has to be zero" << std::endl;
                             std::cout <<  std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
                                           std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
                                           std::setw(4) << J << std::setw(4) << T3 << std::setw(16) << counter <<
                                           std::setw(12) << std::setprecision(6) << vset << std::endl;
                           }
                         }
   
                         if ( std::abs(vset)>1e-8)
                         {
                          SetME_iso_no2b(i1, i2, i3, T12, i4, i5, i6, T45, J, T3, vset );
                         }
   
                       }// for T3
                     }// for T45
                   }// for T12
                 }// for J
               }// for i6
             }// for i5
           }// for i4
         }// for i3
       }// for i2
     }// for i1
 
   }// end of parallel outer block
 
   std::cout << "Done reading" << std::endl;
   IMSRGProfiler::timer["ThreeBodyStorage_no2b::ReadFile"] += omp_get_wtime() - t_start;
}



// Make our concrete instantiations of the ThreeBodyStorage_no2b class. We use either single-precision or half-precision type

template class ThreeBodyStorage_no2b<ME_single_type>;
template class ThreeBodyStorage_no2b<ME_half_type>;




