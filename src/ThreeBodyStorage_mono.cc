
#include "ThreeBodyStorage_mono.hh"
#include "AngMom.hh"


#include <omp.h>



//OrbitIsospin::OrbitIsospin(int idx, int n, int l, int j)
//  : idx(idx), n(n), l(l), j(j), e(2*n+l)
//{}


ThreeBodyChannelMono::ThreeBodyChannelMono()
  : j1(-1), j2(-1), j3(-1), l1(-1), l2(-1), l3(-1), twoT(-1), Ndim(-1)
//  : j1(-1), j2(-1), j3(-1), l1(-1), l2(-1), l3(-1), T3(-1), Ndim(-1)
//  : J2(-1), P2(-1), J1(-1), P1(-1), T3(-1), Ndim(-1)
{}

//ThreeBodyChannelMono::ThreeBodyChannelMono(int J2, int P2, int J1, int P1, int T3,  ThreeBodySpaceMono& thr)
//  : J2(J2), P2(P2), J1(J1), P1(P1), T3(T3)
ThreeBodyChannelMono::ThreeBodyChannelMono(int j1, int j2, int j3, int l1, int l2, int l3, int twoT,  ThreeBodySpaceMono& thr)
  : j1(j1), j2(j2), j3(j3), l1(l1), l2(l2), l3(l3), twoT(twoT)
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
    if(oa.j != j1) continue;
    if(oa.l != l1) continue;
    for (int ib=0; ib < Norbs; ib++)
    {
      OrbitIsospin & ob = thr.GetIsospinOrbit(ib);
      if(ob.e > emax) continue;
      if(oa.e + ob.e > e2max) continue;
      if(ob.j != j2) continue;
      if(ob.l != l2) continue;
//      if( not AngMom::Triangle( oa.j, ob.j, 2*J2) ) continue;
//      if(std::abs(oa.j-ob.j) > 2*J2) continue;
//      if(        (oa.j+ob.j) < 2*J2) continue;
//      if((oa.l+ob.l)%2 != P2) continue;
      for (int ic=0; ic < Norbs; ic++)
      {
        OrbitIsospin & oc = thr.GetIsospinOrbit(ic);
//        if(oc.j != J1) continue;
//        if(oc.l%2 != P1) continue;
        if(oc.e > emax) continue;
        if(oa.e + oc.e > e2max) continue;
        if(ob.e + oc.e > e2max) continue;
        if(oa.e + ob.e + oc.e > e3max) continue;
        if(oc.j != j3) continue;
        if(oc.l != l3) continue;
        for (int T12: {0,1})
        {
//          if( std::abs(2*T2-1) > T3 ) continue;
//          if(         (2*T2+1) < T3 ) continue;
          if( not AngMom::Triangle( 2*T12, 1, twoT) ) continue;
//          if(((J2+T2)%2 == 0) and (ia == ib) ) continue;
//          int ph = 1 - 2*(((oa.j + ob.j)/2 + J2 + T2)%2);
//          int ph = 1; // TODO check this !!!
          int key_ab = Hash_abct(ia, ib, ic, T12);
//          int key_ba = Hash_abct(ib, ia, ic, T12);
//          std::cout << "key " << ia << " " << ib << " "  << ic << " " << T12 << "  is  " << key_ab << " , " << key_ba << "   Ndim = " << Ndim << std::endl;
//          std::cout << "    ( ja,jb,jc ,la,lb,lc = " << oa.j << " " << ob.j << " " << oc.j << "  " << oa.l << " " << ob.l << " " << oc.l << std::endl;
          abct2n[key_ab] = Ndim;
//          abct2n[key_ba] = Ndim;
          iphase[key_ab] = 1;
//          iphase[key_ba] = ph;
          Ndim += 1;
        }// for T2
      }// for ic
    }// for ib
  }// for ia
}

int ThreeBodyChannelMono::Hash_abct(int a, int b, int c, int Tab) const
{
  return ( a + (b<<10) + (c<<20) + (Tab<<30) );
}

void ThreeBodyChannelMono::UnHash_abct(int key, int & a, int & b, int & c, int & Tab) const
{
  int Lowest_ten_bits = 0x3FF;
  a = ( (key>> 0) & Lowest_ten_bits);
  b = ( (key>>10) & Lowest_ten_bits);
  c = ( (key>>20) & Lowest_ten_bits);
  Tab = ( (key>>30) & Lowest_ten_bits);
}



//ThreeBodySpaceMono::~ThreeBodySpaceMono()
//{}

ThreeBodySpaceMono::ThreeBodySpaceMono()
  :  ThreeBodyChannels() , NChannels(-1)
{}


ThreeBodySpaceMono::ThreeBodySpaceMono( int emax, int e2max, int e3max, int lmax)
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
//  std::cout << "iOrbits are" << std::endl;
//  for (size_t ii=0;ii<iOrbits.size(); ii++)
//  {
//     auto& oi = iOrbits[ii];
//     std::cout << ii << " :  " << oi.idx << " " << oi.n << " " << oi.l << " " << oi.j << "  " << oi.e << std::endl;
//  }

//  std::cout << "To begin with, Hashing 1,1,1,0,0,0,1 gives " << Hash_Channel(1,1,1,0,0,0,1) 
//            << "   while hashing 5,3,7,2,1,3,1 gives " << Hash_Channel(5,3,7,2,1,3,1) << std::endl;
//  std::cout << " break it down: " << ((1+1)/2-0) << " + " << ((1+1)/2-0) << " <<1 + " << ((1+1)/2-0) << " <<2  + " << (1/2) << " <<3 + " << ( 0 +Lmax*0 + Lmax*Lmax*0) << "<<4 " << std::endl;
//  std::cout << "   should be " <<  ((1+1)/2-0 ) << " + " << ( ((1+1)/2-0)<<1) << " + " << (((1+1)/2-0)<<2) << " + " << ((1/2)<<3) << " + " << (0<<4) << std::endl;
//  std::cout << "   should be " <<  ((1+1)/2-0 )  +  (((1+1)/2-0)<<1)  +  (((1+1)/2-0)<<2)  + ((1/2)<<3) << " + " << (0<<4) << std::endl;
//  std::cout << "   but it's " <<  (  ((1+1)/2-0 )  +  (((1+1)/2-0)<<1)  +  (((1+1)/2-0)<<2)  + ((1/2)<<3)  +  (0<<4)  ) << std::endl;

// Allocate the ThreeBody channels

  int count = 0;
  for (int l1=0; l1<=lmax; l1+=1)
  {
    int twoj1Min = std::abs(2*l1-1);
    int twoj1Max = 2*l1+1;
    for (int j1=twoj1Min; j1<=twoj1Max; j1+=2)
    {
      for (int l2=0; l2<=lmax; l2+=1)
      {
        int twoj2Min = std::abs(2*l2-1);
        int twoj2Max = 2*l2+1;
        for (int j2=twoj2Min; j2<=twoj2Max; j2+=2)
        {
          for (int l3=0; l3<=lmax; l3+=1)
          {
            int twoj3Min = std::abs(2*l3-1);
            int twoj3Max = 2*l3+1;
            for (int j3=twoj3Min; j3<=twoj3Max; j3+=2)
            {
              for ( int twoT : {1,3} )
              {
                ThreeBodyChannelMono channel(j1,j2,j3, l1,l2,l3, twoT, *this);
                if(channel.Ndim < 1) continue;

                ThreeBodyChannels.push_back(channel);
//                int key = Hash_Channel(J2,P2,J1,P1,T3);
                int key = Hash_Channel(j1,j2,j3,l1,l2,l3,twoT);
//                std::cout << "***ThreeBodyChannel " << j1 << " " << j2 << " " << j3 << " " << l1 << " " << l2 << " " << l3 << " " << twoT << "  key = " << key << "  count = " << count << std::endl;
                idcs2ch[key] = count;
                count += 1;
              }
            }
          }
        }
      }
    }


  }
//  std::cout << "#####   idcs2ch  ######" << std::endl;
//  for (auto& iter : idcs2ch )
//  {
//     int J1,J2,J3,L1,L2,L3,TT;
//     UnHash_Channel(iter.first,  J1,J2,J3,L1,L2,L3,TT);
//     std::cout << " " << iter.first << "  =>  " << iter.second << "  ;   " << J1 << " " << J2 << " " << J3 << " " << L1 << " " << L2 << " " << L3 << " " << TT << std::endl;
//  }


//  int J2max = std::min(2*lmax+1, e2max+1);
//  int count = 0;
//  for (int J2=0; J2<=J2max; ++J2)
//  {
//    for (int P2: {0,1})
//    {
//      for (int J1=1; J1<=2*lmax+1; J1+=2)
//      {
//        for (int P1: {0,1})
//        {
//          for (int T3 : {1,3})
//          {
////            ThreeBodyChannelMono channel(J2, P2, J1, P1, T3, *this);
//            ThreeBodyChannelMono channel(J2, P2, J1, P1, T3, *this);
//            if(channel.Ndim < 1) continue;
//
//            ThreeBodyChannels.push_back(channel);
//            int key = Hash_Channel(J2,P2,J1,P1,T3);
//            idcs2ch[key] = count;
//            count += 1;
//          }
//        }
//      }
//    }
//  }
  NChannels = ThreeBodyChannels.size();
  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
}



//int ThreeBodySpaceMono::Hash_Channel(int J2, int P2, int J1, int P1, int T3) const
int ThreeBodySpaceMono::Hash_Channel(int j1, int j2, int j3, int l1, int l2, int l3, int twoT) const
{
  // 1st bit: 1 if j1>l1 0 if j1<l1, 2nd,3rd bits same for j2,j3.
  // 4th bit: 0 for twoT=1, 1 for twoT=3
  // all other bits give l1,l2,l3

//  std::cout << "@$@$ DEBUG @$@$ " << __func__ << "   " << j1 << " " << j2 << " " << j3 << " " << l1 << " " << l2 << " " << l3 << " " << twoT << "  |  "
//            << ((j1+1)/2-l1)  << "  " << (((j2+1)/2-l2)<<1)  << " " << (((j3+1)/2-l3)<<2) << " " << ((twoT/2)<<3) <<" "  << (( l1 +Lmax*l2 + Lmax*Lmax*l3)<<4) << "   " << "   ::  "
//            << ((j1+1)/2-l1) +           (((j2+1)/2-l2)<<1) +          (((j3+1)/2-l3)<<2)  + ((twoT/2)<<3) + (( l1 +Lmax*l2 + Lmax*Lmax*l3)<<4) << "  Lmax= "<<Lmax << std::endl;


  return ((j1+1)/2-l1) + (((j2+1)/2-l2)<<1) + (((j3+1)/2-l3)<<2)  + ((twoT/2)<<3) + (( l1 +(Lmax+1)*l2 + (Lmax+1)*(Lmax+1)*l3)<<4);  

}

void ThreeBodySpaceMono::UnHash_Channel(int key, int& j1, int& j2, int& j3, int& l1, int& l2, int& l3, int& twoT) const
{
  l1   = ((key>>4) % (Lmax+1) );
  l2   = (((key>>4)/(Lmax+1))%(Lmax+1));
  l3   = ((key>>4) / ((Lmax+1)*(Lmax+1)));
  twoT = ((key>>3)&0x1) * 2 + 1 ;
  j1   = ( key    &0x1) * 2 + 2*l1 -1;
  j2   = ((key>>1)&0x1) * 2 + 2*l2 -1;
  j3   = ((key>>2)&0x1) * 2 + 2*l3 -1;
}

// Take Jab,Pab,Jc,Pc,T and return the channel number for subsequent lookup in the vector ThreeBodyChannels
//int ThreeBodySpaceMono::GetChannelIndex(int Jab, int Pab, int Jc, int Pc, int T)const
int ThreeBodySpaceMono::GetChannelIndex(int j1, int j2, int j3, int l1, int l2, int l3, int twoT)const
{
//  auto key = Hash_Channel(Jab, Pab, Jc, Pc, T);
  auto key = Hash_Channel(j1,j2,j3,l1,l2,l3,twoT);
  auto iter = idcs2ch.find(key);
  int ch =-1;
  if ( iter!=idcs2ch.end() ) ch = iter->second;
//  std::cout << "@$@$@ DEBUG @$@$ " << __func__ << "   " << j1 << " " << j2 << " " << j3 << " " << l1 << " " << l2 << " " << l3 << " " << twoT << "   -> key= " << key << "  => ch= " << ch << std::endl;
  return ch;
}


//////////////////////////////////////////////////////////////////////////////////////////
//////  Begin implementation of ThreeBodyStorage_mono methods. This is the main event.
////////////////////////////////////////////////////////////////////////////////////////

template<class StoreType>
ThreeBodyStorage_mono<StoreType>::ThreeBodyStorage_mono( const ThreeBodyStorage_mono<StoreType>& TBS_in )
//: MatEl(TBS_in.MatEl), threebodyspace(TBS_in.threebodyspace),    ThreeBodyStorage( TBS_in )
:  ThreeBodyStorage( TBS_in ),  MatEl(TBS_in.MatEl), threebodyspace(TBS_in.threebodyspace)
{}


template<class StoreType>
std::unique_ptr<ThreeBodyStorage> ThreeBodyStorage_mono<StoreType>::Clone() const { return std::unique_ptr<ThreeBodyStorage>( new ThreeBodyStorage_mono<StoreType>( *this)); };
//std::shared_ptr<ThreeBodyStorage> ThreeBodyStorage_mono<StoreType>::Clone() const { return std::shared_ptr<ThreeBodyStorage>( new ThreeBodyStorage_mono<StoreType>( *this)); };


template<class StoreType>
void ThreeBodyStorage_mono<StoreType>::Allocate()
{
  // Create the three body space
  threebodyspace = ThreeBodySpaceMono(emax,E2max,E3max, std::min(emax,lmax) );

  // Allocate the vectors for the matrix element storage
  for (int ch=0; ch<threebodyspace.NChannels; ch++)
  {
    ThreeBodyChannelMono& ch_mono = threebodyspace.ThreeBodyChannels.at(ch);
    size_t n = ch_mono.Ndim;

    MatEl[ch] = std::vector<StoreType>( n*(n+1)/2, StoreType(0.0));
  }
  is_allocated = true;
}





template<class StoreType>
void ThreeBodyStorage_mono<StoreType>::Multiply(const double rhs) 
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
void ThreeBodyStorage_mono<StoreType>::Add(const ThreeBodyStorage& rhs) 
{
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_mono*)(&rhs))->MatEl;

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
void ThreeBodyStorage_mono<StoreType>::Subtract(const ThreeBodyStorage& rhs) 
{
    if (rhs.GetStorageMode() == this->GetStorageMode() )
    {
       auto& rhsMatEl = ((ThreeBodyStorage_mono*)(&rhs))->MatEl;

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
void ThreeBodyStorage_mono<StoreType>::SetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int twoT, ThreeBodyStorage::ME_type V)
{
   OrbitIsospin & oa = threebodyspace.iOrbits[a];
   OrbitIsospin & ob = threebodyspace.iOrbits[b];
   OrbitIsospin & oc = threebodyspace.iOrbits[c];
   OrbitIsospin & od = threebodyspace.iOrbits[d];
   OrbitIsospin & oe = threebodyspace.iOrbits[e];
   OrbitIsospin & of = threebodyspace.iOrbits[f];
 
//   if (Tab != Tde) return;
//   if ( oa.j2 != od.j2 ) return;
//   if ( ob.j2 != oe.j2 ) return;
//   if ( oc.j2 != of.j2 ) return;
//   if ( oa.l != od.l ) return;
//   if ( ob.l != oe.l ) return;
//   if ( oc.l != of.l ) return;

//   int P1 = oc.l%2;
//   if(P1 != of.l%2) return;
// 
//   int J1 = oc.j;
//   if(J1 != of.j) return;
// 
//   int P2 = (oa.l + ob.l)%2;
//   if(P2 != (od.l + oe.l)%2) return;
 
//   int ch = threebodyspace.GetChannelIndex(J2,P2,J1,P1,twoT);
   int ch = threebodyspace.GetChannelIndex(oa.j,ob.j,oc.j,oa.l,ob.l,oc.l,twoT);
   if (ch==-1) return;
   ThreeBodyChannelMono & ch_mono = threebodyspace.ThreeBodyChannels[ch];
//   std::cout << "channel is " << ch << "  -> " << oa.j << " " << ob.j << " " << oc.j << " " << oa.l << " " << ob.l << " " << oc.l << " " << twoT << std::endl;
   
//   std::cout << "abcTab = " << a << " " << b << " " << c << " " << Tab << std::endl;
//   int key_bra = ch_mono.GetIndex(a,b,c,Tab);
//   int key_ket = ch_mono.GetIndex(d,e,f,Tde);
   int key_bra = ch_mono.GetIndex(oa.idx,ob.idx,oc.idx,Tab);
   int key_ket = ch_mono.GetIndex(od.idx,oe.idx,of.idx,Tde);
//   std::cout << "key_bra, key_ket = " << key_bra << " " << key_ket << std::endl;
 
//   int ph = ch_mono.GetPhase(key_bra) * ch_mono.GetPhase(key_ket);
   int ph = 1;
//   std::cout << "ph = " << ph << std::endl;
//   std::cout << "size of the abct2n for this channel (ch=" << ch << ") is " << ch_mono.abct2n.size() << std::endl;
//   for ( auto& itabc : ch_mono.abct2n )
//   {
//      std::cout << "     " << itabc.first << " => " << itabc.second << std::endl;
//   }
   int ibra = ch_mono.abct2n.at(key_bra);
   int iket = ch_mono.abct2n.at(key_ket);
//   std::cout << "   ibra,iket = " << ibra << " " << iket << std::endl;
   auto index = threebodyspace.idx1d(ibra,iket);
 
   MatEl[ch][index] = StoreType(V*ph);
}


// Return a three-body matrix element where a,b are coupled to J2 and Tab, and d,e are coupled to J2 and Tde
// Tab,c are coupled to T3, and Tde,c are coupled to T3
// There is no total J quantum number, because it has been summed over with a weight 2J+1
template<class StoreType>
//ThreeBodyStorage::ME_type ThreeBodyStorage_mono<StoreType>::GetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int twoT) const
ThreeBodyStorage::ME_type ThreeBodyStorage_mono<StoreType>::GetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int twoT) const
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


//  int ch = threebodyspace.GetChannelIndex(J2,P2,J1,P1,twoT);
  int ch = threebodyspace.GetChannelIndex(oa.j,ob.j,oc.j,oa.l,ob.l,oc.l,twoT);
  if (ch==-1) return vout;
  const ThreeBodyChannelMono & ch_mono = threebodyspace.ThreeBodyChannels[ch];
  int key_bra = ch_mono.GetIndex(oa.idx,ob.idx,oc.idx,Tab);
  int key_ket = ch_mono.GetIndex(od.idx,oe.idx,of.idx,Tde);
  if ( oa.idx < od.idx ) { std::swap(key_bra,key_ket); };
  if(ch_mono.iphase.find(key_bra) == ch_mono.iphase.end()) return vout;  // return 0.
  if(ch_mono.iphase.find(key_ket) == ch_mono.iphase.end()) return vout; // return 0.
  int ph = ch_mono.GetPhase(key_bra) * ch_mono.GetPhase(key_ket);
  int ibra = ch_mono.abct2n.at(key_bra);
  int iket = ch_mono.abct2n.at(key_ket);

  auto index = threebodyspace.idx1d(ibra,iket);
  vout = ph * MatEl.at(ch)[ index ];

  return vout;
}

//  Get a 3b matrix element in pn mode, summed over total J with weight 2J+1
template<class StoreType>
//ThreeBodyStorage::ME_type ThreeBodyStorage_mono<StoreType>::GetME_pn_mono(int a, int b, int c, int d, int e, int f, int J2) const
ThreeBodyStorage::ME_type ThreeBodyStorage_mono<StoreType>::GetME_pn_mono(int a, int b, int c, int d, int e, int f) const
{


  Orbit & oa = modelspace->GetOrbit(a);
  Orbit & ob = modelspace->GetOrbit(b);
  Orbit & oc = modelspace->GetOrbit(c);
  Orbit & od = modelspace->GetOrbit(d);
  Orbit & oe = modelspace->GetOrbit(e);
  Orbit & of = modelspace->GetOrbit(f);

  int i1 = threebodyspace.nlj2idx.at({oa.n, oa.l, oa.j2});
  int i2 = threebodyspace.nlj2idx.at({ob.n, ob.l, ob.j2});
  int i3 = threebodyspace.nlj2idx.at({oc.n, oc.l, oc.j2});
  int i4 = threebodyspace.nlj2idx.at({od.n, od.l, od.j2});
  int i5 = threebodyspace.nlj2idx.at({oe.n, oe.l, oe.j2});
  int i6 = threebodyspace.nlj2idx.at({of.n, of.l, of.j2});


//  int iA = i1;
//  int iB = i2;
//  int iC = i3;
//  int iD = i4;
//  int iE = i5;
//  int iF = i6;

//  std::cout << "      before resuffling, isospin ordering is " << iA << " " << iB << " " << iC << " " << iD << " " << iE << " " << iF << std::endl;
//  std::cout << "    bra side: " << oa.n << " " << oa.l << " " << oa.j2 << " , " << ob.n << " " << ob.l << " " << ob.j2 << " , " << oc.n << " " << oc.l << " " << oc.j2 << std::endl;



//  if ( iA > iB ) { std::swap(iA,iB); std::swap(iD,iE) ; std::swap(a,b);  std::swap(d,e); };
//  if ( iB > iC ) { std::swap(iB,iC); std::swap(iE,iF) ; std::swap(b,c);  std::swap(e,f); };
//  if ( iA > iB ) { std::swap(iA,iB); std::swap(iD,iE) ; std::swap(a,b);  std::swap(d,e); };



//  oa = modelspace->GetOrbit(a);
//  ob = modelspace->GetOrbit(b);
//  oc = modelspace->GetOrbit(c);
//  od = modelspace->GetOrbit(d);
//  oe = modelspace->GetOrbit(e);
//  of = modelspace->GetOrbit(f);


  ThreeBodyStorage::ME_type vout = 0;

//  int P1 = oc.l%2;
//  if(P1 != of.l%2) return vout; // return 0.
//
//  int J1 = oc.j2;
//  if(J1 != of.j2) return  vout; // return 0.
//
//  int P2 = (oa.l + ob.l)%2;
//  if(P2 != (od.l + oe.l)%2) return   vout; // return 0.
  int tza = oa.tz2;
  int tzb = ob.tz2;
  int tzc = oc.tz2;
  int tzd = od.tz2;
  int tze = oe.tz2;
  int tzf = of.tz2;


  if (  i1 > i4 ) { std::swap(i1,i4); std::swap(i2,i5); std::swap(i3,i6); std::swap(tza,tzd); std::swap(tzb,tze); std::swap(tzc,tzf);};
//  std::cout << "abc " << a << " " << b <<  " " << c << "   tz's  " << tza << " " << tzb << " " << tzc << std::endl;


  int Z3 = tza+tzb+tzc; // total isospin projection, times 2
  if(Z3 != tzd+tze+tzf) return   vout; // return 0.
//  int Z3 = oa.tz2 + ob.tz2 + oc.tz2; // total isospin projection, times 2
//  if(Z3 != od.tz2 + oe.tz2 + of.tz2) return   vout; // return 0.

  // convert to the indices used interally
//  i1 = threebodyspace.nlj2idx.at({oa.n, oa.l, oa.j2});
//  i2 = threebodyspace.nlj2idx.at({ob.n, ob.l, ob.j2});
//  i3 = threebodyspace.nlj2idx.at({oc.n, oc.l, oc.j2});
//  i4 = threebodyspace.nlj2idx.at({od.n, od.l, od.j2});
//  i5 = threebodyspace.nlj2idx.at({oe.n, oe.l, oe.j2});
//  i6 = threebodyspace.nlj2idx.at({of.n, of.l, of.j2});




//  std::cout << "  -----------------   " << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 << " " << i6  << "  Z3 = " << Z3 << std::endl;

  for (int T12 : {0,1})
  {
    if(std::abs(tza + tzb) > 2*T12) continue;
    for (int T45 : {0,1})
    {
      if(std::abs(tzd + tze) > 2*T45) continue;
//      for (int twoT=std::max( std::abs(2*T12-1), std::abs(2*T45-1) );
      for (int twoT=std::abs(Z3);  twoT<= std::min( 2*T12+1, 2*T45+1 ); twoT+=2)
      {
        vout += GetME_iso_mono(i1, i2, i3, T12, i4, i5, i6, T45, twoT) *
             AngMom::CG(0.5, tza*0.5, 0.5, tzb*0.5, T12, (tza+tzb)*0.5) *
             AngMom::CG(0.5, tzd*0.5, 0.5, tze*0.5, T45, (tzd+tze)*0.5) *
             AngMom::CG(T12, (tza+tzb)*0.5, 0.5, tzc*0.5, twoT*0.5, Z3*0.5) *
             AngMom::CG(T45, (tzd+tze)*0.5, 0.5, tzf*0.5, twoT*0.5, Z3*0.5);

//             std::cout << "         " << T12 << " " << T45 << " " << twoT << "     " << GetME_iso_mono(i1, i2, i3, T12, i4, i5, i6, T45, twoT) << " * "
//                       <<  AngMom::CG(0.5, tza*0.5, 0.5, tzb*0.5, T12, (tza+tzb)*0.5) *
//             AngMom::CG(0.5, tzd*0.5, 0.5, tze*0.5, T45, (tzd+tze)*0.5) *
//             AngMom::CG(T12, (tza+tzb)*0.5, 0.5, tzc*0.5, twoT*0.5, Z3*0.5) *
//             AngMom::CG(T45, (tzd+tze)*0.5, 0.5, tzf*0.5, twoT*0.5, Z3*0.5) << std::endl;
      }
    }
  }
  return vout;
}







template<class StoreType>
double ThreeBodyStorage_mono<StoreType>::Norm() const 
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
void ThreeBodyStorage_mono<StoreType>::Erase()  // Set all elements to zero
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
void ThreeBodyStorage_mono<StoreType>::Deallocate() 
{
    std::map<int, std::vector<StoreType>>().swap( MatEl);
}

template<class StoreType>
size_t ThreeBodyStorage_mono<StoreType>::size() const
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
size_t ThreeBodyStorage_mono<StoreType>::CountME(int Emax_file, int E2max_file, int E3max_file, int Lmax_file, std::vector<OrbitIsospin>& file_Orbits) const
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
//    if(e1 > emax) continue;
    if(e1 > Emax_file) continue;
//    for (int i2=0; i2 <= i1; i2++) {
    for (int i2=0; i2 < Norbs; i2++) {
      OrbitIsospin & o2 = file_Orbits[i2];
      int j2 = o2.j;
      int l2 = o2.l;
      int e2 = o2.e;
//      if(e2 > Emax_file) continue;
      if(e1 + e2 > E2max_file) continue;
      for (int i3=0; i3 < Norbs; i3++) {
        OrbitIsospin & o3 = file_Orbits[i3];
        int j3 = o3.j;
        int l3 = o3.l;
        int e3 = o3.e;
//        if(e3 > Emax_file) continue;
        if(e2 + e3 > E2max_file) continue;
        if(e1 + e3 > E2max_file) continue;
        if(e1 + e2 + e3 > E3max_file) continue;

        for (int i4=0; i4 <= i1; i4++) {
          OrbitIsospin & o4 = file_Orbits[i4];
          int j4 = o4.j;
          int l4 = o4.l;
          int e4 = o4.e;
          if ( j1 != j4 ) continue;
          if ( l1 != l4 ) continue;
//          if(e4 > Emax_file) continue;
//          for (int i5=0; i5 <= i4; i5++) {
          for (int i5=0; i5 < Norbs; i5++) {
            OrbitIsospin & o5 = file_Orbits[i5];
            int j5 = o5.j;
            int l5 = o5.l;
            int e5 = o5.e;
            if ( j2 != j5 ) continue;
            if ( l2 != l5 ) continue;
//            if(e5 > Emax_file) continue;
            if(e4 + e5 > E2max_file) continue;
            for (int i6=0; i6 < Norbs; i6++) {
              OrbitIsospin & o6 = file_Orbits[i6];
              int j6 = o6.j;
              int l6 = o6.l;
              int e6 = o6.e;
              if(j6 != j3) continue;
              if(l6 != l3) continue;
//              if(e6 > Emax_file) continue;
              if(e4 + e6 > E2max_file) continue;
              if(e5 + e6 > E2max_file) continue;
              if(e4 + e5 + e6 > E3max_file) continue;
              if( (l1+l2+l3)%2 != (l4+l5+l6)%2 ) continue;


//              int Jmin = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
//              int Jmax =std::min( j1+j2, j4+j5 )/2;
//              if (Jmin>Jmax) continue;
//              size_t JT_block_size = (Jmax+1-Jmin) * ISOSPIN_BLOCK_DIMENSION; // ISOSPIN_BLOCK_DIMENSION inherited form base class

//              counter += JT_block_size;
//              counter += ISOSPIN_BLOCK_DIMENSION; // This should be 5, for  (t2_bra,t2_ket,2T) = { (0,0,1), (0,1,1), (1,0,1), (1,1,1), (1,1,3) } , in that order
              counter += 5;

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
void ThreeBodyStorage_mono<StoreType>::ReadFile( std::vector<std::string>& StringInputs, std::vector<int>& IntInputs )
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
   std::cout << "n elem to read = " << n_elem_to_read << std::endl;
 
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
 
   std::cout << "filemode is " << filemode << std::endl;
   if (filemode=="me3j" or filemode=="gz")
   {
     int linesize= 496;
     char line[linesize];
     if (filemode=="gz")
     {
       zipstream.getline(line,linesize);
     }
     else
     {
       infile.getline(line,linesize);
     }
     std::cout << "the header was : " << line << std::endl;
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
//       std::cout << "   iF1 = " << iF1 << std::endl;
       OrbitIsospin & o1 = file_Orbits[iF1];
       int j1 = o1.j;
       int l1 = o1.l;
       int e1 = o1.e;
       if(e1 > emax) continue;
       if(e1 > Emax_file) continue;
//       for (int iF2=0; iF2 <= iF1; iF2++) 
       for (int iF2=0; iF2 < Norbs_file; iF2++) 
       {
         OrbitIsospin & o2 = file_Orbits[iF2];
         int j2 = o2.j;
         int l2 = o2.l;
         int e2 = o2.e;
//         if(e2 > Emax_file) continue;
         if(e1 + e2 > E2max_file) continue;
//         for (int iF3=0; iF3 < Norbs_file; iF3++)
         for (int iF3=0; iF3 < Norbs_file; iF3++)
         {
           OrbitIsospin & o3 = file_Orbits[iF3];
           int j3 = o3.j;
           int l3 = o3.l;
           int e3 = o3.e;
//           if(e3 > Emax_file) continue;
           if(e2 + e3 > E2max_file) continue;
           if(e1 + e3 > E2max_file) continue;
           if(e1 + e2 + e3 > E3max_file) continue;
   
           for (int iF4=0; iF4 <= iF1; iF4++)
           {
             OrbitIsospin & o4 = file_Orbits[iF4];
             int j4 = o4.j;
             int l4 = o4.l;
             int e4 = o4.e;
//             if(e4 > Emax_file) continue;
             if ( j4 != j1 ) continue;
             if ( l4 != l1 ) continue;
//             int i5max = (iF4==iF1) ? iF2 : iF4;
//             for (int iF5=0; iF5 <= i5max; iF5++)
             for (int iF5=0; iF5 < Norbs_file; iF5++)
             {
               OrbitIsospin & o5 = file_Orbits[iF5];
               int j5 = o5.j;
               int l5 = o5.l;
               int e5 = o5.e;
//               if ( e5 > Emax_file) continue;
               if ( j5 != j2 ) continue;
               if ( l5 != l2 ) continue;
               if(e4 + e5 > E2max_file) continue;
//               int i6max = (iF4==iF1 and iF5==iF2) ? iF3 : iF5;
//               for (int iF6=0; iF6 <= i6max; iF6++)
               for (int iF6=0; iF6 <Norbs_file; iF6++)
               {
                 OrbitIsospin & o6 = file_Orbits[iF6];
//                 std::cout << "        iF6 = " << iF6 << std::endl;
                 int j6 = o6.j;
                 int l6 = o6.l;
                 int e6 = o6.e;
                 if(j6 != j3) continue;
                 if(l6 != l3) continue;
//                 if(e6 > Emax_file) continue;
                 if(e4 + e6 > E2max_file) continue;
                 if(e5 + e6 > E2max_file) continue;
                 if(e4 + e5 + e6 > E3max_file) continue;
                 if( (l1+l2+l3)%2 != (l4+l5+l6)%2 ) continue;
   
                 orbit_counter +=1;
//                 std::cout << "        orbit counter = " << orbit_counter << std::endl;

//                 int Jmin = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
//                 int Jmax =std::min( j1+j2, j4+j5 )/2;
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
//                 for (int J = Jmin; J <= Jmax; J++)
//                 {
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
//                             std::cout << "   FILLING THE BUFFER. " << std::endl;
                             size_t n_this_pass = std::min( buffer_size, n_elem_to_read-total_counter );
                             if (filemode == "bin")        infile.read((char*)&vbuf[0], n_this_pass*nwords);
                             else if (filemode == "gz")    for (size_t iread=0;iread<n_this_pass;iread++) zipstream >> vbuf[iread];
                             else if (filemode == "me3j")  for (size_t iread=0;iread<n_this_pass;iread++) infile >> vbuf[iread];
//                             std::cout << "  FIRST 3 elements are " << vbuf[0] << " " << vbuf[1] << " " << vbuf[2] << std::endl;
                             }
                           #pragma omp barrier
                           // back to work, little threads...
   
                         }
                         counter += 1;
                         total_counter += 1;
   
                         if (orbit_counter % num_threads != this_thread) continue; // This is how we divide the jobs among the threads
//                         std::cout << "        now do some more checks... e values " << e1 << " " << e2 << " " << e3 << " " << e4 << " " << e5 << " " << e6 << std::endl;
   
                         if( (e1 > emax) or (e2 > emax) or (e3 > emax) or (e4 > emax) or (e5 > emax) or (e6 > emax)) continue;
                         if( std::max({e1,e2,e3,e4,e5,e6}) > emax) continue;
                         if( ((e1+e2) > E2max) or ((e1+e3) > E2max) or ((e2+e3) > E2max) or ((e4+e5) > E2max) or ((e4+e6) > E2max) or ((e5+e6) > E2max)) continue;
                         if( (e1+e2+e3 > E3max) or ( e4+e5+e6 > E3max) ) continue;
//                         std::cout << "        passed e checks. total_counter = " << total_counter << "   n to read = " << n_elem_to_read << std::endl;
                         if( total_counter > n_elem_to_read) break;
   
                         int i1 = threebodyspace.nlj2idx.at({o1.n,o1.l,o1.j});
                         int i2 = threebodyspace.nlj2idx.at({o2.n,o2.l,o2.j});
                         int i3 = threebodyspace.nlj2idx.at({o3.n,o3.l,o3.j});
                         int i4 = threebodyspace.nlj2idx.at({o4.n,o4.l,o4.j});
                         int i5 = threebodyspace.nlj2idx.at({o5.n,o5.l,o5.j});
                         int i6 = threebodyspace.nlj2idx.at({o6.n,o6.l,o6.j});

   
                         StoreType vset = vbuf[counter-1];
//                         std::cout << "         vset = " << vset << std::endl;
   
//                         if ( (i1==0) or (i2==0) or (i3==0) )//and (i2==2) and (i3==0) )
//                         {
//                             std::cout <<  std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
//                                           std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
//                                           std::setw(4) << T3 << std::setw(16) << counter <<
//                                           std::setw(12) << std::setprecision(6) << vset << std::endl;
//                         }
                         // For monopole matrix elements, we don't have automatic zeros in there
//                         if( (i1==i2 and (J+T12)%2 ==0 ) or ( i4==i5 and (J+T45)%2 ==0 ) ) {
//                           if( abs(vset) > 1.e-6 ){
//                             std::cout << "Warning: something wrong, this three-body matrix element has to be zero" << std::endl;
//                             std::cout <<  std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
//                                           std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
//                                           std::setw(4) << T3 << std::setw(16) << counter <<
//                                           std::setw(12) << std::setprecision(6) << vset << std::endl;
//                           }
//                         }
   
                         if ( std::abs(vset)>1e-8)
                         {
//                          std::cout << "calling SetME_iso_mono with " << i1 << " " << i1 << " " <<i3 << " " << T12 << "  | " << i4 << " " << i5 << " " << i6 << " " << T45 << " ,  " << T3 << "   ->  " << vset << std::endl;
//                          SetME_iso_mono(i1, i2, i3, T12, i4, i5, i6, T45, J, T3, vset );
                          SetME_iso_mono(i1, i2, i3, T12, i4, i5, i6, T45, T3, vset );
//                          double vtest = GetME_iso_mono(i1,i2,i3,T12,i4,i5,i6,T45,T3);
//                          std::cout <<  "   =======  check Get: " << vtest << std::endl;
                         }
   
                       }// for T3
                     }// for T45
                   }// for T12
//                 }// for J
               }// for i6
             }// for i5
           }// for i4
         }// for i3
       }// for i2
     }// for i1
 
   }// end of parallel outer block
 
//   std::cout << "Done reading" << std::endl;
   IMSRGProfiler::timer["ThreeBodyStorage_mono::ReadFile"] += omp_get_wtime() - t_start;
}



// Make our concrete instantiations of the ThreeBodyStorage_mono class. We use either single-precision or half-precision type

template class ThreeBodyStorage_mono<ME_single_type>;
template class ThreeBodyStorage_mono<ME_half_type>;




