#include "AngMom.hh"
#include "ThreeBodyMENO2B.hh"
#include "ReadWrite.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

OrbitIsospin::~OrbitIsospin()
{}

OrbitIsospin::OrbitIsospin(int idx, int n, int l, int j)
  : idx(idx), n(n), l(l), j(j), e(2*n+l)
{}

ThreeBodyChannelNO2B::~ThreeBodyChannelNO2B()
{}

ThreeBodyChannelNO2B::ThreeBodyChannelNO2B()
  : J2(-1), P2(-1), J1(-1), P1(-1), T3(-1), Ndim(-1)
{}

ThreeBodyChannelNO2B::ThreeBodyChannelNO2B(int J2, int P2, int J1, int P1, int T3, ThreeBodyMENO2B* thr)
  : J2(J2), P2(P2), J1(J1), P1(P1), T3(T3), thr(thr)
{
  Ndim = 0;
  int emax = std::max(thr->Emax, thr->Emax_file);
  int e2max = std::max(thr->E2max, thr->E2max_file);
  int e3max = std::max(thr->E3max, thr->E3max_file);
  int Norbs = thr->iOrbits.size();
  for (int ia=0; ia < Norbs; ia++){
    OrbitIsospin & oa = thr->iOrbits[ia];
    if(oa.e > emax) continue;
    for (int ib=ia; ib < Norbs; ib++){
      OrbitIsospin & ob = thr->iOrbits[ib];
      if(ob.e > emax) continue;
      if(oa.e + ob.e > e2max) continue;
      if(std::abs(oa.j-ob.j) > 2*J2) continue;
      if(        (oa.j+ob.j) < 2*J2) continue;
      if((oa.l+ob.l)%2 != P2) continue;
      for (int ic=0; ic < Norbs; ic++){
        OrbitIsospin & oc = thr->iOrbits[ic];
        if(oc.j != J1) continue;
        if(oc.l%2 != P1) continue;
        if(oc.e > emax) continue;
        if(oa.e + oc.e > e2max) continue;
        if(ob.e + oc.e > e2max) continue;
        if(oa.e + ob.e + oc.e > e3max) continue;
        for (int T2: {0,1}){
          if( std::abs(2*T2-1) > T3 ) continue;
          if(         (2*T2+1) < T3 ) continue;
          if((J2+T2)%2 == 0 and ia == ib) continue;
          int ph = 1 - 2*(((oa.j + ob.j)/2 + J2 + T2)%2);
          int key_ab = Hash_abct(ia, ib, ic, T2);
          int key_ba = Hash_abct(ib, ia, ic, T2);
          abct2n[key_ab] = Ndim;
          abct2n[key_ba] = Ndim;
          iphase[key_ab] = 1;
          iphase[key_ba] = ph;
          Ndim += 1;
        }
      }
    }
  }

}

int ThreeBodyChannelNO2B::Hash_abct(int a, int b, int c, int Tab)
{
  return ( a + (b<<10) + (c<<20) + (Tab<<30) );
}

void ThreeBodyChannelNO2B::UnHash_abct(int key, int & a, int & b, int & c, int & Tab)
{
  int Lowest_ten_bits = 0x3FF;
  a = ( (key>> 0) & Lowest_ten_bits);
  b = ( (key>>10) & Lowest_ten_bits);
  c = ( (key>>20) & Lowest_ten_bits);
  Tab = ( (key>>30) & Lowest_ten_bits);
}

ThreeBodySpaceNO2B::~ThreeBodySpaceNO2B()
{}

ThreeBodySpaceNO2B::ThreeBodySpaceNO2B()
  : NChannels(-1), ThreeBodyChannels()
{}

ThreeBodySpaceNO2B::ThreeBodySpaceNO2B(ThreeBodyMENO2B* thr)
{
  int lmax = std::max(thr->Lmax, thr->Lmax_file);
  int emax = std::max(thr->Emax, thr->Emax_file);
  int e2max = std::max(thr->E2max, thr->E2max_file);
  int e3max = std::max(thr->E3max, thr->E3max_file);
  int J2max = std::min(2*lmax+1, e2max+1);
  int cnt = 0;
  for (int J2=0; J2<=J2max; ++J2){
    for (int P2: {0,1}){

      for (int J1=1; J1<=2*lmax+1; J1+=2){
        for (int P1: {0,1}){
          for (int T3 : {1,3}){
            ThreeBodyChannelNO2B channel(J2, P2, J1, P1, T3, thr);
            if(channel.Ndim < 1) continue;
            std::cout <<
              "   J2 = " << std::setw(4) << J2 <<
              ",  P2 = " << std::setw(4) << P2 <<
              ",  J1 = " << std::setw(4) << J1 <<
              ",  P1 = " << std::setw(4) << P1 <<
              ",  T3 = " << std::setw(4) << T3 <<
              ",  Ndim = " << std::setw(10) << channel.Ndim <<
              std::endl;
            ThreeBodyChannels.push_back(channel);
            int key = Hash_Channel(J2,P2,J1,P1,T3);
            int j2, p2, j1, p1, t3;
            UnHash_Channel(key,j2,p2,j1,p1,t3);
            idcs2ch[GetChannelIndex(J2,P2,J1,P1,T3)] = cnt;
            cnt += 1;
          }
        }
      }
    }
  }
  NChannels = ThreeBodyChannels.size();
}

int ThreeBodySpaceNO2B::Hash_Channel(int J2, int P2, int J1, int P1, int T3)
{
  return ( J2 + (J1<<8) + (P2<<16) + (P1<<17) + (T3<<18) );
}

void ThreeBodySpaceNO2B::UnHash_Channel(int key, int& J2, int& P2, int& J1, int& P1, int& T3)
{
  J2 = ( (key>> 0) & 255 );
  J1 = ( (key>> 8) & 255 );
  P2 = ( (key>>16) & 1 );
  P1 = ( (key>>17) & 1 );
  T3 = ( (key>>18) & 255 );
}


ThreeBodyMENO2B::~ThreeBodyMENO2B()
{}

ThreeBodyMENO2B::ThreeBodyMENO2B()
  : Emax(0), E2max(0), E3max(0), Lmax(0),
  Emax_file(0), E2max_file(0), E3max_file(0), Lmax_file(0)
{}

ThreeBodyMENO2B::ThreeBodyMENO2B(const ThreeBodyMENO2B& tbme)
  : modelspace(tbme.modelspace), threebodyspace(tbme.threebodyspace),
  MatEl( tbme.MatEl ), iOrbits( tbme.iOrbits ),  nlj2idx(tbme.nlj2idx),
  Emax(tbme.Emax), E2max(tbme.E2max),
  E3max(tbme.E3max), Lmax(tbme.Lmax),
  Emax_file(tbme.Emax_file), E2max_file(tbme.E2max_file),
  E3max_file(tbme.E3max_file), Lmax_file(tbme.Lmax_file),
  initialized(tbme.initialized)
{}

ThreeBodyMENO2B& ThreeBodyMENO2B::operator*=(const double rhs)
{
  for ( auto& itmat : MatEl )
  {
    for ( auto& it : itmat.second )
    {
      it *= rhs;
    }
  }
  return *this;
}

ThreeBodyMENO2B& ThreeBodyMENO2B::operator+=(const ThreeBodyMENO2B& rhs)
{
  for ( auto& itmat : rhs.MatEl )
  {
    auto ch = itmat.first;
    for ( size_t i=0; i<itmat.second.size(); i++)
    {
      MatEl[ch][i] += itmat.second[i];
    }
  }
  return *this;
}

ThreeBodyMENO2B& ThreeBodyMENO2B::operator-=(const ThreeBodyMENO2B& rhs)
{
  for ( auto& itmat : rhs.MatEl )
  {
    auto ch = itmat.first;
    for ( size_t i=0; i<itmat.second.size(); i++)
    {
      MatEl[ch][i] -= itmat.second[i];
    }
  }
  return *this;
}

void ThreeBodyMENO2B::Allocate(ModelSpace & ms,
    int emax_file, int e2max_file, int e3max_file, int lmax_file,
    std::string filename)
{
  modelspace = &ms;
  Emax = modelspace->GetEmax();
  E2max = modelspace->GetE2max();
  E3max = modelspace->GetE3max();
  Lmax = std::min(modelspace->GetLmax(), Emax);
  Emax_file = emax_file;
  E2max_file = e2max_file;
  E3max_file = e3max_file;
  Lmax_file = lmax_file;
  FileName = filename;
  initialized = true;
  int idx = 0;
  for (int e=0; e<=std::max(modelspace->GetEmax(),Emax_file); ++e) {
    int lmin = e%2;
    for (int l=lmin; l<=std::min(e,std::max(modelspace->GetLmax(),lmax_file)); l+=2) {
      int n = (e-l)/2;
      int twojMin = std::abs(2*l-1);
      int twojMax = 2*l+1;
      for (int twoj=twojMin; twoj<=twojMax; twoj+=2) {
        OrbitIsospin orb(idx,n,l,twoj);
        iOrbits.push_back(orb);
        nlj2idx[{n,l,twoj}]=idx;
        idx += 1;
      }
    }
  }
  threebodyspace = ThreeBodySpaceNO2B(this);
  for (int ch=0; ch<threebodyspace.NChannels; ch++){
    ThreeBodyChannelNO2B ch_no2b=threebodyspace.ThreeBodyChannels[ch];
    size_t n = ch_no2b.Ndim;
    std::vector<ThreeBME_type> vch(n*(n+1)/2, 0.0);
    MatEl[ch] = vch;
  }
}

void ThreeBodyMENO2B::SetThBME(int a, int b, int c, int Tab,
    int d, int e, int f, int Tde, int J2, int T3, ThreeBME_type V)
{
  OrbitIsospin & oa = iOrbits[a];
  OrbitIsospin & ob = iOrbits[b];
  OrbitIsospin & oc = iOrbits[c];
  OrbitIsospin & od = iOrbits[d];
  OrbitIsospin & oe = iOrbits[e];
  OrbitIsospin & of = iOrbits[f];

  int P1 = oc.l%2;
  if(P1 != of.l%2) return;

  int J1 = oc.j;
  if(J1 != of.j) return;

  int P2 = (oa.l + ob.l)%2;
  if(P2 != (od.l + oe.l)%2) return;

  int ch = threebodyspace.idcs2ch[threebodyspace.GetChannelIndex(J2,P2,J1,P1,T3)];
  ThreeBodyChannelNO2B & ch_no2b = threebodyspace.ThreeBodyChannels[ch];
  int ibra = ch_no2b.GetIndex(a,b,c,Tab);
  int iket = ch_no2b.GetIndex(d,e,f,Tde);
  if(ch_no2b.iphase.find(ibra) == ch_no2b.iphase.end()) return;
  if(ch_no2b.iphase.find(iket) == ch_no2b.iphase.end()) return;
  int ph = ch_no2b.iphase[ibra] * ch_no2b.iphase[iket];
  int bra = ch_no2b.abct2n[ibra];
  int ket = ch_no2b.abct2n[iket];
  auto& Vch = MatEl[ch];
  Vch[idx1d(bra,ket)] = V * ph;
}

ThreeBME_type ThreeBodyMENO2B::GetThBME(int a, int b, int c, int Tab,
    int d, int e, int f, int Tde, int J2, int T3)
{
  OrbitIsospin & oa = iOrbits[a];
  OrbitIsospin & ob = iOrbits[b];
  OrbitIsospin & oc = iOrbits[c];
  OrbitIsospin & od = iOrbits[d];
  OrbitIsospin & oe = iOrbits[e];
  OrbitIsospin & of = iOrbits[f];

  int P1 = oc.l%2;
  if(P1 != of.l%2) return 0.0;

  int J1 = oc.j;
  if(J1 != of.j) return 0.0;

  int P2 = (oa.l + ob.l)%2;
  if(P2 != (od.l + oe.l)%2) return 0.0;

  int ch = threebodyspace.idcs2ch[threebodyspace.GetChannelIndex(J2,P2,J1,P1,T3)];
  ThreeBodyChannelNO2B & ch_no2b = threebodyspace.ThreeBodyChannels[ch];
  int ibra = ch_no2b.GetIndex(a,b,c,Tab);
  int iket = ch_no2b.GetIndex(d,e,f,Tde);
  if(ch_no2b.iphase.find(ibra) == ch_no2b.iphase.end()) return 0.0;
  if(ch_no2b.iphase.find(iket) == ch_no2b.iphase.end()) return 0.0;
  int ph = ch_no2b.iphase[ibra] * ch_no2b.iphase[iket];
  int bra = ch_no2b.abct2n[ibra];
  int ket = ch_no2b.abct2n[iket];
  auto& Vch = MatEl[ch];
  return Vch[idx1d(bra,ket)] * ph;
}

ThreeBME_type ThreeBodyMENO2B::GetThBME(int a, int b, int c, int d, int e, int f,  int J2)
{
  Orbit & oa = modelspace->GetOrbit(a);
  Orbit & ob = modelspace->GetOrbit(b);
  Orbit & oc = modelspace->GetOrbit(c);
  Orbit & od = modelspace->GetOrbit(d);
  Orbit & oe = modelspace->GetOrbit(e);
  Orbit & of = modelspace->GetOrbit(f);

  int P1 = oc.l%2;
  if(P1 != of.l%2) return 0.0;

  int J1 = oc.j2;
  if(J1 != of.j2) return 0.0;

  int P2 = (oa.l + ob.l)%2;
  if(P2 != (od.l + oe.l)%2) return 0.0;

  int Z3 = oa.tz2 + ob.tz2 + oc.tz2;
  if(Z3 != od.tz2 + oe.tz2 + of.tz2) return 0.0;

  int i1 = nlj2idx[{oa.n, oa.l, oa.j2}];
  int i2 = nlj2idx[{ob.n, ob.l, ob.j2}];
  int i3 = nlj2idx[{oc.n, oc.l, oc.j2}];
  int i4 = nlj2idx[{od.n, od.l, od.j2}];
  int i5 = nlj2idx[{oe.n, oe.l, oe.j2}];
  int i6 = nlj2idx[{of.n, of.l, of.j2}];

  double v=0.0;
  for (int T12: {0,1}){
    if(std::abs(oa.tz2 + ob.tz2) > 2*T12) continue;
    for (int T45: {0,1}){
      if(std::abs(od.tz2 + oe.tz2) > 2*T45) continue;
      for (int T=std::max( std::abs(2*T12-1), std::abs(2*T45-1) );
          T<= std::min( 2*T12+1, 2*T45+1 ); T+=2){
        v += GetThBME(i1, i2, i3, T12, i4, i5, i6, T45, J2, T) *
          AngMom::CG(0.5, oa.tz2*0.5, 0.5, ob.tz2*0.5, T12, (oa.tz2+ob.tz2)*0.5) *
          AngMom::CG(0.5, od.tz2*0.5, 0.5, oe.tz2*0.5, T45, (od.tz2+oe.tz2)*0.5) *
          AngMom::CG(T12, (oa.tz2+ob.tz2)*0.5, 0.5, oc.tz2*0.5, T*0.5, Z3*0.5) *
          AngMom::CG(T45, (od.tz2+oe.tz2)*0.5, 0.5, of.tz2*0.5, T*0.5, Z3*0.5);
      }
    }
  }
  return v;
}

void ThreeBodyMENO2B::ReadFile()
{
  long long unsigned int n_elms = CountME();
  if(FileName.find(".gz") != std::string::npos){
    std::ifstream infile(FileName, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream zipstream;
    zipstream.push(boost::iostreams::gzip_decompressor());
    zipstream.push(infile);
    char line[512];
    zipstream.getline(line, 512);
    ReadStream(zipstream, n_elms);
    return;
  }
  if(FileName.find("stream.bin") != std::string::npos){
    std::ifstream infile(FileName, std::ios::binary);
    infile.seekg(0, infile.end);
    size_t n_elem = infile.tellg();
    n_elem /= sizeof(float);
    std::vector<float> v(n_elem);
    infile.read((char*)&v[0], n_elem*sizeof(float));
    VectorStream vecstream(v);
    ReadStream(vecstream, n_elem);
    return;
  }
  if(FileName.find(".me3j") != std::string::npos){
    std::ifstream infile(FileName);
    char line[512];
    infile.getline(line, 512);
    ReadStream(infile, n_elms);
    return;
  }
}

long long unsigned int ThreeBodyMENO2B::CountME()
{
  long long unsigned int counter=0;
  int Norbs = iOrbits.size();
  for (int i1=0; i1 < Norbs; i1++) {
    OrbitIsospin & o1 = iOrbits[i1];
    int j1 = o1.j;
    int l1 = o1.l;
    int e1 = o1.e;
    if(e1 > Emax_file) continue;
    for (int i2=0; i2 <= i1; i2++) {
      OrbitIsospin & o2 = iOrbits[i2];
      int j2 = o2.j;
      int l2 = o2.l;
      int e2 = o2.e;
      if(e2 > Emax_file) continue;
      if(e1 + e2 > E2max_file) continue;
      for (int i3=0; i3 < Norbs; i3++) {
        OrbitIsospin & o3 = iOrbits[i3];
        int j3 = o3.j;
        int l3 = o3.l;
        int e3 = o3.e;
        if(e3 > Emax_file) continue;
        if(e2 + e3 > E2max_file) continue;
        if(e1 + e3 > E2max_file) continue;
        if(e1 + e2 + e3 > E3max_file) continue;

        for (int i4=0; i4 < Norbs; i4++) {
          OrbitIsospin & o4 = iOrbits[i4];
          int j4 = o4.j;
          int l4 = o4.l;
          int e4 = o4.e;
          if(e4 > Emax_file) continue;
          for (int i5=0; i5 <= i4; i5++) {
            OrbitIsospin & o5 = iOrbits[i5];
            int j5 = o5.j;
            int l5 = o5.l;
            int e5 = o5.e;
            if(e5 > Emax_file) continue;
            if(e4 + e5 > E2max_file) continue;
            for (int i6=0; i6 < Norbs; i6++) {
              OrbitIsospin & o6 = iOrbits[i6];
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

              for (int J = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
                  J <= std::min( j1+j2, j4+j5 )/2; J++) {
                for (int T12: {0,1}){
                  for (int T45: {0,1}){
                    for (int T = std::max( std::abs(2*T12-1), std::abs(2*T45-1) );
                        T <= std::min( 2*T12+1, 2*T45+1 ); T+=2) {
                      counter += 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return counter;
}

  template<class T>
void ThreeBodyMENO2B::ReadStream(T & infile, long long unsigned int n_elms)
{
  int buffer_size = 10000000;
  int counter = 0;
  long long unsigned int total_counter = 0;
  int Emax = modelspace->GetEmax();
  int E2max = modelspace->GetE2max();
  int E3max = modelspace->GetE3max();
  int Norbs = iOrbits.size();
  std::vector<ThreeBME_type> v(buffer_size,0.0);

  for (int i1=0; i1 < Norbs; i1++) {
    OrbitIsospin & o1 = iOrbits[i1];
    int j1 = o1.j;
    int l1 = o1.l;
    int e1 = o1.e;
    if(e1 > Emax) continue;
    if(e1 > Emax_file) continue;
    for (int i2=0; i2 <= i1; i2++) {
      OrbitIsospin & o2 = iOrbits[i2];
      int j2 = o2.j;
      int l2 = o2.l;
      int e2 = o2.e;
      if(e2 > Emax_file) continue;
      if(e1 + e2 > E2max_file) continue;
      for (int i3=0; i3 < Norbs; i3++) {
        OrbitIsospin & o3 = iOrbits[i3];
        int j3 = o3.j;
        int l3 = o3.l;
        int e3 = o3.e;
        if(e3 > Emax_file) continue;
        if(e2 + e3 > E2max_file) continue;
        if(e1 + e3 > E2max_file) continue;
        if(e1 + e2 + e3 > E3max_file) continue;

        for (int i4=0; i4 < Norbs; i4++) {
          OrbitIsospin & o4 = iOrbits[i4];
          int j4 = o4.j;
          int l4 = o4.l;
          int e4 = o4.e;
          if(e4 > Emax_file) continue;
          for (int i5=0; i5 <= i4; i5++) {
            OrbitIsospin & o5 = iOrbits[i5];
            int j5 = o5.j;
            int l5 = o5.l;
            int e5 = o5.e;
            if(e5 > Emax_file) continue;
            if(e4 + e5 > E2max_file) continue;
            for (int i6=0; i6 < Norbs; i6++) {
              OrbitIsospin & o6 = iOrbits[i6];
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
              for (int J = std::max( std::abs(j1-j2), std::abs(j4-j5) )/2;
                  J <= std::min( j1+j2, j4+j5 )/2; J++) {
                for (int T12: {0,1}){
                  for (int T45: {0,1}){
                    for (int T3 = std::max( std::abs(2*T12-1), std::abs(2*T45-1) );
                        T3 <= std::min( 2*T12+1, 2*T45+1 ); T3+=2) {
                      if(counter == buffer_size) counter = 0;
                      if(counter == 0){
                        if( n_elms - total_counter >= buffer_size){
                          for (int iread=0; iread<buffer_size; iread++) infile >> v[iread];
                        }
                        else{
                          for (int iread=0; iread<n_elms-total_counter; iread++) infile >> v[iread];
                        }
                      }
                      counter += 1;
                      total_counter += 1;

                      if( e1 > Emax) continue;
                      if( e2 > Emax) continue;
                      if( e3 > Emax) continue;
                      if( e4 > Emax) continue;
                      if( e5 > Emax) continue;
                      if( e6 > Emax) continue;

                      if( e1+e2 > E2max ) continue;
                      if( e1+e3 > E2max ) continue;
                      if( e2+e3 > E2max ) continue;
                      if( e4+e5 > E2max ) continue;
                      if( e4+e6 > E2max ) continue;
                      if( e5+e6 > E2max ) continue;

                      if( e1+e2+e3 > E3max ) continue;
                      if( e4+e5+e6 > E3max ) continue;

                      if( i1==i2 and (J+T12)%2 ==0 ) {
                        if( abs(v[counter-1]) > 1.e-6 ){
                          std::cout << "Warning: something wrong, this three-body matrix element has to be zero" << std::endl;
                          std::cout <<
                            std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
                            std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
                            std::setw(4) << J << std::setw(4) << T3 << std::setw(16) << counter <<
                            std::setw(12) << std::setprecision(6) << v[counter-1] << std::endl;
                        }
                      }

                      if( i4==i5 and (J+T45)%2 ==0 ) {
                        if( abs(v[counter-1]) > 1.e-6 ){
                          std::cout << "Warning: something wrong, this three-body matrix element has to be zero" << std::endl;
                          std::cout <<
                            std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
                            std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
                            std::setw(4) << J << std::setw(4) << T3 << std::setw(16) << counter <<
                            std::setw(12) << std::setprecision(6) << v[counter-1] << std::endl;
                        }
                      }
                      //std::cout <<
                      //  std::setw(4) << i1 << std::setw(4) << i2 << std::setw(4) << i3 << std::setw(4) << T12 <<
                      //  std::setw(4) << i4 << std::setw(4) << i5 << std::setw(4) << i6 << std::setw(4) << T45 <<
                      //  std::setw(4) << J << std::setw(4) << T3 << std::setw(16) << counter <<
                      //  std::setw(12) << std::setprecision(6) << v[counter-1] << std::endl;
                      SetThBME(i1, i2, i3, T12, i4, i5, i6, T45, J, T3, v[counter-1]);

                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

}
