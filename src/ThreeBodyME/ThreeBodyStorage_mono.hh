#ifndef ThreeBodyStorage_mono_h
#define ThreeBodyStorage_mono_h 1

#include "ThreeBodyStorage.hh"
#include "ModelSpace.hh"
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/// Use half-precision floats. For documentation, see http://half.sourceforge.net
 #ifndef NO_x86
// #include <x86intrin.h>
 #endif
 #include "half.hpp"
 typedef half_float::half ME_half_type;
 typedef float ME_single_type;
 typedef double ME_double_type;


class ThreeBodySpaceMono; // forward declaration

/// This should probably be declared in ThreeBodyStorage.hh, not in a specialization.
//class OrbitIsospin
//{
//  public:
//    int idx, n, l, j, e;
////    OrbitIsospin(int idx, int n, int l, int j);
//    OrbitIsospin(int idx, int n, int l, int j) : idx(idx), n(n), l(l), j(j), e(2*n+l) {};
//};

/// A Three body channel for use in the HF approximation.
/// For HF, a matrix element <abc|V|def> enters only if (a,d) (b,e) and (c,f) all differ only in the radial quantum number
class ThreeBodyChannelMono
{
  public:
    int j1;
    int j2;
    int j3;
    int l1;
    int l2;
    int l3;
    int twoT;
    int Ndim;
    std::unordered_map<int, int> abct2n; // maps key hashed from (a,b,c,Tab)  to a single sequential index which is used to store matrix elements
    std::unordered_map<int, int> iphase;


    ThreeBodyChannelMono();
    ThreeBodyChannelMono(int j1, int j2, int j3, int l1, int l2, int l3, int twoT, ThreeBodySpaceMono& thr);
//    ThreeBodyChannelMono(int J2, int P2, int J1, int P1, int T3, ThreeBodySpaceMono& thr);
    int GetIndex(int a, int b, int c, int Tab) const {return Hash_abct(a, b, c, Tab);};
    int GetPhase(int i) const {return iphase.at(i);};
  private:
    int Hash_abct(int, int, int, int) const;
    void UnHash_abct(int, int&, int&, int&, int&) const;
};


class ThreeBodySpaceMono
{
  public:
    int Emax;
    int E2max;
    int E3max;
    int Lmax;
    std::vector<ThreeBodyChannelMono> ThreeBodyChannels;
    std::unordered_map<int, int> idcs2ch;
    int NChannels;
    ThreeBodySpaceMono();
    ThreeBodySpaceMono( int Emax, int E2max, int E3max, int Lmax);
//    int GetChannelIndex(int Jab, int Pab, int Jc, int Pc, int T)const ;
    int GetChannelIndex(int j1, int j2, int j3, int l1, int l2, int l3, int twoT)const ;


    OrbitIsospin& GetIsospinOrbit(size_t i ) { return iOrbits.at(i);}; // move  to ThreeBodySpaceNO2B
    const OrbitIsospin& GetIsospinOrbit(size_t i ) const { return (OrbitIsospin&) iOrbits.at(i);}; // move  to ThreeBodySpaceNO2B
    size_t GetNumberIsospinOrbits() const { return iOrbits.size();}; // move  to ThreeBodySpaceNO2B

    std::vector<OrbitIsospin> iOrbits; // eventually these three can probably get moved to ThreeBodySpaceNO2B
    std::map<std::array<int,3>, int> nlj2idx;
    size_t idx1d(size_t bra, size_t ket) const { return std::max(bra+1,ket+1) * (std::max(bra+1,ket+1)-1)/2 + std::min(bra+1,ket+1)-1;};

  private:
//    int Hash_Channel(int, int, int, int, int) const;
//    void UnHash_Channel(int, int&, int&, int&, int&, int&) const;
    int Hash_Channel(int j1, int j2, int j3, int l1, int l2, int l3, int twoT) const;
    void UnHash_Channel(int key, int& j1, int& j2, int& j3, int& l1, int& l2, int& l3, int& twoT) const;
};



template <class StoreType>
class ThreeBodyStorage_mono : public ThreeBodyStorage
{
  private: 
    std::map<int, std::vector<StoreType>> MatEl;
    ThreeBodySpaceMono threebodyspace;

  public:


    using ThreeBodyStorage::ThreeBodyStorage;  // inherit constructors from the base class
    ThreeBodyStorage_mono( const ThreeBodyStorage_mono& ); // also implement at copy constructor

//    std::shared_ptr<ThreeBodyStorage> Clone() const override;
    std::unique_ptr<ThreeBodyStorage> Clone() const override;

    std::string GetStorageMode() const override {return "mono";};

    void Multiply(const double) override;
    void Add(const ThreeBodyStorage&)  override;
    void Subtract(const ThreeBodyStorage&)  override;

    void Allocate() override;

    void SetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int twoT, ThreeBodyStorage::ME_type V) override;
    ThreeBodyStorage::ME_type GetME_iso_mono(int a, int b, int c, int Tab, int d, int e, int f, int Tde,  int twoT) const override;
    ThreeBodyStorage::ME_type GetME_pn_mono(int a, int b, int c, int d, int e, int f) const override ; 



  double Norm() const override;
  void Erase() override; // Set all elements to zero
  void Deallocate() override;
  size_t size() const override;


/////// FILE READING METHODS THAT MAYBE SHOULD GO IN ReadWrite  ??
    void ReadFile( std::vector<std::string>& StringInputs, std::vector<int>& IntInputs ) override ;

    size_t CountME(int Emax_file, int E2max_file, int E3max_file, int Lmax_file, std::vector<OrbitIsospin>& file_Orbits) const;
//    template<class T> void ReadStream(T & infile, long long unsigned int n_elms);
//    template<class T> void ReadStream(T & infile, size_t n_elms);
//    void ReadBinaryStream( std::vector<ThreeBMENO2B_Store_type>& v, size_t nelms);
//    template <class T> void ReadBinaryStream( std::vector<T>& v, size_t nelms);


};




#endif
