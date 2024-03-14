#ifndef ThreeBodyStorage_no2b_h
#define ThreeBodyStorage_no2b_h 1

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
 #include "half.hpp"
 typedef half_float::half ME_half_type;
 typedef float ME_single_type;
 typedef double ME_double_type;


class ThreeBodySpaceNO2B; // forward declaration

//class OrbitIsospin
//{
//  public:
//    int idx, n, l, j, e;
//    OrbitIsospin(int idx, int n, int l, int j) : idx(idx), n(n), l(l), j(j), e(2*n+l) {};
//};

/// A Three body channel for use in the NO2B approximation.
/// For a state |abc> we have ja,jb coupled to J2 with parity P2
/// J1 = jc, and P1 = parity of orbit c.
/// T3 is twice the total 3-body isospin.
/// Orbit c is the one that will get summed over the reference
/// to get the normal ordered 2-body contribution.
///     --- SRS interpretation of Takayuki's code, so grain of salt.
class ThreeBodyChannelNO2B
{
  public:
    int J2;
    int P2;
    int J1;
    int P1;
    int T3;
    int Ndim;
    std::unordered_map<int, int> abct2n; // maps key hashed from (a,b,c,Tab)  to a single sequential index which is used to store matrix elements
    std::unordered_map<int, int> iphase;


    ThreeBodyChannelNO2B();
    ThreeBodyChannelNO2B(int J2, int P2, int J1, int P1, int T3, ThreeBodySpaceNO2B& thr);
    int GetIndex(int a, int b, int c, int Tab) const {return Hash_abct(a, b, c, Tab);};
    int GetPhase(int i) const {return iphase.at(i);};
  private:
    int Hash_abct(int, int, int, int) const;
    void UnHash_abct(int, int&, int&, int&, int&) const;
};


class ThreeBodySpaceNO2B
{
  public:
    int Emax;
    int E2max;
    int E3max;
    int Lmax;
    std::vector<ThreeBodyChannelNO2B> ThreeBodyChannels;
    std::unordered_map<int, int> idcs2ch;
    int NChannels;
    ThreeBodySpaceNO2B();
    ThreeBodySpaceNO2B( int Emax, int E2max, int E3max, int Lmax);
    int GetChannelIndex(int Jab, int Pab, int Jc, int Pc, int T)const ;


    OrbitIsospin& GetIsospinOrbit(size_t i ) { return iOrbits.at(i);}; // move  to ThreeBodySpaceNO2B
    const OrbitIsospin& GetIsospinOrbit(size_t i ) const { return (OrbitIsospin&) iOrbits.at(i);}; // move  to ThreeBodySpaceNO2B
    size_t GetNumberIsospinOrbits() const { return iOrbits.size();}; // move  to ThreeBodySpaceNO2B

    std::vector<OrbitIsospin> iOrbits; // eventually these three can probably get moved to ThreeBodySpaceNO2B
    std::map<std::array<int,3>, int> nlj2idx;
    size_t idx1d(size_t bra, size_t ket) const { return std::max(bra+1,ket+1) * (std::max(bra+1,ket+1)-1)/2 + std::min(bra+1,ket+1)-1;};

  private:
    int Hash_Channel(int, int, int, int, int) const;
    void UnHash_Channel(int, int&, int&, int&, int&, int&) const;
};



template <class StoreType>
class ThreeBodyStorage_no2b : public ThreeBodyStorage
{
  private: 
    std::map<int, std::vector<StoreType>> MatEl;
//  protected:
    ThreeBodySpaceNO2B threebodyspace;

  public:

//    int Emax_file;  // Don't think we need any of these _file variables
//    int E2max_file;
//    int E3max_file;
//    int Lmax_file;
//    bool initialized=false;
//    std::string FileName;

//    ~ThreeBodyMENO2B();
//    ThreeBodyMENO2B();
    //ThreeBodyMENO2B(ModelSpace & ms, int emax_file, int e2max_file, int e3max_file, int lmax_file, std::string filename);

    using ThreeBodyStorage::ThreeBodyStorage;  // inherit constructors from the base class
    ThreeBodyStorage_no2b( const ThreeBodyStorage_no2b& ); // also implement at copy constructor

//    std::shared_ptr<ThreeBodyStorage> Clone() const override;
    std::unique_ptr<ThreeBodyStorage> Clone() const override;

    std::string GetStorageMode() const override {return "no2b";};

    void Multiply(const double) override;
    void Add(const ThreeBodyStorage&)  override;
    void Subtract(const ThreeBodyStorage&)  override;

    void Allocate() override;

    void SetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int twoT, ThreeBodyStorage::ME_type V) override;
    ThreeBodyStorage::ME_type GetME_iso_no2b(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int twoT) const override;
    ThreeBodyStorage::ME_type GetME_pn_no2b(int a, int b, int c, int d, int e, int f, int J2) const override ; 



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
