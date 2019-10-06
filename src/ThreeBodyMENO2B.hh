#ifndef ThreeBodyMENO2B_h
#define ThreeBodyMENO2B_h 1

#include <fstream>
#include <array>
#include <vector>
#include <map>
#include "ModelSpace.hh"

typedef float ThreeBME_type;
class ThreeBodyMENO2B;

class OrbitIsospin
{
  public:
    int idx, n, l, j, e;
    OrbitIsospin(int idx, int n, int l, int j);
    ~OrbitIsospin();
};

class ThreeBodyChannelNO2B
{
  public:
    ThreeBodyMENO2B* thr;
    int J2;
    int P2;
    int J1;
    int P1;
    int T3;
    int Ndim;
    std::map<std::array<int,4>, int> abct2n;
    std::map<std::array<int,4>, int> iphase;
    ThreeBodyChannelNO2B();
    ThreeBodyChannelNO2B(int, int, int, int, int, ThreeBodyMENO2B*);
    ~ThreeBodyChannelNO2B();
};


class ThreeBodySpaceNO2B
{
  public:
    std::vector<ThreeBodyChannelNO2B> ThreeBodyChannels;
    std::map<std::array<int,5>, int> idcs2ch;
    int NChannels;
    ThreeBodySpaceNO2B();
    ThreeBodySpaceNO2B(ThreeBodyMENO2B*);
    ~ThreeBodySpaceNO2B();
};


class ThreeBodyMENO2B
{
  public:
    ModelSpace * modelspace;
    ThreeBodySpaceNO2B threebodyspace;
    std::map<int, std::vector<ThreeBME_type>> MatEl;
    std::vector<OrbitIsospin> iOrbits;
    std::map<std::array<int,3>, int> nlj2idx;
    std::map<std::array<int,4>, int> nljz2idx;
    std::map<std::array<int,6>, double> cgs;
    int Emax;
    int E2max;
    int E3max;
    int Lmax;
    int Emax_file;
    int E2max_file;
    int E3max_file;
    int Lmax_file;
    bool initialized=false;
    std::string FileName;

    ~ThreeBodyMENO2B();
    ThreeBodyMENO2B();
    ThreeBodyMENO2B(ModelSpace & ms, int emax_file, int e2max_file, int e3max_file, int lmax_file, std::string filename);
    size_t idx1d(size_t bra, size_t ket) { return std::max(bra,ket) * (std::max(bra,ket)-1)/2 + std::min(bra,ket);};
    void SetThBME(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int T3, ThreeBME_type V);
    ThreeBME_type GetThBME(int a, int b, int c, int Tab, int d, int e, int f, int Tde, int J2, int T3);
    ThreeBME_type GetThBME(int a, int b, int c, int d, int e, int f, int J2);
    void ReadFile();
    long long unsigned int CountME();
    template<class T> void ReadStream(T & infile, long long unsigned int n_elms);
};
#endif
