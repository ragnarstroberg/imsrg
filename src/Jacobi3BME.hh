#ifndef Jacobi3BME_h
#define Jacobi3BME_h 1

#include "ModelSpace.hh"
#include <vector>
#include <array>

struct jacobi1_state{ int n; int l; int s; int j; int t; };
struct jacobi2_state{ int n; int l; int j2; };

struct labstate_nav{ int na; int nb; int nc; int la; int lb; int lc; int ja; int jb; int jc; int bj12; int jtot; };
struct relativestate_nav{ int n12; int n3; int ncm; int l12; int l3; int lcm; int j12; int j3; int I3; int s12; int jtot; };
struct bookkeeping_nav{ int matstart; int nrelstart; int nrelnum; int nspstart; int nspnum; };

class HartreeFock; // forward declaration

class Jacobi3BME
{
 public:

  int Nmax;  // maximum number of oscillator quanta in a jacobi state
  int twoJmin;
  int twoJmax;
  int twoTmin;
  int twoTmax;

  int emax;  // cuts on lab frame basis, 2n+l <= emax
  int E2max;
  int E3max; // e1+e2+e3 <= E3max

  std::vector<jacobi1_state> jacobi_1; // quantum numbers of the relative wf of the first 2 particles
  std::vector<jacobi2_state> jacobi_2; // quantum number of the second jacobi coordinate, 3rd -> cm of other 2

  std::vector<std::vector<std::array<size_t,2>>> NAS_jacobi_states;  // for each T,J,p,N block, this maps the NAS index to a pair of jacobi states (actually, to their indices for lookup in jacobi_1 and jacobi_2)



  std::vector<double> meAS;         // antisymmetrized matrix elements
  std::vector<double> meNAS;        // non-antisymmetrized matrix elements
  std::vector<size_t> start_locAS;  // starting element for a given T,J,p channel
  std::vector<size_t> start_locNAS; // starting element for a given T,J,p channel
  std::vector<size_t> dimensionAS;  // antisymmetrized dimension for each T,J,p,N
  std::vector<size_t> dimensionNAS; // non-antisymmetrized dimension for each T,J,p,N

  std::vector<double> cfpvec;      // the cfp's (coefficients of fractional parentage), i.e. the overlaps of the AS and NAS basis states
  std::vector<size_t> cfp_start_loc;  // starting element for a given T,J,p


  Jacobi3BME();
  Jacobi3BME( int nmax, int twojmin, int twojmax, int twotmin, int twotmax );
 
  void Allocate();

  size_t HashTJN( int twoT, int twoJ, int N);
  size_t HashTJNN( int twoT, int twoJ, int Nbra, int Nket);
  void   SetDimensionAS( int twoT, int twoJ, int parity, int N, size_t dim) {dimensionAS.at(HashTJN(twoT,twoJ,N))=dim;};
  size_t GetDimensionAS( int twoT, int twoJ, int parity, int N ) {return dimensionAS.at(HashTJN(twoT,twoJ,N));};
  void   SetDimensionNAS( int twoT, int twoJ, int parity, int N, size_t dim) {dimensionNAS.at(HashTJN(twoT,twoJ,N))=dim;};
  size_t GetDimensionNAS( int twoT, int twoJ, int parity, int N ) {return dimensionNAS.at(HashTJN(twoT,twoJ,N));}; 

  void SetEmax(int e)   {emax=e;};
  void SetE2max(int e2) {E2max=e2;};
  void SetE3max(int e3) {E3max=e3;};

  void GetJacobiStates( int twoT, int twoJ, int parity, int E12, int iNAS, jacobi1_state& jac1, jacobi2_state& jac2);
  void CalculateJacobiLabels();

  size_t GetCFPStartLocation(int t, int j, int N) { return cfp_start_loc[HashTJN(t,j,N)]; };
  double& AccessCFP( int t, int j, int p, int N, int iAS, int iNAS ) { return cfpvec[ GetCFPStartLocation(t,j,N) + iAS*GetDimensionNAS(t,j,p,N) + iNAS] ;  }; 


//  double& ElementAS(size_t ibra,size_t iket, int Nbra, int Nket, int twoJ,int twoT, int p);
//  double& ElementNAS(size_t ibra,size_t iket, int Nbra, int Nket, int twoJ,int twoT, int p); 
  void SetMatElAS( size_t ibra,size_t iket, int Nbra, int Nket, int twoT,int twoJ, int p, double me);
  void SetMatElNAS( size_t ibra,size_t iket, int Nbra, int Nket, int twoT,int twoJ, int p, double me);
  double GetMatElAS( size_t ibra,size_t iket, int Nbra, int Nket, int twoT,int twoJ, int p);
  double GetMatElNAS( size_t ibra,size_t iket, int Nbra, int Nket, int twoT,int twoJ, int p);
  size_t GetStartLocAS(int twoT, int twoJ, int N1, int N2) { return start_locAS.at( HashTJNN(twoT,twoJ,N1,N2) ); };
  size_t GetStartLocNAS(int twoT, int twoJ, int N1, int N2) { return start_locNAS.at( HashTJNN(twoT,twoJ,N1,N2) ); };
  void ComputeNAS_MatrixElements( );
  double GetLabMatEl( Ket3& bra, Ket3& ket, int Jab, int Jcd, int twoJ, int Tab, int Tcd, int twoT);

//  double GetV3mon( size_t a, size_t b, size_t c, size_t d, size_t e, size_t f ); //< Get a single 3-body monopole term, for use in a Hartree-Fock calculation

  std::string TcoeffHash(uint64_t na, uint64_t nb, uint64_t nc, uint64_t Jab, uint64_t twoJ, uint64_t jac1, uint64_t jac2, uint64_t twoJ12, uint64_t Lcm );
//  void TcoeffUnHash(std::string& key, uint64_t& na, uint64_t& nb, uint64_t& nc, uint64_t& Jab, uint64_t& twoJ, uint64_t& jac1, uint64_t& jac2, uint64_t& twoJ12, uint64_t& Lcm );
  void TcoeffUnHash(std::string& key, int& na, int& nb, int& nc, int& Jab, int& twoJ, int& jac1, int& jac2, int& twoJ12, int& Lcm );
//  void GetV3mon_all( std::vector<uint64_t>& keys, std::vector<double>& v3mon, ModelSpace& modelspace ); //< Get all the monopoles in one go, which should be more efficient
  void GetV3mon_all( HartreeFock& hf ); //< Get all the monopoles in one go, which should be more efficient

  double Tcoeff_wrapper( Ket3& ket, int Jab, int twoJ, jacobi1_state& jac1, jacobi2_state& jac2, int twoJ12, int Ncm, int Lcm);

  void ReadTcoeffNavratil( std::string fname, std::vector<double>& tcoeff, std::vector<labstate_nav>& labst, std::vector<relativestate_nav>& relst, std::vector<bookkeeping_nav>& bookkeeping ); //< for unit testing

  void TestReadTcoeffNavratil( std::string fname ); //< for unit testing


};






#endif
