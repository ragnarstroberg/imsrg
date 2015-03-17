
#ifndef Operator_h
#define Operator_h 1

#include "ModelSpace.hh"
#include <armadillo>
#include <string>
#include <vector>
#include <map>

//#define JMAX 30

using namespace std;

typedef uint64_t orbindx3_t;

///
/// The Operator class provides a generic operator up to three-body, scalar or tensor.
/// The class contains lots of methods and overloaded operators so that the resulting
/// code that uses the operators can look as close as possible to the math that is
/// written down.

class Operator
{
 public:
  //Fields
/// The zero body piece of the operator.
  double ZeroBody; 
/// The one body piece of the operator, stored in a single NxN armadillo matrix, where N is
/// the number of single-particle orbits.
  arma::mat OneBody;
/// The two-body piece of the operator, stored in a vector of maps of of armadillo matrices.
/// The index of the vector indicates the J-coupled two-body channel of the ket state, while the
/// map key is the two-body channel of the bra state. This is done to allow for tensor operators
/// which connect different two-body channels without having to store all possible combinations.
/// In the case of a scalar operator, there is only one map key for the bra state, corresponding
/// to that of the ket state.
/// The normalized J-coupled TBME's are stored in the matrices. However, when the TBME's are
/// accessed by GetTBME(), they are returned as
/// \f$ \tilde{\Gamma}_{ijkl} \equiv \sqrt{(1+\delta_{ij})(1+\delta_{kl})} \Gamma_{ijkl} \f$
/// because the flow equations are derived in terms of \f$ \tilde{\Gamma} \f$.
/// For efficiency, only matrix elements with \f$ i\leq j \f$ and \f$ k\leq l \f$
/// are stored.
/// When performing sums that can be cast as matrix multiplication, we have something
/// of the form
/// \f[
/// \tilde{Z}_{ijkl} \sim \frac{1}{2} \sum_{ab}\tilde{X}_{ijab} \tilde{Y}_{abkl}
/// \f]
/// which may be rewritten as a restricted sum, or matrix multiplication
/// \f[
/// Z_{ijkl} \sim \sum_{a\leq b} X_{ijab} Y_{abkl} = \left( X\cdot Y \right)_{ijkl}
/// \f]
  vector<map<int,arma::mat> > TwoBody;  
/// The three-body piece of the operator, stored in a map of vectors of vectors of doubles.
/// The 3BMEs are stored in unnormalized JT coupled form
/// \f$ \langle (ab)J_{ab}t_{ab};cJT | V | (de)J_{de}t_{de};f JT \rangle \f$.
/// To minimize the number of stored matrix elements, only elements with
/// \f$ a\geq b \geq c, a\geq d\geq e \geq f \f$ are stored.
/// The other combinations are obtained on the fly by GetThreeBodyME().
/// The storage format is ThreeBody[orbit_index][Jab_index][JT_index].
//  map< orbindx3_t, map< int, map<int, map<int, vector<double> > > > >ThreeBody;
//  map< orbindx3_t, map< array<int,3>,array<double,5> > >ThreeBody; // << This seems to work.
//  map< array<int,6>, map< array<int,3>,array<double,5> > >ThreeBody; // << This seems to work.
  map< array<int,9>,array<double,5> >ThreeBody; // 
//  map< orbindx3_t, vector< vector<double> > > ThreeBody;

  int rank_J; ///< Spherical tensor rank of the operator
  int rank_T; ///< Isotensor rank of the operator
  int parity; ///< Parity of the operator, 0=even 1=odd
  int particle_rank; ///< Maximum particle rank. Should be 2 or 3.

/// A list of which two-body channels can be connected by this operator
  vector<vector<int> > TwoBodyTensorChannels;

  int E2max; // I don't do anything with this yet...
  int E3max;

  ModelSpace * modelspace;
  bool hermitian;
  bool antihermitian;
  int nChannels;
  static double bch_transform_threshold;
  static double bch_product_threshold;
  static map<string, double> timer;

  void PrintTimes();


  //Constructors
  // In the future, consider using C++11 rvalues / move constructor to avoid copies in certain cases
  Operator(); /// Default constructor
  Operator(ModelSpace&);
  Operator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  Operator( const Operator& rhs);

  void AllocateTwoBody();
  void AllocateThreeBody();

  //Overloaded operators
  Operator& operator=( const Operator& rhs);
  Operator& operator+=( const Operator& rhs);
  Operator operator+( const Operator& rhs) const;
  Operator& operator-=( const Operator& rhs);
  Operator operator-( const Operator& rhs) const;
  Operator operator-( ) const;
  Operator& operator*=( const double rhs);
  Operator operator*( const double rhs) const;
  Operator& operator/=( const double rhs);
  Operator operator/( const double rhs) const;

  //Methods
  // One body setter/getters
  double GetOneBody(int i,int j) {return OneBody(i,j);};
  void SetOneBody(int i, int j, double val) { OneBody(i,j) = val;};

  //TwoBody setter/getters
  double GetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d) const;
  void   SetTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch_bra, int ch_ket, Ket &bra, Ket &ket) const;
  void   SetTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, Ket &bra, Ket& ket, double tbme);
  double GetTBME(int ch_bra, int ch_ket, int ibra, int iket) const;
  void   SetTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme);
  void   AddToTBME(int ch_bra, int ch_ket, int ibra, int iket, double tbme);
  double GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket) const;
  void   SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme);
  void   AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, Ket& bra, Ket& ket, double tbme);
  double GetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d) const;
  void   SetTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int j_bra, int p_bra, int t_bra, int j_ket, int p_ket, int t_ket, int a, int b, int c, int d, double tbme);

  // Scalar setters/getters for backwards compatibility
  double GetTBME(int ch, int a, int b, int c, int d) const;
  void   SetTBME(int ch, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int ch, int a, int b, int c, int d, double tbme);
  double GetTBME(int ch, Ket &bra, Ket &ket) const;
  void   SetTBME(int ch, Ket &bra, Ket& ket, double tbme);
  void   AddToTBME(int ch, Ket &bra, Ket& ket, double tbme);
  double GetTBME(int ch, int ibra, int iket) const;
  void   SetTBME(int ch, int ibra, int iket, double tbme);
  void   AddToTBME(int ch, int ibra, int iket, double tbme);
  double GetTBME(int j, int p, int t, Ket& bra, Ket& ket) const;
  void   SetTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  void   AddToTBME(int j, int p, int t, Ket& bra, Ket& ket, double tbme);
  double GetTBME(int j, int p, int t, int a, int b, int c, int d) const;
  void   SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);
  void   AddToTBME(int j, int p, int t, int a, int b, int c, int d, double tbme);


  double GetTBMEmonopole(int a, int b, int c, int d) const;
  double GetTBMEmonopole(Ket & bra, Ket & ket) const;


//// Three body setter getters
  double AddToThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);
  void   SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n, double V);
  double GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int i, int j, int k, int l, int m, int n);
  double GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n);


///// Some other three body methods

//  void SortThreeBodyOrbits(int& a,int& b,int& c);
  int SortThreeBodyOrbits(int a_in, int b_in, int c_in, int& a,int& b,int& c);
  double ThreeBodyRecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J);
//  int GetRecouplingCase(int a_in, int b_in, int c_in, int a, int b, int c);
  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};




  // Other setter-getters
  ModelSpace * GetModelSpace() const {return modelspace;};

  void Erase();
  void EraseZeroBody(){ZeroBody = 0;}; // set zero-body term to zero
  void EraseOneBody(){OneBody.zeros();}; // set all one-body terms to zero
  void EraseTwoBody(); // set all two-body terms to zero
  void EraseThreeBody(); // set all two-body terms to zero

  void SetHermitian() {hermitian=true;antihermitian=false;};
  void SetAntiHermitian() {antihermitian=true;hermitian=false;};
  void SetNonHermitian() {antihermitian=false;hermitian=false;};
  bool IsHermitian()const {return hermitian;};
  bool IsAntiHermitian()const {return antihermitian;};
  bool IsNonHermitian()const {return not (hermitian or antihermitian);};
  int GetParticleRank()const {return particle_rank;};
  int GetJRank()const {return rank_J;};
  int GetTRank()const {return rank_T;};
  int GetParity()const {return parity;};
  void SetParticleRank(int pr) {particle_rank = pr;};



  void ScaleZeroBody(double x);
  void ScaleOneBody(double x);
  void ScaleTwoBody(double x);

  // The actually interesting methods
  Operator DoNormalOrdering(); ///< Calls DoNormalOrdering2() or DoNormalOrdering3(), depending on the rank of the operator.
  Operator DoNormalOrdering2(); ///< Returns the normal ordered two-body operator
  Operator DoNormalOrdering3(); ///< Returns the normal ordered three-body operator

  Operator Commutator(  Operator& opright) ;
/// X.BCH_Product(Y) returns \f$Z\f$ such that \f$ e^{Z} = e^{X}e^{Y}\f$
/// by employing the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}([X,[X,Y]]+[Y,[Y,X]]) + \ldots \f]
  Operator BCH_Product(  Operator& )  ; 
/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// by employing the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [X,Y] + \frac{1}{2!}[X,[X,Y]] + \frac{1}{3!}[X,[X,[X,Y]]] + \ldots \f]
  Operator BCH_Transform(  Operator& ) ; 

/// Calculates the total kinetic energy (center of mass + relative) for each orbit in the model space,
/// assuming a harmonic oscillator basis.
/// \f[ t_{ij} = \frac{1}{2}\epsilon_{i} \delta_{ij} \f]
/// \f[ t_{ij} = \frac{1}{2} \sqrt{n_>(n_>+\ell+\frac{1}{2})} \delta_{j_ij_j}\delta_{\ell_i\ell_j} \delta_{n_>,n_<+1} \f]
  void CalculateKineticEnergy();
  void Eye(); ///< set to identity operator

  Operator CommutatorScalarScalar( Operator& opright) ;
  Operator CommutatorScalarTensor( Operator& opright) ;
  Operator CommutatorTensorTensor( Operator& opright) ;

/// Obtain the Frobenius norm of the operator, which here is 
/// defined as 
/// \f[ \|X\| = \sqrt{\|X_{(0)}\|^2 +\|X_{(1)}\|^2 +\|X_{(2)}\|^2 +\|X_{(3)}\|^2 } \f]
/// and
/// \f[ \|X_{(1)}\|^2 = \sum\limits_{ij} X_{ij}^2 \f]
  double Norm() const;
  double OneBodyNorm() const;
  double TwoBodyNorm() const;



  void PrintOneBody() const {OneBody.print();};
  void PrintTwoBody(int ch) const {TwoBody.at(ch).at(ch).print();};


  //Methods
  void Copy(const Operator& rhs);

  static void Set_BCH_Transform_Threshold(double x){bch_transform_threshold=x;};
  static void Set_BCH_Product_Threshold(double x){bch_product_threshold=x;};

  
  void DoPandyaTransformation(Operator&) ;
  void CalculateCrossCoupled(vector<arma::mat>&, vector<arma::mat>&) ; 

  void comm110ss( Operator& opright, Operator& opout) ; 
  void comm220ss( Operator& opright, Operator& opout) ;
  void comm111ss( Operator& opright, Operator& opout) ;
  void comm121ss( Operator& opright, Operator& opout) ;
  void comm221ss( Operator& opright, Operator& opout) ;
  void comm122ss( Operator& opright, Operator& opout) ;
  void comm222_pp_hhss( Operator& opright, Operator& opout) ;
  void comm222_phss( Operator& opright, Operator& opout) ;
  void comm222_pp_hh_221ss( Operator& opright, Operator& opout) ;

// make st and tt commutators

  void comm111st( Operator& opright, Operator& opout) ;
  void comm121st( Operator& opright, Operator& opout) ;
  void comm221st( Operator& opright, Operator& opout) ;
  void comm122st( Operator& opright, Operator& opout) ;
  void comm222_pp_hhst( Operator& opright, Operator& opout) ;
  void comm222_pp_hh_221st( Operator& opright, Operator& opout) ;
  void comm222_phst( Operator& opright, Operator& opout) ;
  void comm222_phst_pandya( Operator& opright, Operator& opout) ;

/*
  void comm111tt(const Operator& opright, Operator& opout) const;
  void comm121tt(const Operator& opright, Operator& opout) const;
  void comm221tt(const Operator& opright, Operator& opout) const;
  void comm122tt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hhtt(const Operator& opright, Operator& opout) const;
  void comm222_phtt(const Operator& opright, Operator& opout) const;
  void comm222_pp_hh_221tt(const Operator& opright, Operator& opout) const;
*/

};

/// Non member function, multiply by scalar from left side
Operator operator*(const double lhs, const Operator& rhs);
Operator operator*(const double lhs, const Operator&& rhs);



#endif

