#ifndef DaggerOperator_h
#define DaggerOperator_h 1


#include "ModelSpace.hh"
#include "Operator.hh"
#include <armadillo>


struct ThreeLegME
{
  ModelSpace* modelspace;
  std::map<size_t, arma::mat> MatEl;


  ThreeLegME() {};
  ThreeLegME(ModelSpace& ms): modelspace(&ms) {Allocate();};

  void Allocate();
  arma::mat& GetMatrix(size_t ch){ return MatEl.at(ch);};
  const arma::mat& GetMatrix(size_t ch)const { return MatEl.at(ch);};

  double GetME(size_t ch, size_t i, size_t j, size_t k)const ;
  double GetME_norm(size_t ch, size_t i, size_t j, size_t k)const ;
  double GetME_J(int J, size_t i, size_t j, size_t k)const ;

  void SetME(size_t ch, size_t i, size_t j, size_t k, double me);
  void AddToME(size_t ch, size_t i, size_t j, size_t k, double me);
  void AddToME_J(int J, size_t i, size_t j, size_t k, double me);

  ThreeLegME& operator*=(const double rhs);
  ThreeLegME operator*(const double rhs) const;
  ThreeLegME& operator+=(const ThreeLegME& rhs);
  ThreeLegME operator+(const ThreeLegME& rhs) const;
  ThreeLegME& operator-=(const ThreeLegME& rhs);
  ThreeLegME operator-(const ThreeLegME& rhs) const;
  ThreeLegME operator-()const ;
  

};



class DaggerOperator : public Operator
{
 public:

  ThreeLegME ThreeLeg;

  DaggerOperator(); ///< Default constructor
  DaggerOperator(ModelSpace& ms); ///< Construct a 2-body scalar operator
  DaggerOperator(ModelSpace& ms, size_t Qorb); ///< Construct a 2-body scalar operator

  //Overloaded operators

  DaggerOperator& operator=(const DaggerOperator& rhs) ;
//  DaggerOperator& operator=(DaggerOperator&& rhs) ;

//  DaggerOperator& operator=( const DaggerOperator& rhs);
  DaggerOperator& operator+=( const DaggerOperator& rhs);
  DaggerOperator operator+( const DaggerOperator& rhs) const;
  DaggerOperator& operator-=( const DaggerOperator& rhs);
  DaggerOperator operator-( const DaggerOperator& rhs) const;
  DaggerOperator operator-( ) const;
  DaggerOperator& operator*=( const double rhs);
  DaggerOperator operator*( const double rhs) const;
  DaggerOperator& operator/=( const double rhs);
  DaggerOperator operator/( const double rhs) const;

  double ThreeLegNorm();

  void Erase();

};



#endif
