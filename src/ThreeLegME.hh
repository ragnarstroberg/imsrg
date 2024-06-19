#ifndef ThreeLegME_h
#define ThreeLegME_h 1

#include "ModelSpace.hh"
#include <map>
#include <armadillo>

struct ThreeLegME
{
  ModelSpace* modelspace;
  std::map<size_t, arma::mat> MatEl;


  ThreeLegME() : modelspace(NULL) {};
  ThreeLegME(ModelSpace* ms): modelspace(ms) {Allocate();};

  void Allocate();
  void Deallocate();
  arma::mat& GetMatrix(size_t ch){ return MatEl.at(ch);};
  const arma::mat& GetMatrix(size_t ch)const { return MatEl.at(ch);};

  double GetME(size_t ch, size_t i, size_t j, size_t k)const ;
  double GetME_norm(size_t ch, size_t i, size_t j, size_t k)const ;
  double GetME_J(int J, size_t i, size_t j, size_t k)const ;

  void SetME(size_t ch, size_t i, size_t j, size_t k, double me);
  void AddToME(size_t ch, size_t i, size_t j, size_t k, double me);
  void AddToME_J(int J, size_t i, size_t j, size_t k, double me);

  ThreeLegME& operator*=(const double rhs);
  ThreeLegME  operator*(const double rhs) const;
  ThreeLegME& operator+=(const ThreeLegME& rhs);
  ThreeLegME  operator+(const ThreeLegME& rhs) const;
  ThreeLegME& operator-=(const ThreeLegME& rhs);
  ThreeLegME  operator-(const ThreeLegME& rhs) const;
  ThreeLegME  operator-()const ;

  double Norm() const;
  void Erase(); // set all matrix elements to zero

};


#endif
