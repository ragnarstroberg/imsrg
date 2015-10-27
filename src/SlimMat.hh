#ifndef SlimMat_h
#define SlimMat_h

#include <armadillo>
#include <vector>

using namespace std;

template <class Type>
class SlimMat
{
 public:
  size_t Dim;
  bool symmetric;
  vector<Type> Data;

  SlimMat();
  SlimMat(bool sym);
  SlimMat(int dim);
  SlimMat(const arma::Mat<Type>& m);
  SlimMat(int dim, bool sym);
  SlimMat(const arma::Mat<Type>& m, bool sym);
  
  arma::Mat<Type> GetMat();
  Type operator()(size_t i,size_t j) const; ///< read-only access
  SlimMat<Type>& operator=(const arma::Mat<Type>&);
  void SetSymmetric(){symmetric=true;};
  void SetAntisymmetric(){symmetric=false;};
  bool IsSymmetric(){return symmetric;};
  Type* memptr();
  size_t Size();

};

#endif
