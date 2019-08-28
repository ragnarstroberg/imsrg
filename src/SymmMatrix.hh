#ifndef SymmMatrix_h
#define SymmMatrix_h


#include <vector>
#include <armadillo>


template <typename T>
class SymmMatrix
{
  std::vector<T> matrix_data;

  size_t dimension;
  int herm; // either +1 or -1


 public:

  SymmMatrix();
  SymmMatrix(size_t dim);
  SymmMatrix(size_t dim, int herm);
  SymmMatrix( arma::mat ); // at some point, we can worry about float to double casts...

  T operator()(size_t a, size_t b); // element access
  T operator()(size_t a, size_t b) const; // maybe want this?

  SymmMatrix<T>& operator=( SymmMatrix<T> rhs);
  SymmMatrix<T>& operator+=( SymmMatrix<T> rhs);
  SymmMatrix<T> operator+( SymmMatrix<T> rhs);
  SymmMatrix<T>& operator-=( SymmMatrix<T> rhs);
  SymmMatrix<T> operator-( SymmMatrix<T> rhs);
  SymmMatrix<T>& operator*=( double rhs);
  SymmMatrix<T> operator*( double rhs);

  size_t Size() const {return dimension;};
  T Access(size_t a, size_t b) const;
  void Put(size_t i, size_t j, T val);

  arma::Mat<T> FullMatrix() const;
  arma::rowvec& Row(size_t irow);
  arma::colvec& Col(size_t icol);

  bool IsSymmetric() const ;
  bool IsAntisymmetric() const;



  size_t Index2to1(size_t i, size_t j) const; // roll 2 indices into 1
  void Index1to2(size_t I, size_t& i, size_t& j) const ; // roll 2 indices into 1

};



#endif
