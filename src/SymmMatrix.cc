
#include "SymmMatrix.hh"
#include <armadillo>





template <typename T>
SymmMatrix<T>::SymmMatrix()
 : herm(1)
{}

template <typename T>
SymmMatrix<T>::SymmMatrix(size_t dim)
 : dimension(dim), herm(1), matrix_data( (dim+1)*(dim+2)/2, 0 )
{}

template <typename T>
SymmMatrix<T>::SymmMatrix(size_t dim, int hrm)
 : dimension(dim), herm(hrm), matrix_data( (dim+1)*(dim+2)/2, 0 )
{}


template <typename T>
SymmMatrix<T>::SymmMatrix( arma::mat mtx )
 : herm(1)
{
  dimension = mtx.n_rows;
  if (dimension != mtx.n_cols)
  {
    std::cout << "WHY DID YOU GIVE ME A NON-SQUARE MATRIX????" << std::endl;
    exit(0);
  }
  for ( size_t i=0; i<dimension; i++)
  {
    for (size_t j=0; j<i; j++)
    {
      
    }
  }
}




template <typename T>
T SymmMatrix<T>::Get(size_t i, size_t j) const
{
  if (i>=dimension or j>=dimension)
  {
    std::cout << "i,j = " << i << " " << j << "  dimension = " << dimension << std::endl;
    throw std::domain_error( "Trouble in SymmMatrix" );
  }
  if (i>=j)
  {
    return matrix_data[ Index2to1(i,j) ];
  }
  else
  {
    return herm * matrix_data[ Index2to1(j,i) ];
  }

}


template <typename T>
void SymmMatrix<T>::Put(size_t i, size_t j, T val)
{
  if (i>=dimension or j>=dimension)
  {
    throw std::domain_error( "Trouble in SymmMatrix" );
  }
  if ((i==j) and herm==-1) return;
  if (i>=j)
  {
    matrix_data[ Index2to1(i,j) ] = val;
  }
  else
  {
    matrix_data[ Index2to1(j,i) ] = herm * val;
  }
}




///        j
///   [ 0       ]
/// i [ 1 4     ]
///   [ 2 5 7   ]
///   [ 3 6 8 9 ]
///
///  If we stored the full matrix, we would have
/// I = dim * j + i
/// but we are missing a triangle of size j(j+1)/2
/// before that point in the matrix
/// so we lower the index by that much.
/// I = dim * j + i - j(j+1)/2
template <typename T>
size_t SymmMatrix<T>::Index2to1(size_t i, size_t j) const
{
  return (2*dimension - j - 1)*j/2 + i ;

}


template <typename T>
void SymmMatrix<T>::Index1to2(size_t I, size_t& i, size_t& j) const
{
  j=0;
  while ( Index2to1(0,j+1) <I) j++;
  i = I - Index2to1(0,j);
}



template <typename T>
arma::Mat<T> SymmMatrix<T>::FullMatrix() const
{
   arma::Mat<T> full_mat(dimension,dimension,arma::fill::zeros);
   for (size_t i=0;i<dimension;i++)
   {
     for (size_t j=0;j<i; j++)
     {
       T value = matrix_data[ Index2to1(i,j) ];
       full_mat(i,j) = value;
       full_mat(j,i) = herm*  value;
     }
     if (herm>0) full_mat(i,i) = matrix_data[ Index2to1(i,i) ];
   }
   return full_mat;
}

template <typename T>
void SymmMatrix<T>::zeros()
{
  matrix_data.assign( matrix_data.size(), 0.);
}

template <typename T>
double SymmMatrix<T>::Norm()
{
  double norm = 0;
  for ( auto v : matrix_data ) norm += v*v;
  return sqrt(norm);
}



template <typename T>
T SymmMatrix<T>::operator()(size_t i, size_t j)
{
  return Get(i,j);
}

template <typename T>
SymmMatrix<T>& SymmMatrix<T>::operator=( const SymmMatrix<T>& rhs)
{
  dimension = rhs.dimension;
  herm = rhs.herm;
  matrix_data = rhs.matrix_data;
  return *this;
}

//template <typename T>
//SymmMatrix<T>& SymmMatrix<T>::operator=( const SymmMatrix<T> rhs)
//{
//  dimension = rhs.dimension;
//  herm = rhs.herm;
//  matrix_data = rhs.matrix_data;
//  return *this;
//}

template <typename T>
SymmMatrix<T>& SymmMatrix<T>::operator*=( double rhs)
{
  for ( auto& v : matrix_data) v*= rhs;
  return *this;
}

template <typename T>
SymmMatrix<T> SymmMatrix<T>::operator*( double rhs)
{
  SymmMatrix<T> lhs (*this);
  lhs *= rhs;
  return lhs;
}

template <typename T>
SymmMatrix<T>& SymmMatrix<T>::operator+=( const SymmMatrix<T>& rhs)
{
  for (size_t i=0;i<matrix_data.size(); i++)
  {
    matrix_data[i] += rhs.matrix_data[i];
  }
  return *this;
}

template <typename T>
SymmMatrix<T> SymmMatrix<T>::operator+( const SymmMatrix<T>& rhs)
{
  SymmMatrix<T> lhs (*this);
  lhs += rhs;
  return lhs;
}

template <typename T>
SymmMatrix<T>& SymmMatrix<T>::operator-=( const SymmMatrix<T>& rhs)
{
  for (size_t i=0;i<matrix_data.size(); i++)
  {
    matrix_data[i] -= rhs.matrix_data[i];
  }
  return *this;
}

template <typename T>
SymmMatrix<T> SymmMatrix<T>::operator-( const SymmMatrix<T>& rhs)
{
  SymmMatrix<T> lhs (*this);
  lhs -= rhs;
  return lhs;
}


template class SymmMatrix<float>;
template class SymmMatrix<double>;
