#include "SlimMat.hh"


template< class Type>
SlimMat<Type>::SlimMat()
: Dim(0),symmetric(true)
{}

template< class Type>
SlimMat<Type>::SlimMat(bool sym)
: Dim(0),symmetric(sym)
{}

template< class Type>
SlimMat<Type>::SlimMat(int dim)
: Dim(dim),symmetric(true)
{}

template< class Type>
SlimMat<Type>::SlimMat(int dim, bool sym)
: Dim(dim),symmetric(sym),Data(Dim*(Dim+1)/2)
{
}

template< class Type>
SlimMat<Type>::SlimMat(const arma::Mat<Type>& m)
: Dim(m.size()), symmetric(true), Data(Dim*(Dim+1)/2)
{
  Data.resize(0);
  for (int i=0;i<Dim;i++)
     Data.insert( Data.end(), m.begin_col(i)+i, m.end_col(i) );
}

template< class Type>
SlimMat<Type>::SlimMat(const arma::Mat<Type>& m, bool sym)
: Dim(m.size()), symmetric(sym),Data(Dim*(Dim+1)/2)
{
  Data.resize(0);
  for (int i=0;i<Dim;i++)
     Data.insert(Data.end(),m.begin_col(i)+i,m.end_col(i));
}

template< class Type>
SlimMat<Type>& SlimMat<Type>::operator=(const arma::Mat<Type>& m)
{
  Data.reserve(Dim*(Dim+1)/2);
  Data.resize(0);
  for (int i=0;i<Dim;i++)
     Data.insert(Data.end(),m.begin_col(i)+i,m.end_col(i));
  return *this;
}

// Get the full armadillo matrix back out
template< class Type>
arma::Mat<Type> SlimMat<Type>::GetMat()
{
  arma::Mat<Type> m(Dim,Dim,arma::fill::zeros);
  for (size_t i=0;i<Dim;i++)
  {
     size_t n = (Dim-1)*i-i*(i-1)/2;
     std::copy( Data.begin()+n+i, Data.begin()+n+Dim, m.memptr()+i*(Dim+1) );
  }
  if (symmetric)
    m = arma::symmatl(m);
  else
    m -= m.t();
  return m;
}


// Element access via i,j.
// No bounds checking, because we like to live dangerously.
template< class Type>
Type SlimMat<Type>::operator()(size_t i, size_t j) const
{
  if (i>=j)
    return Data[ (Dim-1)*i+j-i*(i-1)/2 ];
  else
    return symmetric ?  Data[(Dim-1)*j+i-j*(j-1)/2] : -Data[(Dim-1)*j+i-j*(j-1)/2] ;
}


template< class Type>
Type* SlimMat<Type>::memptr(){ return Data.data();}

template< class Type>
size_t SlimMat<Type>::Size(){return Data.size();}


// explicit instantiations
template class SlimMat<float>;
template class SlimMat<double>;


