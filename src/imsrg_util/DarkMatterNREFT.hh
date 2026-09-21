#ifndef DarkMatterNREFT_hh
#define DarkMatterNREFT_hh

#include "ModelSpace.hh"
#include "Operator.hh"



// Operators for Dark matter scattering as defined in
// Gazda, Catena, and Forssen, Phys Rev D 95 103011 (2017)
// The implementation is essentially copied from Daniel Gazda's
// fortran implementation and translated to C++. Any errors
// are almost certainly Ragnar's fault.
//                                         - SRS  June 21 2018
namespace DM_NREFT
{

  double jho( int a, int la, int nb, int lb, int L, double y);
  double jdmho( int a, int la, int nb, int lb, int L, double y);
  double jdpho( int a, int la, int nb, int lb, int L, double y);

  double PhiF(  int la, int j2a, int j2b );
  double PhiS1( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS2( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS3( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );
  double PhiS4( int na, int la, int j2a, int nb, int lb, int j2b, int J, double y );

  Operator Ms(      ModelSpace&, std::string IsoSV, int J, int L, double q ); // helper function
  Operator Mg(      ModelSpace&, std::string IsoSV, int J, int L, double q ); // helper function

  // These are the operators
  Operator M(       ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigma(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigmap(  ModelSpace&, std::string IsoSV, int J, double q );
  Operator Sigmapp( ModelSpace&, std::string IsoSV, int J, double q );
  Operator Delta(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Deltap(  ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phip(    ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phitp(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Phipp(   ModelSpace&, std::string IsoSV, int J, double q );
  Operator Omega(   ModelSpace&, std::string IsoSV, int J, double q ); // Warning: neither Symmetrize nor AntiSymmetrize 
  Operator Omegat(  ModelSpace&, std::string IsoSV, int J, double q );


}


#endif
