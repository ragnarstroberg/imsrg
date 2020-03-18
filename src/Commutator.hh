
#include "Operator.hh"
//#include "DaggerOperator.hh"
#include "TwoBodyME.hh"
#include "ThreeLegME.hh"
#include "armadillo"
#include <map>
#include <deque>
#include <array>
#include <string>



namespace Commutator{

  extern bool use_goose_tank_correction;
  extern bool use_brueckner_bch;
  extern bool use_imsrg3;
  extern bool only_2b_omega;
  extern double bch_transform_threshold;
  extern double bch_product_threshold;

  void Set_BCH_Transform_Threshold(double x);
  void Set_BCH_Product_Threshold(double x);
  void SetUseBruecknerBCH(bool tf);
  void SetUseGooseTank(bool tf);
  void SetUseIMSRG3(bool tf);
  void SetOnly2bOmega(bool tf);


  Operator Commutator(const Operator& X, const Operator& Y) ; 
  Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) ;
  Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) ;
  Operator CommutatorScalarDagger( const Operator& X, const Operator& Y) ;


  Operator BCH_Product(  Operator& X, Operator& Y )  ; 
  Operator BCH_Transform( const Operator& Op, const Operator& Omega ) ; 
  Operator Standard_BCH_Transform( const Operator& Op, const Operator& Omega ) ; 
  Operator Brueckner_BCH_Transform( const Operator& Op, const Operator& Omega ) ;

  double EstimateBCHError( Operator& Omega, Operator H);

  std::deque<arma::mat> InitializePandya(Operator& Z, size_t nch, std::string orientation);
  void DoPandyaTransformation(const Operator& Z, std::deque<arma::mat>&, std::string orientation) ;
  void DoPandyaTransformation_SingleChannel(const Operator& Z, arma::mat& X, int ch_cc, std::string orientation) ;
//  void AddInversePandyaTransformation(Operator& Z, const std::deque<arma::mat>&);
  void AddInversePandyaTransformation(const std::deque<arma::mat>& Zbar, Operator& Z);   // Changed from the above declaration. Not sure how this was compiling...
  void AddInversePandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, int ch_cc);


  void comm110ss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void comm220ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm111ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm121ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm221ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm122ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_phss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z) ;


  void ConstructScalarMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh);
//  void ConstructScalarMpp_Mhh_GooseTank(const Operator& X, const Operator& Y, Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh) ;


// IMSRG(3) commutators. Still a work in progress...
  void comm330ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm331ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm231ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm231ss_slow( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.

  void comm132ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  size_t Hash_comm232_key( std::array<size_t,5>& kljJJ );
  void comm232ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm232ss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm232ss_slow( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z ) ; // implemented and tested.
//  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z ) ;      // implemented and tested.
  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z ) ;      // implemented and tested.
  void comm332_pphhss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;      // implemented and tested.
  
  void comm133ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm223ss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm223ss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm233_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) ;     // implemented and tested.
  void comm233_pp_hhss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;     // implemented and tested.
  void comm233_phss( const Operator& X, const Operator& Y, Operator& Z ) ;        // implemented and tested.
  void comm233_phss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;        // implemented and tested.

  void comm333_ppp_hhhss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented. not tested.
  void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented. not tested.
  void comm333_pph_hhpss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented. not tested.

  bool check_2b_channel_Tz_parity( const Operator& Op, Orbit& o1, Orbit&o2, Orbit& o3, Orbit& o4 );

  Operator GooseTankUpdate( const Operator& Omega, const Operator& Nested);



// scalar-tensor commutators

  void DoTensorPandyaTransformation(const Operator& Z, std::map<std::array<index_t,2>,arma::mat>&) ;
  void DoTensorPandyaTransformation_SingleChannel(const Operator& Z, arma::mat& X, int ch_bra_cc, int ch_ket_cc) ;
  void AddInverseTensorPandyaTransformation(Operator& Z, const std::map<std::array<index_t,2>,arma::mat>&);
  void AddInverseTensorPandyaTransformation_SingleChannel(Operator& Z, arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc);

  void comm111st( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm121st( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm122st( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_pp_hh_221st( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_phst( const Operator& X, const Operator& Y, Operator& Z) ;




  // commutator terms involving a dagger operator. 211 means [two legs, one leg] => one leg
  // sd means scalar-dagger
  void comm211sd( const Operator& X, const Operator& Y, Operator& Z) ; 
  void comm231sd( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm431sd( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm413_233sd( const Operator& X, const Operator& Y, Operator& Z) ; 
  void comm433sd_pphh( const Operator& X, const Operator& Y, Operator& Z) ; 
  void comm433sd_ph( const Operator& X, const Operator& Y, Operator& Z) ; 
  void comm433sd_ph_dumbway( const Operator& X, const Operator& Y, Operator& Z) ; // Do it with loops, not matmult. Easier to check, but much slower.

  void comm433_pp_hh_431sd( const Operator& X, const Operator& Y, Operator& Z ) ; 
//  void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, TwoBodyME& Mpp, TwoBodyME& Mhh);
  void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, const Operator& Z, ThreeLegME& Mpp, ThreeLegME& Mhh);
  void DoPandyaTransformation_SingleChannel_Dagger(const Operator& Z, arma::mat& X, int ch_cc) ;
  void AddInversePandyaTransformation_Dagger(const std::deque<arma::mat>& Zbar, Operator& Z );






}
