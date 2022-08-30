///////////////////////////////////////////////////////////////////////////////////
//    Commutator.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef Commutator_h
#define Commutator_h 1

#include "Operator.hh"
//#include "DaggerOperator.hh"
#include "TwoBodyME.hh"
#include "ThreeLegME.hh"
#include "armadillo"
#include <map>
#include <unordered_set>
#include <deque>
#include <array>
#include <string>



namespace Commutator{

  extern bool use_goose_tank_correction;
  extern bool use_brueckner_bch;
  extern bool use_imsrg3;
  extern bool use_imsrg3_n7;
  extern bool perturbative_triples;
  extern bool bch_skip_ieq1;
  extern bool only_2b_omega;
  extern bool imsrg3_no_qqq;
  extern bool imsrg3_valence_2b;
  extern bool discard_0b_from_3b;
  extern bool discard_1b_from_3b;
  extern bool discard_2b_from_3b;
  extern double bch_transform_threshold;
  extern double bch_product_threshold;
  extern double threebody_threshold;
  extern double imsrg3_dE6max;

  extern std::map<std::string,bool> comm_term_on; // This allows turning on/off individual IMSRG(3) commutator terms for testing.

  void Set_BCH_Transform_Threshold(double x);
  void Set_BCH_Product_Threshold(double x);
  void SetThreebodyThreshold(double x);
  void SetUseBruecknerBCH(bool tf);
  void SetUseGooseTank(bool tf);
  void SetUseIMSRG3(bool tf);
  void SetUseIMSRG3N7(bool tf);
  void SetOnly2bOmega(bool tf);
  void SetBCHSkipiEq1(bool tf);
  void SetIMSRG3Noqqq(bool tf);
  void SetIMSRG3valence2b(bool tf);

  void TurnOffTerm( std::string term ) ;
  void TurnOnTerm( std::string term ) ;

  void Discard0bFrom3b( bool tf);
  void Discard1bFrom3b( bool tf);
  void Discard2bFrom3b( bool tf);

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
  void DoPandyaTransformation_SingleChannel_XandY(const Operator& X, const Operator& Y, arma::mat& X2_CC_ph, arma::mat& Y2_CC_ph, int ch_cc);
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
  void comm223ss_new( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm223ss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented and tested.
  void comm233_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) ;     // implemented and tested.
  void comm233_pp_hhss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;     // implemented and tested.
  void comm233_phss( const Operator& X, const Operator& Y, Operator& Z ) ;        // implemented and tested.
  void comm233_phss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;        // implemented and tested.

  void comm333_ppp_hhhss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented. andtested.
  void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) ;           // implemented. and tested.
  void comm333_pph_hhpss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;     

  bool check_2b_channel_Tz_parity( const Operator& Op, Orbit& o1, Orbit&o2, Orbit& o3, Orbit& o4 );

  Operator GooseTankUpdate( const Operator& Omega, const Operator& Nested);



// scalar-tensor commutators

  void DoTensorPandyaTransformation(const Operator& Z, std::map<std::array<index_t,2>,arma::mat>&) ;
  void DoTensorPandyaTransformation_SingleChannel(const Operator& Z, arma::mat& X, int ch_bra_cc, int ch_ket_cc) ;
  void AddInverseTensorPandyaTransformation(Operator& Z, const std::map<std::array<index_t,2>,arma::mat>&);

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



  void prod110ss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void prod111ss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void prod112ss( const Operator& X, const Operator& Y, Operator& Z ) ; 

  void comm232_new_Determine1BChannels(
    const Operator& Z,
    const Operator& Y,
    std::map<std::array<int,3>,std::vector<size_t>>& local_one_body_channels,
    std::map<std::array<int,3>,std::vector<size_t>>& external_local_one_body_channels
  );

  void comm232_new_Populate1BChannel(
    const Operator& Z,
    const Operator& Y,
    const std::map<int,double>& e_fermi,
    const std::map<std::array<int,3>,std::vector<size_t>>& local_one_body_channels,
    const std::map<std::array<int,3>,std::vector<size_t>>& external_local_one_body_channels,
    const std::vector<std::array<int,3>>& obc_keys,
    std::map<std::array<int,3>,std::vector<std::array<size_t,4>>>& klj_list,
    std::map<std::array<int,3>,arma::mat>& ZMAT_list
  );

  void comm232_new_Determine3BStatesIn1BChannel(
    const Operator& Y,
    const Operator& Z,
    const std::map<std::array<int,3>,std::vector<std::array<size_t,4>>>& klj_list,
    std::array<int, 3> obc_key,
    std::vector<size_t>& abc_list,
    std::vector<double>& abc_occ_list
  );

  void comm232_new_GenerateRequiredRecouplings(
    const Operator& Z,
    const Operator& Y,
    const std::vector<size_t>& abc_list,
    const std::vector<std::array<size_t, 4>>& klj_list_i,
    const std::map<int, double>& e_fermi,
    size_t dim_abc,
    size_t dim_klj,
    std::unordered_set<size_t>& kljJJ_needed);

  void comm232_new_ComputeRequiredRecouplings(
    const Operator& Y,
    const std::unordered_set<size_t>& kljJJ_needed,
    std::vector<double>& recoupling_cache,
    std::unordered_map<size_t, size_t>& recoupling_cache_lookup);

  void comm232_new_FillMatrices(const Operator& Z, const Operator& X, const Operator Y, size_t dim_abc, size_t dim_i, size_t dim_klj, 
    int j2i,
    bool x_has_3,
    bool y_has_3,
    const std::vector<size_t>& abc_list,
    const std::vector<std::array<size_t, 4>>& klj_list_i,
    const std::map<int, double>& e_fermi,
    const std::vector<double>& abc_occ_list,
    const std::vector<size_t>& obc_orbits,
    const std::vector<double>& recoupling_cache,
    const std::unordered_map<size_t, size_t>& recoupling_cache_lookup,
    arma::mat& X2MAT,
    arma::mat& Y2MAT,
    arma::mat& X3MAT,
    arma::mat& Y3MAT
  );

  void comm232_new_Unpack2BResult(
    const Operator& X,
    const Operator& Y,
    size_t nch,
    const std::map<std::array<int,3>,std::vector<size_t>>& external_local_one_body_channels,
    const std::map<std::array<int,3>,std::vector<std::array<size_t,4>>>& klj_list,
    const std::map<std::array<int,3>,arma::mat>& ZMAT_list,
    Operator& Z
    );

}

#endif
