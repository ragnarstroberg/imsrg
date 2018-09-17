#ifndef DaggerOperator_h
#define DaggerOperator_h 1


#include "Operator.hh"


class DaggerOperator: public Operator
{
 public:


  index_t Q_space_orbit; // Modelspace orbit used as a placeholder for where the a+ particle came from.


  //Constructors
 ~DaggerOperator(); // Explicitly define destructor/constructor in order
  DaggerOperator(); // to maintain a counter of how many operators exist
  DaggerOperator(ModelSpace&); ///< Construct a 2-body scalar operator
  DaggerOperator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  DaggerOperator(ModelSpace&, index_t Q_space_orbit);
  DaggerOperator( const DaggerOperator& rhs); ///< Copy constructor
  DaggerOperator( DaggerOperator&&);

  //Assignment operator cannot be inherited.
  DaggerOperator& operator=( const DaggerOperator& rhs);

  index_t GetQSpaceOrbit() const {return Q_space_orbit;};
  void SetQSpaceOrbit(index_t Q) {Q_space_orbit = Q;};


//  DaggerOperator BCH_Product(  Operator& )  ; 
  DaggerOperator BCH_Transform( const Operator& ) ; 
  DaggerOperator Standard_BCH_Transform( const Operator& ) ; 
//  DaggerOperator Brueckner_BCH_Transform( const Operator& ) ;



  friend DaggerOperator Commutator(const Operator& X, const DaggerOperator& Y) ;
  void SetToCommutator( const Operator& X, const DaggerOperator& Y);
  void CommutatorScalarDagger( const Operator& X, const DaggerOperator& Y);

 
  // commutator terms, 211 means [two legs, one leg] => one leg
  // sd means scalar-dagger
  void comm211sd( const Operator& X, const DaggerOperator& Y) ; 
  void comm231sd( const Operator& X, const DaggerOperator& Y) ;
  void comm431sd( const Operator& X, const DaggerOperator& Y) ;
  void comm413_233sd( const Operator& X, const DaggerOperator& Y) ; 
//  void comm413sd( const Operator& X, const DaggerOperator& Y) ; 
//  void comm233sd( const Operator& X, const DaggerOperator& Y) ; 
  void comm433sd_pphh( const Operator& X, const DaggerOperator& Y) ; 
  void comm433sd_ph( const Operator& X, const DaggerOperator& Y) ; 

  void comm433_pp_hh_431sd( const Operator& X, const Operator& Y ) ; 
  void ConstructDaggerMpp_Mhh(const Operator& X, const Operator& Y, TwoBodyME& Mpp, TwoBodyME& Mhh) const;


};

#endif 
