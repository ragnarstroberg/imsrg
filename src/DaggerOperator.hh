#ifndef DaggerOperator_h
#define DaggerOperator_h 1


#include "Operator.hh"


class DaggerOperator: public Operator
{
 public:




  //Constructors
 ~DaggerOperator(); // Explicitly define destructor/constructor in order
  DaggerOperator(); // to maintain a counter of how many operators exist
  DaggerOperator(ModelSpace&); ///< Construct a 2-body scalar operator
  DaggerOperator(ModelSpace&, int Jrank, int Trank, int Parity, int part_rank);
  DaggerOperator( const DaggerOperator& rhs); ///< Copy constructor
  DaggerOperator( DaggerOperator&&);

  //Assignment operator cannot be inherited.
  DaggerOperator& operator=( const DaggerOperator& rhs);



  friend DaggerOperator Commutator(const Operator& X, const DaggerOperator& Y) ;
  void SetToCommutator( const Operator& X, const DaggerOperator& Y);
  void CommutatorScalarDagger( const Operator& X, const DaggerOperator& Y);

 
  // commutator terms, 211 means [two legs, one leg] => one leg
  // sd means scalar-dagger
  void comm211sd( const Operator& X, const DaggerOperator& Y) ; 
  void comm231sd( const Operator& X, const DaggerOperator& Y) ;
  void comm431sd( const Operator& X, const DaggerOperator& Y) ;
  void comm413sd( const Operator& X, const DaggerOperator& Y) ; 
  void comm233sd( const Operator& X, const DaggerOperator& Y) ; 
  void comm433sd_pphh( const Operator& X, const DaggerOperator& Y) ; 
  void comm433sd_ph( const Operator& X, const DaggerOperator& Y) ; 



};

#endif 
