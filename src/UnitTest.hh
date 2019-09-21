#ifndef UnitTest_h
#define UnitTest_h

#include "ModelSpace.hh"
#include "Operator.hh"
#include "Commutator.hh"

class UnitTest
{

 static uint64_t random_seed;

 public:

  ModelSpace* modelspace;

  UnitTest(ModelSpace&);

  void SetRandomSeed( uint64_t s ){ random_seed = s;};

  Operator RandomOp( ModelSpace& modelspace, int jrank, int tz, int parity, int particle_rank, int hermitian);

  Operator RandomDaggerOp( ModelSpace& modelspace, index_t Q);


  double GetMschemeMatrixElement_1b( const Operator& Op, int a, int ma, int b, int mb );
  double GetMschemeMatrixElement_2b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md );
  double GetMschemeMatrixElement_3b( const Operator& Op, int a, int ma, int b, int mb, int c, int mc, int d, int md, int e, int me, int f, int mf );

  double GetMschemeMatrixElement_1leg( const Operator& Op, int a, int ma );
  double GetMschemeMatrixElement_3leg( const Operator& Op, int a, int ma, int b, int mb, int c, int mc );

//  void Test3BodyAntisymmetry();
  void Test3BodyAntisymmetry(Operator& Y);
  void Test3BodyHermiticity(Operator& Y);

  // test strategy: Fill two random operators, calculate a specific commutator term
  // using the J-coupled expression, and in m-scheme (where the formula is simpler)
  // and make sure that they give the same answer
  void TestCommutators();
//  void TestCommutators3();
  void TestCommutators3(Operator& X, Operator& Y); 

  void TestDaggerCommutators(index_t Q);


  bool Test_comm110ss( const Operator& X, const Operator& Y );
  bool Test_comm220ss( const Operator& X, const Operator& Y );
  bool Test_comm111ss( const Operator& X, const Operator& Y );

  bool Test_comm121ss( const Operator& X, const Operator& Y ); 
  bool Test_comm221ss( const Operator& X, const Operator& Y ); 

  bool Test_comm122ss( const Operator& X, const Operator& Y );

  bool Test_comm222_pp_hhss( const Operator& X, const Operator& Y );

  bool Test_comm222_phss( const Operator& X, const Operator& Y ) ;


  bool Test_comm330ss( const Operator& X, const Operator& Y );
  bool Test_comm331ss( const Operator& X, const Operator& Y );
  bool Test_comm231ss( const Operator& X, const Operator& Y );

  bool Test_comm132ss( const Operator& X, const Operator& Y );

  bool Test_comm223ss( const Operator& X, const Operator& Y );


  bool Test_comm211sd( const Operator& X, const Operator& Y );
  bool Test_comm231sd( const Operator& X, const Operator& Y );
  bool Test_comm431sd( const Operator& X, const Operator& Y );
  bool Test_comm233sd( const Operator& X, const Operator& Y );
  bool Test_comm413sd( const Operator& X, const Operator& Y );
  bool Test_comm433_pp_hh_sd( const Operator& X, const Operator& Y );
  bool Test_comm433sd_ph( const Operator& X, const Operator& Y );

};


#endif
