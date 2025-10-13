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

  bool TestNormalOrdering(Operator& Op);

//  void Test3BodySetGet(Operator& Y);

  // test strategy: Fill two random operators, calculate a specific commutator term
  // using the J-coupled expression, and in m-scheme (where the formula is simpler)
  // and make sure that they give the same answer
//  void TestCommutators();
  bool TestCommutators();
  bool TestCommutators_Tensor(Operator& X, Operator& Y);
//  bool TestCommutators_Tensor();
  bool TestCommutators_IsospinChanging();
  bool TestCommutators_ParityChanging();
//  void TestCommutators3();
//  void TestCommutators3(Operator& X, Operator& Y); 
  void TestCommutators3(Operator& X, Operator& Y, std::vector<std::string>& skiplist ); 

  void TestDaggerCommutators(index_t Q);
  void TestDaggerCommutatorsAlln(index_t Q);

  typedef void commutator_func (const Operator&,const Operator&,Operator&) ;
//  bool Test_against_ref_impl(const Operator& X, const Operator& Y, void (*ComRef)(const Operator&,const Operator&), void (*ComOpt)(const Operator&,const Operator&), std::string output_tag="" );
  bool Test_against_ref_impl(const Operator& X, const Operator& Y, commutator_func ComOpt, commutator_func ComRef, std::string output_tag="" );

  bool Mscheme_Test_comm110ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm220ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm111ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm121ss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm221ss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm122ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_pp_hhss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_phss( const Operator& X, const Operator& Y ) ;

  bool Mscheme_Test_comm121st( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm221st( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm222_pp_hhst( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm222_phst( const Operator& X, const Operator& Y ) ;

  bool Test_comm110ss( const Operator& X, const Operator& Y );
  bool Test_comm220ss( const Operator& X, const Operator& Y );
  bool Test_comm111ss( const Operator& X, const Operator& Y );
  bool Test_comm121ss( const Operator& X, const Operator& Y ); 
  bool Test_comm221ss( const Operator& X, const Operator& Y ); 
  bool Test_comm122ss( const Operator& X, const Operator& Y );
  bool Test_comm222_pp_hhss( const Operator& X, const Operator& Y );
  bool Test_comm222_phss( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_pp_hh_221ss( const Operator& X, const Operator& Y );


  // Tensor
  bool Test_comm111st( const Operator& X, const Operator& Y );
  bool Test_comm121st( const Operator& X, const Operator& Y );
  bool Test_comm122st( const Operator& X, const Operator& Y );
  bool Test_comm221st( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_pp_hhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm222_phst( const Operator& X, const Operator& Y ) ;
  // Tensor 3b
  bool Test_comm331st( const Operator& X, const Operator& Y ) ;
  bool Test_comm231st( const Operator& X, const Operator& Y ) ;
  bool Test_comm232st( const Operator& X, const Operator& Y ) ;
  bool Test_comm132st( const Operator& X, const Operator& Y ) ;
  bool Test_comm223st( const Operator& X, const Operator& Y ) ;
  bool Test_comm133st( const Operator& X, const Operator& Y ) ;

  bool Test_comm332_pphhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm332_ppph_hhhpst( const Operator& X, const Operator& Y ) ;
  bool Test_comm233_pp_hhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm233_phst( const Operator& X, const Operator& Y ) ;
  bool Test_comm333_ppp_hhhst( const Operator& X, const Operator& Y ) ;
  bool Test_comm333_pph_hhpst( const Operator& X, const Operator& Y ) ;

  bool Mscheme_Test_comm330ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm331ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm231ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm132ss( const Operator& X, const Operator& Y );
//
  bool Mscheme_Test_comm332_ppph_hhhpss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm332_pphhss( const Operator& X, const Operator& Y ); 
  bool Mscheme_Test_comm133ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm223ss( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm233_pp_hhss( const Operator& X, const Operator& Y );    
  bool Mscheme_Test_comm233_phss( const Operator& X, const Operator& Y );    
  bool Mscheme_Test_comm232ss( const Operator& X, const Operator& Y );

  bool Mscheme_Test_comm333_ppp_hhhss( const Operator& X, const Operator& Y );  
  bool Mscheme_Test_comm333_pph_hhpss( const Operator& X, const Operator& Y );  

  // scalar-tensor commutator with 3b
  bool Mscheme_Test_comm331st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm223st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm231st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm232st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm133st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm132st( const Operator& X, const Operator& Y );
  bool Mscheme_Test_comm332_ppph_hhhpst(const Operator& X, const Operator& Y); 
  bool Mscheme_Test_comm332_pphhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm233_pp_hhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm233_phst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm333_ppp_hhhst(const Operator &X, const Operator &Y);
  bool Mscheme_Test_comm333_pph_hhpst(const Operator &X, const Operator &Y); 


  bool Test_comm330ss( const Operator& X, const Operator& Y );
  bool Test_comm331ss( const Operator& X, const Operator& Y );
  bool Test_comm231ss( const Operator& X, const Operator& Y );

  bool Test_comm132ss( const Operator& X, const Operator& Y );
  bool Test_comm232ss( const Operator& X, const Operator& Y );
  bool Test_comm332_ppph_hhhpss( const Operator& X, const Operator& Y ); 
  bool Test_comm332_pphhss( const Operator& X, const Operator& Y ); 

  bool Test_comm223ss( const Operator& X, const Operator& Y );
  bool Test_comm133ss( const Operator& X, const Operator& Y );

  bool Test_comm233_pp_hhss( const Operator& X, const Operator& Y );    
  bool Test_comm233_phss( const Operator& X, const Operator& Y );      
  bool Test_comm333_ppp_hhhss( const Operator& X, const Operator& Y );  
  bool Test_comm333_pph_hhpss( const Operator& X, const Operator& Y );  


  bool Test_comm211sd(        const Operator& X, const Operator& Y   );
  bool Test_comm231sd(        const Operator& X, const Operator& Y   );
  bool Test_comm431sd(        const Operator& X, const Operator& Y   );
  bool Test_comm233sd(        const Operator& X, const Operator& Y   );
  bool Test_comm413sd(        const Operator& X, const Operator& Y   );
  bool Test_comm433_pp_hh_sd( const Operator& X, const Operator& Y   );
  bool Test_comm433sd_ph(     const Operator& X, const Operator& Y   );

  bool TestRPAEffectiveCharge( const Operator& H, const Operator& OpIn, size_t k, size_t l);

//  bool TestFactorizedDoubleCommutators(ModelSpace& ms);
  bool TestFactorizedDoubleCommutators();

  bool TestPerturbativeTriples();

  bool SanityCheck();

};


#endif
