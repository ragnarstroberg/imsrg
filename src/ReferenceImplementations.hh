
#ifndef ReferenceImplementations_hh
#define ReferenceImplementations_hh

#include "Operator.hh"


namespace ReferenceImplementations
{


  void comm110ss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void comm220ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm111ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm121ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm221ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm122ss( const Operator& X, const Operator& Y, Operator& Z ) ;
  void comm222_pp_hhss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_phss( const Operator& X, const Operator& Y, Operator& Z) ;
  void comm222_pp_hh_221ss( const Operator& X, const Operator& Y, Operator& Z) ;

  void comm330ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm331ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm231ss( const Operator& X, const Operator& Y, Operator& Z ) ;           

  void comm132ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm232ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm332_ppph_hhhpss( const Operator& X, const Operator& Y, Operator& Z ) ; 
  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z ) ;      
  
  void comm133ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm223ss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm233_pp_hhss( const Operator& X, const Operator& Y, Operator& Z ) ;     
  void comm233_phss( const Operator& X, const Operator& Y, Operator& Z ) ;        
  void comm233_phss_debug( const Operator& X, const Operator& Y, Operator& Z ) ;        

  void comm333_ppp_hhhss( const Operator& X, const Operator& Y, Operator& Z ) ;           
  void comm333_pph_hhpss( const Operator& X, const Operator& Y, Operator& Z ) ;          



}// namespace ReferenceImplementations

#endif
