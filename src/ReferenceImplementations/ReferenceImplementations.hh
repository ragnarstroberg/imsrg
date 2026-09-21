
///////////////////////////////////////////////////////////////////////////////////
//    ReferenceImplementations.hh, part of  imsrg++
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

#ifndef ReferenceImplementations_hh
#define ReferenceImplementations_hh 1

#include "Operator.hh"

namespace ReferenceImplementations
{

  void comm110ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm220ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm111ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm121ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm221ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm122ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_phss(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hh_221ss(const Operator &X, const Operator &Y, Operator &Z);

  void comm330ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm331ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm231ss(const Operator &X, const Operator &Y, Operator &Z);

  void comm132ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm232ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm332_ppph_hhhpss(const Operator &X, const Operator &Y, Operator &Z);
  void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z);

  void comm133ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm223ss(const Operator &X, const Operator &Y, Operator &Z);
  void comm233_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm233_phss(const Operator &X, const Operator &Y, Operator &Z);

  void comm333_ppp_hhhss(const Operator &X, const Operator &Y, Operator &Z);
  void comm333_pph_hhpss(const Operator &X, const Operator &Y, Operator &Z);

  // scalar-tensor commutators
  void comm111st(const Operator &X, const Operator &Y, Operator &Z);
  void comm121st(const Operator &X, const Operator &Y, Operator &Z);
  void comm122st(const Operator &X, const Operator &Y, Operator &Z);
  void comm221st(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);
  void comm222_phst(const Operator &X, const Operator &Y, Operator &Z);


  // scalar-tensor with a 3b operator
  void comm331st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm223st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm231st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm232st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm133st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  void comm132st(const Operator &X, const Operator &Y, Operator &Z);            // PASS the unit test (J and T)
  
  void comm332_ppph_hhhpst(const Operator &X, const Operator &Y, Operator &Z);  // PASS the unit test (J and T)
  void comm332_pphhst(const Operator &X, const Operator &Y, Operator &Z);       // PASS the unit test (J and T)
  void comm233_pp_hhst(const Operator &X, const Operator &Y, Operator &Z);      // PASS the unit test (J and T)
  void comm233_phst(const Operator &X, const Operator &Y, Operator &Z);         // PASS the unit test (J and T)
  void comm333_ppp_hhhst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test (J and T)
  void comm333_pph_hhpst(const Operator &X, const Operator &Y, Operator &Z);    // PASS the unit test (J and T)


  /// Two-nested-commutator expressions Z = [X,[X,Y]_3]  where X and Y are 2-body.
  void diagram_CIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIb(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIb(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIc(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIId(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_CIIIb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_DIb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIVa(const Operator &X, const Operator &Y, Operator &Z);
  void diagram_DIVb(const Operator &X, const Operator &Y, Operator &Z);

  void diagram_DIVb_intermediate(const Operator &X, const Operator &Y, Operator &Z);

  // The commutators for [Omega, [Omega, Gamma]]
  void comm223_231_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232_BruteForce(const Operator &Eta, const Operator &Gamma, Operator &Z);

  void comm223_231(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_232(const Operator &Eta, const Operator &Gamma, Operator &Z);

  double TriplesGuess(const Operator &Omega, const Operator &H);
  void comm223_231_fI(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_fII(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_fIIIa(const Operator &Eta, const Operator &Gamma, Operator &Z);
  void comm223_231_fIIIb(const Operator &Eta, const Operator &Gamma, Operator &Z);


  double GetDenom(const Operator& H, const std::vector<index_t>& holes, const std::vector<index_t>& particles);
  double GetMP4_term( const Operator& H, int diagram);
  double GetMP4_F1( const Operator& H) ;
  double GetMP4_F2( const Operator& H) ;
  double GetMP4_F3( const Operator& H) ;
  double GetMP4_F4( const Operator& H) ;
  double GetMP4_F5( const Operator& H) ;
  double GetMP4_F6( const Operator& H) ;
  double GetMP4_F7( const Operator& H) ;
  double GetMP4_F8( const Operator& H) ;
  double GetMP4_F9( const Operator& H) ;
  double GetMP4_F10( const Operator& H) ;
  double GetMP4_F11( const Operator& H) ;
  double GetMP4_F12( const Operator& H) ;
  double GetMP4_F13( const Operator& H) ;
  double GetMP4_F14( const Operator& H) ;
  double GetMP4_F15( const Operator& H) ;
  double GetMP4_F16( const Operator& H) ;
  double GetMP4_F17( const Operator& H) ;
  double GetMP4_F18( const Operator& H) ;
  double GetMP4_F19( const Operator& H) ;
  double GetMP4_F20( const Operator& H) ;
  double GetMP4_F21( const Operator& H) ;
  double GetMP4_F22( const Operator& H) ;
  double GetMP4_F23( const Operator& H) ;
  double GetMP4_F24( const Operator& H) ;
  double GetMP4_F25( const Operator& H) ;
  double GetMP4_F26( const Operator& H) ;
  double GetMP4_F27( const Operator& H) ;
  double GetMP4_F28( const Operator& H) ;
  double GetMP4_F29( const Operator& H) ;
  double GetMP4_F30( const Operator& H) ;
  double GetMP4_F31( const Operator& H) ;
  double GetMP4_F32( const Operator& H) ;
  double GetMP4_F33( const Operator& H) ;
  double GetMP4_F34( const Operator& H) ;
  double GetMP4_F35( const Operator& H) ;
  double GetMP4_F36( const Operator& H) ;
  double GetMP4_F37( const Operator& H) ;
  double GetMP4_F38( const Operator& H) ;
  double GetMP4_F39( const Operator& H) ;


} // namespace ReferenceImplementations

#endif
