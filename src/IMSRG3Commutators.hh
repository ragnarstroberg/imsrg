///////////////////////////////////////////////////////////////////////////////////
//    IMSRGS3Commutators.hh, part of  imsrg++
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

#ifndef IMSRG3Commutators_hh
#define IMSRG3Commutators_hh 1

#include "Operator.hh"
#include <array>

namespace Commutator
{

    extern bool use_imsrg3;
    extern bool use_imsrg3_n7;
    extern bool use_imsrg3_mp4;

    extern bool perturbative_triples;

    extern bool imsrg3_no_qqq;
    extern bool imsrg3_only_vvv;
    extern bool imsrg3_valence_2b;
    extern bool discard_0b_from_3b;
    extern bool discard_1b_from_3b;
    extern bool discard_2b_from_3b;
//    extern bool imsrg3_verbose; // Switch to print out more info for debugging/profiling.
    extern double threebody_threshold;
    extern double imsrg3_dE6max;

    void Discard0bFrom3b(bool tf);
    void Discard1bFrom3b(bool tf);
    void Discard2bFrom3b(bool tf);

    void SetIMSRG3Noqqq(bool tf);
    void SetIMSRG3Onlyvvv(bool tf);
    void SetIMSRG3valence2b(bool tf);
//    void SetIMSRG3Verbose(bool tf);
    void SetSingleThread(bool tf); 
  
  
    // IMSRG(3) commutators. Still a work in progress...
    void comm330ss(const Operator &X, const Operator &Y, Operator &Z);      // implemented and tested.
    void comm331ss(const Operator &X, const Operator &Y, Operator &Z);      // implemented and tested.
    void comm231ss(const Operator &X, const Operator &Y, Operator &Z);      // implemented and tested.

    void comm132ss(const Operator &X, const Operator &Y, Operator &Z); // implemented and tested.
    size_t Hash_comm232_key(std::array<size_t, 5> &kljJJ);
    void comm232ss(const Operator &X, const Operator &Y, Operator &Z);                   // implemented and tested.
    void comm232ss_srs_optimized(const Operator &X, const Operator &Y, Operator &Z);     // implemented and tested.
    void comm232ss_srs_optimized_old(const Operator &X, const Operator &Y, Operator &Z); // implemented and tested.
    void comm232ss_mh_optimized(const Operator &X, const Operator &Y, Operator &Z);      // implemented and tested.
    void comm232ss_expand_full(const Operator &X, const Operator &Y, Operator &Z);       // implemented and tested.
    void comm232ss_expand_reduced(const Operator &X, const Operator &Y, Operator &Z);    // implemented and tested.
    void comm232ss_debug(const Operator &X, const Operator &Y, Operator &Z);             // implemented and tested.
    void comm232ss_slow(const Operator &X, const Operator &Y, Operator &Z);              // implemented and tested.
    void comm332_ppph_hhhpss(const Operator &X, const Operator &Y, Operator &Z);         // implemented and tested.
                                                                                         //  void comm332_pphhss( const Operator& X, const Operator& Y, Operator& Z ) ;      // implemented and tested.
    void comm332_pphhss(const Operator &X, const Operator &Y, Operator &Z);              // implemented and tested.
    void comm332_pphhss_debug(const Operator &X, const Operator &Y, Operator &Z);        // implemented and tested.

    void comm133ss(const Operator &X, const Operator &Y, Operator &Z);             // implemented and tested.
    void comm223ss(const Operator &X, const Operator &Y, Operator &Z);             // implemented and tested.
    void comm223ss_new(const Operator &X, const Operator &Y, Operator &Z);         // implemented and tested.
    void comm223ss_debug(const Operator &X, const Operator &Y, Operator &Z);       // implemented and tested.
    void comm233_pp_hhss(const Operator &X, const Operator &Y, Operator &Z);       // implemented and tested.
    void comm233_pp_hhss_debug(const Operator &X, const Operator &Y, Operator &Z); // implemented and tested.
    void comm233_phss(const Operator &X, const Operator &Y, Operator &Z);          // implemented and tested.
    void comm233_phss_debug(const Operator &X, const Operator &Y, Operator &Z);    // implemented and tested.

    void comm333_ppp_hhhss(const Operator &X, const Operator &Y, Operator &Z); // implemented. andtested.
    void comm333_pph_hhpss(const Operator &X, const Operator &Y, Operator &Z); // implemented. and tested.
    void comm333_pph_hhpss_debug(const Operator &X, const Operator &Y, Operator &Z);


}

#endif
