
#include "Commutator.hh"
#include "Commutator232.hh"
#include "ModelSpace.hh"
#include "Operator.hh"
#include "TwoBodyME.hh"
#include "ThreeBodyME.hh"
#include "armadillo"
#include "PhysicalConstants.hh" // for SQRT2
#include "AngMom.hh"
#include <cstddef>
#include <map>
#include <deque>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <omp.h>

/// Commutator expressions for second-quantized operators
namespace Commutator
{

  bool use_imsrg3 = false;
  bool use_imsrg3_n7 = false;
  bool use_imsrg3_mp4 = false;

  bool single_thread = false;
  bool verbose = false;

  std::map<std::string, bool> comm_term_on = {
      {"comm110ss", true},
      {"comm220ss", true},
      {"comm111ss", true},
      {"comm121ss", true},
      {"comm221ss", true},
      {"comm122ss", true},
      {"comm222_pp_hhss", true},
      {"comm222_phss", true},
      //      {"comm222_pp_hh_221ss", true},
      {"comm111st", true},
      {"comm121st", true},
      {"comm221st", true},
      {"comm122st", true},
      {"comm222_pp_hhst", true},
      {"comm222_phst", true},
      ///////  by default IMSRG(3) terms are turned off.
      {"comm330ss", false},
      {"comm331ss", false},
      {"comm231ss", false},
      {"comm132ss", false},
      {"comm232ss", false},
      {"comm332_ppph_hhhpss", false},
      {"comm332_pphhss", false},
      {"comm133ss", false},
      {"comm223ss", false},
      {"comm233_pp_hhss", false},
      {"comm233_phss", false},
      {"comm333_ppp_hhhss", false},
      {"comm333_pph_hhpss", false},
      ////////  tensor commutators in IMSRG(3)
      ////////  by default IMSRG(3) terms are turned off.
      {"comm331st", false},
      {"comm231st", false},
      {"comm232st", false},
      {"comm132st", false},
      {"comm223st", false},
      {"comm133st", false},

  };

  void TurnOffTerm(std::string term) { comm_term_on[term] = false; }
  void TurnOnTerm(std::string term) { comm_term_on[term] = true; }

  void PrintSettings()
  {
    std::cout << "use_imsrg3 : " << use_imsrg3 << std::endl;
    std::cout << "use_imsrg3_n7 : " << use_imsrg3_n7 << std::endl;
    std::cout << "use_imsrg3_mp4 : " << use_imsrg3_mp4 << std::endl;
    std::cout << "single_thread : " << single_thread << std::endl;
    for (auto &it : comm_term_on)
    {
      std::cout << it.first << " : " << it.second << std::endl;
    }
  }

  void SetUseIMSRG3(bool tf)
  {
    for (std::string term : {
             "comm330ss", "comm331ss", "comm231ss", "comm132ss", "comm232ss",
             "comm332_ppph_hhhpss", "comm332_pphhss", "comm133ss", "comm223ss",
             "comm233_pp_hhss", "comm233_phss", "comm333_ppp_hhhss", "comm333_pph_hhpss"})
    {
      comm_term_on[term] = tf;
    }
    use_imsrg3 = tf;
  }

  void SetUseIMSRG3N7(bool tf)
  {
    use_imsrg3_n7 = tf;
    if (use_imsrg3_n7)
    {
      for (std::string term : {
               "comm332_ppph_hhhpss", "comm332_pphhss",
               "comm233_pp_hhss", "comm233_phss", "comm333_ppp_hhhss", "comm333_pph_hhpss"})
      {
        comm_term_on[term] = false;
      }
      for (std::string term : {
               "comm330ss", "comm331ss", "comm231ss", "comm132ss", "comm232ss",
               "comm133ss", "comm223ss"})
      {
        comm_term_on[term] = true;
      }
    }
  }

  void SetUseIMSRG3N7_Tensor(bool tf)
  {
    use_imsrg3_n7 = tf;
    use_imsrg3 = tf;
    for (std::string term : {
            "comm331st", "comm231st", "comm132st", "comm232st",
            "comm133st", "comm223st"})
    {
      comm_term_on[term] = tf;
    }
  }

  void SetUseIMSRG3_MP4(bool tf)
  {
    if (tf)
    {
      for (std::string term : {
               "comm331ss", "comm231ss", "comm132ss",
               "comm332_ppph_hhhpss", "comm332_pphhss",
               "comm233_pp_hhss", "comm233_phss", "comm333_ppp_hhhss", "comm333_pph_hhpss"})
      {
        comm_term_on[term] = false;
      }
      for (std::string term : {
               "comm330ss", "comm232ss",
               "comm133ss", "comm223ss"})
      {
        comm_term_on[term] = true;
      }
    }
  }

  void SetIMSRG3Noqqq(bool tf)
  {
    imsrg3_no_qqq = tf;
  }

  void SetIMSRG3valence2b(bool tf)
  {
    imsrg3_valence_2b = tf;
  }

  void SetSingleThread(bool tf)
  {
    single_thread = tf;
  }



  void SetVerbose(bool tf)
  {
    verbose = tf;
  }



  // Operator Operator::Commutator( Operator& opright)
  /// Returns \f$ Z = [X,Y] \f$
  Operator Commutator(const Operator &X, const Operator &Y)
  {
    int xrank = X.rank_J;
    int yrank = Y.rank_J;
    int xlegs = X.GetNumberLegs();
    int ylegs = Y.GetNumberLegs();

    X.modelspace->PreCalculateSixJ();

    if (xrank == 0)
    {
      if ((xlegs % 2 == 0) and (ylegs % 2 == 1))
      {
        return CommutatorScalarDagger(X, Y);
      }
      if ((xlegs % 2 == 1) and (ylegs % 2 == 0))
      {
        return -CommutatorScalarDagger(Y, X);
      }
      if ((xlegs % 2 == 1) and (ylegs % 2 == 1))
      {
        std::cout << "TROUBLE!!!! Called commutator with two Dagger operators. This isn't implemented. Dying..." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (yrank == 0)
      {
        if ((not X.IsReduced()) and (not Y.IsReduced()))
        {
          return CommutatorScalarScalar(X, Y); // [S,S]
        }
        else
        {
          if (not X.IsReduced())
          {
            Operator Ynred = Y;
            Ynred.MakeNotReduced();
            Operator Z = CommutatorScalarScalar(X, Ynred);
            Z.MakeReduced();
            return Z;
          }
          if (not Y.IsReduced())
          {
            Operator Xnred = X;
            Xnred.MakeNotReduced();
            Operator Z = CommutatorScalarScalar(Xnred, Y);
            Z.MakeReduced();
            return Z;
          }
          else if (X.GetTRank() == 0 and Y.GetTRank() == 0) // if X and Y are parity-changing, then Z is parity conserving
          {
            Operator Xnred = X;
            Xnred.MakeNotReduced();
            Operator Ynred = Y;
            Ynred.MakeNotReduced();
            Operator Z = CommutatorScalarScalar(Xnred, Ynred);
            return Z;
          }
          else
          {
            std::cout << " TROUBLE IN " << __FILE__ << " line " << __LINE__ << " Calling scalar commutators with two reduced isospin changing operators..." << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
      }
      else
      {
        if ((not X.IsReduced()) and (Y.IsReduced()))
        {
          return CommutatorScalarTensor(X, Y); // [S,T]
        }
        else if (Y.IsReduced())
        {
          Operator Xnred = X;
          Xnred.MakeNotReduced();
          return CommutatorScalarTensor(Xnred, Y);
        }
        else
        {
          std::cout << " TROUBLE IN " << __FILE__ << " line " << __LINE__ << " not sure what to do with this." << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
    else if (yrank == 0)
    {
      return -CommutatorScalarTensor(Y, X); // [T,S]
    }
    else
    {
      std::cout << "In Tensor-Tensor because X.rank_J = " << X.rank_J << "  X.rank_T = " << X.rank_T << "  X.parity = " << X.parity << "   ";
      std::cout << "  Y.rank_J = " << Y.rank_J << "  Y.rank_T = " << Y.rank_T << "  Y.parity = " << Y.parity << std::endl;
      std::cout << " Tensor-Tensor commutator not yet implemented." << std::endl;
    }
    return 0 * Y;
  }

  /// Commutator where \f$ X \f$ and \f$Y\f$ are scalar operators.
  /// Should be called through Commutator()
  Operator CommutatorScalarScalar(const Operator &X, const Operator &Y)
  {
    double t_css = omp_get_wtime();
    int z_Jrank = X.GetJRank() + Y.GetJRank(); // I sure hope this is zero.
    int z_Trank = X.GetTRank() + Y.GetTRank();
    int z_parity = (X.GetParity() + Y.GetParity()) % 2;
    int z_particlerank = std::max(X.GetParticleRank(), Y.GetParticleRank());
    if (use_imsrg3)
      z_particlerank = std::max(z_particlerank, 3);
    ModelSpace &ms = *(Y.GetModelSpace());
    Operator Z(ms, z_Jrank, z_Trank, z_parity, z_particlerank);

    if ((X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()))
      Z.SetAntiHermitian();
    else if ((X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()))
      Z.SetHermitian();
    else
      Z.SetNonHermitian();

    bool save_single_thread = single_thread;

    if (Z.GetParticleRank() > 2)
    {
      Z.ThreeBody.SwitchToPN_and_discard();
    }

    // Here is where we start calling the IMSRG(2) commutator expressions.
    if (comm_term_on["comm110ss"])
      comm110ss(X, Y, Z);
    if (comm_term_on["comm220ss"])
      comm220ss(X, Y, Z);

    if (comm_term_on["comm111ss"])
      comm111ss(X, Y, Z);
    if (comm_term_on["comm121ss"])
      comm121ss(X, Y, Z);
    if (comm_term_on["comm122ss"])
      comm122ss(X, Y, Z);

    // The 222_pp_hh and 221ss terms can share a common intermediate
    // so if we're computing both, we do them together
    if (comm_term_on["comm222_pp_hhss"] and comm_term_on["comm221ss"])
    {
      comm222_pp_hh_221ss(X, Y, Z);
    }
    else // otherwise, just do the one that is turned on.
    {
      if (comm_term_on["comm222_pp_hhss"])
        comm222_pp_hhss(X, Y, Z);
      if (comm_term_on["comm221ss"])
        comm221ss(X, Y, Z);
    }
    //    if (comm_term_on["comm222_pp_hh_221ss"])
    //      comm222_pp_hh_221ss(X, Y, Z);
    if (comm_term_on["comm222_phss"])
      comm222_phss(X, Y, Z);

    if (use_imsrg3 and ((X.Norm() > threebody_threshold) and (Y.Norm() > threebody_threshold)))
    {
      if (Z.modelspace->scalar3b_transform_first_pass)
        SetSingleThread(true);

      // This one is so important we always include it
      // important for suppressing off-diagonal H3
      if (comm_term_on["comm133ss"])
        comm133ss(X, Y, Z); // scales as n^7, but really more like n^6

      // This gets the perturbative energy from the induced 3 body
      if (comm_term_on["comm330ss"])
        comm330ss(X, Y, Z); // scales as n^6

      // This one is essential. If it's not here, then there are no induced 3 body terms
      if (comm_term_on["comm223ss"])
        comm223ss(X, Y, Z); // scales as n^7

      // one of the two most important IMSRG(3) terms
      if (comm_term_on["comm232ss"])
        comm232ss(X, Y, Z); // this is the slowest n^7 term

      // Maybe not so important, but I think relatively cheap
      if (comm_term_on["comm331ss"])
        comm331ss(X, Y, Z); // scales as n^7

      // Demonstrated that this can have some effect
      if (comm_term_on["comm231ss"])
        comm231ss(X, Y, Z); // scales as n^6

      // no demonstrated effect yet, but it's cheap
      if (comm_term_on["comm132ss"])
        comm132ss(X, Y, Z); // scales as n^6

      // Not too bad, though naively n^8
      if (comm_term_on["comm233_pp_hhss"])
        comm233_pp_hhss(X, Y, Z);

      // This one is super slow too. It involves 9js
      // mat mult makes everything better!
      if (comm_term_on["comm233_phss"])
        comm233_phss(X, Y, Z);

      // not too bad, though naively n^8
      if (comm_term_on["comm332_ppph_hhhpss"])
        comm332_ppph_hhhpss(X, Y, Z);

      // naively n^8, but reasonably fast when implemented as a mat mult
      if (comm_term_on["comm332_pphhss"])
        comm332_pphhss(X, Y, Z);

      // naively n^9 but pretty fast as a mat mult
      if (comm_term_on["comm333_ppp_hhhss"])
        comm333_ppp_hhhss(X, Y, Z);

      // This one works, but it's incredibly slow.  naively n^9.
      // Much improvement by going to mat mult
      if (comm_term_on["comm333_pph_hhpss"])
        comm333_pph_hhpss(X, Y, Z);

      X.profiler.counter["N_ScalarCommutators_3b"] += 1;

      // after going through once, we've stored all the 6js (and maybe 9js), so we can run in OMP loops from now on
      X.modelspace->scalar3b_transform_first_pass = false;

      if (discard_2b_from_3b or discard_1b_from_3b or discard_0b_from_3b)
      {
        Operator Zcopy(Z);
        Zcopy.EraseZeroBody();
        Zcopy.EraseOneBody();
        Zcopy.EraseTwoBody();
        Zcopy = Zcopy.DoNormalOrdering(); // Now just have the NO 0,1,2b parts of the 3b piece
        if (not discard_0b_from_3b)
          Zcopy.EraseZeroBody();
        if (not discard_1b_from_3b)
          Zcopy.EraseOneBody();
        if (not discard_2b_from_3b)
          Zcopy.EraseTwoBody();
        Z.EraseThreeBody();
        Zcopy.ThreeBody.SwitchToPN_and_discard();
        Z -= Zcopy; // Remove the NOnB part of the 3b operator.
      }

    } // if imsrg3 and above threshold

    Z.modelspace->scalar_transform_first_pass = false;
    SetSingleThread(save_single_thread);

    X.profiler.timer[__func__] += omp_get_wtime() - t_css;
    X.profiler.counter["N_ScalarCommutators"] += 1;
    return Z;
  }

  /// Commutator \f$[X,Y]\f$ where \f$ X \f$ is a scalar operator and \f$Y\f$ is a tensor operator.
  /// Should be called through Commutator()
  Operator CommutatorScalarTensor(const Operator &X, const Operator &Y)
  {
    X.profiler.counter["N_TensorCommutators"] += 1;
    double t_cst = omp_get_wtime();
    int ZJ = Y.GetJRank();
    int Zparity = (X.GetParity() + Y.GetParity()) % 2;
    int ZTrank = Y.GetTRank();
    int Zpart = Y.GetParticleRank();
    Operator Z(*(Y.modelspace), ZJ, ZTrank, Zparity, Zpart);

    if ((X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()))
      Z.SetAntiHermitian();
    else if ((X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()))
      Z.SetHermitian();
    else
      Z.SetNonHermitian();
    
    if (Z.GetParticleRank() > 2)
    {
      Z.ThreeBody.SwitchToPN_and_discard();
    }


    // If it's the first time we're calling this, then we go single-threaded because there will be some sixj/ninej symbols
    // that we need to compute and store. After the first pass, they're all stored so we can go parallel.
    bool save_single_thread = single_thread;
    if (Z.modelspace->tensor_transform_first_pass[Z.GetJRank() * 4 + X.GetParity() + 2 * Y.GetParity()])
      SetSingleThread(true);

    if (comm_term_on["comm111st"])
      comm111st(X, Y, Z);
    if (comm_term_on["comm121st"])
      comm121st(X, Y, Z);

    if (comm_term_on["comm122st"])
      comm122st(X, Y, Z);
    // The 222_pp_hh and 221 terms can share a common intermediate
    // so if we're computing both, we do them together
    if (comm_term_on["comm222_pp_hhst"] and comm_term_on["comm221st"])
    {
      comm222_pp_hh_221st(X, Y, Z);
    }
    //    else  // otherwise, just do the one that is turned on.
    //    {
    //      if (comm_term_on["comm222_pp_hhst"] )
    //         comm222_pp_hhst(X, Y, Z);
    //      if ( comm_term_on["comm221st"])
    //         comm221st(X, Y, Z);
    //    }
    if (comm_term_on["comm222_phst"])
      comm222_phst(X, Y, Z);

    //    // We shouldn't be here, because Commutator calls CommutatorScalar Scalar if Z had J Rank = 0.
    //    if (use_imsrg3 and X.GetJRank() == 0 and Y.GetJRank() == 0 and Z.GetJRank() == 0)
    //    {
    //      if (Z.GetParticleRank() < 3)
    //      {
    //        Z.ThreeBody.SwitchToPN_and_discard();
    //      }
    //      std::cout << "tensor comm223ss" << std::endl;
    //      comm223ss(X, Y, Z);
    //      std::cout << "tensor comm232ss" << std::endl;
    //      comm232ss_slow(X, Y, Z);
    //      std::cout << "tensor comm231ss" << std::endl;
    //      comm231ss_slow(X, Y, Z);
    //    }

    if (use_imsrg3 and ((X.Norm() > threebody_threshold) and (Y.Norm() > threebody_threshold)))
    {
      if (comm_term_on["comm331st"])
        comm331st(X, Y, Z);
      if (comm_term_on["comm231st"])
        comm231st(X, Y, Z);
      if (comm_term_on["comm132st"])
        comm132st(X, Y, Z);
      if (comm_term_on["comm232st"])
        comm232st(X, Y, Z);
      if (comm_term_on["comm133st"])
        comm133st(X, Y, Z);
      if (comm_term_on["comm223st"])
        comm223st(X, Y, Z);
    } // if imsrg3 and above threshold

    // This is a better place to put this.
    Z.modelspace->tensor_transform_first_pass.at(Z.GetJRank() * 4 + X.GetParity() + 2 * Y.GetParity()) = false;
    SetSingleThread(save_single_thread);

    if (Z.IsHermitian())
      Z.Symmetrize();
    else if (Z.IsAntiHermitian())
      Z.AntiSymmetrize();
    X.profiler.timer["CommutatorScalarTensor"] += omp_get_wtime() - t_cst;
    return Z;
  }

  /// Commutator where \f$ X \f$ is a scalar and \f$Y\f$ is a dagger (i.e. it creates an additional particle).
  /// Should be called through Commutator()
  Operator CommutatorScalarDagger(const Operator &X, const Operator &Y)
  {
    X.profiler.counter["N_DaggerCommutators"] += 1;
    double t_css = omp_get_wtime();
    Operator Z = Y;
    Z.EraseZeroBody();
    Z.EraseOneBody();
    Z.EraseThreeLeg();

    comm211sd(X, Y, Z);
    comm231sd(X, Y, Z);
    comm413_233sd(X, Y, Z);
    comm433_pp_hh_431sd(X, Y, Z);
    comm433sd_ph(X, Y, Z);

    X.profiler.timer["CommutatorScalarDagger"] += omp_get_wtime() - t_css;
    return Z;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Below is the implementation of the commutators in the various channels
  ///////////////////////////////////////////////////////////////////////////////////////////

  //*****************************************************************************************
  //                ____Y    __         ____X
  //          X ___(_)             Y___(_)
  //
  //  [X1,Y1](0) = Sum_ab (2j_a+1) x_ab y_ba  (n_a-n_b)
  //             = Sum_a  (2j_a+1)  (xy-yx)_aa n_a
  //
  // -- AGREES WITH NATHAN'S RESULTS
  /// \f[
  ///  [X_{1)},Y_{(1)}]_{(0)} = \sum_{a} n_a (2j_a+1) \left(X_{(1)}Y_{(1)}-Y_{(1)}X_{(1)}\right)_{aa}
  /// \f]
  void comm110ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    if (X.IsHermitian() and Y.IsHermitian())
      return; // I think this is the case
    if (X.IsAntiHermitian() and Y.IsAntiHermitian())
      return; // I think this is the case
    if (Z.GetJRank() > 0 or Z.GetTRank() > 0 or Z.GetParity() != 0)
      return;

    arma::mat xyyx = X.OneBody * Y.OneBody - Y.OneBody * X.OneBody;
    for (auto &a : Z.modelspace->holes)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Z.ZeroBody += (oa.j2 + 1) * oa.occ * xyyx(a, a);
    }
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //         __Y__       __X__
  //        ()_ _()  -  ()_ _()
  //           X           Y
  //
  //  [ X^(2), Y^(2) ]^(0) = 1/2 Sum_abcd  Sum_J (2J+1) x_abcd y_cdab (n_a n_b nbar_c nbar_d)
  //                       = 1/2 Sum_J (2J+1) Sum_abcd x_abcd y_cdab (n_a n_b nbar_c nbar_d)
  //                       = 1/2 Sum_J (2J+1) Sum_ab  (X*P_pp*Y)_abab  P_hh
  //
  //  -- AGREES WITH NATHAN'S RESULTS (within < 1%)
  /// \f[
  /// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) \tilde{X}_{abcd}^{J} \tilde{Y}_{cdab}^{J}
  /// \f]
  /// may be rewritten as
  /// \f[
  /// [X_{(2)},Y_{(2)}]_{(0)} = 2 \sum_{J} (2J+1) Tr(X_{hh'pp'}^{J} Y_{pp'hh'}^{J})
  /// \f] where we obtain a factor of four from converting two unrestricted sums to restricted sums, i.e. \f$\sum_{ab} \rightarrow \sum_{a\leq b} \f$,
  /// and using the normalized TBME.
  void comm220ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    if (X.GetParticleRank() < 2 or Y.GetParticleRank() < 2)
      return;
    if (X.IsHermitian() and Y.IsHermitian())
      return; // I think this is the case
    if (X.IsAntiHermitian() and Y.IsAntiHermitian())
      return; // I think this is the case
    if (Z.GetJRank() > 0 or Z.GetTRank() > 0 or Z.GetParity() != 0)
      return;

    for (int ch = 0; ch < Y.nChannels; ++ch)
    {
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto hh = tbc.GetKetIndex_hh();
      auto ph = tbc.GetKetIndex_ph();
      auto pp = tbc.GetKetIndex_pp();
      arma::uvec nbar_indices = arma::join_cols(hh, ph);
      nbar_indices = arma::join_cols(nbar_indices, pp);
      if (hh.size() == 0)
        continue;
      auto nn = tbc.Ket_occ_hh;
      arma::vec nbarnbar = arma::join_cols(tbc.Ket_unocc_hh, tbc.Ket_unocc_ph);
      auto &X2 = X.TwoBody.GetMatrix(ch).submat(hh, nbar_indices);
      arma::mat Y2 = Y.TwoBody.GetMatrix(ch).submat(nbar_indices, hh);
      Y2.head_rows(nbarnbar.size()).each_col() %= nbarnbar;
      Z.ZeroBody += 2 * (2 * tbc.J + 1) * arma::sum(arma::diagvec(X2 * Y2) % nn); // This could be made more efficient, but who cares?
    }
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //
  //        |____. Y          |___.X
  //        |        _        |
  //  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
  //        |                 |
  //
  // -- AGREES WITH NATHAN'S RESULTS
  /// \f[
  /// [X_{(1)},Y_{(1)}]_{(1)} = X_{(1)}Y_{(1)} - Y_{(1)}X_{(1)}
  /// \f]
  void comm111ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    Z.OneBody += X.OneBody * Y.OneBody - Y.OneBody * X.OneBody;
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //                                       |
  //      i |              i |             |
  //        |    ___.Y       |__X__        |
  //        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab
  //      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)
  //                                       |
  //---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
  //                                                  * sum_b y_ab x_biaj - yba x_aibj
  //
  //                     (note: I think this should actually be)
  //                                                = sum_ab (n_a nbar_b) sum_J (2J+1)/(2j_i+1)
  //                                                      * y_ab xbiaj - yba x_aibj
  //
  // -- AGREES WITH NATHAN'S RESULTS
  /// Returns \f$ [X_{(1)},Y_{(2)}] - [Y_{(1)},X_{(2)}] \f$, where
  /// \f[
  /// [X_{(1)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1}\sum_{ab} (n_a \bar{n}_b) \sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj})
  /// \f]
  void comm121ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    index_t norbits = Z.modelspace->all_orbits.size();
    int hZ = Z.IsHermitian() ? 1 : -1;

#pragma omp parallel for
    for (index_t indexi = 0; indexi < norbits; ++indexi)
    {
      auto i = indexi;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      index_t jmin = Z.IsNonHermitian() ? 0 : i;
      for (auto j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        if (j < jmin)
          continue; // only calculate upper triangle
        double zij = 0;
        for (auto &a : Z.modelspace->holes) // C++11 syntax
        {
          Orbit &oa = Z.modelspace->GetOrbit(a);

          if (Y.particle_rank > 1)
          {
            for (auto b : X.GetOneBodyChannel(oa.l, oa.j2, oa.tz2))
            {

              Orbit &ob = Z.modelspace->GetOrbit(b);
              double nanb = oa.occ * (1 - ob.occ);
              if (std::abs(nanb) < ModelSpace::OCC_CUT)
                continue;

              double ybiaj = Y.TwoBody.GetTBMEmonopole(b, i, a, j);
              double yaibj = Y.TwoBody.GetTBMEmonopole(a, i, b, j);
              zij += (ob.j2 + 1) * nanb * X.OneBody(a, b) * ybiaj;
              zij -= (oa.j2 + 1) * nanb * X.OneBody(b, a) * yaibj;
            }
          }
          if (X.particle_rank > 1)
          {
            for (auto b : Y.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
            {
              // comm211 part
              Orbit &ob = Z.modelspace->GetOrbit(b);
              double nanb = oa.occ * (1 - ob.occ);
              if (std::abs(nanb) < ModelSpace::OCC_CUT)
                continue;
              zij -= (ob.j2 + 1) * nanb * Y.OneBody(a, b) * X.TwoBody.GetTBMEmonopole(b, i, a, j);
              zij += (oa.j2 + 1) * nanb * Y.OneBody(b, a) * X.TwoBody.GetTBMEmonopole(a, i, b, j);
            }
          }
        }
        Z.OneBody(i, j) += zij;
        if (jmin == i and j != i)
          Z.OneBody(j, i) += hZ * zij;
      }
    }

    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //
  //      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1)
  //        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
  //        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
  //        |   (  )   _     |   (  )
  //        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
  //      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
  //
  //
  // -- AGREES WITH NATHAN'S RESULTS
  //   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b)
  // \f[
  // [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{c}
  // \left( \mathcal{P}_{pp} (X \mathcal{P}_{hh} Y^{J}
  // - Y^{J} \mathcal{P}_{hh} X^{J}) \mathcal{P}_{pp}
  //  - \mathcal{P}_{hh} (X^{J} \mathcal{P}_{pp} Y^{J}
  //  -  Y^{J} \mathcal{P}_{pp} X^{J}) \mathcal{P}_{hh} \right)_{cicj}
  // \f]
  /// \f[
  /// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2(2j_i+1)}\sum_{J}(2J+1)\sum_{abc} (\bar{n}_a\bar{n}_bn_c + n_an_b\bar{n}_c)
  ///  (X^{J}_{ciab} Y^{J}_{abcj} - Y^{J}_{ciab}X^{J}_{abcj})
  /// \f]
  /// This may be rewritten as
  /// \f[
  /// [X_{(2)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1} \sum_{c} \sum_{J} (2J+1) \left( n_c \mathcal{M}^{J}_{pp,icjc} + \bar{n}_c\mathcal{M}^{J}_{hh,icjc} \right)
  /// \f]
  /// With the intermediate matrix \f[ \mathcal{M}^{J}_{pp} \equiv \frac{1}{2} (X^{J}\mathcal{P}_{pp} Y^{J} - Y^{J}\mathcal{P}_{pp}X^{J}) \f]
  /// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$
  // void Operator::comm221ss( const Operator& X, const Operator& Y)
  void comm221ss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();

    int hZ = Z.IsHermitian() ? 1 : -1;

    TwoBodyME Mpp(Z.modelspace, Z.GetJRank(), Z.GetTRank(), Z.GetParity());
    TwoBodyME Mhh(Z.modelspace, Z.GetJRank(), Z.GetTRank(), Z.GetParity());
    ConstructScalarMpp_Mhh(X, Y, Z, Mpp, Mhh);

    int norbits = Z.modelspace->all_orbits.size();
    std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
    //   for (int i=0;i<norbits;++i)
    for (int indexi = 0; indexi < norbits; ++indexi)
    {
      //      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      //      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      for (int j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        if (j < jmin)
          continue;
        Orbit &oj = Z.modelspace->GetOrbit(j);
        double zij = 0;
        for (auto &c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double nc = oc.occ;
          double nbarc = 1.0 - nc;
          int Jmin = std::max(std::abs(oc.j2 - oi.j2), std::abs(oc.j2 - oj.j2)) / 2;
          int Jmax = (oc.j2 + std::min(oi.j2, oj.j2)) / 2;

          if (std::abs(nc) > 1e-9)
          {
            for (int J = Jmin; J <= Jmax; J++)
            {
              zij += (2 * J + 1) * nc * Mpp.GetTBME_J(J, c, i, c, j);
            }
          }
          if (std::abs(nbarc) > 1e-9)
          {
            for (int J = Jmin; J <= Jmax; J++)
            {
              zij += (2 * J + 1) * nbarc * Mhh.GetTBME_J(J, c, i, c, j);
            }
          }
        }

        Z.OneBody(i, j) += zij / (oi.j2 + 1.0);
        if (jmin == i and i != j)
          Z.OneBody(j, i) += hZ * zij / (oi.j2 + 1.0);
      } // for j
    }

    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //
  //    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
  //    |     |___.Y          |__X___|
  //    |     |         _     |      |
  //    |_____|               |      |_____.Y
  //    |  X  |               |      |
  //
  // -- AGREES WITH NATHAN'S RESULTS
  /// Returns \f$ [X_{(1)},Y_{(2)}]_{(2)} - [Y_{(1)},X_{(2)}]_{(2)} \f$, where
  /// \f[
  /// [X_{(1)},Y_{(2)}]^{J}_{ijkl} = \sum_{a} ( X_{ia}Y^{J}_{ajkl} + X_{ja}Y^{J}_{iakl} - X_{ak} Y^{J}_{ijal} - X_{al} Y^{J}_{ijka} )
  /// \f]
  /// here, all TBME are unnormalized, i.e. they should have a tilde.
  // This is still too slow...
  // TODO: MODIFY THIS TO ACCOMMODATE J=0, T!=0 operators
  void comm122ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    if (not(X.GetParity() == 0 and Y.GetParity() == 0 and Z.GetParity() == 0 and X.GetTRank() == 0 and Y.GetTRank() == 0 and Z.GetTRank() == 0))
    {
      comm122ss_slower(X, Y, Z);
      return;
    }
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    int hZ = Z.IsHermitian() ? 1 : -1;

    int n_nonzero = Z.modelspace->SortedTwoBodyChannels.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < n_nonzero; ++ich)
    {
      int ch = Z.modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      auto &X2 = X.TwoBody.GetMatrix(ch, ch);
      auto &Y2 = Y.TwoBody.GetMatrix(ch, ch);
      auto &Z2 = Z.TwoBody.GetMatrix(ch, ch);
      arma::mat W2(size(Z2), arma::fill::zeros); // temporary intermediate matrix

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0; indx_ij < npq; ++indx_ij)
      {
        Ket &bra = tbc.GetKet(indx_ij);
        int i = bra.p;
        int j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);
        int flipphaseij = -Z.modelspace->phase((oi.j2 + oj.j2) / 2 - tbc.J);

        // make lists of the indices we want, then do matrix multiplication.
        // there may be a more efficient way to find these
        std::vector<index_t> ind1_ia, ind1_ja, ind2_aj, ind2_ai;
        std::vector<double> factor_ia, factor_ja;
        //         for (int a : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
        for (int a : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
        {
          size_t ind2 = tbc.GetLocalIndex(std::min(a, j), std::max(a, j));
          if (ind2 < 0 or ind2 >= tbc.GetNumberKets())
            continue;
          ind1_ia.push_back(a);
          ind2_aj.push_back(ind2);
          factor_ia.push_back(a > j ? flipphaseij : (a == j ? PhysConst::SQRT2 : 1));
        }
        if (i != j)
        {
          //           for (int a : Z.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
          for (int a : Z.GetOneBodyChannel(oj.l, oj.j2, oj.tz2))
          {
            size_t ind2 = tbc.GetLocalIndex(std::min(a, i), std::max(a, i));
            if (ind2 < 0 or ind2 >= tbc.GetNumberKets())
              continue;
            ind1_ja.push_back(a);
            ind2_ai.push_back(ind2);
            factor_ja.push_back(i > a ? flipphaseij : (i == a ? PhysConst::SQRT2 : 1));
          }
        }

        arma::uvec u_ind1_ia(ind1_ia);
        arma::uvec u_ind1_ja(ind1_ja);
        arma::uvec u_ind2_aj(ind2_aj);
        arma::uvec u_ind2_ai(ind2_ai);
        arma::vec v_factor_ia(factor_ia);
        arma::vec v_factor_ja(factor_ja);
        if (i == j)
        {
          v_factor_ia /= PhysConst::SQRT2;
          v_factor_ja /= PhysConst::SQRT2;
        }

        // This is fairly obfuscated, but hopefully faster for bigger calculations
        // SRS THERE MUST BE A BETTER WAY...
        if (X.particle_rank > 1 and Y.particle_rank > 1)
        {
          W2.col(indx_ij) = join_horiz(Y2.cols(join_vert(u_ind2_aj, u_ind2_ai)), X2.cols(join_vert(u_ind2_aj, u_ind2_ai))) * join_vert(join_vert(X1.unsafe_col(i).rows(u_ind1_ia) % v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja) % v_factor_ja),
                                                                                                                                       -join_vert(Y1.unsafe_col(i).rows(u_ind1_ia) % v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja) % v_factor_ja));
        }
        else if (X.particle_rank < 2 and Y.particle_rank > 1)
        {
          W2.col(indx_ij) = Y2.cols(join_vert(u_ind2_aj, u_ind2_ai)) * join_vert(X1.unsafe_col(i).rows(u_ind1_ia) % v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja) % v_factor_ja);
        }
        else if (X.particle_rank > 1 and Y.particle_rank < 2)
        {
          W2.col(indx_ij) = -X2.cols(join_vert(u_ind2_aj, u_ind2_ai)) * join_vert(Y1.unsafe_col(i).rows(u_ind1_ia) % v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja) % v_factor_ja);
        }

        if (i == j)
          W2.col(indx_ij) *= 2;
      }
      Z2 -= W2 + hZ * W2.t();
    }
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  // A more straightforward implementation, which can handle parity and isospin changing operators
  // It's not dramatically slower, but slow enough that we should only use it when we need it.
  void comm122ss_slower(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X1 = X.OneBody;
    auto &Y1 = Y.OneBody;
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z2 = Z.TwoBody;

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = *(bra.op);
        Orbit &oj = *(bra.oq);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;

          Orbit &ok = *(ket.op);
          Orbit &ol = *(ket.oq);

          double zijkl = 0;
          // X1 Y2
          for (size_t a : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
          {
            zijkl += X1(i, a) * Y2.GetTBME_J(J, J, a, j, k, l);
          }
          for (size_t a : X.OneBodyChannels.at({oj.l, oj.j2, oj.tz2}))
          {
            zijkl += X1(j, a) * Y2.GetTBME_J(J, J, i, a, k, l);
          }
          // X2 Y1
          for (size_t a : Y.OneBodyChannels.at({ok.l, ok.j2, ok.tz2}))
          {
            zijkl += X2.GetTBME_J(J, J, i, j, a, l) * Y1(a, k);
          }
          for (size_t a : Y.OneBodyChannels.at({ol.l, ol.j2, ol.tz2}))
          {
            zijkl += X2.GetTBME_J(J, J, i, j, k, a) * Y1(a, l);
          }

          // Y1 X2
          for (size_t a : Y.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
          {
            zijkl -= Y1(i, a) * X2.GetTBME_J(J, J, a, j, k, l);
          }
          for (size_t a : Y.OneBodyChannels.at({oj.l, oj.j2, oj.tz2}))
          {
            zijkl -= Y1(j, a) * X2.GetTBME_J(J, J, i, a, k, l);
          }
          // Y2 X1
          for (size_t a : X.OneBodyChannels.at({ok.l, ok.j2, ok.tz2}))
          {
            zijkl -= Y2.GetTBME_J(J, J, i, j, a, l) * X1(a, k);
          }
          for (size_t a : X.OneBodyChannels.at({ol.l, ol.j2, ol.tz2}))
          {
            zijkl -= Y2.GetTBME_J(J, J, i, j, k, a) * X1(a, l);
          }

          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //*****************************************************************************************
  //
  //  |     |      |     |
  //  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
  //  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
  //  |__X__|      |__Y__|
  //  |     |      |     |
  //
  // -- AGREES WITH NATHAN'S RESULTS
  //   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b)
  /// Calculates the part of the commutator \f$ [X_{(2)},Y_{(2)}]_{(2)} \f$ which involves particle-particle
  /// or hole-hole intermediate states.
  /// \f[
  /// [X_{(2)},Y_{(2)}]^{J}_{ijkl} = \frac{1}{2} \sum_{ab} (\bar{n}_a\bar{n}_b - n_an_b) (X^{J}_{ijab}Y^{J}_{ablk} - Y^{J}_{ijab}X^{J}_{abkl})
  /// \f]
  /// This may be written as
  /// \f[
  /// [X_{(2)},Y_{(2)}]^{J} = \mathcal{M}^{J}_{pp} - \mathcal{M}^{J}_{hh}
  /// \f]
  /// With the intermediate matrices
  /// \f[
  /// \mathcal{M}^{J}_{pp} \equiv \frac{1}{2}(X^{J} \mathcal{P}_{pp} Y^{J} - Y^{J} \mathcal{P}_{pp} X^{J})
  /// \f]
  /// and likewise for \f$ \mathcal{M}^{J}_{hh} \f$.
  // void Operator::comm222_pp_hhss( Operator& opright, Operator& opout )
  // void Operator::comm222_pp_hhss( const Operator& X, const Operator& Y )
  void comm222_pp_hhss(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();

    TwoBodyME Mpp = Z.TwoBody;
    TwoBodyME Mhh = Z.TwoBody;
    Mpp.Erase();
    Mhh.Erase();

    ConstructScalarMpp_Mhh(X, Y, Z, Mpp, Mhh);

    Z.TwoBody += Mpp;
    Z.TwoBody -= Mhh;
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  /// Construct the intermediate matrices \f$ \mathcal{M}^J_{pp}\f$ and \f$ \mathcal{M}^{J}_{hh} \f$
  /// for use in Commutator::comm222_pp_hh_221ss().
  ///  \f$ \mathcal{M}^{J}_{pp} \equiv \frac{1}{2} (X^{J} \mathcal{P}_{pp} Y^{J} - Y^{J}\mathcal{P}_{pp} X^{J} ) \f$
  /// where \f$ \mathcal{P}_{pp} \f$ is a projector onto particle-particle states.
  void ConstructScalarMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z, TwoBodyME &Mpp, TwoBodyME &Mhh)
  {

    double t_start = omp_get_wtime();

    int hX = X.IsHermitian() ? +1 : -1;
    int hY = Y.IsHermitian() ? +1 : -1;

    std::vector<size_t> ch_bra_list, ch_ket_list;
    auto ch_iter = Z.TwoBody.MatEl;
    // TODO: We'll need to be more careful about this special case.
    if (Z.GetParticleRank() < 2 and Y.GetParticleRank() > 1)
      ch_iter = Y.TwoBody.MatEl;

    for (auto &iter : ch_iter)
    {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
      // Optimization sugggested by Antoine Belley.
      // The form of Hod will make Eta zero for some channels which would not be zero by symmetry.
      // Check if they're zero and if so skip them.
      if (X.GetTRank() == 0 and X.GetParity() == 0 and Y.GetTRank() == 0 and Y.GetParity() == 0)
      {
        if (X.TwoBody.GetMatrix(ch_bra, ch_ket).is_zero() or Y.TwoBody.GetMatrix(ch_bra, ch_ket).is_zero())
          continue;
      }

      ch_bra_list.push_back(ch_bra);
      ch_ket_list.push_back(ch_ket);
    }
    int nch = ch_bra_list.size();
#ifndef OPENBLAS_NOUSEOMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int ich = 0; ich < nch; ++ich)
    {
      int ch_bra = ch_bra_list[ich];
      int ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

      // Figure out what the intermediate channel should be
      int ch_ab_XY = ch_bra;                           // For Hamiltonian-type X,Y,Z  ij, ab, kl all belong to the same channel.
      int ch_ab_YX = ch_ket;                           // For Hamiltonian-type X,Y,Z  ij, ab, kl all belong to the same channel.
      if ((X.GetTRank() != 0) or (X.GetParity() != 0)) // if X changes parity or Tz, take more care.
      {
        int parity_ab_XY = (tbc_bra.parity + X.GetParity()) % 2;
        int parity_ab_YX = (tbc_ket.parity + X.GetParity()) % 2;
        // If both X and Y are isospin-changing operators, then there could be multiple intermediate channels.
        // I'm not sure why that would come up, so for now let's just throw an error if that's the case.
        int Tz_ab_XY = tbc_bra.Tz; //
        int Tz_ab_YX = tbc_ket.Tz; //
        if (X.GetTRank() != 0)
        {
          Tz_ab_XY = tbc_ket.Tz;
          Tz_ab_YX = tbc_bra.Tz;
          if (Y.GetTRank() != 0)
          {
            std::cout << "=======  Uh Oh. Taking the commutator of two isospin-changing operators."
                      << " That's not yet implemented. How did I get here?  LINE " << __LINE__ << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
        ch_ab_XY = Z.modelspace->GetTwoBodyChannelIndex(tbc_bra.J, parity_ab_XY, Tz_ab_XY);
        ch_ab_YX = Z.modelspace->GetTwoBodyChannelIndex(tbc_bra.J, parity_ab_YX, Tz_ab_YX);
      }

      TwoBodyChannel &tbc_ab_XY = Z.modelspace->GetTwoBodyChannel(ch_ab_XY);
      TwoBodyChannel &tbc_ab_YX = Z.modelspace->GetTwoBodyChannel(ch_ab_YX);

      // If X or Y change parity or isospin, then we need to worry about the fact that we only store
      // ch_bra <= ch_ket. If we need the other ordering we get it by Hermiticity.
      auto &X_ijab = (ch_bra <= ch_ab_XY) ? X.TwoBody.GetMatrix(ch_bra, ch_ab_XY) : X.TwoBody.GetMatrix(ch_ab_XY, ch_bra).t() * hX;
      auto &X_abkl = (ch_ab_YX <= ch_ket) ? X.TwoBody.GetMatrix(ch_ab_YX, ch_ket) : X.TwoBody.GetMatrix(ch_ket, ch_ab_YX).t() * hX;
      auto &Y_ijab = (ch_bra <= ch_ab_YX) ? Y.TwoBody.GetMatrix(ch_bra, ch_ab_YX) : Y.TwoBody.GetMatrix(ch_ab_YX, ch_bra).t() * hY;
      auto &Y_abkl = (ch_ab_XY <= ch_ket) ? Y.TwoBody.GetMatrix(ch_ab_XY, ch_ket) : Y.TwoBody.GetMatrix(ch_ket, ch_ab_XY).t() * hY;

      auto &Matrixpp = Mpp.GetMatrix(ch_bra, ch_ket);
      auto &Matrixhh = Mhh.GetMatrix(ch_bra, ch_ket);

      auto &bras_pp = tbc_ab_XY.GetKetIndex_pp();
      auto &bras_hh = tbc_ab_XY.GetKetIndex_hh();
      auto &bras_ph = tbc_ab_XY.GetKetIndex_ph();
      auto &kets_pp = tbc_ab_YX.GetKetIndex_pp();
      auto &kets_hh = tbc_ab_YX.GetKetIndex_hh();
      auto &kets_ph = tbc_ab_YX.GetKetIndex_ph();
      auto &nanb_bra = tbc_ab_XY.Ket_occ_hh;
      auto &nanb_ket = tbc_ab_YX.Ket_occ_hh;
      auto &nbarnbar_hh_bra = tbc_ab_XY.Ket_unocc_hh;
      auto &nbarnbar_ph_bra = tbc_ab_XY.Ket_unocc_ph;
      auto &nbarnbar_hh_ket = tbc_ab_YX.Ket_unocc_hh;
      auto &nbarnbar_ph_ket = tbc_ab_YX.Ket_unocc_ph;

      // explicitly checking sizes because for size=0, we get DGEMM errors
      if (bras_pp.size() > 0)
        Matrixpp += X_ijab.cols(bras_pp) * Y_abkl.rows(bras_pp);
      if (bras_hh.size() > 0)
        Matrixhh += X_ijab.cols(bras_hh) * arma::diagmat(nanb_bra) * Y_abkl.rows(bras_hh);
      if (bras_hh.size() > 0)
        Matrixpp += X_ijab.cols(bras_hh) * arma::diagmat(nbarnbar_hh_bra) * Y_abkl.rows(bras_hh);
      if (bras_ph.size() > 0)
        Matrixpp += X_ijab.cols(bras_ph) * arma::diagmat(nbarnbar_ph_bra) * Y_abkl.rows(bras_ph);

      if (Z.IsHermitian() and ch_bra == ch_ket)
      {
        Matrixpp += Matrixpp.t();
        Matrixhh += Matrixhh.t();
      }
      else if (Z.IsAntiHermitian() and ch_bra == ch_ket) // i.e. LHS and RHS are both hermitian or ant-hermitian
      {
        Matrixpp -= Matrixpp.t();
        Matrixhh -= Matrixhh.t();
      }
      else
      {
        if (kets_pp.size() > 0)
          Matrixpp -= Y_ijab.cols(kets_pp) * X_abkl.rows(kets_pp);
        if (kets_hh.size() > 0)
          Matrixhh -= Y_ijab.cols(kets_hh) * arma::diagmat(nanb_ket) * X_abkl.rows(kets_hh);
        if (kets_hh.size() > 0)
          Matrixpp -= Y_ijab.cols(kets_hh) * arma::diagmat(nbarnbar_hh_ket) * X_abkl.rows(kets_hh);
        if (kets_ph.size() > 0)
          Matrixpp -= Y_ijab.cols(kets_ph) * arma::diagmat(nbarnbar_ph_ket) * X_abkl.rows(kets_ph);
      }

    } // for ch

    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  /// Since comm222_pp_hhss() and comm221ss() both require the construction of
  /// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
  /// only calculate the intermediates once.
  void comm222_pp_hh_221ss(const Operator &X, const Operator &Y, Operator &Z)
  {

    double t_start = omp_get_wtime();
    int hZ = Z.IsHermitian() ? 1 : -1;
    if (X.GetParticleRank() < 2 or Y.GetParticleRank() < 2)
      return;

    TwoBodyME Mpp = Z.TwoBody;
    TwoBodyME Mhh = Z.TwoBody;
    Mpp.Erase();
    Mhh.Erase();
    double t_internal = omp_get_wtime();

    ConstructScalarMpp_Mhh(X, Y, Z, Mpp, Mhh);

    if (Commutator::verbose)
    {
       X.profiler.timer["_ConstructScalarMpp_Mhh"] += omp_get_wtime() - t_internal;
       t_internal = omp_get_wtime();
    }

    Z.TwoBody += Mpp;
    Z.TwoBody -= Mhh;

    if (Commutator::verbose)
    {
       X.profiler.timer["_pphh TwoBody bit"] += omp_get_wtime() - t_internal;
       t_internal = omp_get_wtime();
    }

    //   int norbits = Z.modelspace->GetNumberOrbits();
    int norbits = Z.modelspace->all_orbits.size();
    // The one body part
    std::vector<index_t> allorb_vec(Z.modelspace->all_orbits.begin(), Z.modelspace->all_orbits.end());
#pragma omp parallel for schedule(dynamic, 1)
    //   for (int i=0;i<norbits;++i)
    for (int indexi = 0; indexi < norbits; ++indexi)
    {
      //      auto i = Z.modelspace->all_orbits[indexi];
      auto i = allorb_vec[indexi];
      Orbit &oi = Z.modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      //      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      for (int j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2))
      {
        if (j < jmin)
          continue;

        Orbit &oj = Z.modelspace->GetOrbit(j);
        double zij = 0;
        for (auto &c : Z.modelspace->all_orbits)
        {
          Orbit &oc = Z.modelspace->GetOrbit(c);
          double nc = oc.occ;
          double nbarc = 1.0 - nc;
          int Jmin = std::max(std::abs(oc.j2 - oi.j2), std::abs(oc.j2 - oj.j2)) / 2;
          int Jmax = (oc.j2 + std::min(oi.j2, oj.j2)) / 2;

          if (std::abs(nc) > 1e-9)
          {
            for (int J = Jmin; J <= Jmax; J++)
            {
              zij += (2 * J + 1) * nc * Mpp.GetTBME_J(J, c, i, c, j);
            }
          }
          if (std::abs(nbarc) > 1e-9)
          {
            for (int J = Jmin; J <= Jmax; J++)
            {
              zij += (2 * J + 1) * nbarc * Mhh.GetTBME_J(J, c, i, c, j);
            }
          }
        }

        Z.OneBody(i, j) += zij / (oi.j2 + 1.0);
        if (jmin == i and i != j)
          Z.OneBody(j, i) += hZ * zij / (oi.j2 + 1.0);
      } // for j
    } // for i

    if (Commutator::verbose)
    {
       X.profiler.timer["_pphh One Body bit"] += omp_get_wtime() - t_internal;
    }
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  }

  //**************************************************************************
  //
  //  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
  //                        { k l J'}
  // SCALAR VARIETY
  /// The scalar Pandya transformation is defined as
  /// \f[
  ///  \bar{X}^{J}_{i\bar{j}k\bar{l}} = - \sum_{J'} (2J'+1)
  ///  \left\{ \begin{array}{lll}
  ///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
  ///  X^{J}_{ilkj}
  /// \f]
  /// where the overbar indicates time-reversed orbits.
  // void Operator::DoPandyaTransformation_SingleChannel(arma::mat& TwoBody_CC_ph, int ch_cc, std::string orientation="normal") const
  void DoPandyaTransformation_SingleChannel(const Operator &Z, arma::mat &TwoBody_CC_ph, int ch_cc, std::string orientation = "normal")
  {
    int herm = Z.IsHermitian() ? 1 : -1;
    TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph());
    int nph_kets = kets_ph.n_rows;
    int J_cc = tbc_cc.J;

    if (orientation == "normal")
      TwoBody_CC_ph.zeros(2 * nph_kets, nKets_cc);
    else if (orientation == "transpose")
      TwoBody_CC_ph.zeros(nKets_cc, 2 * nph_kets);
    else
    {
      std::cout << __PRETTY_FUNCTION__ << " =>  Unknown orientation input  " << orientation << ". Don't know what to do with this." << std::endl;
      return;
    }

    // loop over cross-coupled ph bras <ab| in this channel
    // (this is the side that gets summed over in the matrix multiplication)
    for (int ibra = 0; ibra < nph_kets; ++ibra)
    {
      Ket &bra_cc = tbc_cc.GetKet(kets_ph[ibra]);
      // we want to evaluate a<=b and a>=b, so to avoid code duplication, we turn this into a loop over the two orderings
      std::vector<size_t> ab_switcheroo = {bra_cc.p, bra_cc.q};
      for (int ab_case = 0; ab_case <= 1; ab_case++)
      {
        int a = ab_switcheroo[ab_case]; // this little bit gives us a,b if ab_case=0 and b,a if ab_case=1
        int b = ab_switcheroo[1 - ab_case];
        size_t bra_shift = ab_case * nph_kets; // if we switch a<->b, we offset the bra index by nph_kets

        Orbit &oa = Z.modelspace->GetOrbit(a);
        Orbit &ob = Z.modelspace->GetOrbit(b);
        double ja = oa.j2 * 0.5;
        double jb = ob.j2 * 0.5;
        double na_nb_factor = oa.occ - ob.occ;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc = 0; iket_cc < nKets_cc; ++iket_cc)
        {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc % nKets_cc);
          int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
          int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
          Orbit &oc = Z.modelspace->GetOrbit(c);
          Orbit &od = Z.modelspace->GetOrbit(d);
          double jc = oc.j2 * 0.5;
          double jd = od.j2 * 0.5;

          int jmin = std::max(std::abs(ja - jd), std::abs(jc - jb));
          int jmax = std::min(ja + jd, jc + jb);
          double Xbar = 0;
          for (int J_std = jmin; J_std <= jmax; ++J_std)
          {
            double sixj = Z.modelspace->GetSixJ(ja, jb, J_cc, jc, jd, J_std);
            if (std::abs(sixj) < 1e-8)
              continue;
            double tbme = Z.TwoBody.GetTBME_J(J_std, a, d, c, b);
            Xbar -= (2 * J_std + 1) * sixj * tbme;
          }
          if (orientation == "normal")
          {
            TwoBody_CC_ph(ibra + bra_shift, iket_cc) = Xbar;
          }
          else // "transpose"
          {
            TwoBody_CC_ph(iket_cc, ibra + bra_shift) = herm * Xbar * na_nb_factor; // we slap the (na-nb) on the transposed one.
          }
        }
      }
    }
  }

  void DoPandyaTransformation_SingleChannel_XandY(const Operator &X, const Operator &Y, arma::mat &X2_CC_ph, arma::mat &Y2_CC_ph, int ch_cc)
  {
    //   int hX = X.IsHermitian() ? 1 : -1;
    int hY = Y.IsHermitian() ? 1 : -1;
    TwoBodyChannel_CC &tbc_cc = X.modelspace->GetTwoBodyChannel_CC(ch_cc);
    int nKets_cc = tbc_cc.GetNumberKets();
    arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph());
    int nph_kets = kets_ph.n_rows;
    int J_cc = tbc_cc.J;

    X2_CC_ph.zeros(nKets_cc, 2 * nph_kets);
    Y2_CC_ph.zeros(2 * nph_kets, nKets_cc);

    // TODO: For a set of one-body channels a,b the J coupling stuff should all be the same, so
    //       there is probably some gain to be had by looping over radial quantum numbers na,nb
    //       and nc,nd  inside the Jph loop.

    // loop over cross-coupled ph bras <ab| in this channel
    // (this is the side that gets summed over in the matrix multiplication)
    for (int ibra = 0; ibra < nph_kets; ++ibra)
    {
      Ket &bra_cc = tbc_cc.GetKet(kets_ph[ibra]);
      // we want to evaluate a<=b and a>=b, so to avoid code duplication, we turn this into a loop over the two orderings
      std::vector<size_t> ab_switcheroo = {bra_cc.p, bra_cc.q};
      for (int ab_case = 0; ab_case <= 1; ab_case++)
      {
        int a = ab_switcheroo[ab_case]; // this little bit gives us a,b if ab_case=0 and b,a if ab_case=1
        int b = ab_switcheroo[1 - ab_case];
        size_t bra_shift = ab_case * nph_kets; // if we switch a<->b, we offset the bra index by nph_kets

        Orbit &oa = X.modelspace->GetOrbit(a);
        Orbit &ob = X.modelspace->GetOrbit(b);
        double ja = oa.j2 * 0.5;
        double jb = ob.j2 * 0.5;
        int jjai = oa.j2;
        int jjbi = ob.j2;
        double na_nb_factor = oa.occ - ob.occ;

        // loop over cross-coupled kets |cd> in this channel
        for (int iket_cc = 0; iket_cc < nKets_cc; ++iket_cc)
        {
          Ket &ket_cc = tbc_cc.GetKet(iket_cc % nKets_cc);
          int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
          int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
          Orbit &oc = X.modelspace->GetOrbit(c);
          Orbit &od = X.modelspace->GetOrbit(d);

          // Check the isospin projection. If this isn't conserved in the usual channel,
          // then all the xcbad and yadcb will be zero and we don't need to bother computing SixJs.
          if ((std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != 2 * X.GetTRank()) and (std::abs(oa.tz2 + od.tz2 - ob.tz2 - oc.tz2) != 2 * Y.GetTRank()))
            continue;

          double jc = oc.j2 * 0.5;
          double jd = od.j2 * 0.5;
          int jjci = oc.j2;
          int jjdi = od.j2;

          double normfactor = 1.0;
          if (a == d)
            normfactor *= PhysConst::SQRT2;
          if (c == b)
            normfactor *= PhysConst::SQRT2;

          int jmin = std::max(std::abs(ja - jd), std::abs(jc - jb));
          int jmax = std::min(ja + jd, jc + jb);
          double Xbar = 0;
          double Ybar = 0;
          int dJ_std = 1;
          if (a == d or c == b)
          {
            dJ_std = 2;
            jmin += jmin % 2;
          }
          for (int J_std = jmin; J_std <= jmax; J_std += dJ_std)
          {
            double sixj = X.modelspace->GetCachedSixJ(jjai, jjbi, J_cc, jjci, jjdi, J_std);
            if (std::abs(sixj) < 1e-8)
              continue;

            // Since we want the same element of two different operators, we use GetTBME_J_norm_twoOps
            // which does all the phase/index lookup stuff once and then just accesses the two matrices.
            double xcbad = 0;
            double yadcb = 0;
            X.TwoBody.GetTBME_J_norm_twoOps(Y.TwoBody, J_std, J_std, c, b, a, d, xcbad, yadcb);
            Xbar -= (2 * J_std + 1) * sixj * xcbad;
            Ybar -= (2 * J_std + 1) * sixj * yadcb * hY;
          }
          X2_CC_ph(iket_cc, ibra + bra_shift) = Xbar * normfactor * na_nb_factor;
          Y2_CC_ph(ibra + bra_shift, iket_cc) = Ybar * normfactor;

        } // for iket_cc
      } // for ab_case
    } // for ibra
  }

  void DoPandyaTransformation(const Operator &Z, std::deque<arma::mat> &TwoBody_CC_ph, std::string orientation = "normal")
  {
    // loop over cross-coupled channels
    int n_nonzero = Z.modelspace->SortedTwoBodyChannels_CC.size();
#pragma omp parallel for schedule(dynamic, 1) if (not Z.modelspace->scalar_transform_first_pass)
    for (int ich = 0; ich < n_nonzero; ++ich)
    {
      int ch_cc = Z.modelspace->SortedTwoBodyChannels_CC[ich];
      DoPandyaTransformation_SingleChannel(Z, TwoBody_CC_ph[ch_cc], ch_cc, orientation);
    }
  }

  // Take Zbar, and perform a Pandya transform on it and put the result in Z (technically there's a minus sign
  // missing in what is done here, but that's on purpose).
  void AddInversePandyaTransformation(const std::deque<arma::mat> &Zbar, Operator &Z)
  {
    // Do the inverse Pandya transform
    int nch = Z.modelspace->GetNumberTwoBodyChannels();
    int hZ = Z.IsHermitian() ? 1 : -1;

    // Collapse two outer loops into one for better load balancing
    std::vector<int> ch_vec;
    std::vector<int> ibra_vec;
    for (int ch = 0; ch < nch; ++ch)
    {
      int nKets = Z.modelspace->GetTwoBodyChannel(ch).GetNumberKets();
      for (int ibra = 0; ibra < nKets; ++ibra)
      {
        ch_vec.push_back(ch);
        ibra_vec.push_back(ibra);
      }
    }
    size_t nch_and_ibra = ch_vec.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (size_t ichbra = 0; ichbra < nch_and_ibra; ichbra++)
    {
      int ch = ch_vec[ichbra];
      int ibra = ibra_vec[ichbra];

      TwoBodyChannel &tbc = Z.modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();
      auto &ZMat = Z.TwoBody.GetMatrix(ch, ch);

      Ket &bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      Orbit &oi = Z.modelspace->GetOrbit(i);
      Orbit &oj = Z.modelspace->GetOrbit(j);
      double ji = oi.j2 / 2.;
      double jj = oj.j2 / 2.;
      int jji = oi.j2;
      int jjj = oj.j2;
      int ketmin = Z.IsHermitian() ? ibra : ibra + 1;
      for (int iket = ketmin; iket < nKets; ++iket)
      {
        Ket &ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit &ok = Z.modelspace->GetOrbit(k);
        Orbit &ol = Z.modelspace->GetOrbit(l);
        double jk = ok.j2 / 2.;
        double jl = ol.j2 / 2.;
        int jjk = ok.j2;
        int jjl = ol.j2;

        double commij = 0;
        double commji = 0;

        int parity_cc = (oi.l + ol.l) % 2;
        int Tz_cc = std::abs(oi.tz2 - ol.tz2) / 2;
        int Jpmin = std::max(std::abs(int(ji - jl)), std::abs(int(jk - jj)));
        int Jpmax = std::min(int(ji + jl), int(jk + jj));
        for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
        {
          double sixj = Z.modelspace->GetCachedSixJ(jji, jjj, J, jjk, jjl, Jprime);
          if (std::abs(sixj) < 1e-8)
            continue;
          int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
          TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
          int nkets_cc = tbc_cc.GetNumberKets();
          int indx_il = tbc_cc.GetLocalIndex(std::min(i, l), std::max(i, l)) + (i > l ? nkets_cc : 0);
          int indx_kj = tbc_cc.GetLocalIndex(std::min(j, k), std::max(j, k)) + (k > j ? nkets_cc : 0);
          double me1 = Zbar[ch_cc](indx_il, indx_kj); // do we need to use at() or is it safe to just use []?
          commij -= (2 * Jprime + 1) * sixj * me1;
        }

        if (k == l)
        {
          commji = commij;
        }
        else if (i == j)
        {
          commji = Z.modelspace->phase(ji + jj + jk + jl) * commij;
        }
        else
        {
          // now loop over the cross coupled TBME's
          parity_cc = (oi.l + ok.l) % 2;
          Tz_cc = std::abs(oi.tz2 - ok.tz2) / 2;
          Jpmin = std::max(std::abs(int(jj - jl)), std::abs(int(jk - ji)));
          Jpmax = std::min(int(jj + jl), int(jk + ji));
          for (int Jprime = Jpmin; Jprime <= Jpmax; ++Jprime)
          {
            double sixj = Z.modelspace->GetCachedSixJ(jjj, jji, J, jjk, jjl, Jprime);
            if (std::abs(sixj) < 1e-8)
              continue;
            int ch_cc = Z.modelspace->GetTwoBodyChannelIndex(Jprime, parity_cc, Tz_cc);
            TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch_cc);
            int nkets_cc = tbc_cc.GetNumberKets();
            int indx_ik = tbc_cc.GetLocalIndex(std::min(i, k), std::max(i, k)) + (i > k ? nkets_cc : 0);
            int indx_lj = tbc_cc.GetLocalIndex(std::min(l, j), std::max(l, j)) + (l > j ? nkets_cc : 0);
            // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
            // the phase we get from that flip combines with the phase from Pij, to give the phase included below
            //                 double me1 = Zbar.at(ch_cc)(indx_ik, indx_lj) ;
            double me1 = Zbar[ch_cc](indx_ik, indx_lj);
            commji -= (2 * Jprime + 1) * sixj * me1;
          }
        }

        double norm = bra.delta_pq() == ket.delta_pq() ? 1 + bra.delta_pq() : PhysConst::SQRT2;
        double zijkl = -(commij - Z.modelspace->phase(jk + jl - J) * commji) / norm;

        ZMat(ibra, iket) += zijkl;
        if (ibra != iket)
          ZMat(iket, ibra) += hZ * zijkl;
      } // for iket
    } // for ichbra
  }

  //*****************************************************************************************
  //
  //  THIS IS THE BIG UGLY ONE.
  //
  //   |          |      |          |
  //   |     __Y__|      |     __X__|
  //   |    /\    |      |    /\    |
  //   |   (  )   |  _   |   (  )   |
  //   |____\/    |      |____\/    |
  //   |  X       |      |  Y       |
  //
  //
  // -- This appears to agree with Nathan's results
  //
  /// Calculates the part of \f$ [X_{(2)},Y_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ Z^{J}_{ijkl} \f$
  /// \f[
  /// Z^{J}_{ijkl} = \sum_{ab}(n_a\bar{n}_b-\bar{n}_an_b)\sum_{J'} (2J'+1)
  /// \left[
  ///  \left\{ \begin{array}{lll}
  ///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
  /// \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} -
  ///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
  ///  -(-1)^{j_i+j_j-J}
  ///  \left\{ \begin{array}{lll}
  ///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
  /// \left( \bar{X}^{J'}_{j\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{i}} -
  ///   \bar{Y}^{J'}_{j\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{i}} \right)
  /// \right]
  /// \f]
  /// This is implemented by defining an intermediate matrix
  /// \f[
  /// \bar{Z}^{J}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
  /// \left[ \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} -
  ///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
  /// -\left( \bar{X}^{J'}_{i\bar{l}b\bar{a}}\bar{Y}^{J'}_{b\bar{a}k\bar{j}} -
  ///    \bar{Y}^{J'}_{i\bar{l}b\bar{a}}\bar{X}^{J'}_{b\bar{a}k\bar{j}} \right)\right]
  /// \f]
  /// The Pandya-transformed matrix elements are obtained with DoPandyaTransformation().
  /// The matrices \f$ \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} \f$
  /// and \f$ \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \f$
  /// are related by a Hermitian conjugation, which saves two matrix multiplications.
  /// The commutator is then given by
  /// \f[
  /// Z^{J}_{ijkl} = \sum_{J'} (2J'+1)
  /// \left[
  ///  \left\{ \begin{array}{lll}
  ///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
  ///  \bar{Z}^{J'}_{i\bar{l}k\bar{j}}
  ///  -(-1)^{j_i+j_j-J}
  ///  \left\{ \begin{array}{lll}
  ///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
  ///  \bar{Z}^{J'}_{j\bar{l}k\bar{i}}
  ///  \right]
  ///  \f]
  ///
  void comm222_phss(const Operator &X, const Operator &Y, Operator &Z)
  {
    if (X.GetParticleRank() < 2 or Y.GetParticleRank() < 2)
      return;

    // We allow for the possibility that X or Y is a scalar under rotations, but changes parity or isospin.
    // For the moment, we just call the scalar tensor routine anyway, but this means we need to be careful
    // about when we're using reduced matrix elements. The scalar-tensor routines Z = [X,Y] assume that
    // Z and Y are reduced, while X is not reduced.
    if (not(X.GetParity() == 0 and Y.GetParity() == 0 and Z.GetParity() == 0 and X.GetTRank() == 0 and Y.GetTRank() == 0 and Z.GetTRank() == 0))
    {
      if (X.GetParity() == 0 and X.GetTRank() == 0)
      {
        Operator Yred = Y;
        Yred.MakeReduced();
        Z.MakeReduced();
        comm222_phst(X, Yred, Z);
      }
      else if (Y.GetParity() == 0 and Y.GetTRank() == 0)
      {
        Operator Xred = -X; // because [X,Y] = -[Y,-X].
        Xred.MakeReduced();
        Z.MakeReduced();
        comm222_phst(Y, Xred, Z);
      }
      else if (X.GetTRank() == 0 and Y.GetTRank() == 0) // We can treat two parity-changing operators, but not two isospin-changing operators
      {
        Operator Xnred = X;
        Operator Yred = Y;
        if (Xnred.IsReduced())
          Xnred.MakeNotReduced();
        Yred.MakeReduced();
        Z.MakeReduced();
        comm222_phst(Xnred, Yred, Z);
      }
      else
      {
        std::cout << " Uh Oh. " << __func__ << " line " << __LINE__ << "  not supported with these two operators: "
                  << "  parity " << X.GetParity() << " " << Y.GetParity()
                  << "  isospin " << X.GetTRank() << " " << Y.GetTRank()
                  << "  I quit." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      Z.MakeNotReduced();
      return;
    }

    int hy = Y.IsHermitian() ? 1 : -1;
    // Create Pandya-transformed hp and ph matrix elements
    double t_start = omp_get_wtime();
    double t_start_full = t_start;

    // Construct the intermediate matrix Z_bar
    size_t nch = Z.modelspace->GetNumberTwoBodyChannels_CC();

    std::deque<arma::mat> Z_bar(nch);
    for (size_t ch = 0; ch < nch; ch++)
    {
      size_t nKets_cc = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros(nKets_cc, 2 * nKets_cc);
    }

#ifndef OPENBLAS_NOUSEOMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t ch = 0; ch < nch; ++ch)
    {
      const TwoBodyChannel_CC &tbc_cc = Z.modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      size_t nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

      arma::mat Y_bar_ph;
      arma::mat Xt_bar_ph;

      // Y has dimension (2*nph , nKets_CC)
      // X has dimension (nKets_CC, 2*nph )
      // This function gives  <ij`| X |ab`> *na(1-nb)  and   <ab`| Y |ij`>,  for both orderings ab` and ba`
      DoPandyaTransformation_SingleChannel_XandY(X, Y, Xt_bar_ph, Y_bar_ph, ch);

      //      auto& Zbar_ch = Z_bar.at(ch);
      auto &Zbar_ch = Z_bar[ch];

      if (Y_bar_ph.size() < 1 or Xt_bar_ph.size() < 1)
        continue;

      // get the phases for taking the transpose of a Pandya-transformed operator
      arma::mat PhaseMatZ(nKets_cc, nKets_cc, arma::fill::ones);
      for (index_t iket = 0; iket < nKets_cc; iket++)
      {
        const Ket &ket = tbc_cc.GetKet(iket);
        if (Z.modelspace->phase((ket.op->j2 + ket.oq->j2) / 2) < 0)
        {
          PhaseMatZ.col(iket) *= -1;
          PhaseMatZ.row(iket) *= -1;
        }
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph());
      auto PhaseMatY = PhaseMatZ.rows(phkets) * hy;

      //                                           [      |     ]
      //     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_y
      //                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
      //                                           [  Yph | Y'hp]   having kets |kj> with k>j.
      //
      //
      //      so <il|Zbar|kj> =  <il|Xbar|hp><hp|Ybar|kj> + <il|Xbar|ph><ph|Ybar|kj>
      //
      arma::mat Y_bar_ph_flip = arma::join_vert(Y_bar_ph.tail_rows(nph_kets) % PhaseMatY, Y_bar_ph.head_rows(nph_kets) % PhaseMatY);
      Zbar_ch = Xt_bar_ph * arma::join_horiz(Y_bar_ph, Y_bar_ph_flip);

      // If Z is hermitian, then XY is anti-hermitian, and so XY - YX = XY + (XY)^T
      if (Z.IsHermitian())
      {
        Zbar_ch.head_cols(nKets_cc) += Zbar_ch.head_cols(nKets_cc).t();
      }
      else // Z is antihermitian, so XY is hermitian, XY - YX = XY - XY^T
      {
        Zbar_ch.head_cols(nKets_cc) -= Zbar_ch.head_cols(nKets_cc).t();
      }

      // By taking the transpose, we get <il|Zbar|kj> with i>l and k<j, and we want the opposite
      // By the symmetries of the Pandya-transformed matrix element, that means we pick
      // up a factor hZ * phase(i+j+k+l). The hZ cancels the hXhY we have for the "head" part of the matrix
      // so we end up adding in either case.
      Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t() % PhaseMatZ;
    }

    X.profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

    // Perform inverse Pandya transform on Z_bar to get Z
    t_start = omp_get_wtime();

    // Actually, the Pandya transform has a minus sign in the definition,
    // and the ph commutator has an overall minus sign, so we're technically subtracting
    // the inverse Pandya transformation. Also, the inverse Pandya transformation
    // is just the regular Pandya transformation. The distinction in the code
    // is because some other commutator-specific things are done at the same time.
    AddInversePandyaTransformation(Z_bar, Z);

    //   Z.modelspace->scalar_transform_first_pass = false;
    X.profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;
    X.profiler.timer[__func__] += omp_get_wtime() - t_start_full;
  }

  // A much slower, but more straighforward implementation which can handle parity and isospin changing operators.
  void comm222_phss_slower(const Operator &X, const Operator &Y, Operator &Z)
  {
    double t_start = omp_get_wtime();
    auto &X2 = X.TwoBody;
    auto &Y2 = Y.TwoBody;
    auto &Z2 = Z.TwoBody;

    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto &iter : Z.TwoBody.MatEl)
    {
      ch_bra_list.push_back(iter.first[0]);
      ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

    //   int nch = Z.modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++)
    {
      size_t ch_bra = ch_bra_list[ich];
      size_t ch_ket = ch_ket_list[ich];
      TwoBodyChannel &tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);

      int J = tbc_bra.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      //     for (int ibra=0; ibra<nkets; ibra++)
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);
        size_t i = bra.p;
        size_t j = bra.q;
        Orbit &oi = Z.modelspace->GetOrbit(i);
        Orbit &oj = Z.modelspace->GetOrbit(j);

        int ketmin = 0;
        if (ch_bra == ch_ket)
          ketmin = ibra;
        for (int iket = ketmin; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);
          size_t k = ket.p;
          size_t l = ket.q;
          Orbit &ok = Z.modelspace->GetOrbit(k);
          Orbit &ol = Z.modelspace->GetOrbit(l);

          double zijkl = 0;
          int Jpmin = std::min(std::max(std::abs(oi.j2 - ol.j2), std::abs(oj.j2 - ok.j2)),
                               std::max(std::abs(oj.j2 - ol.j2), std::abs(oi.j2 - ok.j2))) /
                      2;
          int Jpmax = std::max(std::min(oi.j2 + ol.j2, oj.j2 + ok.j2),
                               std::min(oj.j2 + ol.j2, oi.j2 + ok.j2)) /
                      2;

          for (int Jp = Jpmin; Jp <= Jpmax; Jp++)
          {

            double zbar_ilkj = 0;
            double zbar_jlki = 0;
            double sixj_ijkl = AngMom::SixJ(oi.j2 * 0.5, oj.j2 * 0.5, J, ok.j2 * 0.5, ol.j2 * 0.5, Jp);
            double sixj_jikl = AngMom::SixJ(oj.j2 * 0.5, oi.j2 * 0.5, J, ok.j2 * 0.5, ol.j2 * 0.5, Jp);

            for (size_t a : Z.modelspace->all_orbits)
            {
              Orbit &oa = Z.modelspace->GetOrbit(a);
              for (size_t b : Z.modelspace->all_orbits)
              {
                Orbit &ob = Z.modelspace->GetOrbit(b);

                // "Direct" term
                double xbar_ilab = 0;
                double ybar_ilab = 0;

                int Jppmin = std::max(std::abs(oi.j2 - ob.j2), std::abs(oa.j2 - ol.j2)) / 2;
                int Jppmax = std::min(oi.j2 + ob.j2, oa.j2 + ol.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_ilab = AngMom::SixJ(oi.j2 * 0.5, ol.j2 * 0.5, Jp, oa.j2 * 0.5, ob.j2 * 0.5, Jpp);
                  double xibal = X2.GetTBME_J(Jpp, i, b, a, l);
                  double yibal = Y2.GetTBME_J(Jpp, i, b, a, l);
                  xbar_ilab -= (2 * Jpp + 1) * sixj_ilab * xibal;
                  ybar_ilab -= (2 * Jpp + 1) * sixj_ilab * yibal;
                }

                double xbar_abkj = 0;
                double ybar_abkj = 0;
                Jppmin = std::max(std::abs(oa.j2 - oj.j2), std::abs(ok.j2 - ob.j2)) / 2;
                Jppmax = std::min(oa.j2 + oj.j2, ok.j2 + ob.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_abkj = AngMom::SixJ(oa.j2 * 0.5, ob.j2 * 0.5, Jp, ok.j2 * 0.5, oj.j2 * 0.5, Jpp);
                  double xajkb = X2.GetTBME_J(Jpp, a, j, k, b);
                  double yajkb = Y2.GetTBME_J(Jpp, a, j, k, b);
                  xbar_abkj -= (2 * Jpp + 1) * sixj_abkj * xajkb;
                  ybar_abkj -= (2 * Jpp + 1) * sixj_abkj * yajkb;
                }

                zbar_ilkj += (oa.occ - ob.occ) * (xbar_ilab * ybar_abkj - ybar_ilab * xbar_abkj);

                // "Exchange" term obtained by exchanging i<->j.
                double xbar_jlab = 0;
                double ybar_jlab = 0;
                Jppmin = std::max(std::abs(oj.j2 - ob.j2), std::abs(oa.j2 - ol.j2)) / 2;
                Jppmax = std::min(oj.j2 + ob.j2, oa.j2 + ol.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_jlab = AngMom::SixJ(oj.j2 * 0.5, ol.j2 * 0.5, Jp, oa.j2 * 0.5, ob.j2 * 0.5, Jpp);
                  double xjbal = X2.GetTBME_J(Jpp, j, b, a, l);
                  double yjbal = Y2.GetTBME_J(Jpp, j, b, a, l);
                  xbar_jlab -= (2 * Jpp + 1) * sixj_jlab * xjbal;
                  ybar_jlab -= (2 * Jpp + 1) * sixj_jlab * yjbal;
                }

                double xbar_abki = 0;
                double ybar_abki = 0;
                Jppmin = std::max(std::abs(oa.j2 - oi.j2), std::abs(ok.j2 - ob.j2)) / 2;
                Jppmax = std::min(oa.j2 + oi.j2, ok.j2 + ob.j2) / 2;
                for (int Jpp = Jppmin; Jpp <= Jppmax; Jpp++)
                {
                  double sixj_abki = AngMom::SixJ(oa.j2 * 0.5, ob.j2 * 0.5, Jp, ok.j2 * 0.5, oi.j2 * 0.5, Jpp);
                  double xaikb = X2.GetTBME_J(Jpp, a, i, k, b);
                  double yaikb = Y2.GetTBME_J(Jpp, a, i, k, b);
                  xbar_abki -= (2 * Jpp + 1) * sixj_abki * xaikb;
                  ybar_abki -= (2 * Jpp + 1) * sixj_abki * yaikb;
                }

                zbar_jlki += (oa.occ - ob.occ) * (xbar_jlab * ybar_abki - ybar_jlab * xbar_abki);

              } // b
            } // a
            int phase_ij = AngMom::phase((oi.j2 + oj.j2 - 2 * J) / 2);
            zijkl += (2 * Jp + 1) * sixj_ijkl * zbar_ilkj;            // Direct
            zijkl -= (2 * Jp + 1) * sixj_jikl * zbar_jlki * phase_ij; // Exchange, with phase

          } // Jp

          // Need to normalize here, because AddToTBME expects a normalized TBME.
          if (i == j)
            zijkl /= PhysConst::SQRT2;
          if (k == l)
            zijkl /= PhysConst::SQRT2;

          Z2.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);

        } // iket
      } // ibra
    } // ch

    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
  } // comm222_phss

  /// Product of two operators. Similar to commutator, and hopefully useful for some debugging
  //
  //
  void prod110ss(const Operator &X, const Operator &Y, Operator &Z)
  {
    //     if (X.IsHermitian() and Y.IsHermitian()) return ; // I think this is the case
    //     if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return ; // I think this is the case
    if (Z.GetJRank() > 0 or Z.GetTRank() > 0 or Z.GetParity() != 0)
      return;

    arma::mat xyyx = X.OneBody * Y.OneBody - Y.OneBody * X.OneBody;
    for (auto &a : Z.modelspace->holes)
    {
      Orbit &oa = Z.modelspace->GetOrbit(a);
      Z.ZeroBody += (oa.j2 + 1) * oa.occ * xyyx(a, a);
    }
  }

} // namespace Commutator
