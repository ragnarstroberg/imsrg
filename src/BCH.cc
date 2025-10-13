
#include "BCH.hh"
#include "Commutator.hh"

namespace BCH
{

  bool use_goose_tank_correction = false;
  bool use_brueckner_bch = false;
  bool bch_skip_ieq1 = false;
  bool use_factorized_correction = false;
  bool use_factorized_correction_BCH_Product = false; 
  bool use_factorized_correct_ZBterm = false;                // correct zero body term at the end of IMSRG2
  bool only_2b_omega = false;
  bool TurnOn_comm223_231 = true;
  bool TurnOn_comm223_232 = true;

  double bch_transform_threshold = 1e-9;
  double bch_product_threshold = 1e-4;






  void SetUseBruecknerBCH(bool tf)
  {
    use_brueckner_bch = tf;
  }

  void SetUseGooseTank(bool tf)
  {
    use_goose_tank_correction = tf;
  }

  void SetBCHSkipiEq1(bool tf)
  {
    bch_skip_ieq1 = tf;
  }

  void SetUseFactorizedCorrection(bool tf)
  {
    use_factorized_correction = tf;
    use_factorized_correction_BCH_Product = tf;
  }

  void SetUseFactorizedCorrectionBCH_product(bool tf)
  {
    use_factorized_correction_BCH_Product = tf;
  }

  void SetComm223_231(bool tf)
  {
    TurnOn_comm223_231 = tf;
  }

  void SetComm223_232(bool tf)
  {
    TurnOn_comm223_232 = tf;
  }

  void SetUseFactorized_Correct_ZBTerm(bool tf)
  {
    use_factorized_correct_ZBterm = tf;
  }

  void SetOnly2bOmega(bool tf)
  {
    only_2b_omega = tf;
  }


  void Set_BCH_Transform_Threshold(double x)
  {
    bch_transform_threshold = x;
  }

  void Set_BCH_Product_Threshold(double x)
  {
    bch_product_threshold = x;
  }



  //*****************************************************************************************
  /// BCH_Transform(X,Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
  /// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
  /// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
  /// with all commutators truncated at the two-body level.
  Operator BCH_Transform(const Operator &OpIn, const Operator &Omega)
  {
    return use_brueckner_bch ? Brueckner_BCH_Transform(OpIn, Omega) : Standard_BCH_Transform(OpIn, Omega);
  }

  /// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
  /// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
  /// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
  /// with all commutators truncated at the two-body level.
  Operator Standard_BCH_Transform(const Operator &OpIn, const Operator &Omega)
  {

    double t_start = omp_get_wtime();
    int max_iter = 40;
    int warn_iter = 12;
    double nx = OpIn.Norm();
    double ny = Omega.Norm();
    Operator OpOut = OpIn;
    if ((OpOut.GetNumberLegs() % 2 == 0) and OpOut.GetParticleRank() < 2)
    {
      OpOut.SetParticleRank(2);
    }
    if (OpOut.GetJRank() == 0 and (OpOut.GetTRank() != 0 or OpOut.GetParity() != 0) and OpOut.IsReduced() == false)
    {
      OpOut.MakeReduced();
    }
      //   if (use_imsrg3 and not OpOut.ThreeBody.is_allocated )
      //   This was a bug, intruduced Feb 2022, fixed May 2022.
      //   Only call SetMode and SetParticleRank if ThreeBody is not already allocated.
      //   Otherwise, we are erasing the 3N inherited from OpIn.
      //   if (use_imsrg3  )
      if (Commutator::use_imsrg3 and not(OpOut.GetParticleRank() >= 3 and OpOut.ThreeBody.IsAllocated() and OpOut.ThreeBody.Is_PN_Mode())) // IMSRG(3) commutators work in PN mode only so far.
      {
        std::cout << __func__ << "  Allocating Three Body with pn mode" << std::endl;
        OpOut.ThreeBody.SetMode("pn");
        OpOut.SetParticleRank(3);
      }
    double factorial_denom = 1.0;

    Operator chi, chi2; // auxiliary one-body operator used to recover 4th-order quadruples.
    if (use_goose_tank_correction)
    {
      chi = OpIn;
      chi.SetParticleRank(1);
      chi.Erase();
    }

    if (nx > bch_transform_threshold)
    {
      //     Operator OpNested = OpIn;
      Operator OpNested = OpOut;
//      Operator OpNested_last = OpOut*0;
//      Operator OpNested_last_last = OpOut*0;
//      Operator OpNested_last_last_last = OpOut*0;

      double epsilon = nx * exp(-2 * ny) * bch_transform_threshold / (2 * ny); // this should probably be explained somewhere...
      for (int i = 1; i <= max_iter; ++i)
      {
        // Put this before the bch_skip step, so we still copy and allocate chi2.
        if (use_factorized_correction)
        {
          chi2 = chi;     // keep nested commutator from two steps ago
          chi = OpNested; // save nested commutator from previous step
        }

        // Specifically for the perturbative triples, we need 1/(i+1)! rather than 1/i!
        // This is because we have Wbar = [Omega,H]_3b + 1/2![Omega,[Omega,H]]_3b + 1/3![Omega,[Omega,[Omega,H]]]_3b + ...
        //                              = [Omega, Htilde]_3b,  where Htilde = H + 1/2![Omega,H] + 1/3![Omega,[Omega,H]] + ...
        if (bch_skip_ieq1 and i == 1)
          continue;

        if (use_goose_tank_correction)
        {
          auto chi_last = chi.OneBody;
          chi = GooseTankUpdate(Omega, OpNested);
          OpNested.OneBody += chi_last; // add the chi from the previous step to OpNested.
        }



//        OpNested_last_last_last = OpNested_last_last;
//        OpNested_last_last = OpNested_last;
//        OpNested_last = OpNested;

        OpNested = Commutator::Commutator(Omega, OpNested); // the ith nested commutator

        int i_min_factorized = bch_skip_ieq1 ? 2 : 1;
        if (i > i_min_factorized and use_factorized_correction)
        {
            Operator Op_DoubleNested = OpNested;
            Op_DoubleNested.Erase();

            if ( TurnOn_comm223_231 )
              Commutator::FactorizedDoubleCommutator::comm223_231(Omega, chi2, Op_DoubleNested);
            if ( TurnOn_comm223_232 )
              Commutator::FactorizedDoubleCommutator::comm223_232(Omega, chi2, Op_DoubleNested);
            OpNested += Op_DoubleNested;
        }


        factorial_denom /= i;

        OpOut += factorial_denom * OpNested;

        if (OpOut.rank_J > 0)
        {
          auto id5 = OpOut.modelspace->GetOrbitIndex(0,2,5,+1);
          std::cout << "Tensor BCH, i=" << i << "  Norm = " << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.OneBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.TwoBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.ThreeBody.Norm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << std::endl;
//                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << "    d5d5 = " << OpNested.OneBody(id5,id5) << " sum " << OpOut.OneBody(id5,id5) << std::endl;
        }
        epsilon *= i + 1;
        if (OpNested.Norm() < epsilon)
          break;
        if (i == warn_iter)
          std::cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << std::endl;
        else if (i == max_iter)
          std::cout << "Warning: BCH_Transform didn't coverge after " << max_iter << " nested commutators" << std::endl;
      }
    }
    //   std::cout << "Done with BCH_Transform, 3-body norm of OpOut = " << OpOut.ThreeBodyNorm() << std::endl;
//    OpIn.profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
    OpIn.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return OpOut;
  }

  //  Update the auxiliary one-body operator chi, using Omega and the ith nested commutator
  //  This has not been tested for tensor commutators, but it *should* work.
  //
  Operator GooseTankUpdate(const Operator &Omega, const Operator &OpNested)
  {
    double t_start = omp_get_wtime();
    Operator goosetank_chi = Operator(*(OpNested.modelspace), OpNested.rank_J, OpNested.rank_T, OpNested.parity, 1);
    goosetank_chi.EraseOneBody();
    if (goosetank_chi.rank_J == 0)
    {
      Commutator::comm221ss(Omega, OpNested, goosetank_chi); // update chi.
    }
    else
    {
      Commutator::comm222_pp_hh_221st(Omega, OpNested, goosetank_chi); // update chi.
    }
    goosetank_chi.Symmetrize();                    // the commutator call only does half the matrix, so we symmetrize
                                                   //   int norbits = OpNested.modelspace->GetNumberOrbits();
                                                   //   for (int i=0;i<norbits;++i)  // enforce n_in_j + nbar_i nbar_j
    for (auto i : OpNested.modelspace->all_orbits) // enforce n_in_j + nbar_i nbar_j
    {
      Orbit &oi = OpNested.modelspace->GetOrbit(i);
      //     for (int j=0;j<norbits;++j)
      for (auto j : OpNested.modelspace->all_orbits)
      {
        Orbit &oj = OpNested.modelspace->GetOrbit(j);
        goosetank_chi.OneBody(i, j) *= oi.occ * oj.occ + (1.0 - oi.occ) * (1.0 - oj.occ);
      }
    }
    OpNested.profiler.timer[__func__] += omp_get_wtime() - t_start;
    return goosetank_chi;
  }

  /// Variation of the BCH transformation procedure
  /// requested by a one Dr. T.D. Morris
  /// \f[ e^{\Omega_1 + \Omega_2} X e^{-\Omega_1 - \Omega_2}
  ///    \rightarrow
  ///  e^{\Omega_2} e^{\Omega_1}  X e^{-\Omega_1} e^{-\Omega_2} \f]
  Operator Brueckner_BCH_Transform(const Operator &OpIn, const Operator &Omega)
  {
    Operator Omega1 = Omega;
    Operator Omega2 = Omega;
    Omega1.SetParticleRank(1);
    Omega1.EraseTwoBody();
    Omega2.EraseOneBody();
    Operator OpOut = Standard_BCH_Transform(OpIn, Omega1);
    OpOut = Standard_BCH_Transform(OpOut, Omega2);
    return OpOut;
  }

  //*****************************************************************************************
  // Baker-Campbell-Hausdorff formula
  //  returns Z, where
  //  exp(Z) = exp(X) * exp(Y).
  //  Z = X + Y + 1/2[X, Y]
  //     + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]]
  //     - 1/24 [Y,[X,[X,Y]]]
  //     - 1/720 [Y,[Y,[Y,[Y,X]]]] - 1/720 [X,[X,[X,[X,Y]]]]
  //     + ...
  //*****************************************************************************************
  /// X.BCH_Product(Y) returns \f$Z\f$ such that \f$ e^{Z} = e^{X}e^{Y}\f$
  /// by employing the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
  /// \f[ Z = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}([X,[X,Y]]+[Y,[Y,X]]) + \ldots \f]
  //*****************************************************************************************
  Operator BCH_Product(Operator &X, Operator &Y)
  {
    double tstart = omp_get_wtime();
    double nx = X.Norm();
    double ny = Y.Norm();
    std::vector<double> bernoulli = {1.0, -0.5, 1. / 6, 0.0, -1. / 30, 0.0, 1. / 42, 0, -1. / 30};
    std::vector<double> factorial = {1.0, 1.0, 2.0, 6.0, 24., 120., 720., 5040., 40320.};

    Operator Z = X + Y;
  

    // if we set the option only_2b_omega
    // then temporarily switch to imsrg2 for updating Omega
    bool _save_imsrg3 = Commutator::use_imsrg3;
    if (only_2b_omega)
    {
      Z.ThreeBody.Erase();
      Commutator::use_imsrg3 = false;
    }
    Operator Nested = Commutator::Commutator(Y, X); // [Y,X]

    double nxy = Nested.Norm();
    // We assume X is small, but just in case, we check if we should include the [X,[X,Y]] term.
    if (nxy * nx > bch_product_threshold)
    {
      Z += (1. / 12) * Commutator::Commutator(Nested, X);
      if (use_factorized_correction_BCH_Product)
      {
        Operator Op_2b1b = Nested;
        Op_2b1b.Erase();
        Commutator::FactorizedDoubleCommutator::comm223_231(X, Y, Op_2b1b);
        Commutator::FactorizedDoubleCommutator::comm223_232(X, Y, Op_2b1b);
        Z += (1. / 12) * Op_2b1b ;
      }
    }

    Operator chi, chi2; // intermidate operator for factorization code  B.C. He
    if (use_factorized_correction_BCH_Product)
    {
      chi2 = X; // keep nested commutator from two steps ago
      chi = X;  // save nested commutator from previous step
    }
    size_t k = 1;
    // k=1 adds 1/2[X,Y],  k=2 adds 1/12 [Y,[Y,X]], k=4 adds -1/720 [Y,[Y,[Y,[Y,X]]]], and so on.
    //   while( Nested.Norm() > bch_product_threshold and k<9)
    while (nxy > bch_product_threshold)
    {
      if ((k < 2) or (k % 2 == 0))
      {
        Z += (bernoulli[k] / factorial[k]) * Nested;
      }

      k++;
      if (k >= bernoulli.size())
        break; // don't evaluate the commutator if we're not going to use it
      if (2 * ny * nxy < bch_product_threshold)
        break;

      if (use_factorized_correction_BCH_Product)
      {
        chi2 = chi;   // keep nested commutator from two steps ago
        chi = Nested; // save nested commutator from previous step
      }

      Nested = Commutator::Commutator(Y, Nested);

      if (k >= 2 and use_factorized_correction_BCH_Product)
      {
        Operator Op_2b1b = Nested;
        Op_2b1b.Erase();
        Commutator::FactorizedDoubleCommutator::comm223_231(Y, chi2, Op_2b1b);
        Commutator::FactorizedDoubleCommutator::comm223_232(Y, chi2, Op_2b1b);
        Nested += Op_2b1b;
      }

      nxy = Nested.Norm();
    }

    Commutator::use_imsrg3 = _save_imsrg3; // set it back to how it was.

    X.profiler.timer[__func__] += omp_get_wtime() - tstart;
    return Z;
  }

  double EstimateBCHError(Operator &Omega, Operator H)
  {
     return 0; // This is not ready for prime time.
//    double normOmega = Omega.Norm();
//    double normOmega2 = Omega.TwoBodyNorm();
//    double normH2 = H.TwoBodyNorm();
//    Operator H223 = H;
//    H223.Erase();
//    comm222_phss(Omega, H, H223);
//    comm222_pp_hh_221ss(Omega, H, H223);
//    if (H223.GetParticleRank() < 3)
//      H223.SetParticleRank(3);
//    //  if ( not H223.ThreeBody.IsAllocated() )  H223.ThreeBody.SwitchToPN_and_discard();
//    if (not H223.ThreeBody.IsAllocated())
//      H223.ThreeBody.Allocate();
//
//    comm223ss(Omega, H, H223);
//    comm220ss(Omega, H, H223);
//
//    double Norm220 = H223.ZeroBody;
//
//    double Norm223 = H223.ThreeBodyNorm();
//    Operator H2223 = H223;
//    H2223.Erase();
//    comm223ss(Omega, H223, H2223);
//    double Norm2223 = H2223.ThreeBodyNorm();
//
//    // now 4 nested commutators
//    H223.Erase();
//    H2223.Erase();
//    H223 = Commutator::Commutator(Omega, H);     // 1 nested
//    H2223 = Commutator(Omega, H223); // 2 nested
//    H223 = Commutator(Omega, H2223); // 3 nested
//    H2223 = Commutator(Omega, H223);
//    double Norm4nested = H2223.Norm();
//
//    double est_err = 2. / 3 * normOmega * normOmega * Norm223 + 1. / 6 * normOmega * normOmega * Norm2223 + exp(2. * normOmega) / 24 * Norm4nested;
//    std::cout << "Contributions to err " << 2. / 3 * normOmega * normOmega * Norm223 << "  " << 1. / 6 * normOmega * normOmega * Norm2223
//              << "  " << exp(2. * normOmega) / 24 * Norm4nested << "    220 = " << Norm220 << "  ||H2|| ,||Omega2|| = " << normH2 << " " << normOmega2 << std::endl;
//
//    return est_err;
  }
  std::tuple<Operator, Operator> BCH_ProductPV(Operator &X, Operator &XPV, Operator &Y, Operator &YPV)
  {
    double tstart = omp_get_wtime();
    double nx = X.Norm();
    double ny = Y.Norm();
    double nxPV = XPV.Norm();
    double nyPV = YPV.Norm();
    std::vector<double> bernoulli = {1.0, -0.5, 1. / 6, 0.0, -1. / 30, 0.0, 1. / 42, 0, -1. / 30};
    std::vector<double> factorial = {1.0, 1.0, 2.0, 6.0, 24., 120., 720., 5040., 40320.};

    Operator Z = X + Y;
    Operator ZPV = XPV+YPV;

    
    Operator Nested = Commutator::Commutator(Y, X) + Commutator::Commutator(YPV, XPV); // [Y,X]
    Operator NestedPV = Commutator::Commutator(YPV, X) + Commutator::Commutator(Y, XPV); // [Y,X]

    double nxy = Nested.Norm();
    double nxyPV = NestedPV.Norm();
    // We assume X is small, but just in case, we check if we should include the [X,[X,Y]] term.
    if (sqrt(nxy * nx* nxy * nx + nxyPV * nxyPV * nxPV * nxPV) > bch_product_threshold)
    {
      Z += (1. / 12) * (Commutator::Commutator(Nested, X)+Commutator::Commutator(NestedPV,XPV));
      ZPV += (1. / 12) * (Commutator::Commutator(NestedPV, X)+Commutator::Commutator(Nested,XPV));
    }

    
    size_t k = 1;
    // k=1 adds 1/2[X,Y],  k=2 adds 1/12 [Y,[Y,X]], k=4 adds -1/720 [Y,[Y,[Y,[Y,X]]]], and so on.
    //   while( Nested.Norm() > bch_product_threshold and k<9)
    while (sqrt(nxy*nxy+nxyPV+nxyPV) > bch_product_threshold)
    {
      if ((k < 2) or (k % 2 == 0))
      {
        Z += (bernoulli[k] / factorial[k]) * Nested;
        ZPV += (bernoulli[k] / factorial[k]) * NestedPV;
      }

      k++;
      if (k >= bernoulli.size())
        break; // don't evaluate the commutator if we're not going to use it
      if (2 * sqrt(ny * nxy * ny * nxy + nyPV * nyPV * nxyPV * nxyPV) < bch_product_threshold)
        break;

      Operator Optmp = Commutator::Commutator(Y, Nested) + Commutator::Commutator(YPV, NestedPV);
      Operator OptmpPV = Commutator::Commutator(YPV, Nested) + Commutator::Commutator(Y, NestedPV);
      Nested = Optmp;
      NestedPV = OptmpPV;

      nxy = Nested.Norm();
      nxyPV = NestedPV.Norm();
    }

    X.profiler.timer["BCH_Product"] += omp_get_wtime() - tstart;
    return {Z, ZPV};
  }

  std::tuple<Operator, Operator> BCH_TransformPV(const Operator &Op, const Operator &OpPV, const Operator &Omega, const Operator &OmegaPV)
  {
    return  Standard_BCH_TransformPV(Op, OpPV, Omega, OmegaPV);
  }

  std::tuple<Operator, Operator> BCH_TransformHPV(const Operator &Op, const Operator &OpPV, const Operator &Omega, const Operator &OmegaPV)
  {
    return  Standard_BCH_TransformHPV(Op, OpPV, Omega, OmegaPV);
  }

  /// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
  /// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
  /// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
  /// with all commutators truncated at the two-body level.
  std::tuple<Operator, Operator> Standard_BCH_TransformPV(const Operator &OpIn, const Operator &OpInPV, const Operator &Omega, const Operator &OmegaPV)
  {

    double t_start = omp_get_wtime();
    int max_iter = 40;
    int warn_iter = 12;
    double nx = OpIn.Norm();
    double ny = Omega.Norm();
    double nxPV = OpInPV.Norm();
    double nyPV = OmegaPV.Norm();
    Operator OpOut = OpIn;
    Operator OpOutPV = OpInPV;
    if ((OpOut.GetNumberLegs() % 2 == 0) and OpOut.GetParticleRank() < 2)
    {
      OpOut.SetParticleRank(2);
    }

    if (OpOut.GetJRank() == 0 and (OpOut.GetTRank() != 0 or OpOut.GetParity() != 0) and OpOut.IsReduced() == false)
    {
      OpOut.MakeReduced();
    }

    if (OpOutPV.GetJRank() == 0 and (OpOutPV.GetTRank() != 0 or OpOutPV.GetParity() != 0) and OpOutPV.IsReduced() == false)
    {
      OpOutPV.MakeReduced();
    }

    double factorial_denom = 1.0;

    
    if (sqrt(nx*nx+nxPV*nxPV) > bch_transform_threshold)
    {
      Operator OpNested = OpOut;
      Operator OpNestedPV = OpOutPV;

      double epsilon = sqrt(nx * nx + nxPV * nxPV) * exp(-2 * sqrt(ny * ny + nyPV * nyPV)) * bch_transform_threshold / (2 * sqrt(ny * ny + nyPV * nyPV)); // this should probably be explained somewhere...
      for (int i = 1; i <= max_iter; ++i)
      {

        // Specifically for the perturbative triples, we need 1/(i+1)! rather than 1/i!
        // This is because we have Wbar = [Omega,H]_3b + 1/2![Omega,[Omega,H]]_3b + 1/3![Omega,[Omega,[Omega,H]]]_3b + ...
        //                              = [Omega, Htilde]_3b,  where Htilde = H + 1/2![Omega,H] + 1/3![Omega,[Omega,H]] + ...
        if (bch_skip_ieq1 and i == 1)
          continue;
        Operator Optmp = Commutator::Commutator(OmegaPV, OpNestedPV) + Commutator::Commutator(Omega, OpNested)  ;
        Operator OptmpPV = Commutator::Commutator(OmegaPV, OpNested) + Commutator::Commutator(Omega, OpNestedPV);
        OpNested = Optmp; // the ith nested commutator
        OpNestedPV = OptmpPV;
        factorial_denom /= i;
        OpOut += factorial_denom * OpNested;
        OpOutPV += factorial_denom * OpNestedPV;
        if (OpOut.rank_J > 0)
        {
          std::cout << "Tensor BCH, i=" << i << "  Norm = " << std::setw(12) << std::setprecision(8) << std::fixed << OpOut.OneBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.TwoBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.ThreeBody.Norm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << std::endl;
          std::cout << "Tensor BCH, i=" << i << "  NormPV = " << std::setw(12) << std::setprecision(8) << std::fixed << OpOutPV.OneBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.TwoBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.ThreeBody.Norm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.Norm() << std::endl;
        }
        epsilon *= i + 1;
        if (sqrt(OpNested.Norm()*OpNested.Norm()+OpNestedPV.Norm()*OpNestedPV.Norm()) < epsilon)
          break;
        if (i == warn_iter)
          std::cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << std::endl;
        else if (i == max_iter)
          std::cout << "Warning: BCH_Transform didn't coverge after " << max_iter << " nested commutators" << std::endl;
      }
    }

    OpIn.profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
    return {OpOut,OpOutPV};
  }

  std::tuple<Operator, Operator> Standard_BCH_TransformHPV(const Operator &OpIn, const Operator &OpInPV, const Operator &Omega, const Operator &OmegaPV)
  {

    double t_start = omp_get_wtime();
    int max_iter = 40;
    int warn_iter = 12;
    double nx = OpIn.Norm();
    double ny = Omega.Norm();
    double nxPV = OpInPV.Norm();
    double nyPV = OmegaPV.Norm();
    Operator OpOut = OpIn;
    Operator OpOutPV = OpInPV;
    if ((OpOut.GetNumberLegs() % 2 == 0) and OpOut.GetParticleRank() < 2)
    {
      OpOut.SetParticleRank(2);
    }

    if (OpOut.GetJRank() == 0 and (OpOut.GetTRank() != 0 or OpOut.GetParity() != 0) and OpOut.IsReduced() == false)
    {
      OpOut.MakeReduced();
    }

    if (OpOutPV.GetJRank() == 0 and (OpOutPV.GetTRank() != 0 or OpOutPV.GetParity() != 0) and OpOutPV.IsReduced() == false)
    {
      OpOutPV.MakeReduced();
    }

    double factorial_denom = 1.0;

    if (sqrt(nx * nx + nxPV * nxPV) > bch_transform_threshold)
    {
      Operator OpNested = OpOut;
      Operator OpNestedPV = OpOutPV;

      double epsilon = sqrt(nx * nx + nxPV * nxPV) * exp(-2 * sqrt(ny * ny + nyPV * nyPV)) * bch_transform_threshold / (2 * sqrt(ny * ny + nyPV * nyPV)); // this should probably be explained somewhere...
      for (int i = 1; i <= max_iter; ++i)
      {

        // Specifically for the perturbative triples, we need 1/(i+1)! rather than 1/i!
        // This is because we have Wbar = [Omega,H]_3b + 1/2![Omega,[Omega,H]]_3b + 1/3![Omega,[Omega,[Omega,H]]]_3b + ...
        //                              = [Omega, Htilde]_3b,  where Htilde = H + 1/2![Omega,H] + 1/3![Omega,[Omega,H]] + ...
        if (bch_skip_ieq1 and i == 1)
          continue;
        Operator Optmp = Commutator::Commutator(Omega, OpNested);
        Operator OptmpPV = Commutator::Commutator(OmegaPV, OpNested) + Commutator::Commutator(Omega, OpNestedPV);
        OpNested = Optmp; // the ith nested commutator
        OpNestedPV = OptmpPV;
        factorial_denom /= i;
        OpOut += factorial_denom * OpNested;
        OpOutPV += factorial_denom * OpNestedPV;
        if (OpOut.rank_J > 0)
        {
          std::cout << "Tensor BCH, i=" << i << "  Norm = " << std::setw(12) << std::setprecision(8) << std::fixed << OpOut.OneBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.TwoBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.ThreeBody.Norm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNested.Norm() << std::endl;
          std::cout << "Tensor BCH, i=" << i << "  NormPV = " << std::setw(12) << std::setprecision(8) << std::fixed << OpOutPV.OneBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.TwoBodyNorm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.ThreeBody.Norm() << " "
                    << std::setw(12) << std::setprecision(8) << std::fixed << OpNestedPV.Norm() << std::endl;
        }
        epsilon *= i + 1;
        if (sqrt(OpNested.Norm() * OpNested.Norm() + OpNestedPV.Norm() * OpNestedPV.Norm()) < epsilon)
          break;
        if (i == warn_iter)
          std::cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << std::endl;
        else if (i == max_iter)
          std::cout << "Warning: BCH_Transform didn't coverge after " << max_iter << " nested commutators" << std::endl;
      }
    }

    OpIn.profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
    return {OpOut, OpOutPV};
  }
}// namespace BCH
