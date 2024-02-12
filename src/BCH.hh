#ifndef BCH_hh
#define BCH_hh 1

#include "Operator.hh"

namespace BCH
{

    extern bool use_goose_tank_correction;
    extern bool use_brueckner_bch;
    extern bool bch_skip_ieq1;
    extern bool use_factorized_correction;
    extern bool use_factorized_correction_BCH_Product;
    extern bool use_factorized_correct_ZBterm;
    extern bool only_2b_omega;

    extern double bch_transform_threshold;
    extern double bch_product_threshold;



    Operator BCH_Product(Operator &X, Operator &Y);
    Operator BCH_Transform(const Operator &Op, const Operator &Omega);
    Operator Standard_BCH_Transform(const Operator &Op, const Operator &Omega);
    Operator Brueckner_BCH_Transform(const Operator &Op, const Operator &Omega);

    double EstimateBCHError(Operator &Omega, Operator H);

    Operator GooseTankUpdate(const Operator &Omega, const Operator &Nested);


    void SetUseBruecknerBCH(bool tf);
    void SetUseGooseTank(bool tf);
    void SetBCHSkipiEq1(bool tf);
    void SetUseFactorizedCorrection(bool tf);
    void SetUseFactorizedCorrectionBCH_product(bool tf);
    void SetUseFactorized_Correct_ZBTerm(bool tf);
    void SetOnly2bOmega(bool tf);

    void Set_BCH_Transform_Threshold(double x);
    void Set_BCH_Product_Threshold(double x);

}// namespace BCH

#endif
