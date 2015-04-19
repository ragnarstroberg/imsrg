
#include "Operator.hh"
#include "AngMom.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif

using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

double  Operator::bch_transform_threshold = 1e-6;
double  Operator::bch_product_threshold = 1e-4;
map<string, double> Operator::timer;



/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator() :
   modelspace(NULL), nChannels(0), hermitian(true), antihermitian(false),
    rank_J(0), rank_T(0), parity(0), particle_rank(2)
{}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    hermitian(true), antihermitian(false), modelspace(&ms), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms,Jrank,Trank,p),
    ThreeBody(&ms, ms.GetN3max()),E3max(ms.GetN3max()),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank)
{
  TwoBody.Allocate();
  if (particle_rank >=3) ThreeBody.Allocate();
}



Operator::Operator(ModelSpace& ms) :
    hermitian(true), antihermitian(false), modelspace(&ms), ZeroBody(0) ,
    nChannels(ms.GetNumberTwoBodyChannels()) ,
    rank_J(0), rank_T(0), parity(0), particle_rank(2),
    OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms,0,0,0),
    ThreeBody(&ms,ms.GetN3max()), E3max(ms.GetN3max())
{
  TwoBody.Allocate();
}

Operator::Operator(const Operator& op)
{
   Copy(op);
}


/////////// COPY METHOD //////////////////////////
void Operator::Copy(const Operator& op)
{
   modelspace    = op.modelspace;
   nChannels     = op.nChannels;
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
   rank_J        = op.rank_J;
   rank_T        = op.rank_T;
   parity        = op.parity;
   particle_rank = op.particle_rank;
   E2max         = op.E2max;
   E3max         = op.E3max;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   TwoBody       = op.TwoBody;
   ThreeBody     = op.ThreeBody;
}

/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////
Operator& Operator::operator=(const Operator& rhs)
{
   Copy(rhs);
   return *this;
}

// multiply operator by a scalar
Operator& Operator::operator*=(const double rhs)
{
   ZeroBody *= rhs;
   OneBody *= rhs;
   TwoBody *= rhs;
   return *this;
}

Operator Operator::operator*(const double rhs) const
{
   Operator opout = Operator(*this);
   opout *= rhs;
   return opout;
}

// Add non-member operator so we can multiply an operator
// by a scalar from the lhs, i.e. s*O = O*s
Operator operator*(const double lhs, const Operator& rhs)
{
   return rhs * lhs;
}
Operator operator*(const double lhs, const Operator&& rhs)
{
   return rhs * lhs;
}


// divide operator by a scalar
Operator& Operator::operator/=(const double rhs)
{
   return *this *=(1.0/rhs);
}

Operator Operator::operator/(const double rhs) const
{
   Operator opout = Operator(*this);
   opout *= (1.0/rhs);
   return opout;
}

// Add operators
Operator& Operator::operator+=(const Operator& rhs)
{
   ZeroBody += rhs.ZeroBody;
   OneBody  += rhs.OneBody;
   TwoBody  += rhs.TwoBody;
   return *this;
}

Operator Operator::operator+(const Operator& rhs) const
{
   return ( Operator(*this) += rhs );
}

// Subtract operators
Operator& Operator::operator-=(const Operator& rhs)
{
   ZeroBody -= rhs.ZeroBody;
   OneBody -= rhs.OneBody;
   TwoBody -= rhs.TwoBody;
   return *this;
}

Operator Operator::operator-(const Operator& rhs) const
{
   return ( Operator(*this) -= rhs );
}

Operator Operator::operator-() const
{
   return (*this)*-1.0;
}



//*******************************************************************
/// Get three body matrix element in proton-neutron formalism.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
/*
double Operator::GetThreeBodyME_pn(int Jab_in, int Jde_in, int J2, int a, int b, int c, int d, int e, int f)
{
  return ThreeBody.GetThreeBodyME_pn(Jab_in, Jde_in, J2, a, b, c, d, e, f);
}
*/

//*******************************************************************
/// Get three body matrix element in isospin formalism
/// \f$ V_{abcdef}^{J_{ab}J_{de}Jt_{ab}t_{de}T} \f$
///  (which is how they're stored).
/// The elements are stored with the following restrictions: \f$ a\geq b \geq c\f$,
/// \f$ d\geq e \geq f\f$, \f$ a\geq d\f$. If \f$ a=d\f$ then \f$ b \geq e\f$,
/// and if \f$ b=e \f$ then \f$ c \geq f \f$.
/// Other orderings are obtained by recoupling on the fly.
//*******************************************************************
/*
double Operator::GetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   return ThreeBody.GetThreeBodyME(Jab_in, Jde_in, J2, tab_in, tde_in, T2, a_in, b_in, c_in, d_in, e_in, f_in);
}
*/
//*******************************************************************
/// Set a three body matrix element. Since only a subset of orbit
/// orderings are stored, we need to recouple if the input ordering
/// is different.
//*******************************************************************
/*
void Operator::SetThreeBodyME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V)
{
   ThreeBody.SetThreeBodyME(Jab_in, Jde_in, J2, tab_in, tde_in, T2, a_in, b_in, c_in, d_in, e_in, f_in,V);
}
*/


void Operator::PrintTimes()
{
   cout << "==== TIMES ====" << endl;
   for ( auto it : timer )
   {
     cout << it.first << ":  " << it.second  << endl;
   }
}

////////////////// MAIN INTERFACE METHODS //////////////////////////

Operator Operator::DoNormalOrdering()
{
   if (particle_rank==3)
   {
      return DoNormalOrdering3();
   }
   else
   {
      return DoNormalOrdering2();
   }
}

//*************************************************************
///  Normal ordering of a 2body operator
///  currently this only handles scalar operators
//*************************************************************
Operator Operator::DoNormalOrdering2()
{
   Operator opNO = *this;


   for (auto& k : modelspace->holes) // loop over hole orbits
   {
      opNO.ZeroBody += (modelspace->GetOrbit(k).j2+1) * OneBody(k,k);
   }


   int norbits = modelspace->GetNumberOrbits();

   for ( auto& itmat : TwoBody.MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      auto& matrix = itmat.second;
      
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J_ket = tbc_ket.J;

      // Zero body part
      arma::vec diagonals = matrix.diag();
      auto hh = tbc_ket.GetKetIndex_hh();
      opNO.ZeroBody += arma::sum( diagonals.elem(hh) ) * (2*J_ket+1);

      // One body part
      for (long long unsigned int a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         int bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for ( auto& b : modelspace->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) // OneBodyChannels should be moved to the operator, to accommodate tensors
         {
            if (b < bstart) continue;
            for (auto& h : modelspace->holes)  // C++11 syntax
            {
               opNO.OneBody(a,b) += (2*J_ket+1.0)/(oa.j2+1)  * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
            }
         }
      }
   } // loop over channels

   if (hermitian) opNO.Symmetrize();
   if (antihermitian) opNO.AntiSymmetrize();

   return opNO;
}



//*******************************************************************************
///   Normal ordering of a three body operator. Start by generating the normal ordered
///   two body piece, then use DoNormalOrdering2() to get the rest. (Note that there
///   are some numerical factors).
///   The normal ordered two body piece is 
///   \f[ \Gamma^J_{ijkl} = V^J_{ijkl} + \sum_a n_a  \sum_K \frac{2K+1}{2J+1} V^{(3)JJK}_{ijakla} \f]
///   Right now, this is only set up for scalar operators, but I don't anticipate
///   handling 3body tensor operators in the near future.
//*******************************************************************************
Operator Operator::DoNormalOrdering3()
{
   Operator opNO3 = Operator(*modelspace);
//   cout << "size of ThreeBody = " << ThreeBody.MatEl.size() << endl;
//   #pragma omp parallel for
   for ( auto& itmat : opNO3.TwoBody.MatEl )
   {
      int ch = itmat.first[0]; // assume ch_bra = ch_ket for 3body...
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& Gamma = (arma::mat&) itmat.second;
      for (int ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         for (int iket=ibra; iket<tbc.GetNumberKets(); ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            for (auto& a : modelspace->holes)
            {
               Orbit & oa = modelspace->GetOrbit(a);
               if ( (2*(oi.n+oj.n+oa.n)+oi.l+oj.l+oa.l)>E3max) continue;
               if ( (2*(ok.n+ol.n+oa.n)+ok.l+ol.l+oa.l)>E3max) continue;
               int kmin2 = abs(2*tbc.J-oa.j2);
               int kmax2 = 2*tbc.J+oa.j2;
               for (int K2=kmin2; K2<=kmax2; K2+=2)
               {
                  Gamma(ibra,iket) += (K2+1) * ThreeBody.GetME_pn(tbc.J,tbc.J,K2,i,j,a,k,l,a);
               }
            }
            Gamma(ibra,iket) /= (2*tbc.J+1);
         }
      }
   }
   opNO3.Symmetrize();
   Operator opNO2 = opNO3.DoNormalOrdering2();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);

   // Also normal order the 1 and 2 body pieces
   opNO2 += DoNormalOrdering2();
   return opNO2;

}




void Operator::Erase()
{
  EraseZeroBody();
  EraseOneBody();
  TwoBody.Erase();
  if (particle_rank >=3)
    ThreeBody.Erase();
}

void Operator::EraseOneBody()
{
   OneBody.zeros();
}

void Operator::EraseTwoBody()
{
 TwoBody.Erase();
}


void Operator::ScaleZeroBody(double x)
{
   ZeroBody *= x;
}

void Operator::ScaleOneBody(double x)
{
   OneBody *= x;
}

void Operator::ScaleTwoBody(double x)
{
   TwoBody.Scale(x);
}

void Operator::Eye()
{
   ZeroBody = 1;
   OneBody.eye();
   TwoBody.Eye();
}


//***********************************************
/// Calculates the kinetic energy operator in the 
/// harmonic oscillator basis.
/// \f[ t_{ab} = \frac{1}{2}\hbar\omega
/// \delta_{\ell_a \ell_b} \delta_{j_aj_b} \delta_{t_{za}t_{zb}}
/// \left\{
/// \begin{array}{ll}
/// 2n_a + \ell_a + \frac{3}{2} &: n_a=n_b\\
/// \sqrt{n_{a}(n_{a}+\ell_a + \frac{1}{2})} &: n_a=n_b+1\\
/// \end{array} \right. \f]
//***********************************************
void Operator::CalculateKineticEnergy()
{
   OneBody.zeros();
   int norbits = modelspace->GetNumberOrbits();
   int A = modelspace->GetTargetMass();
   double hw = modelspace->GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace->GetOrbit(a);
      OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for (int b=a+1;b<norbits;++b)  // make this better once OneBodyChannel is implemented
      {
         Orbit & ob = modelspace->GetOrbit(b);
         if (oa.l == ob.l and oa.j2 == ob.j2 and oa.tz2 == ob.tz2)
         {
            if (oa.n == ob.n+1)
               OneBody(a,b) = 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
            else if (oa.n == ob.n-1)
               OneBody(a,b) = 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
            OneBody(b,a) = OneBody(a,b);
         }
      }
   }
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
/// This function is designed for use with comm222_phss() and so it takes in
/// two arrays of matrices, one for hp terms and one for ph terms.
void Operator::DoPandyaTransformation(vector<arma::mat>& TwoBody_CC_hp, vector<arma::mat>& TwoBody_CC_ph)
{
   // loop over cross-coupled channels
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      arma::uvec kets_ph = tbc_cc.GetKetIndex_ph();
      int nph_kets = kets_ph.n_rows;
      int J_cc = tbc_cc.J;

      TwoBody_CC_hp[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);
      TwoBody_CC_ph[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);

      // loop over cross-coupled ph bras <ac| in this channel
      for (int ibra=0; ibra<nph_kets; ++ibra)
      {
         Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
         int a = bra_cc.p;
         int b = bra_cc.q;
         Orbit & oa = modelspace->GetOrbit(a);
         Orbit & ob = modelspace->GetOrbit(b);
         double ja = oa.j2*0.5;
         double jb = ob.j2*0.5;

         // loop over cross-coupled kets |bd> in this channel
         // we go to 2*nKets to include |bd> and |db>
         for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
         {
            Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
            int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
            int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
            Orbit & oc = modelspace->GetOrbit(c);
            Orbit & od = modelspace->GetOrbit(d);
            double jc = oc.j2*0.5;
            double jd = od.j2*0.5;


            int jmin = max(abs(ja-jd),abs(jc-jb));
            int jmax = min(ja+jd,jc+jb);
            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
               if (abs(sixj) < 1e-8) continue;
               double tbme = TwoBody.GetTBME_J(J_std,a,d,c,b);
               sm -= (2*J_std+1) * sixj * tbme ;
            }
            TwoBody_CC_hp[ch_cc](ibra,iket_cc) = sm;


            // Exchange (a <-> b) to account for the (n_a - n_b) term
            // Get Tz,parity and range of J for <bd || ca > coupling
            jmin = max(abs(jb-jd),abs(jc-ja));
            jmax = min(jb+jd,jc+ja);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(jb,ja,J_cc,jc,jd,J_std);
               if (abs(sixj) < 1e-8) continue;
               double tbme = TwoBody.GetTBME_J(J_std,b,d,c,a);
               sm -= (2*J_std+1) * sixj * tbme ;
            }
            TwoBody_CC_ph[ch_cc](ibra,iket_cc) = sm;

         }
      }
   }
}







//*****************************************************************************************
//  I think this explanation needs to be revisited...
//   Transform to a cross-coupled basis, for use in the 2body commutator relation
//    using a Pandya tranformation 
//                                               ____________________________________________________________________________
//                                              |                                                                            |
//   |a    |b  <ab|_J'    \a   /c  <ac|_J       |  <ac|V|bd>_J = sum_J' (2J'+1) (-1)^(ja+jd+J+J'+1) { ja jc J } <ab|V|cd>_J' |
//   |     |               \  /                 |                                                   { jd jb J'}              |
//   |_____|                \/_____             |____________________________________________________________________________|
//   |  V  |      ===>         V  /\ 
//   |     |                     /  \
//   |c    |d  |cd>_J'          /b   \d   |bd>_J
//
//   STANDARD                 CROSS  
//   COUPLING                 COUPLED  
//                                      
void Operator::CalculateCrossCoupled(vector<arma::mat> &TwoBody_CC_left, vector<arma::mat> &TwoBody_CC_right)
{
   // loop over cross-coupled channels
   #pragma omp parallel for schedule(dynamic,1)
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      arma::uvec kets_ph = tbc_cc.GetKetIndex_ph();
      int nph_kets = kets_ph.n_rows;
      int J_cc = tbc_cc.J;

//   These matrices don't actually need to be square, since we only care about
//   the particle-hole kets, which we sum over via matrix multiplication:
//   [ ]                  [     ]
//   |L|  x  [  R   ]  =  |  N  |
//   [ ]                  [     ]
//
      TwoBody_CC_left[ch_cc]  = arma::mat(2*nKets_cc, nph_kets,   arma::fill::zeros);
      TwoBody_CC_right[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);

      // loop over cross-coupled ph bras <ac| in this channel
      for (int i_ph=0; i_ph<nph_kets; ++i_ph)
      {
         Ket & bra_cc = tbc_cc.GetKet( kets_ph[i_ph] );
         int a = bra_cc.p;
         int b = bra_cc.q;
         Orbit & oa = modelspace->GetOrbit(a);
         Orbit & ob = modelspace->GetOrbit(b);
         double ja = oa.j2/2.0;
         double jb = ob.j2/2.0;

         // loop over cross-coupled kets |bd> in this channel
         // we go to 2*nKets to include |bd> and |db>
         for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
         {
            Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
            int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
            int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
            Orbit & oc = modelspace->GetOrbit(c);
            Orbit & od = modelspace->GetOrbit(d);
            double jc = oc.j2/2.0;
            double jd = od.j2/2.0;

            int phase_ad = modelspace->phase(ja+jb);

            // Get Tz,parity and range of J for <ca || bd > coupling
            int Tz_std = (oa.tz2 + oc.tz2)/2;
            int parity_std = (oa.l + oc.l)%2;
            int jmin = max(abs(ja-jc),abs(jb-jd));
            int jmax = min(ja+jc,jb+jd);
            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jb,J_cc,jd,jc,J_std);
               if (abs(sixj)<1e-8) continue;
               double tbme = TwoBody.GetTBME(J_std,parity_std,Tz_std,a,c,d,b);
               sm += (2*J_std+1) * sixj * tbme ; 
            }
            TwoBody_CC_left[ch_cc](iket_cc,i_ph) = sm * phase_ad;


            // Get Tz,parity and range of J for <bc || da > coupling
            Tz_std = (oa.tz2 + od.tz2)/2;
            parity_std = (oa.l + od.l)%2;
            jmin = max(abs(ja-jd),abs(jc-jb));
            jmax = min(ja+jd,jb+jc);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
               if (abs(sixj) < 1e-8) continue;
               double tbme = TwoBody.GetTBME(J_std,parity_std,Tz_std,b,c,d,a);
               sm -= (2*J_std+1) * sixj * tbme ;
            }
            TwoBody_CC_right[ch_cc](i_ph,iket_cc) = sm;

         }
      }
   }
}




//*****************************************************************************************
/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [X,Y] + \frac{1}{2!}[X,[X,Y]] + \frac{1}{3!}[X,[X,[X,Y]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Operator::BCH_Transform(  Operator &Omega)
{
   double t = omp_get_wtime();
   int max_iter = 20;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   Operator OpNested = *this;
   for (int i=1; i<max_iter; ++i)
   {
      OpNested = Omega.Commutator(OpNested) / i;

      OpOut += OpNested;

      if (OpNested.Norm() < (nx+ny)*bch_transform_threshold)
      {
        timer["BCH_Transform"] += omp_get_wtime() - t;
        return OpOut;
      }
      if (i == warn_iter)
      {
         cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << endl;
      }

   }
   cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << endl;
   timer["BCH_Transform"] += omp_get_wtime() - t;
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
Operator Operator::BCH_Product(  Operator &Y)
{
   double t = omp_get_wtime();
   Operator& X = *this;

//   Operator Z = X + Y; 
//   Operator XY = X.Commutator(Y);
   Operator Z = 0.5*(X.Commutator(Y));
//   Z += XY*(1./2);    // [X,Y]
   double nx = X.Norm();
   double ny = Y.Norm();
//   double nc1 = XY.Norm();
   double nxy = Z.Norm();

   if ( nxy < (nx+ny)*bch_product_threshold )
   {
     Z += X;
     Z += Y;
     timer["BCH_Product"] += omp_get_wtime() - t;
     return Z;
   }
//   Operator XYDiff = Y-X;
//   Z += (1./6)* Z.Commutator( XYDiff );
   Y -= X;
   Z += (1./6)* Z.Commutator( Y );
   Z += Y;
   Z += 2*X;

   timer["BCH_Product"] += omp_get_wtime() - t;
   return Z;

//   Operator YYX = XY.Commutator(Y); // [[X,Y],Y] = [Y,[Y,X]]
//   double nc2 = YYX.Norm();
//   Z += YYX * (1/12.); 
//
//   if ( nc2/12 < (nx+ny)*bch_product_threshold ) return Z;
//
//   Operator XXY = X.Commutator(XY); // [X,[X,Y]]
//   double nc3 = XXY.Norm();
//   Z += XXY * (1/12.);      // [X,[X,Y]]

//   if ( nc3/12 < (nx+ny)*bch_product_threshold ) return Z;
//
//   cout << "Warning: BCH product expansion not converged after 3 nested commutators!" << endl;
//
//   Operator YXXY = Y.Commutator(XXY); // [Y,[X,[X,Y]]]
//   double nc4 = YXXY.Norm();
//   Z += YXXY*(-1./24);
//
//   Operator YYYX = Y.Commutator(YYX) ; // [Y,[Y,[Y,X]]]
//   Operator YYYYX = Y.Commutator(YYYX) ; // [Y,[Y,[Y,[Y,X]]]]
//   double nc5 = YYYYX.Norm();
//   Operator XXXY =  X.Commutator(XXY) ; // [X,[X,[X,[X,Y]]]]
//   Operator XXXXY = X.Commutator(XXXY) ; // [X,[X,[X,[X,Y]]]]
//   double nc6 = XXXXY.Norm();
//   Z += (YYYYX + XXXXY)*(-1./720);
//
//   if ( nc6/720. < (nx+ny)*bch_product_threshold ) return Z;
//   cout << "Warning: BCH product expansion not converged after 5 nested commutators!" << endl;

//   return Z;
}

/// Obtain the Frobenius norm of the operator, which here is 
/// defined as 
/// \f[ \|X\| = \sqrt{\|X_{(1)}\|^2 +\|X_{(2)}\|^2 } \f]
/// and
/// \f[ \|X_{(1)}\|^2 = \sum\limits_{ij} X_{ij}^2 \f]
double Operator::Norm() const
{
   double nrm = 0;
   double n1 = OneBodyNorm();
   double n2 = TwoBody.Norm();
   return sqrt(n1*n1+n2*n2);
}

double Operator::OneBodyNorm() const
{
   return arma::norm(OneBody,"fro");
}



double Operator::TwoBodyNorm() const
{
  return TwoBody.Norm();
}

void Operator::Symmetrize()
{
   OneBody = arma::symmatu(OneBody);
   TwoBody.Symmetrize();
}

void Operator::AntiSymmetrize()
{
   int norb = modelspace->GetNumberOrbits();
   for (int i=0;i<norb;++i)
   {
      for(int j=i+1;j<norb;++j)
      {
        OneBody(j,i) = -OneBody(i,j);
      }
   }
   TwoBody.AntiSymmetrize();
}

Operator Operator::Commutator( Operator& opright)
{
   timer["N_Commutators"] += 1;
   if (rank_J==0)
   {
      if (opright.rank_J==0)
      {
         return CommutatorScalarScalar(opright); // [S,S]
      }
      else
      {
         return CommutatorScalarTensor(opright); // [S,T]
      }
   }
   else if(opright.rank_J==0)
   {
      return (-1)*opright.CommutatorScalarTensor(*this); // [T,S]
   }
   else
   {
      cout << "In Tensor-Tensor because rank_J = " << rank_J << "  and opright.rank_J = " << opright.rank_J << endl;
      cout << " Tensor-Tensor commutator not yet implemented." << endl;
      return *this;
   }
}


Operator Operator::CommutatorScalarScalar( Operator& opright) 
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   if ( not out.IsAntiHermitian() )
   {
      comm110ss(opright, out);
      if (particle_rank>1 and opright.particle_rank>1)
        comm220ss(opright, out) ;
   }

    double t = omp_get_wtime();
   comm111ss(opright, out);
    timer["comm111ss"] += omp_get_wtime() - t;

    t = omp_get_wtime();
//   comm111st(opright, out);  // << equivalent in scalar case
   comm121ss(opright, out);
//   comm121st(opright, out);  // << equivalent in scalar case
    timer["comm121ss"] += omp_get_wtime() - t;

    t = omp_get_wtime();
   comm122ss(opright, out); //  This is the slow one for some reason.
    timer["comm122ss"] += omp_get_wtime() - t;

   if (particle_rank>1 and opright.particle_rank>1)
   {
    t = omp_get_wtime();
    comm222_pp_hh_221ss(opright, out);
    timer["comm222_pp_hh_221ss"] += omp_get_wtime() - t;
     
////   comm222_pp_hh_221st(opright, out); // << equivalent in scalar case

    t = omp_get_wtime();
    comm222_phss(opright, out);
//   comm222_phst_pandya(opright, out);
////   comm222_phst(opright, out);
    timer["comm222_phss"] += omp_get_wtime() - t;
   }


   if ( out.IsHermitian() )
   {
      out.Symmetrize();
   }
   else if (out.IsAntiHermitian() )
   {
      out.AntiSymmetrize();
   }


   return out;
}


// Calculate [S,T]
Operator Operator::CommutatorScalarTensor( Operator& opright) 
{
   Operator out = opright; // This ensures the commutator has the same tensor rank as opright
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   comm111st(opright, out);
   comm121st(opright, out);

   comm122st(opright, out);
   comm222_pp_hh_221st(opright, out);
   comm222_phst(opright, out);

   return out;
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
void Operator::comm110ss( Operator& opright, Operator& out) 
{
  if (IsHermitian() and opright.IsHermitian()) return ; // I think this is the case
  if (IsAntiHermitian() and opright.IsAntiHermitian()) return ; // I think this is the case

   int norbits = modelspace->GetNumberOrbits();
   arma::mat xyyx = OneBody*opright.OneBody - opright.OneBody*OneBody;
   for ( auto& a : modelspace->holes) 
   {
      out.ZeroBody += (modelspace->GetOrbit(a).j2+1) * xyyx(a,a);
   }
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
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) X_{abcd}^{J} Y_{cdab}^{J}
/// \f]
/// may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{ab} (\mathcal{P}_{hh} X_{(2)}^{J} \mathcal{P}_{pp} Y_{(2)}^{J})_{abab}
/// \f] where \f$ \mathcal{P}_{hh} \f$ is a projector onto hole-hole two body states.
void Operator::comm220ss(  Operator& opright, Operator& out) 
{
   if (IsHermitian() and opright.IsHermitian()) return; // I think this is the case
   if (IsAntiHermitian() and opright.IsAntiHermitian()) return; // I think this is the case

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto hh = tbc.GetKetIndex_hh();
      auto pp = tbc.GetKetIndex_pp();
      out.ZeroBody += 2 * (2*tbc.J+1) * arma::trace( TwoBody.GetMatrix(ch).submat(hh,pp) * opright.TwoBody.GetMatrix(ch).submat(pp,hh) );
//      out.ZeroBody += 2 * (2*tbc.J+1) * arma::trace( TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(),tbc.GetKetIndex_pp()) * opright.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(),tbc.GetKetIndex_hh()) );

   }
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
void Operator::comm111ss( Operator & opright, Operator& out) 
{
   out.OneBody += OneBody*opright.OneBody - opright.OneBody*OneBody;
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
//                                                      * y_ab xbiag - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
/// Returns \f$ [X_{(1)},Y_{(2)}] - [Y_{(1)},X_{(2)}] \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]_{ij} = \frac{1}{2j_i+1}\sum_{ab} (n_a \bar{n}_b) \sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj})
/// \f]
void Operator::comm121ss( Operator& opright, Operator& out) 
{
   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for 
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      int jmin = out.IsNonHermitian() ? 0 : i;
      for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (auto& b : modelspace->particles)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                out.OneBody(i,j) += (ob.j2+1) *  OneBody(a,b) * opright.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                out.OneBody(i,j) -= (oa.j2+1) *  OneBody(b,a) * opright.TwoBody.GetTBMEmonopole(a,i,b,j) ;

                // comm211 part
                out.OneBody(i,j) -= (ob.j2+1) *  opright.OneBody(a,b) * TwoBody.GetTBMEmonopole(b,i,a,j) ;
                out.OneBody(i,j) += (oa.j2+1) *  opright.OneBody(b,a) * TwoBody.GetTBMEmonopole(a,i,b,j) ;
             }
          }
      }
   }
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
void Operator::comm221ss( Operator& opright, Operator& out) 
{

   int norbits = modelspace->GetNumberOrbits();
   TwoBodyME Mpp = opright.TwoBody;
   TwoBodyME Mhh = opright.TwoBody;

   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = (arma::mat&) TwoBody.GetMatrix(ch,ch);
      auto& RHS = (arma::mat&) opright.TwoBody.GetMatrix(ch,ch);
      arma::mat& Matrixpp =  Mpp.GetMatrix(ch,ch);
      arma::mat& Matrixhh =  Mpp.GetMatrix(ch,ch);
      
      Matrixpp = ( LHS.rows(tbc.GetKetIndex_pp()) * RHS.cols(tbc.GetKetIndex_pp()));
      Matrixpp -= Matrixpp.t();
      Matrixhh = ( LHS.rows(tbc.GetKetIndex_hh()) * RHS.cols(tbc.GetKetIndex_hh()));
      Matrixhh -= Matrixhh.t();

      // If commutator is hermitian or antihermitian, we only
      // need to do half the sum. Add this.
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
         {
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,i,c,j,c);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,i,c,j,c);
            }
            cijJ *= (2*tbc.J+1.0)/(oi.j2 +1.0);
            #pragma omp critical
            out.OneBody(i,j) +=  cijJ;
         } // for j
      } // for i
   } //for ch
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
// Right now, this is the slowest one...
/// Returns \f$ [X_{(1)},Y_{(2)}]_{(2)} - [Y_{(1)},X_{(2)}]_{(2)} \f$, where
/// \f[
/// [X_{(1)},Y_{(2)}]^{J}_{ijkl} = \sum_{a} ( X_{ia}Y^{J}_{ajkl} + X_{ja}Y^{J}_{iakl} - X_{ak} Y^{J}_{ijal} - X_{al} Y^{J}_{ijka} )
/// \f]
void Operator::comm122ss( Operator& opright, Operator& opout ) 
{
   auto& L1 = OneBody;
   auto& R1 = opright.OneBody;

   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& L2 = TwoBody.GetMatrix(ch,ch);
      auto& R2 = opright.TwoBody.GetMatrix(ch,ch);
      arma::mat& OUT = opout.TwoBody.GetMatrix(ch,ch);


      int npq = tbc.GetNumberKets();
      int norbits = modelspace->GetNumberOrbits();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         double pre_ij = i==j ? SQRT2 : 1;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         arma::Row<double> L2_ij = L2.row(indx_ij); // trying this to better use the cache. not sure if it helps.
         arma::Row<double> R2_ij = R2.row(indx_ij);
         int klmin = opout.IsNonHermitian() ? 0 : indx_ij;
         for (int indx_kl=klmin;indx_kl<npq; ++indx_kl)
         {
            Ket & ket = tbc.GetKet(indx_kl);
            int k = ket.p;
            int l = ket.q;
            double pre_kl = k==l ? SQRT2 : 1;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            arma::vec L2_kl = L2.unsafe_col(indx_kl);
            arma::vec R2_kl = R2.unsafe_col(indx_kl);

            double cijkl = 0;


            for (int a : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
                 int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
                 if (indx_aj < 0) continue;
                 double pre_aj = a>j ? tbc.GetKet(indx_aj).Phase(tbc.J) : (a==j ? SQRT2 : 1);
                 cijkl += pre_kl * pre_aj  * ( L1(i,a) * R2_kl(indx_aj) - R1(i,a) * L2_kl(indx_aj) );
            }

            for (int a : modelspace->OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
            {
                 int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
                 if (indx_ia < 0) continue;
                 double pre_ia = i>a ? tbc.GetKet(indx_ia).Phase(tbc.J) : (i==a ? SQRT2 : 1);
                 cijkl += pre_kl * pre_ia * ( L1(j,a) * R2_kl(indx_ia) - R1(j,a) * L2_kl(indx_ia) );
             }

            for (int a : modelspace->OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
                int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
                if (indx_al < 0) continue;
                double pre_al = a>l ? tbc.GetKet(indx_al).Phase(tbc.J) : (a==l ? SQRT2 : 1);
                cijkl += pre_ij * pre_al * ( R1(a,k) * L2_ij(indx_al) - L1(a,k) * R2_ij(indx_al) );
            }

            for (int a : modelspace->OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
            {
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));
               if (indx_ka < 0) continue;
               double pre_ka = k>a ? tbc.GetKet(indx_ka).Phase(tbc.J) : (k==a ? SQRT2 : 1);
               cijkl += pre_ij * pre_ka * ( R1(a,l) * L2_ij(indx_ka) - L1(a,l) * R2_ij(indx_ka) );
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            OUT(indx_ij,indx_kl) += cijkl / norm;
         }
      }
   }

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
void Operator::comm222_pp_hhss( Operator& opright, Operator& opout ) 
{
//   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& LHS = (arma::mat&) TwoBody.GetMatrix(ch,ch);
      arma::mat& RHS = (arma::mat&) opright.TwoBody.GetMatrix(ch,ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.GetMatrix(ch,ch);

      arma::mat Mpp = (LHS.rows(tbc.GetKetIndex_pp()) * RHS.cols(tbc.GetKetIndex_pp()));
      arma::mat Mhh = (LHS.rows(tbc.GetKetIndex_hh()) * RHS.cols(tbc.GetKetIndex_hh()));
      OUT += Mpp - Mpp.t() - Mhh + Mhh.t();
   }
}







/// Since comm222_pp_hhss() and comm221ss() both require the ruction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
void Operator::comm222_pp_hh_221ss( Operator& opright, Operator& opout )  
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();

   TwoBodyME Mpp = opout.TwoBody;
   TwoBodyME Mhh = opout.TwoBody;

   double t = omp_get_wtime();
   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = TwoBody.GetMatrix(ch,ch);
      auto& RHS = opright.TwoBody.GetMatrix(ch,ch);
      arma::mat& OUT =  opout.TwoBody.GetMatrix(ch,ch);

      arma::mat & Matrixpp = Mpp.GetMatrix(ch,ch);
      arma::mat & Matrixhh = Mhh.GetMatrix(ch,ch);

      arma::uvec kets_pp = tbc.GetKetIndex_pp();
      arma::uvec kets_hh = tbc.GetKetIndex_hh();
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * RHS.rows(kets_hh);

      if (opout.IsHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (opout.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * LHS.rows(kets_hh);
      }

      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
   timer["pphh TwoBody bit"] += omp_get_wtime() - t;

   t = omp_get_wtime();
   // The one body part
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      int jmin = opout.IsNonHermitian() ? 0 : i;
      for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         if (j<jmin) continue;
         double cijJ = 0;
         for (int ch=0;ch<nChannels;++ch)
         {
            TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
            double Jfactor = (2*tbc.J+1.0);
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               cijJ += Mpp.GetTBME(ch,c,i,c,j) * Jfactor;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ += Mhh.GetTBME(ch,c,i,c,j) * Jfactor;
            }
         }
         opout.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   timer["pphh One Body bit"] += omp_get_wtime() - t;
}






void Operator::InversePandyaTransformation(vector<arma::mat>& W, vector<arma::mat>& opout, bool hermitian)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();

      opout[ch] = arma::mat(nKets,nKets,arma::fill::zeros);
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         int ketmin = hermitian ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            int parity_cc = (oi.l+ol.l)%2;
            int Tz_cc = abs(oi.tz2+ol.tz2)/2;
            int jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            int jmax = min(int(ji+jl),int(jk+jj));

            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (abs(sixj)<1e-8) continue;
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc = modelspace->GetTwoBodyChannel_CC(ch_cc);
               int indx_il = tbc.GetLocalIndex(min(i,l),max(i,l));
               int indx_kj = tbc.GetLocalIndex(min(j,k),max(j,k));
               if (i>l) indx_il += tbc.GetNumberKets();
               if (k>j) indx_kj += tbc.GetNumberKets();
               double me1 = W[ch_cc](indx_il,indx_kj);
               commij += (2*Jprime+1) * sixj * me1;
            }

            // now loop over the cross coupled TBME's
            parity_cc = (oi.l+ok.l)%2;
            Tz_cc = abs(oi.tz2+ok.tz2)/2;
            jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
               if (abs(sixj)<1e-8) continue;
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc = modelspace->GetTwoBodyChannel_CC(ch_cc);
               int indx_jl = tbc.GetLocalIndex(min(j,l),max(j,l));
               int indx_ki = tbc.GetLocalIndex(min(i,k),max(i,k));
               if (j>l) indx_jl += tbc.GetNumberKets();
               if (k>i) indx_ki += tbc.GetNumberKets();
               double me1 = W[ch_cc](indx_jl,indx_ki);
               commji += (2*Jprime+1) *  sixj * me1;
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            opout[ch](ibra,iket) = (commij - modelspace->phase(ji+jj-J)*commji) / norm;;
         }
      }
   }
 
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
// This version "works" with Pandya transform
/// Calculates the part of \f$ [X_{(2)},Y_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ Z^{J}_{ijkl} \f$
/// \f[
/// Z^{J}_{ijkl} = \sum_{ab}(n_a-n_b)\sum_{J'} (2J'+1)
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
/// \bar{W}^{J}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J'}_{i\bar{l}a\bar{b}}\bar{Y}^{J'}_{a\bar{b}k\bar{j}} - 
///   \bar{Y}^{J'}_{i\bar{l}a\bar{b}}\bar{X}^{J'}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J'}_{i\bar{l}b\bar{a}}\bar{Y}^{J'}_{b\bar{a}k\bar{j}} - 
///    \bar{Y}^{J'}_{i\bar{l}b\bar{a}}\bar{X}^{J'}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoPandyaTransformation().
/// The commutator is then given by
/// \f[
/// Z^{J}_{ijkl} = \sum_{J'} (2J'+1)
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_j  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{W}^{J'}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_i  &  J \\
///  j_k  &  j_l  &  J' \\
///  \end{array} \right\}
///  \bar{W}^{J'}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
void Operator::comm222_phss( Operator& opright, Operator& opout ) 
{

   // Create Pandya-transformed hp and ph matrix elements
   vector<arma::mat> X_bar_hp (nChannels );
   vector<arma::mat> X_bar_ph (nChannels );

   vector<arma::mat> Y_bar_hp (nChannels );
   vector<arma::mat> Y_bar_ph (nChannels );

   double t = omp_get_wtime();
   DoPandyaTransformation(X_bar_hp, X_bar_ph );
   opright.DoPandyaTransformation(Y_bar_hp, Y_bar_ph );
   timer["DoPandyaTransformation"] += omp_get_wtime() - t;

   t = omp_get_wtime();
   // Construct the intermediate matrix W_bar
   vector<arma::mat> W_bar (nChannels );
   int hx = IsHermitian() ? 1 : -1;
   int hy = opright.IsHermitian() ? 1 : -1;
   // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
   for (int ch=0;ch<nChannels;++ch)
   {
      W_bar[ch] =  hx*( X_bar_hp[ch].t() * Y_bar_hp[ch] - X_bar_ph[ch].t() * Y_bar_ph[ch]) ;
      if (hx*hy<0)
        W_bar[ch] += W_bar[ch].t();
      else
        W_bar[ch] -= W_bar[ch].t();
   }
   timer["Build W_bar"] += omp_get_wtime() - t;

   vector<arma::mat> W (nChannels );
   t = omp_get_wtime();
   InversePandyaTransformation(W_bar, W, opout.IsHermitian());
   timer["InversePandyaTransformation"] += omp_get_wtime() - t;
   for (int ch=0; ch<nChannels; ++ch)
   {
     opout.TwoBody.GetMatrix(ch,ch) += W[ch];
   }

}



/*
// This version works with Cross Coupling
void Operator::comm222_phss( Operator& opright, Operator& opout ) 
{

   // Update Cross-coupled matrix elements
   vector<arma::mat> X_TwoBody_CC_left (nChannels );
   vector<arma::mat> X_TwoBody_CC_right (nChannels );

   vector<arma::mat> Y_TwoBody_CC_left (nChannels );
   vector<arma::mat> Y_TwoBody_CC_right (nChannels );

   double t = omp_get_wtime();
   CalculateCrossCoupled(X_TwoBody_CC_left, X_TwoBody_CC_right );
   opright.CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   t = omp_get_wtime() - t;
   timer["CalculateCrossCoupled"] += t;

   // Construct the intermediate matrices N1 and N2
   vector<arma::mat> N1 (nChannels, arma::mat() );
 // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] = X_TwoBody_CC_left[ch] * Y_TwoBody_CC_right[ch] - Y_TwoBody_CC_left[ch] * X_TwoBody_CC_right[ch];
      N1[ch] += N1[ch].t(); // maybe this works?  Yes !
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.GetMatrix(ch,ch);
      int J = tbc.J;

      int nKets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2*0.5;
         double jj = oj.j2*0.5;
         int ketmin = opout.IsHermitian() ? ibra : ( opout.IsAntiHermitian() ? ibra+1 : 0);
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2*0.5;
            double jl = ol.j2*0.5;


            double commij = 0;
            double commji = 0;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi.l+ok.l)%2;
            int Tz_cc = abs(oi.tz2+ok.tz2)/2;
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
               if (abs(sixj)<1e-8) continue;
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc = modelspace->GetTwoBodyChannel_CC(ch_cc);
               int indx_ik = tbc.GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = tbc.GetLocalIndex(min(j,l),max(j,l));
               if (i>k) indx_ik += tbc.GetNumberKets();
               if (j>l) indx_jl += tbc.GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               commji += (2*Jprime+1) * sixj * (me1);
            }

            parity_cc = (oi.l+ol.l)%2;
            Tz_cc = abs(oi.tz2+ol.tz2)/2;
            jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (abs(sixj)<1e-8) continue;
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc = modelspace->GetTwoBodyChannel_CC(ch_cc);
               int indx_il = tbc.GetLocalIndex(min(i,l),max(i,l));
               int indx_jk = tbc.GetLocalIndex(min(j,k),max(j,k));
               if (i>l)  indx_il += tbc.GetNumberKets();
               if (j>k)  indx_jk += tbc.GetNumberKets();
               double me1 = N1[ch_cc](indx_il,indx_jk);
               commij += (2*Jprime+1) * sixj * (me1);
            }

            double comm = commij*modelspace->phase(ji+jl) + commji*modelspace->phase(ji+jl+J);

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            comm /= norm;
            OUT(ibra,iket) += comm;
         }
      }
   }
}
*/




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////   BEGIN SCALAR-TENSOR COMMUTATORS      //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// This is no different from the scalar-scalar version
void Operator::comm111st( Operator & opright, Operator& out) 
{
   out.OneBody += OneBody*opright.OneBody - opright.OneBody*OneBody;
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
// X is scalar, Y is tensor
// 
void Operator::comm121st( Operator& opright, Operator& out) 
{
   int norbits = modelspace->GetNumberOrbits();
   int nch = modelspace->GetNumberTwoBodyChannels();
   int Lambda = opright.rank_J;
   #pragma omp parallel for
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          Orbit &oj = modelspace->GetOrbit(j);
          double jj = oj.j2/2.0;
          if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue; // At some point, make a OneBodyChannel class...

          double zij = 0;
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             double ja = oa.j2/2.0;
             for (auto& b : modelspace->particles)  // C++11 syntax
             {
                Orbit &ob = modelspace->GetOrbit(b);
                double jb = ob.j2/2.0;
                for ( auto& itmat : TwoBody.MatEl )
                {
                  int ch_bra = itmat.first[0];
                  int ch_ket = itmat.first[1];
//                for (int ch_bra=0;ch_bra<nch;++ch_bra)
//                {
                   TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
                   TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
                   double J = tbc_bra.J;

                   // part with opleft two body and opright one body, i.e. [X(2),Y(1)]  ==  -[Y(1),X(2)]
//                   {
                      double sixj_ab = modelspace->GetSixJ(jj, ja, J, jb, ji, Lambda);
                      double sixj_ba = modelspace->GetSixJ(jj, jb, J, ja, ji, Lambda);
                      int phase_ab = modelspace->phase(ji+jb+J);
                      int phase_ba = modelspace->phase(ji+ja+J);
                      double hatfactor_ab = (2*J+1)*sqrt((oa.j2+1));
                      double hatfactor_ba = (2*J+1)*sqrt((ob.j2+1));
   
//                      zij -= opright.OneBody(a,b)*GetTBME(ch_bra,ch_ket,b,i,a,j) * phase_ab * sixj_ab * hatfactor_ab;
//                      zij += opright.OneBody(b,a)*GetTBME(ch_bra,ch_ket,a,i,b,j) * phase_ba * sixj_ba * hatfactor_ba;
                      zij -= opright.OneBody(a,b)*TwoBody.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase_ab * sixj_ab * hatfactor_ab;
                      zij += opright.OneBody(b,a)*TwoBody.GetTBME(ch_bra,ch_ket,a,i,b,j) * phase_ba * sixj_ba * hatfactor_ba;
//                   }

                   // part with opleft one body and opright two body, i.e. [X(1),Y(2)]
                   if (oa.j2 != ob.j2) continue;  // X(1) is scalar -> ja==jb
//                   {
//                      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
                      double Jprime = tbc_ket.J;
                      double sixj = modelspace->GetSixJ(Jprime, J, Lambda, ji, jj, ja);
                      int phase = modelspace->phase(jj+ja+J);
                      double hatfactor = (2*J+1)*sqrt((2*Jprime+1));

//                      zij += OneBody(a,b)*opright.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase * sixj * hatfactor;
//                      zij -= OneBody(b,a)*opright.GetTBME(ch_bra,ch_ket,a,i,b,j) * phase * sixj * hatfactor; // Here hatfactor_ab == hatfactor_ba
                      zij += OneBody(a,b)*opright.TwoBody.GetTBME(ch_bra,ch_ket,b,i,a,j) * phase * sixj * hatfactor;
                      zij -= OneBody(b,a)*opright.TwoBody.GetTBME(ch_bra,ch_ket,a,i,b,j) * phase * sixj * hatfactor; // Here hatfactor_ab == hatfactor_ba
//                   }
                }
             }
          }
          out.OneBody(i,j) += zij/sqrt(oi.j2+1);
      }
   }
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
void Operator::comm221st( Operator& opright, Operator& out) 
{

   int norbits = modelspace->GetNumberOrbits();
//   Operator Mpp = opright;
//   Operator Mhh = opright;
   TwoBodyME Mpp = opright.TwoBody;
   TwoBodyME Mhh = opright.TwoBody;

//   #pragma omp parallel for schedule(dynamic,5)
   for ( auto& itmat : opright.TwoBody.MatEl )
   {
//   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
//   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int npq = tbc_bra.GetNumberKets();
//      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);

//         arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_ket);
//         arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
//         arma::mat& Matrixpp = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);
//         arma::mat& Matrixhh = (arma::mat&) Mpp.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& LHS = (arma::mat&) TwoBody.GetMatrix(ch_bra,ch_ket);
         arma::mat& RHS = (arma::mat&) opright.TwoBody.GetMatrix(ch_bra,ch_ket);
         arma::mat& Matrixpp = (arma::mat&) Mpp.GetMatrix(ch_bra,ch_ket);
         arma::mat& Matrixhh = (arma::mat&) Mpp.GetMatrix(ch_bra,ch_ket);

         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);
         Matrixpp = (LHS.rows(tbc_bra.GetKetIndex_pp()) * RHS.cols(tbc_bra.GetKetIndex_pp()));
         Matrixpp -= Matrixpp.t();
         Matrixhh = (LHS.rows(tbc_bra.GetKetIndex_hh()) * RHS.cols(tbc_bra.GetKetIndex_hh()));
         Matrixhh -= Matrixhh.t();
//         Matrixpp = (LHS * tbc_bra.Proj_pp * RHS  -  RHS * tbc_bra.Proj_pp * LHS);
//         Matrixhh = (LHS * tbc_bra.Proj_hh * RHS  -  RHS * tbc_bra.Proj_hh * LHS);
      

      // If commutator is hermitian or antihermitian, we only
      // need to do half the sum. I should add this in later.
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         for (int j=0;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mpp.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               int phase = modelspace->phase((oj.j2+oc.j2)/2+J2);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ +=   Mhh.GetTBME(ch_bra,ch_ket,i,c,j,c)*phase*sixj;
            }
            #pragma omp critical
            out.OneBody(i,j) += cijJ / (oi.j2 +1.0);
         } // for j
      } // for i
//   } //for ch_ket
   } //for ch_bra
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
// Right now, this is the slowest one...
// Agrees with previous code in the scalar-scalar limit
void Operator::comm122st( Operator& opright, Operator& opout ) 
{
   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();
   int Lambda = opout.rank_J;

//   #pragma omp parallel for schedule(dynamic,5)
    for ( auto& itmat : opright.TwoBody.MatEl )
    {
     int ch_bra = itmat.first[0];
     int ch_ket = itmat.first[1];
//   for (int ch_bra=0; ch_bra<nChannels; ++ch_bra)
//   {

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int J = tbc_bra.J;
      int npq = tbc_bra.GetNumberKets();



//      {
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int Jprime = tbc_ket.J;
//         arma::mat& LHS = (arma::mat&) TwoBody.at(ch_bra).at(ch_bra);
//         arma::mat& RHS = (arma::mat&) opright.TwoBody.at(ch_bra).at(ch_ket);
//         arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch_bra).at(ch_ket);
         arma::mat& OUT = (arma::mat&) opout.TwoBody.GetMatrix(ch_bra,ch_ket);

      for (int ibra = 0;ibra<npq; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.0;
         double jj = oj.j2/2.0;
         for (int iket=ibra;iket<npq; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.0;
            double jl = ol.j2/2.0;

            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
              Orbit& oa = modelspace->GetOrbit(a);
              double ja = oa.j2/2.0;

               cijkl += OneBody(i,a) * opright.TwoBody.GetTBME(ch_bra,ch_ket,a,j,k,l);
               cijkl += OneBody(j,a) * opright.TwoBody.GetTBME(ch_bra,ch_ket,i,a,k,l);
               cijkl -= OneBody(a,k) * opright.TwoBody.GetTBME(ch_bra,ch_ket,i,j,a,l);
               cijkl -= OneBody(a,l) * opright.TwoBody.GetTBME(ch_bra,ch_ket,i,j,k,a);


               double sixj_ia = modelspace->GetSixJ(J,Jprime,Lambda,ja,ji,jj);
               int phase_ia = modelspace->phase(ji+jj+J);
               cijkl -= opright.OneBody(i,a) * TwoBody.GetTBME(ch_bra,ch_bra,a,j,k,l) * sixj_ia * phase_ia * sqrt((oi.j2+1)*(2*Jprime+1));

               double sixj_ja = modelspace->GetSixJ(J,Jprime,Lambda,ja,jj,ji);
               int phase_ja = modelspace->phase(ji+jj+J);
               cijkl -= opright.OneBody(j,a) * TwoBody.GetTBME(ch_bra,ch_bra,i,a,k,l) * sixj_ja * phase_ja * sqrt((oj.j2+1)*(2*Jprime+1));

               double sixj_ak = modelspace->GetSixJ(J,Jprime,Lambda,jk,ja,jl);
               int phase_ak = modelspace->phase(ja+jl+J+Lambda);
               cijkl += opright.OneBody(a,k) * TwoBody.GetTBME(ch_bra,ch_bra,i,j,a,l) * sixj_ak * phase_ak * sqrt((oa.j2+1)*(2*Jprime+1));

               double sixj_al = modelspace->GetSixJ(J,Jprime,Lambda,jl,ja,jk);
               int phase_al = modelspace->phase(ja+jk+J+Lambda);
               cijkl += opright.OneBody(a,l) * TwoBody.GetTBME(ch_bra,ch_bra,i,j,k,a) * sixj_al * phase_al * sqrt((oa.j2+1)*(2*Jprime+1));

            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            cijkl /= norm;
            #pragma omp critical
            {
            OUT(ibra,iket) += cijkl;
            if (ibra != iket)
               OUT(iket,ibra) += herm * cijkl;
            }
         }
      }
  //    }
   }
}





//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |              ( note that   1-n_a-n_b  =  nbar_a nbar_b - n_an_b )
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm222_pp_hhst( Operator& opright, Operator& opout ) 
{
//   #pragma omp parallel for schedule(dynamic,5)
   for ( auto & itmat : opright.TwoBody.MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      int npq = tbc_bra.GetNumberKets();
//      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);


         arma::mat& LHS = (arma::mat&) TwoBody.GetMatrix(ch_bra,ch_bra);
         arma::mat& RHS = (arma::mat&) opright.TwoBody.GetMatrix(ch_bra,ch_ket);
         arma::mat& OUT = (arma::mat&) opout.TwoBody.GetMatrix(ch_bra,ch_ket);
         auto brapp = tbc_bra.GetKetIndex_pp();
         auto ketpp = tbc_ket.GetKetIndex_pp();
         auto brahh = tbc_bra.GetKetIndex_hh();
         auto kethh = tbc_ket.GetKetIndex_hh();

         // More thinking needed...
         arma::mat Matrixpp = LHS.rows(ketpp) * RHS.cols(brapp);
         arma::mat Matrixhh = LHS.rows(kethh) * RHS.cols(brahh);


         if (opout.IsHermitian())
         {
            Matrixpp -= Matrixpp.t();
            Matrixhh -= Matrixhh.t();
         }
         else if (opout.IsAntiHermitian())
         {
            Matrixpp += Matrixpp.t();
            Matrixhh += Matrixhh.t();
         }
         else //Figure this out later...
         {
//            Matrixpp -= RHS.t().rows(tbc_ket.GetKetIndex_pp()) * LHS.cols(tbc_bra.GetKetIndex_pp());
//            Matrixhh -= RHS.t().rows(tbc_ket.GetKetIndex_hh()) * LHS.cols(tbc_bra.GetKetIndex_hh());
         }

      
         OUT += Matrixpp - Matrixhh;
  //    }
   }
}







// Since comm222_pp_hh and comm211 both require the ruction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
void Operator::comm222_pp_hh_221st( Operator& opright, Operator& opout )  
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();

//   Operator Mpp = opout;
//   Operator Mhh = opout;
   TwoBodyME Mpp = opout.TwoBody;
   TwoBodyME Mhh = opout.TwoBody;

//   #pragma omp parallel for schedule(dynamic,5)
//   for (int ch_bra=0;ch_bra<nChannels;++ch_bra)
//   {
     for ( auto& itmat : opright.TwoBody.MatEl )
     {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);

//      auto& LHS = TwoBody.at(ch_bra).at(ch_bra);
      auto& LHS = TwoBody.GetMatrix(ch_bra,ch_bra);

//      {
         TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
         if ( ch_ket != ch_bra )
         {
            cout << " AAH! ch_bra = " << ch_bra << "  ch_ket = " << ch_ket << endl;
         }

//         auto& RHS  =  opright.TwoBody.at(ch_bra).at(ch_ket);
         auto& RHS  =  opright.TwoBody.GetMatrix(ch_bra,ch_ket);
         arma::mat& OUT2 =    opout.TwoBody.GetMatrix(ch_bra,ch_ket);

         arma::mat& Matrixpp =  Mpp.GetMatrix(ch_bra,ch_ket);
         arma::mat& Matrixhh =  Mhh.GetMatrix(ch_bra,ch_ket);
     
//         arma::uvec kets_pp = arma::uvec(tbc_bra.KetIndex_pp);
//         arma::uvec kets_hh = arma::uvec(tbc_bra.KetIndex_hh);
         arma::uvec kets_pp = tbc_bra.GetKetIndex_pp();
         arma::uvec kets_hh = tbc_bra.GetKetIndex_hh();
      
         Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
         Matrixhh =  LHS.cols(kets_hh) * RHS.rows(kets_hh);
 
         if (opout.IsHermitian())
         {
            Matrixpp += Matrixpp.t();
            Matrixhh += Matrixhh.t();
         }
         else if (opout.IsAntiHermitian()) // i.e. LHS and RHS are hermitian
         {
            Matrixpp -= Matrixpp.t();
            Matrixhh -= Matrixhh.t();
         }
         else // Figure this out later...
         {
//            Matrixpp -=  RHS.t().cols(kets_pp) * LHS.t().rows(kets_pp);
//            Matrixhh -=  RHS.t().cols(kets_hh) * LHS.t().rows(kets_hh);
         }

         // The two body part
         OUT2 += Matrixpp - Matrixhh;


         int J1 = tbc_bra.J;
         int J2 = tbc_ket.J;
         double hatfactor = (2*J1+1)*sqrt(2*J2+1);

      // The one body part
      for (int i=0;i<norbits;++i)
      {
         Orbit &oi = modelspace->GetOrbit(i);
         int jmin = opout.IsNonHermitian() ? 0 : i;
         for (int j=jmin;j<norbits;++j)
         {
            Orbit &oj = modelspace->GetOrbit(j);
            if (oi.j2 != oj.j2 or oi.l != oj.l or oi.tz2 != oj.tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (auto& c : modelspace->holes)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ -=   Mpp.GetTBME(ch_bra,ch_ket,c,i,j,c)*sixj;
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               Orbit &oc = modelspace->GetOrbit(c);
               double sixj = modelspace->GetSixJ(J1, J2, opright.rank_J, oj.j2/2.0, oi.j2/2.0, oc.j2/2.0);
               cijJ -=   Mhh.GetTBME(ch_bra,ch_ket,c,i,j,c)*sixj;
            }
            #pragma omp critical
            opout.OneBody(i,j) += cijJ *hatfactor / sqrt(oi.j2 +1.0);
            if (not opout.IsNonHermitian())
            {
               opout.OneBody(j,i) = herm * opout.OneBody(i,j);
            }
         } // for j
       } // for i
//     } //for ch_ket
   } //for ch_bra
}




//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.             [X2,Y2](2)_ph = -sum_ab (na-nb) sum_J' (-1)^(J'+ja+jb)
//
//                                                                             v----|    v----|                                    v----|    v----|
//   |          |      |          |              * [ (-1)^(ji+jl) { jk jl J } <bi|X|aj> <ak|Y|bl>   - (-1)^(jj+jl-J') { jk jl J } <bj|X|ai> <ak|Y|bl>
//   |     __Y__|      |     __X__|                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |                                            v----|    v----|                                    v----|    v----|
//   |____\/    |      |____\/    |            -  (-1)^(ji+jk-J') { jk jl J } <bi|X|aj> <al|Y|bk>   + (-1)^(jj+jk)    { jk jl J } <bj|X|ai> <al|Y|bk>  ]
//   |  X       |      |  Y       |                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//           
//            
// -- This appears to agree with Nathan's results
//
// Haven't converted this one yet...
void Operator::comm222_phst( Operator& opright, Operator& opout ) 
{
   int herm = opout.IsHermitian() ? 1 : -1;


   Operator Xcc = Operator(*modelspace, 0, 1, 0, 2);
   Operator Ycc = Operator(Xcc);

//   DoPandyaTransformation(Xcc);
//   opright.DoPandyaTransformation(Ycc);


   // Update Cross-coupled matrix elements
/*
   vector<arma::mat> X_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> X_TwoBody_CC_right (nChannels, arma::mat() );

   vector<arma::mat> Y_TwoBody_CC_left (nChannels, arma::mat() );
   vector<arma::mat> Y_TwoBody_CC_right (nChannels, arma::mat() );


   CalculateCrossCoupled(X_TwoBody_CC_left, X_TwoBody_CC_right );
   //CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   opright.CalculateCrossCoupled(Y_TwoBody_CC_left, Y_TwoBody_CC_right );
   // Construct the intermediate matrices N1 and N2
   vector<arma::mat> N1 (nChannels, arma::mat() );
   vector<arma::mat> N2 (nChannels, arma::mat() );
 // probably better not to parallelize here, since the armadillo matrix muliplication can use OMP
 // I think that N1 == N2.t() ??? Check this. No.
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] = X_TwoBody_CC_left[ch] * Y_TwoBody_CC_right[ch] - Y_TwoBody_CC_left[ch] * X_TwoBody_CC_right[ch];
      N1[ch] += N1[ch].t(); // maybe this works?  Yes !
//      N2[ch] = Y_TwoBody_CC_left[ch] *  X_TwoBody_CC_right[ch];
//      N1[ch] -= N2[ch]; // maybe this works? Yes!
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& OUT = (arma::mat&) opout.TwoBody.at(ch).at(ch);

      int nKets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;

         for (int iket=ibra; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;


            double comm = 0;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi.l+ok.l)%2;
            int Tz_cc = abs(oi.tz2+ok.tz2)/2;
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+tbc.J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
            }

            parity_cc = (oi.l+ol.l)%2;
            Tz_cc = abs(oi.tz2+ol.tz2)/2;
            jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_il = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,l),max(i,l));
               int indx_jk = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,k),max(j,k));
               if (i>l)
                 indx_il += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>k)
                 indx_jk += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               comm -= (2*Jprime+1) * phase * sixj * (me1);
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            comm /= norm;
            #pragma omp critical
            {
               OUT(ibra,iket) += comm;
              if (iket != ibra)
              {
                 OUT(iket,ibra) += herm*comm;
              }
            }
            if (abs(comm) > 1e-7)
            {
               cout << modelspace->phase(ji+jj+jk+jl) << ":  "  << comm << endl;
            }

         }
      }
   }
*/
}



