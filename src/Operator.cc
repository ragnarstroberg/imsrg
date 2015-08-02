
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

//double  Operator::bch_transform_threshold = 1e-6;
double  Operator::bch_transform_threshold = 1e-9;
double  Operator::bch_product_threshold = 1e-4;
map<string, double> Operator::timer;


Operator::~Operator()
{
//   cout << "calling Operator destructor" << endl;
}

/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator()
 :   modelspace(NULL), 
    rank_J(0), rank_T(0), parity(0), particle_rank(2),
    hermitian(true), antihermitian(false), nChannels(0)
{
}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms,Jrank,Trank,p),  ThreeBody(&ms),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank),
    E3max(ms.GetN3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels()) 
{
  SetUpOneBodyChannels();
  if (particle_rank >=3) ThreeBody.Allocate();
}



Operator::Operator(ModelSpace& ms) :
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms),  ThreeBody(&ms),
    rank_J(0), rank_T(0), parity(0), particle_rank(2),
    E3max(ms.GetN3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels())
{
  SetUpOneBodyChannels();
}

Operator::Operator(const Operator& op)
: modelspace(op.modelspace),  ZeroBody(op.ZeroBody),
  OneBody(op.OneBody), TwoBody(op.TwoBody) ,ThreeBody(op.ThreeBody),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank),
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels)
{
//   cout << "Calling copy constructor for Operator" << endl;
}

Operator::Operator(Operator&& op)
: modelspace(op.modelspace), ZeroBody(op.ZeroBody),
  OneBody(move(op.OneBody)), TwoBody(move(op.TwoBody)) , ThreeBody(move(op.ThreeBody)),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank),
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels)
{
//   cout << "Calling move constructor for Operator" << endl;
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
   OneBodyChannels = op.OneBodyChannels;
}

/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////
Operator& Operator::operator=(const Operator& rhs)
{
//   cout << "Using copy assignment" << endl;
   Copy(rhs);
   return *this;
}

Operator& Operator::operator=(Operator&& rhs)
{
//   cout << "Using move assignment" << endl;
   modelspace    = rhs.modelspace;
   nChannels     = rhs.nChannels;
   hermitian     = rhs.hermitian;
   antihermitian = rhs.antihermitian;
   rank_J        = rhs.rank_J;
   rank_T        = rhs.rank_T;
   parity        = rhs.parity;
   particle_rank = rhs.particle_rank;
   E2max         = rhs.E2max;
   E3max         = rhs.E3max;
   ZeroBody      = rhs.ZeroBody;
   OneBody       = move(rhs.OneBody);
   TwoBody       = move(rhs.TwoBody);
   ThreeBody     = move(rhs.ThreeBody);
   OneBodyChannels = move(rhs.OneBodyChannels);
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



void Operator::PrintTimes()
{
   cout << "====================== TIMES =======================" << endl;
   cout.setf(ios::fixed);
   for ( auto it : timer )
   {
     cout << setw(40) << std::left << it.first + ":  " << setw(12) << setprecision(5) << std::right << it.second  << endl;
   }
}


void Operator::SetUpOneBodyChannels()
{
  for ( int i=0; i<modelspace->GetNumberOrbits(); ++i )
  {
    Orbit& oi = modelspace->GetOrbit(i);
//    int lmin = max( oi.l - rank_J - parity, (oi.l+parity)%2);
    int lmin = max( oi.l - rank_J + (rank_J+parity)%2, (oi.l+parity)%2);
    int lmax = min( oi.l + rank_J, modelspace->Nmax);
    for (int l=lmin; l<=lmax; l+=2)
    {
      int j2min = max(max(oi.j2 - 2*rank_J, 2*l-1),1);
      int j2max = min(oi.j2 + 2*rank_J, 2*l+1);
      for (int j2=j2min; j2<=j2max; j2+=2)
      {
        int tz2min = max( oi.tz2 - 2*rank_T, -1);
        int tz2max = min( oi.tz2 + 2*rank_T, 1);
        for (int tz2=tz2min; tz2<=tz2max; tz2+=2)
        {
          OneBodyChannels[ {l, j2, tz2} ].push_back(i);
        }
      }
    }
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
///  set up for scalar or tensor operators, but
///  the tensor part hasn't been tested
//*************************************************************
Operator Operator::DoNormalOrdering2()
{
//   Operator opNO = *this;
   Operator opNO(*this);

   if (opNO.rank_J==0 and opNO.rank_T==0 and opNO.parity==0)
   {
     for (auto& k : modelspace->holes) // loop over hole orbits
     {
        opNO.ZeroBody += (modelspace->GetOrbit(k).j2+1) * OneBody(k,k);
     }
   }

   index_t norbits = modelspace->GetNumberOrbits();

   for ( auto& itmat : TwoBody.MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      auto& matrix = itmat.second;
      
      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J_bra = tbc_bra.J;
      int J_ket = tbc_ket.J;

      // Zero body part
      if (opNO.rank_J==0 and opNO.rank_T==0 and opNO.parity==0)
      {
        arma::vec diagonals = matrix.diag();
        auto hh = tbc_ket.GetKetIndex_hh();
        opNO.ZeroBody += arma::sum( diagonals.elem(hh) ) * (2*J_ket+1);
      }

      // One body part
      for (index_t a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         double ja = oa.j2/2.0;
         index_t bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for ( auto& b : opNO.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) // OneBodyChannels should be moved to the operator, to accommodate tensors
         {
            if (b < bstart) continue;
            Orbit &ob = modelspace->GetOrbit(b);
            double jb = ob.j2/2.0;
            for (auto& h : modelspace->holes)  // C++11 syntax
            {
              if (opNO.rank_J==0)
              {
                 opNO.OneBody(a,b) += (2*J_ket+1.0)/(2*ja+1)  * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
              }
              else
              {
                 Orbit &oh = modelspace->GetOrbit(h);
                 double jh = oh.j2/2.0;
                 opNO.OneBody(a,b) += sqrt((2*J_bra+1.0)*(2*J_ket+1.0)) *modelspace->phase(ja+jh+J_ket+opNO.rank_J)
                                             * modelspace->GetSixJ(J_bra,J_ket,opNO.rank_J,jb,ja,jh) * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
              }
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
            Gamma(ibra,iket) /= (2*tbc.J+1)* sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
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


Operator Operator::UndoNormalOrdering()
{
   Operator opNO = *this;
   cout << "Undoing Normal ordering. Initial ZeroBody = " << opNO.ZeroBody << endl;

   for (auto& k : modelspace->holes) // loop over hole orbits
   {
      cout << "*** " << opNO.ZeroBody << endl;
      opNO.ZeroBody -= (modelspace->GetOrbit(k).j2+1) * OneBody(k,k);
      cout << "k = " << k << "  0body = " << opNO.ZeroBody << endl;
   }

   index_t norbits = modelspace->GetNumberOrbits();

   for ( auto& itmat : TwoBody.MatEl )
   {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      auto& matrix = itmat.second;
      
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J_ket = tbc_ket.J;

      // Zero body part
      arma::vec diagonals = matrix.diag();
      auto hh = tbc_ket.GetKetIndex_hh();
      opNO.ZeroBody += arma::sum( diagonals.elem(hh) ) * (2*J_ket+1);
      cout << "ch_bra = " << ch_bra << "  0body = " << opNO.ZeroBody << endl;

      // One body part
      for (index_t a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         index_t bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for ( auto& b : opNO.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) // OneBodyChannels should be moved to the operator, to accommodate tensors
         {
            if (b < bstart) continue;
            for (auto& h : modelspace->holes)  // C++11 syntax
            {
               opNO.OneBody(a,b) -= (2*J_ket+1.0)/(oa.j2+1)  * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
            }
         }
      }
   } // loop over channels

   if (hermitian) opNO.Symmetrize();
   if (antihermitian) opNO.AntiSymmetrize();

   cout << "Zero-body piece is now " << opNO.ZeroBody << endl;
   return opNO;

}


ModelSpace* Operator::GetModelSpace()
{
   return modelspace;
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

void Operator::SetHermitian()
{
  hermitian = true;
  antihermitian = false;
  TwoBody.SetHermitian();
}

void Operator::SetAntiHermitian()
{
  hermitian = false;
  antihermitian = true;
  TwoBody.SetAntiHermitian();
}

void Operator::SetNonHermitian()
{
  hermitian = false;
  antihermitian = false;
  TwoBody.SetNonHermitian();
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




//*****************************************************************************************
/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Operator::BCH_Transform(  Operator &Omega)
{
   double t = omp_get_wtime();
   int max_iter = 40;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   Operator OpNested = *this;
   if (nx<bch_transform_threshold) return OpOut;
   double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
   for (int i=1; i<max_iter; ++i)
   {
//      OpNested = Omega.Commutator(OpNested);
      OpNested = Commutator(Omega,OpNested);
      OpNested /= i;

      OpOut += OpNested;

      if (OpNested.Norm() < epsilon *(i+1))
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
   double nx = X.Norm();
   double ny = Y.Norm();
   if (nx < 1e-7) return Y;
   if (ny < 1e-7) return X;

//   Operator Z = X.Commutator(Y);
   Operator Z = Commutator(X,Y);
   Z *= 0.5;
   double nxy = Z.Norm();

   if ( nxy < (nx+ny)*bch_product_threshold )
   {
     Z += X;
     Z += Y;
     timer["BCH_Product"] += omp_get_wtime() - t;
     return Z;
   }

   Y -= X;
//   Z += (1./6)* Z.Commutator( Y );
   Z += (1./6)* Commutator( Z, Y );
   Z += Y;
   X *=2;
   Z += X;

   timer["BCH_Product"] += omp_get_wtime() - t;
   return Z;
}

/// Obtain the Frobenius norm of the operator, which here is 
/// defined as 
/// \f[ \|X\| = \sqrt{\|X_{(1)}\|^2 +\|X_{(2)}\|^2 } \f]
/// and
/// \f[ \|X_{(1)}\|^2 = \sum\limits_{ij} X_{ij}^2 \f]
double Operator::Norm() const
{
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
   if (rank_J==0)
     OneBody = arma::symmatu(OneBody);
   else
   {
     int norb = modelspace->GetNumberOrbits();
     for (int i=0;i<norb; ++i)
     {
       Orbit& oi = modelspace->GetOrbit(i);
       for ( int j : OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
       {
         if (j<= i) continue;
         Orbit& oj = modelspace->GetOrbit(j);
         OneBody(j,i) = modelspace->phase((oi.j2-oj.j2)/2) * OneBody(i,j);
       }
     }
   }
   TwoBody.Symmetrize();
}

void Operator::AntiSymmetrize()
{
   if (rank_J==0)
   {
     OneBody = arma::trimatu(OneBody) - arma::trimatu(OneBody).t();
   }
   else
   {
     int norb = modelspace->GetNumberOrbits();
     for (int i=0;i<norb; ++i)
     {
       Orbit& oi = modelspace->GetOrbit(i);
       for ( int j : OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
       {
         if (j<= i) continue;
         if (rank_J==0)
          OneBody(j,i) = -OneBody(i,j);
         else
         {
           Orbit& oj = modelspace->GetOrbit(j);
           OneBody(j,i) = -modelspace->phase((oi.j2-oj.j2)/2) * OneBody(i,j);
         }
       }
     }
   }
   TwoBody.AntiSymmetrize();
}

//Operator Operator::Commutator( Operator& opright)
/// Returns \f$ Z = [X,Y] \f$
/// @relates Operator
Operator Commutator( const Operator& X, const Operator& Y)
{
   X.timer["N_Commutators"] += 1;
   if (X.rank_J==0)
   {
      if (Y.rank_J==0)
      {
//         return X.CommutatorScalarScalar(Y); // [S,S]
         return CommutatorScalarScalar(X,Y); // [S,S]
      }
      else
      {
//         return X.CommutatorScalarTensor(Y); // [S,T]
         return CommutatorScalarTensor(X,Y); // [S,T]
      }
   }
   else if(Y.rank_J==0)
   {
//      return (-1)*Y.CommutatorScalarTensor(X); // [T,S]
      return (-1)*CommutatorScalarTensor(Y,X); // [T,S]
   }
   else
   {
      cout << "In Tensor-Tensor because X.rank_J = " << X.rank_J << "  and Y.rank_J = " << Y.rank_J << endl;
      cout << " Tensor-Tensor commutator not yet implemented." << endl;
      return X;
   }
}


/// Commutator where \f$ X \f$ and \f$Y\f$ are scalar operators.
/// Should be called through Commutator()
//Operator Operator::CommutatorScalarScalar( Operator& opright) 
Operator CommutatorScalarScalar( const Operator& X, const Operator& Y) 
{
   Operator Z = Y;
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

//   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
//   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
//   else out.SetNonHermitian();
   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   if ( not Z.IsAntiHermitian() )
   {
      //X.comm110ss(Y, Z);
      //X.comm220ss(Y, Z) ;
      Z.comm110ss(X, Y);
      Z.comm220ss(X, Y) ;
   }

   double t = omp_get_wtime();
//   X.comm111ss(Y, Z);
   Z.comm111ss(X, Y);
   Z.timer["comm111ss"] += omp_get_wtime() - t;

    t = omp_get_wtime();
//   X.comm121ss(opright, out);
   Z.comm121ss(X,Y);
    Z.timer["comm121ss"] += omp_get_wtime() - t;

    t = omp_get_wtime();
//   X.comm122ss(Y, Z); 
   Z.comm122ss(X,Y); 
    Z.timer["comm122ss"] += omp_get_wtime() - t;

   if (X.particle_rank>1 and Y.particle_rank>1)
   {
    t = omp_get_wtime();
//    X.comm222_pp_hh_221ss(Y, Z);
    Z.comm222_pp_hh_221ss(X, Y);
    Z.timer["comm222_pp_hh_221ss"] += omp_get_wtime() - t;
     
    t = omp_get_wtime();
//    X.comm222_phss(Y, Z);
    Z.comm222_phss(X, Y);
    Z.timer["comm222_phss"] += omp_get_wtime() - t;
   }


   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();


   return Z;
}


/// Commutator \f$[X,Y]\f$ where \f$ X \f$ is a scalar operator and \f$Y\f$ is a tensor operator.
/// Should be called through Commutator()
//Operator Operator::CommutatorScalarTensor( Operator& opright) 
Operator CommutatorScalarTensor( const Operator& X, const Operator& Y) 
{
   Operator Z = Y; // This ensures the commutator has the same tensor rank as Y
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   Z.comm111st(X, Y);
   Z.comm121st(X, Y);

   Z.comm122st(X, Y);
   Z.comm222_pp_hh_221st(X, Y);
   Z.comm222_phst(X, Y);

   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();

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
//void Operator::comm110ss( Operator& Y, Operator& Z) 
void Operator::comm110ss( const Operator& X, const Operator& Y) 
{
  Operator& Z = *this;
  if (X.IsHermitian() and Y.IsHermitian()) return ; // I think this is the case
  if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return ; // I think this is the case

   arma::mat xyyx = X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
   for ( auto& a : modelspace->holes) 
   {
      Z.ZeroBody += (modelspace->GetOrbit(a).j2+1) * xyyx(a,a);
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
/// [X_{(2)},Y_{(2)}]_{(0)} = \frac{1}{2} \sum_{J} (2J+1) \sum_{abcd} (n_a n_b \bar{n}_c \bar{n}_d) \tilde{X}_{abcd}^{J} \tilde{Y}_{cdab}^{J}
/// \f]
/// may be rewritten as
/// \f[
/// [X_{(2)},Y_{(2)}]_{(0)} = 2 \sum_{J} (2J+1) Tr(X_{hh'pp'}^{J} Y_{pp'hh'}^{J})
/// \f] where we obtain a factor of four from converting two unrestricted sums to restricted sums, i.e. \f$\sum_{ab} \rightarrow \sum_{a\leq b} \f$,
/// and using the normalized TBME.
void Operator::comm220ss( const Operator& X, const Operator& Y) 
{
   Operator& Z = *this;
   if (X.IsHermitian() and Y.IsHermitian()) return; // I think this is the case
   if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return; // I think this is the case

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto hh = tbc.GetKetIndex_hh();
      auto pp = tbc.GetKetIndex_pp();
      auto & X2 = X.TwoBody.GetMatrix(ch);
      auto & Y2 = Y.TwoBody.GetMatrix(ch);
      Z.ZeroBody += 2 * (2*tbc.J+1) * arma::trace( X2.submat(hh,pp) * Y2.submat(pp,hh) );
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
//void Operator::comm111ss( Operator & Y, Operator& Z) 
void Operator::comm111ss( const Operator & X, const Operator& Y) 
{
   Operator& Z = *this;
   Z.OneBody += X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
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
//void Operator::comm121ss( Operator& Y, Operator& Z) 
void Operator::comm121ss( const Operator& X, const Operator& Y) 
{
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for 
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (auto& b : modelspace->particles)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                Z.OneBody(i,j) += (ob.j2+1) *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                Z.OneBody(i,j) -= (oa.j2+1) *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;

                // comm211 part
                Z.OneBody(i,j) -= (ob.j2+1) *  Y.OneBody(a,b) * X.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                Z.OneBody(i,j) += (oa.j2+1) *  Y.OneBody(b,a) * X.TwoBody.GetTBMEmonopole(a,i,b,j) ;
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
//void Operator::comm221ss( Operator& Y, Operator& Z) 
void Operator::comm221ss( const Operator& X, const Operator& Y) 
{

   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();
   TwoBodyME Mpp = Y.TwoBody;
   TwoBodyME Mhh = Y.TwoBody;

   #pragma omp parallel for schedule(dynamic,1)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = (arma::mat&) X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = (arma::mat&) Y.TwoBody.GetMatrix(ch,ch);
      arma::mat& Matrixpp =  Mpp.GetMatrix(ch,ch);
      arma::mat& Matrixhh =  Mpp.GetMatrix(ch,ch);
      
      Matrixpp = ( LHS.rows(tbc.GetKetIndex_pp()) * RHS.cols(tbc.GetKetIndex_pp()));
      Matrixpp -= Matrixpp.t();
      Matrixhh = ( LHS.rows(tbc.GetKetIndex_hh()) * RHS.cols(tbc.GetKetIndex_hh()));
      Matrixhh -= Matrixhh.t();
   }

   // If commutator is hermitian or antihermitian, we only
   // need to do half the sum. Add this.
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         double cijJ = 0;
         // Sum c over holes and include the nbar_a * nbar_b terms
         for (auto& c : modelspace->holes)
         {
           Orbit &oc = modelspace->GetOrbit(c);
           int jmin = abs(oi.j2-oc.j2)/2;
           int jmax = (oi.j2+oc.j2)/2;
           for (int J=jmin; J<=jmax; ++J)
             cijJ +=   (2*J+1) * Mpp.GetTBME_J(J,i,c,j,c);
         // Sum c over particles and include the n_a * n_b terms
         }
         for (auto& c : modelspace->particles)
         {
           Orbit &oc = modelspace->GetOrbit(c);
           int jmin = abs(oi.j2-oc.j2)/2;
           int jmax = (oi.j2+oc.j2)/2;
           for (int J=jmin; J<=jmax; ++J)
             cijJ +=  (2*J+1) * Mhh.GetTBME_J(J,i,c,j,c);
         }
         Z.OneBody(i,j) +=  cijJ / (oi.j2+1.0);
      } // for j
   } // for i
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
//void Operator::comm122ss( Operator& Y, Operator& Z ) 
void Operator::comm122ss( const Operator& X, const Operator& Y ) 
{
   Operator& Z = *this;
   auto& X1 = X.OneBody;
   auto& Y1 = Y.OneBody;

   int n_nonzero = modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& X2 = X.TwoBody.GetMatrix(ch,ch);
      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
      auto& Z2 = Z.TwoBody.GetMatrix(ch,ch);


      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         double pre_ij = i==j ? SQRT2 : 1;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         arma::Row<double> X2_ij = X2.row(indx_ij); // trying this to better use the cache. not sure if it helps.
         arma::Row<double> Y2_ij = Y2.row(indx_ij);
         int klmin = Z.IsNonHermitian() ? 0 : indx_ij;
         for (int indx_kl=klmin;indx_kl<npq; ++indx_kl)
         {
            Ket & ket = tbc.GetKet(indx_kl);
            int k = ket.p;
            int l = ket.q;
            double pre_kl = k==l ? SQRT2 : 1;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            arma::vec X2_kl = X2.unsafe_col(indx_kl);
            arma::vec Y2_kl = Y2.unsafe_col(indx_kl);

            double cijkl = 0;


            for (int a : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
                 int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
                 if (indx_aj < 0) continue;
                 double pre_aj = a>j ? tbc.GetKet(indx_aj).Phase(tbc.J) : (a==j ? SQRT2 : 1);
                 cijkl += pre_kl * pre_aj  * ( X1(i,a) * Y2(indx_aj,indx_kl) - Y1(i,a) * X2_kl(indx_aj) );
            }

            for (int a : modelspace->OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
            {
                 int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
                 if (indx_ia < 0) continue;
                 double pre_ia = i>a ? tbc.GetKet(indx_ia).Phase(tbc.J) : (i==a ? SQRT2 : 1);
                 cijkl += pre_kl * pre_ia * ( X1(j,a) * Y2(indx_ia,indx_kl)  - Y1(j,a) * X2_kl(indx_ia) );
             }

            for (int a : modelspace->OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
                int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
                if (indx_al < 0) continue;
                double pre_al = a>l ? tbc.GetKet(indx_al).Phase(tbc.J) : (a==l ? SQRT2 : 1);
                cijkl -= pre_ij * pre_al * ( X1(a,k) * Y2(indx_ij,indx_al) - Y1(a,k) * X2_ij(indx_al) );
            }

            for (int a : modelspace->OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
            {
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));
               if (indx_ka < 0) continue;
               double pre_ka = k>a ? tbc.GetKet(indx_ka).Phase(tbc.J) : (k==a ? SQRT2 : 1);
               cijkl -= pre_ij * pre_ka * ( X1(a,l) * Y2(indx_ij,indx_ka) - Y1(a,l) * X2_ij(indx_ka) );
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Z2(indx_ij,indx_kl) += cijkl / norm;
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
//void Operator::comm222_pp_hhss( Operator& opright, Operator& opout ) 
void Operator::comm222_pp_hhss( const Operator& X, const Operator& Y ) 
{
   Operator& Z = *this;
//   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& LHS = (arma::mat&) X.TwoBody.GetMatrix(ch,ch);
      arma::mat& RHS = (arma::mat&) Y.TwoBody.GetMatrix(ch,ch);
      arma::mat& OUT = (arma::mat&) Z.TwoBody.GetMatrix(ch,ch);

      arma::mat Mpp = (LHS.rows(tbc.GetKetIndex_pp()) * RHS.cols(tbc.GetKetIndex_pp()));
      arma::mat Mhh = (LHS.rows(tbc.GetKetIndex_hh()) * RHS.cols(tbc.GetKetIndex_hh()));
      OUT += Mpp - Mpp.t() - Mhh + Mhh.t();
   }
}







/// Since comm222_pp_hhss() and comm221ss() both require the ruction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
//void Operator::comm222_pp_hh_221ss( Operator& Y, Operator& Z )  
void Operator::comm222_pp_hh_221ss( const Operator& X, const Operator& Y )  
{

//   int herm = Z.IsHermitian() ? 1 : -1;
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();

   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   // Don't use omp, because the matrix multiplication is already
   // parallelized by armadillo.
   int nch = modelspace->SortedTwoBodyChannels.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.TwoBody.GetMatrix(ch,ch);
      auto& OUT = Z.TwoBody.GetMatrix(ch,ch);

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * RHS.rows(kets_hh);

      if (Z.IsHermitian())
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (Z.IsAntiHermitian()) // i.e. LHS and RHS are both hermitian or ant-hermitian
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
      int jmin = Z.IsNonHermitian() ? 0 : i;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
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
               cijJ += Jfactor * Mpp.GetTBME(ch,c,i,c,j);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   timer["pphh One Body bit"] += omp_get_wtime() - t;
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
//void Operator::DoPandyaTransformation(TwoBodyME& TwoBody_CC_hp, TwoBodyME& TwoBody_CC_ph)
void Operator::DoPandyaTransformation(vector<arma::mat>& TwoBody_CC_hp, vector<arma::mat>& TwoBody_CC_ph) const
{
   // loop over cross-coupled channels
   int n_nonzero = modelspace->SortedTwoBodyChannels_CC.size();
   int herm = IsHermitian() ? 1 : -1;
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      arma::uvec& kets_ph = tbc_cc.GetKetIndex_ph();
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

         // loop over cross-coupled kets |bc> in this channel
         // we go to 2*nKets to include |bc> and |cb>
         for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
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
            TwoBody_CC_ph[ch_cc](ibra,iket_cc+nKets_cc) = herm* modelspace->phase(ja+jb+jc+jd) * sm;


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
            TwoBody_CC_hp[ch_cc](ibra,iket_cc+nKets_cc) = herm* modelspace->phase(ja+jb+jc+jd) * sm;

         }
      }
   }
}


//void Operator::InversePandyaTransformation(vector<arma::mat>& W, vector<arma::mat>& opout)
void Operator::AddInversePandyaTransformation(vector<arma::mat>& Zbar)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
//   for (int ch=0;ch<nChannels;++ch)
   int n_nonzeroChannels = modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
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
         int ketmin = IsHermitian() ? ibra : ibra+1;
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
               double me1 = Zbar[ch_cc](indx_il,indx_kj);
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
               double me1 = Zbar[ch_cc](indx_jl,indx_ki);
               commji += (2*Jprime+1) *  sixj * me1;
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            TwoBody.GetMatrix(ch,ch)(ibra,iket) += (commij - modelspace->phase(ji+jj-J)*commji) / norm;
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
//void Operator::comm222_phss( Operator& Y, Operator& Z ) 
void Operator::comm222_phss( const Operator& X, const Operator& Y ) 
{

   Operator& Z = *this;
   // Create Pandya-transformed hp and ph matrix elements
   vector<arma::mat> X_bar_hp (nChannels );
   vector<arma::mat> X_bar_ph (nChannels );
   vector<arma::mat> Y_bar_hp (nChannels );
   vector<arma::mat> Y_bar_ph (nChannels );

   double t = omp_get_wtime();
   X.DoPandyaTransformation(X_bar_hp, X_bar_ph );
   Y.DoPandyaTransformation(Y_bar_hp, Y_bar_ph );
   timer["DoPandyaTransformation"] += omp_get_wtime() - t;

   // Construct the intermediate matrix Z_bar
   t = omp_get_wtime();
   vector<arma::mat> Z_bar (nChannels );

//   for (int ch : modelspace->SortedTwoBodyChannels_CC )
   int nch = modelspace->SortedTwoBodyChannels_CC.size();
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich )
   {
      int ch = modelspace->SortedTwoBodyChannels_CC[ich];
      if ( X.IsHermitian() ) // keep track of minus sign from taking transpose of X
         Z_bar[ch] =  X_bar_hp[ch].t() * Y_bar_hp[ch] - X_bar_ph[ch].t() * Y_bar_ph[ch] ;
      else
         Z_bar[ch] =  X_bar_ph[ch].t() * Y_bar_ph[ch] - X_bar_hp[ch].t() * Y_bar_hp[ch] ;

      if ( Z.IsHermitian() ) // if Z is hermitian, then XY is antihermitian
        Z_bar[ch] += Z_bar[ch].t();
      else
        Z_bar[ch] -= Z_bar[ch].t();
   }
   timer["Build Z_bar"] += omp_get_wtime() - t;

   // Perform inverse Pandya transform on W_bar to get Z
   t = omp_get_wtime();
   Z.AddInversePandyaTransformation(Z_bar);
   timer["InversePandyaTransformation"] += omp_get_wtime() - t;

}






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
void Operator::comm111st( const Operator & X, const Operator& Y)
{
  comm111ss(X,Y);
}
//{
//   Operator& Z = *this;
//   Z.OneBody += X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
//}

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
// X is scalar one-body, Y is tensor two-body
// There must be a better way to do this looping. 
//
//void Operator::comm121st( Operator& Y, Operator& Z) 
void Operator::comm121st( const Operator& X, const Operator& Y) 
{
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();
   int Lambda = Z.GetJRank();
   #pragma omp parallel for // for starters, don't do it parallel
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          Orbit &oj = modelspace->GetOrbit(j);
          double jj = oj.j2/2.0;
          if (j<i) continue; // only calculate upper triangle
          double& Zij = Z.OneBody(i,j);
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             double ja = oa.j2/2.0;
               for (auto& b : modelspace->particles) // is this is slow, it can probably be sped up by looping over OneBodyChannels
               {
                  Orbit &ob = modelspace->GetOrbit(b);
                  double jb = ob.j2/2.0;
                  if (ob.j2 == oa.j2 and ob.l == oa.l and ob.tz2 == oa.tz2)
                  {
                    int J1min = min(abs(ji-ja),abs(jj-ja));
                    int J1max = max(ji,jj) + ja;
                    for (int J1=J1min; J1<=J1max; ++J1)
                    {
                      int phasefactor = modelspace->phase(jj+ja+J1+Lambda);
                      int J2min = max(abs(Lambda - J1),J1min);
                      int J2max = min(Lambda + J1,J1max);
                      for (int J2=J2min; J2<=J2max; ++J2)
                      {
//                        double toscalar = sqrt((2*J1+1)/(2*ji+1));
//                        double prefactor = toscalar * phasefactor * sqrt((2*J1+1)*(2*J2+1)) * modelspace->GetSixJ(J1,J2,Lambda,jj,ji,ja);
                        double prefactor = phasefactor * sqrt((2*J1+1)*(2*J2+1)) * modelspace->GetSixJ(J1,J2,Lambda,jj,ji,ja);
                        if (J1>=abs(ja-ji) and J1<=ja+ji and J2>=abs(ja-jj) and J2<=ja+jj )
                          Zij +=  prefactor * ( X.OneBody(a,b) * Y.TwoBody.GetTBME_J(J1,J2,b,i,a,j) );
                        if (J1>=abs(ja-jj) and J1<=ja+jj and J2>=abs(ja-ji) and J2<=ja+ji )
                          Zij -= prefactor * X.OneBody(b,a) * Y.TwoBody.GetTBME_J(J1,J2,a,i,b,j)  ;
                      }
                    }
                  }

                  // Now, X is scalar two-body and Y is tensor one-body
                  if ( (abs(ja-jb)>Lambda) or (ja+jb<Lambda) ) continue;
                  int J1min = max(abs(ji-ja),abs(jj-jb));
                  int J1max = min(ji+ja,jj+jb);
                  for (int J1=J1min; J1<=J1max; ++J1)
                  {
                    double toscalar = sqrt((2*ja+1)/(2*ji+1));
                    double prefactor = toscalar*  modelspace->phase(ji+jb+J1) * (2*J1+1) * modelspace->GetSixJ(ja,jb,Lambda,ji,jj,J1);
                    Zij += prefactor * X.TwoBody.GetTBME_J(J1,J1,b,i,a,j) * Y.OneBody(a,b);
                  }

                  J1min = max(abs(ji-jb),abs(jj-ja));
                  J1max = min(ji+jb,jj+ja);
                  for (int J1=J1min; J1<=J1max; ++J1)
                  {
//                    double toscalar = sqrt((2*ja+1)/(2*ji+1));
//                    double prefactor = toscalar* modelspace->phase(ji+ja+J1) * (2*J1+1) * modelspace->GetSixJ(jb,ja,Lambda,ji,jj,J1);
                    double prefactor = modelspace->phase(ji+ja+J1) * (2*J1+1) * modelspace->GetSixJ(jb,ja,Lambda,ji,jj,J1);
                    Zij += prefactor * X.TwoBody.GetTBME_J(J1,J1,a,i,b,j) * Y.OneBody(b,a);
                  }
               }
               
               
             }
          }
      }
   
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
//void Operator::comm122st( Operator& Y, Operator& Z ) 
void Operator::comm122st( const Operator& X, const Operator& Y ) 
{
   Operator& Z = *this;
   int Lambda = Z.rank_J;

    vector< array<int,2> > channels;
    for ( auto& itmat : Z.TwoBody.MatEl ) channels.push_back( itmat.first );
    int nmat = channels.size();
   #pragma omp parallel for schedule(dynamic,1)
    for (int ii=0; ii<nmat; ++ii)
    {
     int ch_bra = channels[ii][0];
     int ch_ket = channels[ii][1];

      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nbras = tbc_bra.GetNumberKets();
      int nkets = tbc_ket.GetNumberKets();
      double hatfactor = sqrt((2*J1+1)*(2*J2+1));
      arma::mat& Z2 = Z.TwoBody.GetMatrix(ch_bra,ch_ket);

      for (int ibra = 0;ibra<nbras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.0;
         double jj = oj.j2/2.0;
//         for (int iket=ibra;iket<nkets; ++iket)
         for (int iket=0;iket<nkets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit& ok = modelspace->GetOrbit(k);
            Orbit& ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.0;
            double jl = ol.j2/2.0;

            double cijkl = 0;

            for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
              cijkl += X.OneBody(i,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,a,j,k,l);
            for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
               cijkl += X.OneBody(j,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,a,k,l);
            for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
               cijkl -= X.OneBody(a,k) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,a,l);
            for ( int a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
               cijkl -= X.OneBody(a,l) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,k,a);


            double prefactor = hatfactor * modelspace->phase(ji+jj+J2+Lambda) ;
            for ( int a : Y.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2/2.0;
               cijkl -= prefactor *  modelspace->GetSixJ(J2,J1,Lambda,ji,ja,jj) * Y.OneBody(i,a) * X.TwoBody.GetTBME(ch_bra,ch_bra,a,j,k,l) ;
            }
            prefactor = hatfactor * modelspace->phase(ji+jj-J1+Lambda) ;
            for ( int a : Y.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2/2.0;
               cijkl -= prefactor * modelspace->GetSixJ(J2,J1,Lambda,jj,ja,ji) * Y.OneBody(j,a) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,a,k,l);
            }
            prefactor = hatfactor * modelspace->phase(jk+jl+J2+Lambda) ;
            for ( int a : Y.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2/2.0;
               cijkl += prefactor * modelspace->GetSixJ(J1,J2,Lambda,jk,ja,jl) * Y.OneBody(a,k) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,a,l) ;
            }
            prefactor = hatfactor * modelspace->phase(jk+jl-J1+Lambda) ;
            for ( int a : Y.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2/2.0;
               cijkl += prefactor * modelspace->GetSixJ(J1,J2,Lambda,jl,ja,jk) * Y.OneBody(a,l) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,k,a) ;
            }


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            cijkl /= norm;
            Z2(ibra,iket) += cijkl;
         }
      }
   }
}





// Since comm222_pp_hh and comm211 both require the ruction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
// X is a scalar, Y is a tensor
//void Operator::comm222_pp_hh_221st( Operator& Y, Operator& Z )  
void Operator::comm222_pp_hh_221st( const Operator& X, const Operator& Y )  
{

   Operator& Z = *this;
   int Lambda = Z.GetJRank();

   TwoBodyME Mpp = Z.TwoBody;
   TwoBodyME Mhh = Z.TwoBody;

//   #pragma omp parallel for schedule(dynamic,5)
   for ( auto& itmat : Y.TwoBody.MatEl )
   {
    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];

    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);

    auto& LHS1 = X.TwoBody.GetMatrix(ch_bra,ch_bra);
    auto& LHS2 = X.TwoBody.GetMatrix(ch_ket,ch_ket);

    auto& RHS  =  itmat.second;
    arma::mat& OUT2 =    Z.TwoBody.GetMatrix(ch_bra,ch_ket);

    arma::mat& Matrixpp =  Mpp.GetMatrix(ch_bra,ch_ket);
    arma::mat& Matrixhh =  Mhh.GetMatrix(ch_bra,ch_ket);
   
    arma::uvec& bras_pp = tbc_bra.GetKetIndex_pp();
    arma::uvec& bras_hh = tbc_bra.GetKetIndex_hh();
    arma::uvec& kets_pp = tbc_ket.GetKetIndex_pp();
    arma::uvec& kets_hh = tbc_ket.GetKetIndex_hh();
    
    Matrixpp =  LHS1.cols(bras_pp) * RHS.rows(bras_pp) - RHS.cols(kets_pp)*LHS2.rows(kets_pp);
    Matrixhh =  LHS1.cols(bras_hh) * RHS.rows(bras_hh) - RHS.cols(kets_hh)*LHS2.rows(kets_hh);
 

    // Now, the two body part is easy
    OUT2 += Matrixpp - Matrixhh;

   }// for itmat

      // The one body part takes some additional work

   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for schedule(dynamic,1)
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      double ji = oi.j2/2.0;
      for (int j : Z.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit &oj = modelspace->GetOrbit(j);
         double jj = oj.j2/2.0;
         double cijJ = 0;
         // Sum c over holes and include the nbar_a * nbar_b terms
           for (auto& c : modelspace->holes)
           {
              Orbit &oc = modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = max(int(abs(jc-jj)),abs(Lambda-J1));
               int j2max = min(int(jc+jj),J1+Lambda); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * modelspace->phase(jj + jc + J1 + Lambda) * Mpp.GetTBME_J(J1,J2,c,i,c,j);
               }
              }
           // Sum c over particles and include the n_a * n_b terms
           }
           for (auto& c : modelspace->particles)
           {
              Orbit &oc = modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = max(int(abs(jc-jj)),abs(Lambda-J1));
               int j2max = min(int(jc+jj),J1+Lambda); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * modelspace->phase(jj + jc + J1 + Lambda) * Mhh.GetTBME_J(J1,J2,c,i,c,j);
               }
              }
           }
         #pragma omp critical
         Z.OneBody(i,j) += cijJ ;
      } // for j
    } // for i
}






//**************************************************************************
//
//  X^J_ij`kl` = - sum_J' { i j J } (2J'+1) X^J'_ilkj
//                        { k l J'}
// TENSOR VARIETY
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
//void Operator::DoTensorPandyaTransformation(vector<arma::mat>& TwoBody_CC_hp, vector<arma::mat>& TwoBody_CC_ph)
void Operator::DoTensorPandyaTransformation(map<array<int,2>,arma::mat>& TwoBody_CC_hp, map<array<int,2>,arma::mat>& TwoBody_CC_ph) const
{
   int Lambda = rank_J;
   // loop over cross-coupled channels
//   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
   for ( int ch_bra_cc : modelspace->SortedTwoBodyChannels_CC )
   {
      TwoBodyChannel& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      int Jbra_cc = tbc_bra_cc.J;
      arma::uvec& bras_ph = tbc_bra_cc.GetKetIndex_ph();
      int nph_bras = bras_ph.n_rows;

      for ( int ch_ket_cc : modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        int Jket_cc = tbc_ket_cc.J;
        if ( (Jbra_cc+Jket_cc < rank_J) or abs(Jbra_cc-Jket_cc)>rank_J ) continue;
        if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + parity)%2>0 ) continue;
        if ( abs(tbc_bra_cc.Tz + tbc_ket_cc.Tz)%2 > 0 ) continue; // need an even number of protons, even number of neutrons.

        int nKets_cc = tbc_ket_cc.GetNumberKets();

        // Need to make these maps to account for ch_bra and ch_ket.
        TwoBody_CC_hp[{ch_bra_cc,ch_ket_cc}] = arma::mat(nph_bras,   2*nKets_cc, arma::fill::zeros);
        TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}] = arma::mat(nph_bras,   2*nKets_cc, arma::fill::zeros);

        // loop over ph bras <ad| in this channel
        for (int ibra=0; ibra<nph_bras; ++ibra)
        {
           Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
           int a = bra_cc.p;
           int b = bra_cc.q;
           Orbit & oa = modelspace->GetOrbit(a);
           Orbit & ob = modelspace->GetOrbit(b);
           double ja = oa.j2*0.5;
           double jb = ob.j2*0.5;

           // loop over kets |bc> in this channel
           // we go to 2*nKets to include |bc> and |cb>
           for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
           {
              Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
              int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
              int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
              Orbit & oc = modelspace->GetOrbit(c);
              Orbit & od = modelspace->GetOrbit(d);
              double jc = oc.j2*0.5;
              double jd = od.j2*0.5;


              int j1min = abs(ja-jd);
              int j1max = ja+jd;
              double sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = max(int(abs(jc-jb)),abs(J1-Lambda));
                int j2max = min(int(jc+jb),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = TwoBody.GetTBME_J(J1,J2,a,d,c,b);
                  sm -= hatfactor * modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              TwoBody_CC_hp[{ch_bra_cc,ch_ket_cc}](ibra,iket_cc) = sm;


              // Exchange (a <-> b) to account for the (n_a - n_b) term
              // Get Tz,parity and range of J for <bd || ca > coupling
              j1min = abs(jb-jd);
              j1max = jb+jd;
              sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = max(int(abs(jc-ja)),abs(J1-Lambda));
                int j2max = min(int(jc+ja),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = TwoBody.GetTBME_J(J1,J2,b,d,c,a);
                  sm -= hatfactor * modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
              TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}](ibra,iket_cc) = sm;

           }
        }
    }
   }
}


void Operator::AddInverseTensorPandyaTransformation(map<array<int,2>,arma::mat>& Zbar)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
//   for (int ch=0;ch<nChannels;++ch)
//   int n_nonzeroChannels = modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())
//   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   Operator& Z = *this;
   int Lambda = Z.rank_J;
   for (auto& iter : Z.TwoBody.MatEl)
   {
      int ch_bra = iter.first[0];
      int ch_ket = iter.first[1];
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nBras = tbc_bra.GetNumberKets();
      int nKets = tbc_ket.GetNumberKets();
      arma::mat& Zijkl = iter.second;

      for (int ibra=0; ibra<nBras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
//         int ketmin = hermitian ? ibra : ibra+1;
         int ketmin = 0;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int parity_bra_cc = (oi.l+ol.l)%2;
            int parity_ket_cc = (ok.l+oj.l)%2;
            int Tz_bra_cc = abs(oi.tz2+ol.tz2)/2;
            int Tz_ket_cc = abs(ok.tz2+oj.tz2)/2;
            int j3min = abs(int(ji-jl));
            int j3max = ji+jl;

            for (int J3=j3min; J3<=j3max; ++J3)
            {
              int ch_bra_cc = modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
              TwoBodyChannel_CC& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
              int indx_il = tbc_bra_cc.GetLocalIndex(min(i,l),max(i,l));
              if (i>l) indx_il += tbc_bra_cc.GetNumberKets();
              int j4min = max(abs(int(jk-jj)),abs(J3-Lambda));
              int j4max = min(int(jk+jj),J3+Lambda);
              for (int J4=j4min; J4<=j4max; ++J4)
              {
                 int ch_ket_cc = modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                 TwoBodyChannel_CC& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                 int indx_kj = tbc_ket_cc.GetLocalIndex(min(j,k),max(j,k));
                 if (k>j) indx_kj += tbc_ket_cc.GetNumberKets();

                  double ninej = modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
                  if (abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                  double tbme = Zbar[{ch_bra_cc,ch_ket_cc}](indx_il,indx_kj);
                  commij += hatfactor * modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
              }
            }

            // Transform Z_jlki
            parity_bra_cc = (oj.l+ol.l)%2;
            parity_ket_cc = (ok.l+oi.l)%2;
            Tz_bra_cc = abs(oj.tz2+ol.tz2)/2;
            Tz_ket_cc = abs(ok.tz2+oi.tz2)/2;
            j3min = abs(int(jj-jl));
            j3max = jj+jl;

            for (int J3=j3min; J3<=j3max; ++J3)
            {
              int ch_bra_cc = modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
              TwoBodyChannel_CC& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
              int indx_jl = tbc_bra_cc.GetLocalIndex(min(j,l),max(j,l));
              if (j>l) indx_jl += tbc_bra_cc.GetNumberKets();
              int j4min = max(abs(int(jk-ji)),abs(J3-Lambda));
              int j4max = min(int(jk+ji),J3+Lambda);
              for (int J4=j4min; J4<=j4max; ++J4)
              {
                 int ch_ket_cc = modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                 TwoBodyChannel_CC& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                 int indx_ki = tbc_ket_cc.GetLocalIndex(min(k,i),max(k,i));
                 if (k>i) indx_ki += tbc_ket_cc.GetNumberKets();

                  double ninej = modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                  if (abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                  double tbme = Zbar[{ch_bra_cc,ch_ket_cc}](indx_jl,indx_ki);
                  commji += hatfactor * modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Zijkl(ibra,iket) += (commij - modelspace->phase(ji+jj-J1)*commji) / norm;
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
/// Calculates the part of \f$ [X_{(2)},\mathbb{Y}^{\Lambda}_{(2)}]_{ijkl} \f$ which involves ph intermediate states, here indicated by \f$ \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} \f$
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{abJ_3J_4}(n_a-n_b) \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
///  -(-1)^{j_i+j_j-J_1}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
/// \left( \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// \right]
/// \f]
/// This is implemented by defining an intermediate matrix
/// \f[
/// \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}} \equiv \sum_{ab}(n_a\bar{n}_b)
/// \left[ \left( \bar{X}^{J3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} - 
///   \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_4}_{a\bar{b}k\bar{j}} \right)
/// -\left( \bar{X}^{J_3}_{i\bar{l}b\bar{a}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{b\bar{a}k\bar{j}} - 
///    \bar{\mathbb{Y}}^{J_3J_4\Lambda}_{i\bar{l}b\bar{a}}\bar{X}^{J_4}_{b\bar{a}k\bar{j}} \right)\right]
/// \f]
/// The Pandya-transformed matrix elements are obtained with DoTensorPandyaTransformation().
/// The matrices \f$ \bar{X}^{J_3}_{i\bar{l}a\bar{b}}\bar{\mathbb{Y}}^{J_3J_4\Lambda}_{a\bar{b}k\bar{j}} \f$
/// and \f$ \bar{\mathbb{Y}}^{J_4J_3\Lambda}_{i\bar{l}a\bar{b}}\bar{X}^{J_3}_{a\bar{b}k\bar{j}} \f$
/// are related by a Hermitian conjugation, which saves two matrix multiplications, provided we
/// take into account the phase \f$ (-1)^{J_3-J_4} \f$ from conjugating the spherical tensor.
/// The commutator is then given by
/// \f[
/// \mathbb{Z}^{J_1J_2\Lambda}_{ijkl} = \sum_{J_3J_4} \hat{J_1}\hat{J_2}\hat{J_3}\hat{J_4}
/// \left[
///  \left\{ \begin{array}{lll}
///  j_i  &  j_l  &  J_3 \\
///  j_j  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{i\bar{l}k\bar{j}}
///  -(-1)^{j_i+j_j-J}
///  \left\{ \begin{array}{lll}
///  j_j  &  j_l  &  J_3 \\
///  j_i  &  j_k  &  J_4 \\
///  J_1  &  J_2  &  \Lambda \\
///  \end{array} \right\}
///  \bar{\mathbb{Z}}^{J_3J_4\Lambda}_{j\bar{l}k\bar{i}}
///  \right]
///  \f]
///
//void Operator::comm222_phst( Operator& Y, Operator& Z ) 
void Operator::comm222_phst( const Operator& X, const Operator& Y ) 
{

   Operator& Z = *this;
   // Create Pandya-transformed hp and ph matrix elements
   vector<arma::mat> X_bar_hp (nChannels );
   vector<arma::mat> X_bar_ph (nChannels );

   map<array<int,2>,arma::mat> Y_bar_hp;
   map<array<int,2>,arma::mat> Y_bar_ph;


   double t = omp_get_wtime();
   X.DoPandyaTransformation(X_bar_hp, X_bar_ph );
   Y.DoTensorPandyaTransformation(Y_bar_hp, Y_bar_ph );
   timer["DoTensorPandyaTransformation"] += omp_get_wtime() - t;


   t = omp_get_wtime();
   // Construct the intermediate matrix Z_bar
   map<array<int,2>,arma::mat> Z_bar;

   for (auto& iter : Y_bar_hp )
   {
      int ch_bra_cc = iter.first[0];
      int ch_ket_cc = iter.first[1];
      int Jbra = modelspace->GetTwoBodyChannel_CC(ch_bra_cc).J;
      int Jket = modelspace->GetTwoBodyChannel_CC(ch_ket_cc).J;
      int flipphase = modelspace->phase( Jbra - Jket ) * ( Z.IsHermitian() ? -1 : 1 );
      if (X.IsHermitian())
      {
        Z_bar[{ch_bra_cc,ch_ket_cc}] =  ( X_bar_hp[ch_bra_cc].t() * Y_bar_hp[{ch_bra_cc,ch_ket_cc}] - X_bar_ph[ch_bra_cc].t() * Y_bar_ph[{ch_bra_cc,ch_ket_cc}]) 
                           -flipphase * ( X_bar_hp[ch_ket_cc].t() * Y_bar_hp[{ch_ket_cc,ch_bra_cc}] - X_bar_ph[ch_ket_cc].t() * Y_bar_ph[{ch_ket_cc,ch_bra_cc}]).t() ;
      }
      else
      {
        Z_bar[{ch_bra_cc,ch_ket_cc}] =  ( X_bar_ph[ch_bra_cc].t() * Y_bar_ph[{ch_bra_cc,ch_ket_cc}]  - X_bar_hp[ch_bra_cc].t() * Y_bar_hp[{ch_bra_cc,ch_ket_cc}] )
                           -flipphase * ( X_bar_ph[ch_ket_cc].t() * Y_bar_ph[{ch_ket_cc,ch_bra_cc}]  - X_bar_hp[ch_ket_cc].t() * Y_bar_hp[{ch_ket_cc,ch_bra_cc}]).t() ;
      }
   }
   timer["Build Z_bar"] += omp_get_wtime() - t;

   t = omp_get_wtime();
   Z.AddInverseTensorPandyaTransformation(Z_bar);
   timer["InverseTensorPandyaTransformation"] += omp_get_wtime() - t;

}









