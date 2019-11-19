
#include "Operator.hh"
#include "AngMom.hh"
#include "IMSRGProfiler.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <deque>
#include <gsl/gsl_math.h>
#include <math.h>
#include "omp.h"

#ifndef SQRT2
  #define SQRT2 1.4142135623730950488L
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
bool Operator::use_brueckner_bch = false;
bool Operator::use_goose_tank_correction = false;
bool Operator::use_goose_tank_correction_titus = false;

Operator& Operator::TempOp(size_t n)
{
  static deque<Operator> TempArray;
  if (n >= TempArray.size()) TempArray.resize(n+1,*this);
  return TempArray[n];
}

//vector<arma::mat>& Operator::TempMatVec(size_t n)
//{
//  static deque<vector<arma::mat>> TempMatVecArray;
//  if (n>= TempMatVecArray.size()) TempMatVecArray.resize(max(n,(size_t)5));
//  return TempMatVecArray[n];
//}

//////////////////// DESTRUCTOR //////////////////////////////////////////
Operator::~Operator()
{
  profiler.counter["N_Operators"] --;
}

/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator()
 :   modelspace(NULL), 
    rank_J(0), rank_T(0), parity(0), particle_rank(2),
    hermitian(true), antihermitian(false), nChannels(0)
{
  profiler.counter["N_Operators"] ++;
}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms,Jrank,Trank,p),  ThreeBody(&ms),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank),
    E3max(ms.GetE3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels()) 
{
  SetUpOneBodyChannels();
  if (particle_rank >=3) ThreeBody.Allocate();
  profiler.counter["N_Operators"] ++;
}

Operator::Operator(ModelSpace& ms) :
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms),  ThreeBody(&ms),
    rank_J(0), rank_T(0), parity(0), particle_rank(2),
    E3max(ms.GetE3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels())
{
  SetUpOneBodyChannels();
  profiler.counter["N_Operators"] ++;
}

Operator::Operator(const Operator& op)
: modelspace(op.modelspace),  ZeroBody(op.ZeroBody),
  OneBody(op.OneBody), TwoBody(op.TwoBody) ,ThreeBody(op.ThreeBody),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank),
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels)
{
  profiler.counter["N_Operators"] ++;
}

Operator::Operator(Operator&& op)
: modelspace(op.modelspace), ZeroBody(op.ZeroBody),
  OneBody(move(op.OneBody)), TwoBody(move(op.TwoBody)) , ThreeBody(move(op.ThreeBody)),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank),
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels)
{
  profiler.counter["N_Operators"] ++;
}



/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////

Operator& Operator::operator=(const Operator& rhs) = default;
Operator& Operator::operator=(Operator&& rhs) = default;

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
   if (rhs.GetParticleRank() > 1)
     TwoBody  += rhs.TwoBody;
   return *this;
}

Operator Operator::operator+(const Operator& rhs) const
{
   if (GetParticleRank() >= rhs.GetParticleRank())
     return ( Operator(*this) += rhs );
   else
     return ( Operator(rhs) += *this );
}

Operator& Operator::operator+=(const double& rhs)
{
   ZeroBody += rhs;
   return *this;
}

Operator Operator::operator+(const double& rhs) const
{
   return ( Operator(*this) += rhs );
}

// Subtract operators
Operator& Operator::operator-=(const Operator& rhs)
{
   ZeroBody -= rhs.ZeroBody;
   OneBody -= rhs.OneBody;
   if (rhs.GetParticleRank() > 1)
     TwoBody -= rhs.TwoBody;
   return *this;
}

Operator Operator::operator-(const Operator& rhs) const
{
   return ( Operator(*this) -= rhs );
}

Operator& Operator::operator-=(const double& rhs)
{
   ZeroBody -= rhs;
   return *this;
}

Operator Operator::operator-(const double& rhs) const
{
   return ( Operator(*this) -= rhs );
}

// Negation operator
Operator Operator::operator-() const
{
   return (*this)*-1.0;
}



void Operator::SetUpOneBodyChannels()
{
  for ( int i=0; i<modelspace->GetNumberOrbits(); ++i )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    // The +-1 comes from the spin [LxS](J)
    int lmin = max( oi.l - rank_J-1, 0);
    int lmax = min( oi.l + rank_J+1, modelspace->GetEmax() );
    for (int l=lmin; l<=lmax; l+=1)
    {
      if ((l + oi.l + parity)%2>0) continue;
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
  for (auto& it: OneBodyChannels)  it.second.shrink_to_fit();
}


size_t Operator::Size()
{
   return sizeof(ZeroBody) + OneBody.size()*sizeof(double) + TwoBody.size() + ThreeBody.size();
}


void Operator::SetOneBody(int i, int j, double val)
{
 OneBody(i,j) = val;
 if ( IsNonHermitian() ) return;
 int flip_phase = IsHermitian() ? 1 : -1;
 if (rank_J > 0)
   flip_phase *= modelspace->phase( (modelspace->GetOrbit(i).j2 - modelspace->GetOrbit(j).j2) / 2 );
 OneBody(j,i) = flip_phase * val;

}

void Operator::SetTwoBody(int J1, int p1, int T1, int J2, int p2, int T2, int i, int j, int k, int l, double v)
{
  TwoBody.SetTBME( J1,  p1,  T1,  J2,  p2,  T2,  i,  j,  k,  l, v);
}

double Operator::GetTwoBody(int ch_bra, int ch_ket, int ibra, int iket)
{
  if ( ch_bra <= ch_ket or IsNonHermitian() )
  {
    return TwoBody.GetMatrix(ch_bra, ch_ket)(ibra,iket);
  }
  else
  {
    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int flipphase = modelspace->phase( tbc_bra.J - tbc_ket.J) * ( IsHermitian() ? 1 : -1 ) ;
    return  flipphase * TwoBody.GetMatrix(ch_ket,ch_bra)(iket,ibra);
  }
}


void Operator::WriteBinary(ofstream& ofs)
{
  double tstart = omp_get_wtime();
  ofs.write((char*)&rank_J,sizeof(rank_J));
  ofs.write((char*)&rank_T,sizeof(rank_T));
  ofs.write((char*)&parity,sizeof(parity));
  ofs.write((char*)&particle_rank,sizeof(particle_rank));
  ofs.write((char*)&E2max,sizeof(E2max));
  ofs.write((char*)&E3max,sizeof(E3max));
  ofs.write((char*)&hermitian,sizeof(hermitian));
  ofs.write((char*)&antihermitian,sizeof(antihermitian));
  ofs.write((char*)&nChannels,sizeof(nChannels));
  ofs.write((char*)&ZeroBody,sizeof(ZeroBody));
  ofs.write((char*)OneBody.memptr(),OneBody.size()*sizeof(double));
  if (particle_rank > 1)
    TwoBody.WriteBinary(ofs);
  if (particle_rank > 2)
    ThreeBody.WriteBinary(ofs);
  profiler.timer["Write Binary Op"] += omp_get_wtime() - tstart;
}


void Operator::ReadBinary(ifstream& ifs)
{
  double tstart = omp_get_wtime();
  ifs.read((char*)&rank_J,sizeof(rank_J));
  ifs.read((char*)&rank_T,sizeof(rank_T));
  ifs.read((char*)&parity,sizeof(parity));
  ifs.read((char*)&particle_rank,sizeof(particle_rank));
  ifs.read((char*)&E2max,sizeof(E2max));
  ifs.read((char*)&E3max,sizeof(E3max));
  ifs.read((char*)&hermitian,sizeof(hermitian));
  ifs.read((char*)&antihermitian,sizeof(antihermitian));
  ifs.read((char*)&nChannels,sizeof(nChannels));
  SetUpOneBodyChannels();
  ifs.read((char*)&ZeroBody,sizeof(ZeroBody));
  ifs.read((char*)OneBody.memptr(),OneBody.size()*sizeof(double));
  if (particle_rank > 1)
    TwoBody.ReadBinary(ifs);
  if (particle_rank > 2)
    ThreeBody.ReadBinary(ifs);
  profiler.timer["Read Binary Op"] += omp_get_wtime() - tstart;
}





////////////////// MAIN INTERFACE METHODS //////////////////////////

Operator Operator::DoNormalOrdering()
{
   if (particle_rank==3)
      return DoNormalOrdering3();
   else
      return DoNormalOrdering2();
}

//*************************************************************
///  Normal ordering of a 2body operator
///  set up for scalar or tensor operators, but
///  the tensor part hasn't been tested
//*************************************************************
Operator Operator::DoNormalOrdering2()
{
   Operator opNO(*this);
   bool scalar = (opNO.rank_J==0 and opNO.rank_T==0 and opNO.parity==0);
   if (scalar)
   {
     for (auto& k : modelspace->holes) // loop over hole orbits
     {
        Orbit& ok = modelspace->GetOrbit(k);
        opNO.ZeroBody += (ok.j2+1) * ok.occ * OneBody(k,k);
     }
   }
   cout << "OneBody contribution: " << opNO.ZeroBody << endl;

   index_t norbits = modelspace->GetNumberOrbits();
   if (TwoBody.Norm() > 1e-7)
   {
     for ( auto& itmat : TwoBody.MatEl )
     {
        int ch_bra = itmat.first[0];
        int ch_ket = itmat.first[1];
        auto& matrix = itmat.second;
        
        TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
        int J_bra = tbc_bra.J;
        int J_ket = tbc_ket.J;
        double hatfactor = sqrt((2*J_bra+1.0)*(2*J_ket+1.0));
  
        // Zero body part
        if (scalar)
        {
          arma::vec diagonals = matrix.diag();
          auto hh = tbc_ket.GetKetIndex_hh();
          auto hocc = tbc_ket.Ket_occ_hh;
          opNO.ZeroBody += arma::sum( hocc % diagonals.elem(hh) ) * hatfactor;
        }
  
        // One body part
        for (index_t a=0;a<norbits;++a)
        {
           Orbit &oa = modelspace->GetOrbit(a);
           double ja = oa.j2/2.0;
           index_t bstart = (IsNonHermitian() or ch_bra!=ch_ket )? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
           for ( auto& b : opNO.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
           {
              if (b < bstart) continue;
              Orbit &ob = modelspace->GetOrbit(b);
              double jb = ob.j2/2.0;
              for (auto& h : modelspace->holes)  // C++11 syntax
              {
                Orbit &oh = modelspace->GetOrbit(h);
                if (opNO.rank_J==0)
                {
                   opNO.OneBody(a,b) += hatfactor /(2*ja+1.0) * oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
                }
                else
                {
                   double jh = oh.j2*0.5;
                   if ((ja+jh < J_bra) or (std::abs(ja-jh)>J_bra) or (jb+jh < J_ket) or (std::abs(jb-jh)>J_ket) ) continue;
                   if ((oa.l + oh.l + tbc_bra.parity)%2 >0) continue;
                   if ((ob.l + oh.l + tbc_ket.parity)%2 >0) continue;
                   if ((oa.tz2 + oh.tz2) != tbc_bra.Tz*2) continue;
                   if ((ob.tz2 + oh.tz2) != tbc_ket.Tz*2) continue;
                   double ME = hatfactor  * oh.occ *modelspace->phase(ja+jh-J_ket-opNO.rank_J)
                                   * modelspace->GetSixJ(J_bra,J_ket,opNO.rank_J,jb,ja,jh) * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
                   if (a>b)
                   {
                     int herm = IsHermitian() ? 1 : -1;
                     opNO.OneBody(b,a) += herm * modelspace->phase(ja-jb) * ME;
                   }
                   else
                   {
                     opNO.OneBody(a,b) += ME;
                   }
                }
             }
           }
        }
     } // loop over channels
//     cout << "------------------------------------------" << endl;
   }

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
               int kmin2 = std::abs(2*tbc.J-oa.j2);
               int kmax2 = 2*tbc.J+oa.j2;
               for (int K2=kmin2; K2<=kmax2; K2+=2)
               {
                  Gamma(ibra,iket) += (K2+1) * oa.occ * ThreeBody.GetME_pn(tbc.J,tbc.J,K2,i,j,a,k,l,a); // This is unnormalized, but it should be normalized!!!!
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



/// Convert to a basis normal ordered wrt the vacuum.
/// This doesn't handle 3-body terms. In that case,
/// the 2-body piece is unchanged.
Operator Operator::UndoNormalOrdering() const
{
   Operator opNO = *this;
//   cout << "Undoing Normal ordering. Initial ZeroBody = " << opNO.ZeroBody << endl;

   if (opNO.GetJRank()==0 and opNO.GetTRank()==0 and opNO.GetParity()==0)
   {
     for (auto& k : modelspace->holes) // loop over hole orbits
     {
        Orbit& ok = modelspace->GetOrbit(k);
        opNO.ZeroBody -= (ok.j2+1) * ok.occ * OneBody(k,k);
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
      double hatfactor = sqrt((2*J_bra+1.0)*(2*J_ket+1.0));

      // Zero body part
      if (opNO.GetJRank()==0 and opNO.GetTRank()==0 and opNO.GetParity()==0)
      {
        arma::vec diagonals = matrix.diag();
        auto hh = tbc_ket.GetKetIndex_hh();
        auto hocc = tbc_ket.Ket_occ_hh;
        opNO.ZeroBody +=  arma::sum( hocc % diagonals.elem(hh) ) *hatfactor;
      }

      // One body part
      for (index_t a=0;a<norbits;++a)
      {
         Orbit &oa = modelspace->GetOrbit(a);
         double ja = oa.j2*0.5;
//         index_t bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
           index_t bstart = (IsNonHermitian() or ch_bra!=ch_ket )? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for ( auto& b : opNO.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
         {
            if (b < bstart) continue;
            Orbit &ob = modelspace->GetOrbit(b);
            double jb = ob.j2*0.5;
            for (auto& h : modelspace->holes)  // C++11 syntax
            {
              Orbit& oh = modelspace->GetOrbit(h);

              if (opNO.rank_J==0)
              {
                 opNO.OneBody(a,b) -= hatfactor/(2*ja+1) * oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
              }
              else
              {
                 double jh = oh.j2*0.5;
                 if ((ja+jh < J_bra) or (std::abs(ja-jh)>J_bra) or (jb+jh < J_ket) or (std::abs(jb-jh)>J_ket) ) continue;

                 if ((oa.l + oh.l + tbc_bra.parity)%2 >0) continue;
                 if ((ob.l + oh.l + tbc_ket.parity)%2 >0) continue;
                 if ((oa.tz2 + oh.tz2) != tbc_bra.Tz*2) continue;
                 if ((ob.tz2 + oh.tz2) != tbc_ket.Tz*2) continue;
                 double ME = hatfactor  * oh.occ *modelspace->phase(ja+jh-J_ket-opNO.rank_J)
                                             * modelspace->GetSixJ(J_bra,J_ket,opNO.rank_J,jb,ja,jh) * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
                 if (a>b)
                 {
                   int herm = IsHermitian() ? 1 : -1;
                   opNO.OneBody(b,a) -= herm * modelspace->phase(ja-jb) * ME;
                 }
                 else
                 {
                   opNO.OneBody(a,b) -= ME;
                 }

              }
            }

         }
      }
   } // loop over channels

   if (hermitian) opNO.Symmetrize();
   if (antihermitian) opNO.AntiSymmetrize();

   return opNO;

}

//********************************************
/// Truncate an operator to a smaller emax
/// A corresponding ModelSpace object must be
/// created at the appropriate scope. That's why
/// the new operator is passed as a 
//********************************************
Operator Operator::Truncate(ModelSpace& ms_new)
{
  Operator OpNew(ms_new, rank_J, rank_T, parity, particle_rank);
  
  int new_emax = ms_new.GetEmax();
  if ( new_emax > modelspace->GetEmax() )
  {
    cout << "Error: Cannot truncate an operator with emax = " << modelspace->GetEmax() << " to one with emax = " << new_emax << endl;
    return OpNew;
  }
//  OpNew.rank_J=rank_J; 
//  OpNew.rank_T=rank_T; 
//  OpNew.parity=parity; 
//  OpNew.particle_rank = particle_rank; 
  OpNew.ZeroBody = ZeroBody;
  OpNew.hermitian = hermitian;
  OpNew.antihermitian = antihermitian;
  int norb = ms_new.GetNumberOrbits();
  OpNew.OneBody = OneBody.submat(0,0,norb-1,norb-1);
  for (auto& itmat : OpNew.TwoBody.MatEl )
  {
    int ch = itmat.first[0];
    TwoBodyChannel& tbc_new = ms_new.GetTwoBodyChannel(ch);
    int chold = modelspace->GetTwoBodyChannelIndex(tbc_new.J,tbc_new.parity,tbc_new.Tz);
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(chold);
    auto& Mat_new = itmat.second;
    auto& Mat = TwoBody.GetMatrix(chold,chold);
    int nkets = tbc_new.GetNumberKets();
    arma::uvec ibra_old(nkets);
    for (int ibra=0;ibra<nkets;++ibra)
    {
      ibra_old(ibra) = tbc.GetLocalIndex(tbc_new.GetKetIndex(ibra));
    }
    Mat_new = Mat.submat(ibra_old,ibra_old);
  }
  return OpNew;
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

void Operator::EraseThreeBody()
{
  ThreeBody = ThreeBodyME();
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

void Operator::MakeReduced()
{
  if (rank_J>0)
  {
    cout << "Trying to reduce an operator with J rank = " << rank_J << ". Not good!!!" << endl;
    return;
  }
  for ( int a=0;a<modelspace->GetNumberOrbits();++a )
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
    {
      if (b<a) continue;
      OneBody(a,b) *= sqrt(oa.j2+1);
      OneBody(b,a) = OneBody(a,b);
    }
  }
  for ( auto& itmat : TwoBody.MatEl )
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(itmat.first[0]);
    itmat.second *= sqrt(2*tbc.J+1);
  }
}

void Operator::MakeNotReduced()
{
  if (rank_J>0)
  {
    cout << "Trying to un-reduce an operator with J rank = " << rank_J << ". Not good!!!" << endl;
    return;
  }
  for ( int a=0;a<modelspace->GetNumberOrbits();++a )
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
    {
      if (b<a) continue;
      OneBody(a,b) /= sqrt(oa.j2+1);
      OneBody(b,a) = OneBody(a,b);
    }
  }
  for ( auto& itmat : TwoBody.MatEl )
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(itmat.first[0]);
    itmat.second /= sqrt(2*tbc.J+1);
  }
}



//// this routine then multiplies the TBME <ab|Op|cd> by coeff if a==b, and again if c==d
//void Operator::ChangeNormalization( double coeff )
//{
//  for (auto& it_mat : TwoBody.MatEl )
//  {
//    int ch_bra = it_mat.first[0];
//    int ch_ket = it_mat.first[1];
//    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
//    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
//    int nbras = tbc_bra.GetNumberKets();
//    int nkets = tbc_ket.GetNumberKets();
//    for (int ibra=0; ibra<nbras; ++ibra)
//    {
//      Ket& bra = tbc_bra.GetKet(ibra);
//      if ( bra.p == bra.q ) it_mat.second.row(ibra) *= coeff;
//    }
//    for (int iket=0; iket<nkets; ++iket)
//    {
//      Ket& ket = tbc_ket.GetKet(iket);
//      if ( ket.p == ket.q ) it_mat.second.col(iket) *= coeff;
//    }
//  }
//
//}



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

// This is unused
void Operator::Eye()
{
   ZeroBody = 1;
   OneBody.eye();
   TwoBody.Eye();
}


/// Calculate the second-order perturbation theory correction to the energy
/// \f[
/// E^{(2)} = \sum_{ia} (2 j_a +1) \frac{|f_{ia}|^2}{f_{aa}-f_{ii}}
/// +  \sum_{\substack{i\leq j // a\leq b}}\sum_{J} (2J+1)\frac{|\Gamma_{ijab}^{J}|^2}{f_{aa}+f_{bb}-f_{ii}-f_{jj}}
/// \f]
///
double Operator::GetMP2_Energy()
{
   double t_start = omp_get_wtime();
   double Emp2 = 0;
   int nparticles = modelspace->particles.size();
   #pragma omp parallel for reduction(+:Emp2)
   for ( int ii=0;ii<nparticles;++ii)
   {
     index_t i = modelspace->particles[ii];
     double ei = OneBody(i,i);
     Orbit& oi = modelspace->GetOrbit(i);
     for (auto& a : modelspace->holes)
     {
       Orbit& oa = modelspace->GetOrbit(a);
       double ea = OneBody(a,a);
       if (std::abs(OneBody(i,a))>1e-6)
         Emp2 += (oa.j2+1) * oa.occ * OneBody(i,a)*OneBody(i,a)/(OneBody(a,a)-OneBody(i,i));
       for (index_t j : modelspace->particles)
       {
         if (j<i) continue;
         double ej = OneBody(j,j);
         Orbit& oj = modelspace->GetOrbit(j);
         for ( auto& b: modelspace->holes)
         {
           if (b<a) continue;
           Orbit& ob = modelspace->GetOrbit(b);
           double eb = OneBody(b,b);
           double denom = ea+eb-ei-ej;
           int Jmin = max(std::abs(oi.j2-oj.j2),std::abs(oa.j2-ob.j2))/2;
           int Jmax = min(oi.j2+oj.j2,oa.j2+ob.j2)/2;
           for (int J=Jmin; J<=Jmax; ++J)
           {
             double tbme = TwoBody.GetTBME_J_norm(J,a,b,i,j);
             if (std::abs(tbme)>1e-6)
              Emp2 += (2*J+1)* oa.occ * ob.occ * tbme*tbme/denom; // no factor 1/4 because of the restricted sum
           }
         }
       }
     }
   }
   profiler.timer["GetMP2_Energy"] += omp_get_wtime() - t_start;
   return Emp2;
}


//*************************************************************
/// Calculate the third order perturbation correction to energy
/// \f[
/// \frac{1}{8}\sum_{abijpq}\sum_J (2J+1)\frac{\Gamma_{abij}^J\Gamma_{ijpq}^J\Gamma_{pqab}^J}{(f_a+f_b-f_p-f_q)(f_a+f_b-f_i-f_j)}
/// +\sum_{abcijk}\sum_J (2J+1)^2 \frac{\bar{\Gamma}_{a\bar{i}j\bar{b}}^J\bar{\Gamma}_{j\bar{b}k\bar{c}}^J\bar{\Gamma}_{k\bar{c}a\bar{i}}^J}{(f_a+f_b-f_i-f_j)(f_a+f_c-f_i-f_k)}
/// \f]
///
//*************************************************************
double Operator::GetMP3_Energy()
{
   // So far, the pp and hh parts seem to work. No such luck for the ph.
   double t_start = omp_get_wtime();
   double Emp3 = 0;
   // This can certainly be optimized, but I'll wait until this is the bottleneck.
   int nch = modelspace->GetNumberTwoBodyChannels();

   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Emp3)
   for (int ich=0;ich<nch;++ich)
   {
     TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ich);
     auto& Mat = TwoBody.GetMatrix(ich,ich);
     int J = tbc.J;
     for (auto iket_ab : tbc.GetKetIndex_hh() )
     {
       Ket& ket_ab = tbc.GetKet(iket_ab);
       index_t a = ket_ab.p;
       index_t b = ket_ab.q;

       for (auto iket_ij : tbc.GetKetIndex_pp() )
       {
         Ket& ket_ij = tbc.GetKet(iket_ij);
         index_t i = ket_ij.p;
         index_t j = ket_ij.q;
         double Delta_abij = OneBody(a,a) + OneBody(b,b) - OneBody(i,i) - OneBody(j,j);


       // hh term
         for (auto iket_cd : tbc.GetKetIndex_hh() )
         {
           Ket& ket_cd = tbc.GetKet(iket_cd);
           index_t c = ket_cd.p;
           index_t d = ket_cd.q;
           double Delta_cdij = OneBody(c,c) + OneBody(d,d) - OneBody(i,i) - OneBody(j,j);
           Emp3 += (2*J+1)*Mat(iket_ab,iket_ij) * Mat(iket_ij,iket_cd) * Mat(iket_cd,iket_ab) / (Delta_abij * Delta_cdij);
         }

       // pp term
         for (auto iket_kl : tbc.GetKetIndex_pp() )
         {
           Ket& ket_kl = tbc.GetKet(iket_kl);
           index_t k = ket_kl.p;
           index_t l = ket_kl.q;
           double Delta_abkl = OneBody(a,a) + OneBody(b,b) - OneBody(k,k) - OneBody(l,l);
           Emp3 += (2*J+1)*Mat(iket_ab,iket_ij)*Mat(iket_ij,iket_kl)*Mat(iket_kl,iket_ab) / (Delta_abij * Delta_abkl);
         }

       } // for ij
     } // for ab
   } // for ich
   cout << "done with pp and hh. E(3) = " << Emp3 << endl;




   index_t nparticles = modelspace->particles.size();
   modelspace->PreCalculateSixJ();
   #pragma omp parallel for schedule(dynamic,1)  reduction(+:Emp3)
   for (index_t ii=0;ii<nparticles;ii++)
   {
     auto i = modelspace->particles[ii];
     double ji = 0.5*modelspace->GetOrbit(i).j2;
     for (auto a : modelspace->holes)
     {
      double ja = 0.5*modelspace->GetOrbit(a).j2;
      int J_min = std::abs(ja-ji);
      int J_max = ja+ji;
      for (int J_tot=J_min;J_tot<=J_max;++J_tot)
      {
       double Jfactor = (2*J_tot + 1) ; // I don't yet understand why it's (2J+1)**2, but this is what came from Johannes.
//       double Jfactor = (2*J_tot + 1)*(2*J_tot + 1) ; // I don't yet understand why it's (2J+1)**2, but this is what came from Johannes.
       for (auto b : modelspace->holes)
       {
        double jb = 0.5*modelspace->GetOrbit(b).j2;
        for(auto j : modelspace->particles)
        {
         double jj = 0.5*modelspace->GetOrbit(j).j2;
         double Delta_abij = OneBody(a,a) + OneBody(b,b) - OneBody(i,i) - OneBody(j,j);
         int J1min = std::max(std::abs(ja-jb),std::abs(ji-jj));
         int J1max = std::min(ja+jb,ji+jj);
         double tbme_abij = 0;
         if ( AngMom::Triangle(jj,jb,J_tot) )
         {
          for (int J1=J1min;J1<=J1max;++J1)  //Pandya 1: <ai`| V |jb`>_Jtot
          {
            tbme_abij -= modelspace->GetSixJ(ja,ji,J_tot,jj,jb,J1)  * (2*J1 + 1) *  TwoBody.GetTBME_J(J1,a,b,j,i);
          }
         }
         for (auto c : modelspace->holes )
          {
           double jc = 0.5*modelspace->GetOrbit(c).j2;
           for (auto k : modelspace->particles )
            {
             double jk = 0.5*modelspace->GetOrbit(k).j2;
             if ( not AngMom::Triangle(jc,jk,J_tot) ) continue;
             double Delta_acik = OneBody(a,a) + OneBody(c,c) - OneBody(i,i) - OneBody(k,k);
             int J2min = std::max(std::abs(jc-jj),std::abs(jk-jb));
             int J2max = std::min(jc+jj,jk+jb);
             double tbme_cjkb = 0;
             if ( AngMom::Triangle(jj,jb,J_tot) )
             {
               for (int J2=J2min;J2<=J2max;++J2) // Pandya 2:  <jb` | V | kc`>_Jtot
               {
                 tbme_cjkb -= modelspace->GetSixJ(jj,jb,J_tot,jk,jc,J2) * (2*J2 + 1) *  TwoBody.GetTBME_J(J2,j,c,k,b);
               }
             }
             int J3min = std::max(std::abs(ji-jk),std::abs(ja-jc));
             int J3max = std::min(ji+jk,ja+jc);
             double tbme_ikac = 0;
             for (int J3=J3min;J3<=J3max;++J3) // Pandya 3:   <kc`| V | ai`>_Jtot
             {
               tbme_ikac -= modelspace->GetSixJ(jk,jc,J_tot,ja,ji,J3) * (2*J3 + 1) *  TwoBody.GetTBME_J(J3,k,i,a,c);
             }
             Emp3 +=  Jfactor * tbme_abij * tbme_cjkb * tbme_ikac / (Delta_abij * Delta_acik);
            } // for k
          } // for c
        } // for j
       } // for b
      } // for J_tot
     } // for a
   } // for i


   profiler.timer["GetMP3_Energy"] += omp_get_wtime() - t_start;
   return Emp3;
}


/*
double Operator::GetMP3_Energy()
{
   double t_start = omp_get_wtime();
   double Emp3 = 0;
   // This can certainly be optimized, but I'll wait until this is the bottleneck.
   index_t nholes = modelspace->holes.size();
   index_t norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for reduction(+:Emp3)
   for (index_t a=0;a<nholes;++a)
   {
     Orbit& oa = modelspace->GetOrbit(a);
     for (index_t b : modelspace->holes)
     {
       if (a>b) continue;
       Orbit& ob = modelspace->GetOrbit(b);
       for (index_t i : modelspace->particles)
       {
         Orbit& oi = modelspace->GetOrbit(i);
         for (index_t j : modelspace->particles)
         {
           if (i>j) continue;
           Orbit& oj = modelspace->GetOrbit(j);
           {
             double Delta_abij  = OneBody(a,a)+OneBody(b,b)-OneBody(i,i)-OneBody(j,j);
             for (index_t p=0;p<norbits;++p)
             {
               Orbit& op = modelspace->GetOrbit(p);
               for (index_t q=p;q<norbits;++q)
               {
                 Orbit& oq = modelspace->GetOrbit(q);
                 if (p==a and q==b) continue; 
                 if (p==i and q==j) continue; 
                 double Delta_abpq = OneBody(a,a)+OneBody(b,b)-OneBody(p,p)-OneBody(q,q);
                 double Delta_apiq = OneBody(a,a)+OneBody(p,p)-OneBody(i,i)-OneBody(q,q);
                 int Jmin =  max(std::abs(op.j2-oq.j2),max(std::abs(oi.j2-oj.j2),std::abs(oa.j2-ob.j2)))/2;
                 int Jmax =  min(op.j2+oq.j2,min(oi.j2+oj.j2,oa.j2+ob.j2))/2;
                 for (int J=Jmin;J<=Jmax;++J)
                 {
                   double tbme1 = TwoBody.GetTBME_J_norm(J,a,b,i,j);
                   if (q<nholes or p>=nholes)
                   {
                     double tbme2 = TwoBody.GetTBME_J_norm(J,i,j,p,q);
                     double tbme3 = TwoBody.GetTBME_J_norm(J,p,q,a,b);
                     Emp3 += 2*(2*J+1)*tbme1*tbme2*tbme3/(Delta_abij*Delta_abpq);
                   }
                   else
                   {
                     double tbme4 = TwoBody.GetTBME_J_norm(J,p,j,q,b);
                     double tbme5 = TwoBody.GetTBME_J_norm(J,i,q,a,p);
                     Emp3 += 16*(2*J+1)*tbme1*tbme4*tbme5/(Delta_abij*Delta_apiq);
                   }
                 }
               }
             }
           }
         }
       }
     }
   }

   profiler.timer["GetMP3_Energy"] += omp_get_wtime() - t_start;
   return Emp3;
}
*/

//*************************************************************
/// Evaluate first order perturbative correction to the operator's
/// ground-state expectation value. A HF basis is assumed.
/// \f[
///  \mathcal{O}^{(1)} = 2\sum_{abij} \frac{H_{abij}\mathcal{O}_{ijab}}{\Delta_{abij}}
/// \f]
///
//*************************************************************
double Operator::MP1_Eval(Operator& H)
{
  auto& Op = *this;
  double opval = 0;
  int nch = modelspace->GetNumberTwoBodyChannels();
  for (int ich=0;ich<nch;++ich)
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ich);
    int J = tbc.J;
    auto& Hmat = H.TwoBody.GetMatrix(ich,ich);
    auto& Opmat = Op.TwoBody.GetMatrix(ich,ich);
    for (auto ibra : tbc.GetKetIndex_hh() )
    {
      Ket& bra = tbc.GetKet(ibra);
      for (auto iket : tbc.GetKetIndex_pp() )
      {
        Ket& ket = tbc.GetKet(iket);
        double Delta_abij = H.OneBody(bra.p,bra.p) + H.OneBody(bra.q,bra.q) -H.OneBody(ket.p,ket.p) - H.OneBody(ket.q,ket.q);
        opval += 2*(2*J+1)*Hmat(ibra,iket) * Opmat(iket,ibra) / Delta_abij;
      }
    }
  }
  return opval;
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
/*
[[deprecated]] void Operator::CalculateKineticEnergy()
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
*/



//*****************************************************************************************
/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Operator::BCH_Transform( const Operator &Omega)
{
   return use_brueckner_bch ? Brueckner_BCH_Transform( Omega ) :  Standard_BCH_Transform( Omega );
}

/// X.BCH_Transform(Y) returns \f$ Z = e^{Y} X e^{-Y} \f$.
/// We use the [Baker-Campbell-Hausdorff formula](http://en.wikipedia.org/wiki/Baker-Campbell-Hausdorff_formula)
/// \f[ Z = X + [Y,X] + \frac{1}{2!}[Y,[Y,X]] + \frac{1}{3!}[Y,[Y,[Y,X]]] + \ldots \f]
/// with all commutators truncated at the two-body level.
Operator Operator::Standard_BCH_Transform( const Operator &Omega)
{
   double t_start = omp_get_wtime();
   int max_iter = 40;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   double factorial_denom = 1.0;
   Operator goosetank_chi;  // auxiliary one-body operator used to recover 4th-order quadruples.
   if (use_goose_tank_correction)
   {
     goosetank_chi = *this;
     goosetank_chi.SetParticleRank(1);
     goosetank_chi.Erase();
   }
   if (nx>bch_transform_threshold)
   {
     Operator OpNested = *this;
     double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
     for (int i=1; i<=max_iter; ++i)
     {

        if (use_goose_tank_correction  )
        {
          auto chi_last = goosetank_chi.OneBody;
          goosetank_chi.GooseTankUpdate( Omega, OpNested);
//          OpNested.OneBody += goosetank_chi.OneBody;  // add the chi from the previous step to OpNested.
          OpNested.OneBody += chi_last;  // add the chi from the previous step to OpNested.
        }
        

        OpNested = Commutator(Omega,OpNested); // the ith nested commutator
//        Operator tmp1 = Commutator(Omega,OpNested); // the ith nested commutator
        factorial_denom /= i;
//        tmp1 /= i;
//        OpNested = tmp1;
        OpOut += factorial_denom * OpNested;
  
        if (this->rank_J > 0)
        {
            cout << "Tensor BCH, i=" << i << "  Norm = " << setw(12) << setprecision(8) << fixed << OpNested.OneBodyNorm() << " " 
                                                         << setw(12) << setprecision(8) << fixed << OpNested.TwoBodyNorm() << " "
                                                         << setw(12) << setprecision(8) << fixed << OpNested.Norm() << endl;
        }
        epsilon *= i+1;
        if (OpNested.Norm() < epsilon)  break;
//        if (OpNested.Norm() < epsilon *(i+1))  break;
        if (i == warn_iter)  cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << endl;
        else if (i == max_iter)   cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << endl;
     }
   }
   profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
   return OpOut;
}


//  Update the auxiliary one-body operator chi, using Omega and the ith nested commutator
//  This has not been tested for tensor commutators, but it *should* work.
//  Of course, there's not as clean a motivation in terms of perturbation theory for the tensors...
// 
void Operator::GooseTankUpdate( const Operator& Omega, const Operator& OpNested)
{
   double t_start = omp_get_wtime();
   auto& goosetank_chi = *this;
   goosetank_chi.EraseOneBody();
   if (this->rank_J==0 )
   {
     goosetank_chi.comm221ss( Omega, OpNested );  // update chi.
   }
   else
   {
     goosetank_chi.comm222_pp_hh_221st( Omega, OpNested );  // update chi.
   }
   goosetank_chi.Symmetrize(); // the commutator call only does half the matrix, so we symmetrize
   int norbits = modelspace->GetNumberOrbits();
   for (int i=0;i<norbits;++i)  // enforce n_in_j + nbar_i nbar_j
   {
     Orbit &oi = modelspace->GetOrbit(i);
     for (int j=0;j<norbits;++j)
     {
      Orbit &oj = modelspace->GetOrbit(j);
      goosetank_chi.OneBody(i,j) *=  oi.occ*oj.occ + (1.0-oi.occ)*(1.0-oj.occ) ;
      }
   }
   profiler.timer["GooseTankUpdate"] += omp_get_wtime() - t_start;
}


/*
Operator Operator::Standard_BCH_Transform( const Operator &Omega)
{
   double t_start = omp_get_wtime();
   int max_iter = 40;
   int warn_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   if (nx>bch_transform_threshold)
   {
     Operator OpNested = *this;
     double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
     for (int i=1; i<=max_iter; ++i)
     {
        Operator tmp1 = Commutator(Omega,OpNested);
         tmp1 /= i;
        OpNested = tmp1;
        OpOut += OpNested;
  
        if (this->rank_J > 0)
        {
            cout << "Tensor BCH, i=" << i << "  Norm = " << setw(12) << setprecision(8) << fixed << OpNested.OneBodyNorm() << " " 
                                                         << setw(12) << setprecision(8) << fixed << OpNested.TwoBodyNorm() << " "
                                                         << setw(12) << setprecision(8) << fixed << OpNested.Norm() << endl;
        }
//        if (OpNested.Norm() < bch_transform_threshold )  break;
        if (OpNested.Norm() < epsilon *(i+1))  break;
        if (i == warn_iter)  cout << "Warning: BCH_Transform not converged after " << warn_iter << " nested commutators" << endl;
        else if (i == max_iter)   cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << endl;
     }
   }
   profiler.timer["BCH_Transform"] += omp_get_wtime() - t_start;
   return OpOut;
}


*/


/// Variation of the BCH transformation procedure
/// requested by a one Dr. T.D. Morris
/// \f[ e^{\Omega_1 + \Omega_2} X e^{-\Omega_1 - \Omega_2}
///    \rightarrow 
///  e^{\Omega_2} e^{\Omega_1}  X e^{-\Omega_1} e^{-\Omega_2} \f]
Operator Operator::Brueckner_BCH_Transform( const Operator &Omega)
{
   Operator Omega1 = Omega;
   Operator Omega2 = Omega;
   Omega1.SetParticleRank(1);
   Omega1.EraseTwoBody();
   Omega2.EraseOneBody();
   Operator OpOut = this->Standard_BCH_Transform(Omega1);
   OpOut = OpOut.Standard_BCH_Transform(Omega2);
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
   double tstart = omp_get_wtime();
   Operator& X = *this;
   double nx = X.Norm();
   vector<double> bernoulli = {1.0, -0.5, 1./6, 0.0, -1./30,  0.0 ,  1./42,     0,  -1./30};
   vector<double> factorial = {1.0,  1.0,  2.0, 6.0,    24.,  120.,   720., 5040.,  40320.};


//   Operator goosetank_chi;  // auxiliary one-body operator used to recover 4th-order quadruples.





   Operator Z = X + Y;
//   if (use_goose_tank_correction) return Z; // Not sure why this is here
   Operator Nested = Y;
   Nested.SetToCommutator(Y,X);

//   if (use_goose_tank_correction)
//   {
//     goosetank_chi = *this;
//     goosetank_chi.SetParticleRank(1);
//     goosetank_chi.Erase();
//     goosetank_chi.GooseTankUpdate( Y, Nested);
//   }


   double nxy = Nested.Norm();
   // We assume X is small, but just in case, we check if we should include the [X,[X,Y]] term.
   if ( nxy*nx > bch_product_threshold)
   {
     Z += (1./12) * Commutator(Nested,X);
//     cout << "Operator::BCH_Product -- Included X^2 term. " << nx << " " << ny << " " << nxy << endl;
   }
   
   int k = 1;
   while( Nested.Norm() > bch_product_threshold and k<9)
   {
     if (k<2 or k%2==0)
        Z += (bernoulli[k]/factorial[k]) * Nested;

//     if (use_goose_tank_correction  and k<2 ) // for the sake of speed, only apply goose tank for first nested commutator here.
//     {
//          auto chi_last = goosetank_chi.OneBody;
//          goosetank_chi.GooseTankUpdate( Y, Nested);
//          Nested.OneBody += chi_last;  // add the chi from the previous step to OpNested.
//     }

     Nested = Commutator(Y,Nested);
     k++;
   }

   profiler.timer["BCH_Product"] += omp_get_wtime() - tstart;
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
   double nrm = 0;
   for ( int p=0; p<modelspace->GetNumberOrbits(); ++p)
   {
     Orbit& op = modelspace->GetOrbit(p);
     for ( auto q : OneBodyChannels.at({op.l, op.j2, op.tz2}) )
     {
       Orbit& oq = modelspace->GetOrbit(q);
       int degeneracy_factor = (op.j2+1) * ( (min(oq.j2,op.j2+rank_J)-max(-oq.j2,op.j2-rank_J))/2+1 );
       nrm += OneBody(p,q)*OneBody(p,q) * degeneracy_factor * degeneracy_factor;
     }
   }
   return sqrt(nrm);
//   return arma::norm(OneBody,"fro");
}



double Operator::TwoBodyNorm() const
{
  return TwoBody.Norm();
}


void Operator::MakeNormalized(){ ChangeNormalization( 1./SQRT2)  ;}
void Operator::MakeUnNormalized(){ ChangeNormalization( SQRT2)  ;}
void Operator::ChangeNormalization(double factor)
{
  for (auto& itmat : TwoBody.MatEl)
  {
    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];
    auto& TBME = itmat.second;
    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
    for (int ibra=0;ibra<nbras;ibra++)
    {
      Ket& bra = tbc_bra.GetKet(ibra);
      if (bra.p == bra.q)
      {
        TBME.col(ibra) *= factor;
      }
    }
    for (int iket=0;iket<nkets;iket++)
    {
      Ket& ket = tbc_ket.GetKet(iket);
      if (ket.p == ket.q)
      {
        TBME.row(iket) *= factor;
      }
    }
  }
}



double Operator::Trace(int Atrace, int Ztrace) const
{
  double t_start = omp_get_wtime();
  int Ntrace = Atrace - Ztrace;
  Operator OpVac = UndoNormalOrdering();
  double trace = OpVac.ZeroBody;
  int emax = modelspace->GetEmax();
  double M = (emax+1)*(emax+2)*(emax+3)/3; // number of m-scheme orbits for the proton or neutron partitions
  double norm_p = Ztrace / M;
  double norm_n = Ntrace / M;
  double norm_pp = Ztrace*(Ztrace-1) / (M *(M-1) ) ;
  double norm_nn = Ntrace*(Ntrace-1) / (M *(M-1) ) ;
  double norm_pn = norm_p * norm_n ;

//  for (int i=0; i<modelspace->GetNumberOrbits(); ++i)
  for (auto i : modelspace->proton_orbits)
  {
     Orbit& oi = modelspace->GetOrbit(i);
     trace += OpVac.OneBody(i,i) * ( oi.j2 +1) * norm_p;
  }
  for (auto i : modelspace->neutron_orbits)
  {
     Orbit& oi = modelspace->GetOrbit(i);
      trace += OpVac.OneBody(i,i) * ( oi.j2 +1) * norm_n;
  }

  for ( int ch=0; ch<modelspace->GetNumberTwoBodyChannels(); ++ch )
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    double trVijij = arma::trace( OpVac.TwoBody.GetMatrix(ch,ch) );
    switch ( tbc.Tz )
    {
      case ( -1 ) :  trace += (tbc.J*2+1)* trVijij * norm_pp; break;
      case (  0 ) :  trace += (tbc.J*2+1)* trVijij * norm_pn; break;
      case (  1 ) :  trace += (tbc.J*2+1)* trVijij * norm_nn; break;
      default: cout << "AAAHHH blew the switch statement. tbc.Tz = " << tbc.Tz << endl;
    }
  }
  profiler.timer["Operator::Trace"] += omp_get_wtime() - t_start;
  return trace;
}


void Operator::ScaleFermiDirac(Operator& H, double T, double Efermi)
{
  int norb = modelspace->GetNumberOrbits();
  
  for (int i=0; i<norb; ++i)
  {
    Orbit& oi = modelspace->GetOrbit(i);
    double ei = OneBody(i,i);
    for (auto j : OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
    {
      double ej = OneBody(j,j);
      OneBody(i,j) *= 1./(1 + exp( (ei+ej-2*Efermi)/(2*T))  );
    }
  }
  for (auto itmat : TwoBody.MatEl )
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel( itmat.first[0] );
    for (int ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      double ei = OneBody(i,i);
      double ej = OneBody(j,j);
      for (int iket=0; iket<tbc.GetNumberKets(); ++iket)
      {
        Ket& ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        double ek = OneBody(k,k);
        double el = OneBody(l,l);
        itmat.second.row(ibra) *= 1./(1 + exp( (ei+ej+ek+el-4*Efermi)/(2*T))  );
      }
    }
  }
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
  int jrank = max(X.rank_J,Y.rank_J);
  int trank = max(X.rank_T,Y.rank_T);
  int parity = (X.parity+Y.parity)%2;
  int particlerank = max(X.particle_rank,Y.particle_rank);
  Operator Z(*(X.modelspace),jrank,trank,parity,particlerank);
  Z.SetToCommutator(X,Y);
  return Z;
}

void Operator::SetToCommutator( const Operator& X, const Operator& Y)
{
//   profiler.counter["N_Commutators"] += 1;
   double t_start = omp_get_wtime();
   Operator& Z = *this;
   modelspace->PreCalculateSixJ();
   int xrank = X.rank_J + X.rank_T + X.parity;
   int yrank = Y.rank_J + Y.rank_T + Y.parity;
   if (xrank==0)
   {
      if (yrank==0)
      {
         Z.CommutatorScalarScalar(X,Y); // [S,S]
      }
      else
      {
         Z.CommutatorScalarTensor(X,Y); // [S,T]
      }
   }
   else if(yrank==0)
   {
      Z.CommutatorScalarTensor(Y,X); // [T,S]
      Z *= -1;
   }
   else
   {
      cout << "In Tensor-Tensor because X.rank_J = " << X.rank_J << "  X.rank_T = " << X.rank_T << "  X.parity = " << X.parity << "   ";
      cout <<                        "  Y.rank_J = " << Y.rank_J << "  Y.rank_T = " << Y.rank_T << "  Y.parity = " << Y.parity << endl;
      cout << " Tensor-Tensor commutator not yet implemented." << endl;
   }
   profiler.timer["Commutator"] += omp_get_wtime() - t_start;
}


/// Commutator where \f$ X \f$ and \f$Y\f$ are scalar operators.
/// Should be called through Commutator()
void Operator::CommutatorScalarScalar( const Operator& X, const Operator& Y) 
{
   profiler.counter["N_ScalarCommutators"] += 1;
   double t_css = omp_get_wtime();
   Operator& Z = *this;
   Z = X.GetParticleRank()>Y.GetParticleRank() ? X : Y;
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   if ( not Z.IsAntiHermitian() )
   {
      Z.comm110ss(X, Y);
      if (X.particle_rank>1 and Y.particle_rank>1)
        Z.comm220ss(X, Y) ;
   }

   double t_start = omp_get_wtime();
   Z.comm111ss(X, Y);
   profiler.timer["comm111ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
   Z.comm121ss(X,Y);
   profiler.timer["comm121ss"] += omp_get_wtime() - t_start;

    t_start = omp_get_wtime();
   Z.comm122ss(X,Y); 
   profiler.timer["comm122ss"] += omp_get_wtime() - t_start;

   if (X.particle_rank>1 and Y.particle_rank>1)
   {
     t_start = omp_get_wtime();
     Z.comm222_pp_hh_221ss(X, Y);
     profiler.timer["comm222_pp_hh_221ss"] += omp_get_wtime() - t_start;
      
     t_start = omp_get_wtime();
     Z.comm222_phss(X, Y);
     profiler.timer["comm222_phss"] += omp_get_wtime() - t_start;
   }


   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();

   profiler.timer["CommutatorScalarScalar"] += omp_get_wtime() - t_css;

}


/// Commutator \f$[X,Y]\f$ where \f$ X \f$ is a scalar operator and \f$Y\f$ is a tensor operator.
/// Should be called through Commutator()
void Operator::CommutatorScalarTensor( const Operator& X, const Operator& Y) 
{
   profiler.counter["N_TensorCommutators"] += 1;
   double t_cst = omp_get_wtime();
   Operator& Z = *this;
   Z = Y; // This ensures the commutator has the same tensor rank as Y
   Z.EraseZeroBody();
   Z.EraseOneBody();
   Z.EraseTwoBody();

   if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Z.SetAntiHermitian();
   else if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Z.SetHermitian();
   else Z.SetNonHermitian();

   double t_start = omp_get_wtime();
//   cout << "comm111st" << endl;
   Z.comm111st(X, Y);
   profiler.timer["comm111st"] += omp_get_wtime() - t_start;
//   cout << "comm121st" << endl;
   t_start = omp_get_wtime();
   Z.comm121st(X, Y);
   profiler.timer["comm121st"] += omp_get_wtime() - t_start;

//   cout << "comm122st" << endl;
   t_start = omp_get_wtime();
   Z.comm122st(X, Y);
   profiler.timer["comm122st"] += omp_get_wtime() - t_start;
//   cout << "comm222_pp_hh_st" << endl;
   t_start = omp_get_wtime();
   Z.comm222_pp_hh_221st(X, Y);
   profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - t_start;
//   cout << "comm222_phst" << endl;
   t_start = omp_get_wtime();
   Z.comm222_phst(X, Y);
   profiler.timer["comm222_phst"] += omp_get_wtime() - t_start;

//   cout << "symmetrize" << endl;
   if ( Z.IsHermitian() )
      Z.Symmetrize();
   else if (Z.IsAntiHermitian() )
      Z.AntiSymmetrize();
//   cout << "done." << endl;
   profiler.timer["CommutatorScalarTensor"] += omp_get_wtime() - t_cst;
//   return Z;
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
void Operator::comm110ss( const Operator& X, const Operator& Y) 
{
  Operator& Z = *this;
  if (X.IsHermitian() and Y.IsHermitian()) return ; // I think this is the case
  if (X.IsAntiHermitian() and Y.IsAntiHermitian()) return ; // I think this is the case

   arma::mat xyyx = X.OneBody*Y.OneBody - Y.OneBody*X.OneBody;
   for ( auto& a : modelspace->holes) 
   {
      Orbit& oa = modelspace->GetOrbit(a);
      Z.ZeroBody += (oa.j2+1) * oa.occ * xyyx(a,a);
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
      auto ph = tbc.GetKetIndex_ph();
      auto pp = tbc.GetKetIndex_pp();
      arma::uvec nbar_indices = arma::join_cols(hh,ph);
      nbar_indices = arma::join_cols(nbar_indices,pp);
      if (hh.size()==0 ) continue;
      auto nn = tbc.Ket_occ_hh;
      arma::vec nbarnbar = arma::join_cols(tbc.Ket_unocc_hh, tbc.Ket_unocc_ph);
      auto & X2 = X.TwoBody.GetMatrix(ch).submat(hh,nbar_indices);
      arma::mat Y2 = Y.TwoBody.GetMatrix(ch).submat(nbar_indices,hh);
      Y2.head_rows(nbarnbar.size()).each_col() %= nbarnbar;
      Z.ZeroBody += 2 * (2*tbc.J+1) * arma::sum(arma::diagvec( X2 * Y2 ) % nn); // This could be made more efficient, but who cares?
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
void Operator::comm121ss( const Operator& X, const Operator& Y) 
{
   Operator& Z = *this;
   index_t norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for 
   for (index_t i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      index_t jmin = Z.IsNonHermitian() ? 0 : i;
      for (auto j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          if (j<jmin) continue; // only calculate upper triangle
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             for (index_t b=0; b<norbits; ++b)
             {
                Orbit &ob = modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                if (std::abs(nanb)<OCC_CUT) continue;
                if (Y.particle_rank>1)
                {
                  Z.OneBody(i,j) += (ob.j2+1) * nanb *  X.OneBody(a,b) * Y.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) -= (oa.j2+1) * nanb *  X.OneBody(b,a) * Y.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
                // comm211 part
                if (X.particle_rank>1)
                {
                  Z.OneBody(i,j) -= (ob.j2+1) * nanb * Y.OneBody(a,b) * X.TwoBody.GetTBMEmonopole(b,i,a,j) ;
                  Z.OneBody(i,j) += (oa.j2+1) * nanb * Y.OneBody(b,a) * X.TwoBody.GetTBMEmonopole(a,i,b,j) ;
                }
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
void Operator::comm221ss( const Operator& X, const Operator& Y) 
{

   double t_start = omp_get_wtime();
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Y.TwoBody;
   static TwoBodyME Mhh = Y.TwoBody;

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

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


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
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


   } //for ch

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
               Orbit& oc = modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ     * Mpp.GetTBME(ch,c,i,c,j);
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (auto& c : modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   }

   profiler.timer["comm221ss"] += omp_get_wtime() - t_start;

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
   int hZ = Z.IsHermitian() ? 1 : -1;

   int n_nonzero = modelspace->SortedTwoBodyChannels.size();
   #pragma omp parallel for schedule(dynamic,1)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      auto& X2 = X.TwoBody.GetMatrix(ch,ch);
      auto& Y2 = Y.TwoBody.GetMatrix(ch,ch);
      auto& Z2 = Z.TwoBody.GetMatrix(ch,ch);
      arma::mat W2(size(Z2),arma::fill::zeros); // temporary intermediate matrix

      int npq = tbc.GetNumberKets();
      for (int indx_ij = 0;indx_ij<npq; ++indx_ij)
      {
         Ket & bra = tbc.GetKet(indx_ij);
         int i = bra.p;
         int j = bra.q;
         Orbit& oi = modelspace->GetOrbit(i);
         Orbit& oj = modelspace->GetOrbit(j);
         int flipphaseij = -modelspace->phase((oi.j2+oj.j2)/2-tbc.J);

         // make lists of the indices we want, then do matrix multiplication.
         // there may be a more efficient way to find these
         vector<index_t> ind1_ia, ind1_ja,ind2_aj,ind2_ai;
         vector<double> factor_ia,factor_ja;
         for (int a : OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
         {
            int ind2 = tbc.GetLocalIndex( min(a,j), max(a,j));
            if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
            ind1_ia.push_back(a);
            ind2_aj.push_back(ind2);
            factor_ia.push_back( a>j ? flipphaseij : (a==j ? SQRT2 : 1));
         }
         if (i!=j)
         {
           for (int a : OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
           {
              int ind2 = tbc.GetLocalIndex( min(a,i), max(a,i));
              if (ind2<0 or ind2>=tbc.GetNumberKets()) continue;
              ind1_ja.push_back(a);
              ind2_ai.push_back(ind2);
              factor_ja.push_back( i>a ? flipphaseij : (i==a ? SQRT2 : 1));
           }

         }

         arma::uvec u_ind1_ia(ind1_ia);
         arma::uvec u_ind1_ja(ind1_ja);
         arma::uvec u_ind2_aj(ind2_aj);
         arma::uvec u_ind2_ai(ind2_ai);
         arma::vec  v_factor_ia(factor_ia);
         arma::vec  v_factor_ja(factor_ja);
         if (i==j)
         {
           v_factor_ia /= SQRT2;
           v_factor_ja /= SQRT2;
         }

         // This is fairly obfuscated, but hopefully faster for bigger calculations
         if (X.particle_rank>1 and Y.particle_rank>1)
         {
            W2.col(indx_ij) = join_horiz(    Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))  , X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) )
                             * join_vert(  join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ),
                                          -join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja ));
         }
         else if (X.particle_rank<2 and Y.particle_rank>1)
         {
            W2.col(indx_ij) =     Y2.cols(join_vert( u_ind2_aj,u_ind2_ai))    
                             *   join_vert( X1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, X1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }
         else if (X.particle_rank>1 and Y.particle_rank<2)
         {
            W2.col(indx_ij) =      -X2.cols(join_vert(u_ind2_aj,u_ind2_ai) ) 
                                      *    join_vert( Y1.unsafe_col(i).rows(u_ind1_ia)%v_factor_ia, Y1.unsafe_col(j).rows(u_ind1_ja)%v_factor_ja );
         }

      if (i==j) W2.col(indx_ij) *= 2;

      }
      Z2 -= W2 + hZ*W2.t();
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

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

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
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
//      auto& nabar_nbbar = tbc.Ket_unocc_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


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
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
   profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;
}








void Operator::ConstructScalarMpp_Mhh(const Operator& X, const Operator& Y, TwoBodyME& Mpp, TwoBodyME& Mhh) const
{
   int nch = modelspace->SortedTwoBodyChannels.size();
//   Operator& Z = *this;
   bool z_is_hermitian = IsHermitian();
   bool z_is_antihermitian = IsAntiHermitian();
//   bool z_is_hermitian = X.IsHermitian() xor Y.IsHermitian();
//   bool z_is_antihermitian = (X.IsHermitian() == Y.IsHermitian()) and (X.IsAntiHermitian() == Y.IsAntiHermitian());
   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1)
   #endif
   for (int ich=0; ich<nch; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

      auto& LHS = X.TwoBody.GetMatrix(ch,ch);
      auto& RHS = Y.TwoBody.GetMatrix(ch,ch);
//      auto& OUT = Z.TwoBody.GetMatrix(ch,ch);

      auto& Matrixpp = Mpp.GetMatrix(ch,ch);
      auto& Matrixhh = Mhh.GetMatrix(ch,ch);

      auto& kets_pp = tbc.GetKetIndex_pp();
      auto& kets_hh = tbc.GetKetIndex_hh();
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


      if (z_is_hermitian)
      {
         Matrixpp +=  Matrixpp.t();
         Matrixhh +=  Matrixhh.t();
      }
      else if (z_is_antihermitian) // i.e. LHS and RHS are both hermitian or ant-hermitian
      {
         Matrixpp -=  Matrixpp.t();
         Matrixhh -=  Matrixhh.t();
      }
      else
      {
        Matrixpp -=  RHS.cols(kets_pp) * LHS.rows(kets_pp);
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


      // The two body part
//      OUT += Matrixpp - Matrixhh;
   } //for ch

}


/// Since comm222_pp_hhss() and comm221ss() both require the construction of 
/// the intermediate matrices \f$\mathcal{M}_{pp} \f$ and \f$ \mathcal{M}_{hh} \f$, we can combine them and
/// only calculate the intermediates once.
void Operator::comm222_pp_hh_221ss( const Operator& X, const Operator& Y )  
{

//   int herm = Z.IsHermitian() ? 1 : -1;
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();

   static TwoBodyME Mpp = Z.TwoBody;
   static TwoBodyME Mhh = Z.TwoBody;

   double t = omp_get_wtime();
   ConstructScalarMpp_Mhh( X, Y, Mpp, Mhh);
/*
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
      auto& kets_ph = tbc.GetKetIndex_ph();
      auto& nanb = tbc.Ket_occ_hh;
      auto& nbarnbar_hh = tbc.Ket_unocc_hh;
      auto& nbarnbar_ph = tbc.Ket_unocc_ph;
      
      Matrixpp =  LHS.cols(kets_pp) * RHS.rows(kets_pp);
      Matrixhh =  LHS.cols(kets_hh) * arma::diagmat(nanb) *  RHS.rows(kets_hh) ;
      if (kets_hh.size()>0)
        Matrixpp +=  LHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  RHS.rows(kets_hh); 
      if (kets_ph.size()>0)
        Matrixpp += LHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  RHS.rows(kets_ph) ;


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
        Matrixhh -=  RHS.cols(kets_hh) * arma::diagmat(nanb) *  LHS.rows(kets_hh) ;
        Matrixpp -=  RHS.cols(kets_hh) * arma::diagmat(nbarnbar_hh) *  LHS.rows(kets_hh) ;
        if (kets_ph.size()>0)
          Matrixpp -=  RHS.cols(kets_ph) * arma::diagmat(nbarnbar_ph) *  LHS.rows(kets_ph) ;
      }


      // The two body part
      OUT += Matrixpp - Matrixhh;
   } //for ch
*/
//   Z.TwoBody += (Mpp - Mhh);
   Z.TwoBody += Mpp;
   Z.TwoBody -= Mhh;
//   OUT += Matrixpp - Matrixhh;
   profiler.timer["pphh TwoBody bit"] += omp_get_wtime() - t;

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
               Orbit& oc = modelspace->GetOrbit(c);
               cijJ += Jfactor * oc.occ * Mpp.GetTBME(ch,c,i,c,j); 
               cijJ += Jfactor * (1-oc.occ) * Mhh.GetTBME(ch,c,i,c,j); 
            }
            // Sum c over particles and include the n_a * n_b terms
            for (auto& c : modelspace->particles)
            {
               cijJ += Jfactor * Mhh.GetTBME(ch,c,i,c,j);
            }
         }
         Z.OneBody(i,j) += cijJ /(oi.j2+1.0);
      } // for j
   } // for i
   profiler.timer["pphh One Body bit"] += omp_get_wtime() - t;
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
void Operator::DoPandyaTransformation_SingleChannel(arma::mat& TwoBody_CC_ph, int ch_cc, string orientation="normal") const
{
   int herm = IsHermitian() ? 1 : -1;
   TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
   int nKets_cc = tbc_cc.GetNumberKets();
   arma::uvec kets_ph = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
   int nph_kets = kets_ph.n_rows;
   int J_cc = tbc_cc.J;

   if (orientation=="normal") TwoBody_CC_ph.zeros( 2*nph_kets, nKets_cc);
   else if (orientation=="transpose") TwoBody_CC_ph.zeros( nKets_cc, 2*nph_kets);
   else
   {
     cout << __PRETTY_FUNCTION__ << " =>  Unknown orientation input  " << orientation << ". Don't know what to do with this." << endl;
     return;
   }

   // loop over cross-coupled ph bras <ab| in this channel
   // (this is the side that gets summed over in the matrix multiplication)
   for (int ibra=0; ibra<nph_kets; ++ibra)
   {
      Ket & bra_cc = tbc_cc.GetKet( kets_ph[ibra] );
      int a = bra_cc.p;
      int b = bra_cc.q;
      Orbit & oa = modelspace->GetOrbit(a);
      Orbit & ob = modelspace->GetOrbit(b);
      double ja = oa.j2*0.5;
      double jb = ob.j2*0.5;
      double na_nb_factor = oa.occ - ob.occ;

      // loop over cross-coupled kets |cd> in this channel
      for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
      {
         Ket & ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
         int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
         int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
         Orbit & oc = modelspace->GetOrbit(c);
         Orbit & od = modelspace->GetOrbit(d);
         double jc = oc.j2*0.5;
         double jd = od.j2*0.5;


         int jmin = max(std::abs(ja-jd),std::abs(jc-jb));
         int jmax = min(ja+jd,jc+jb);
         double Xbar = 0;
         for (int J_std=jmin; J_std<=jmax; ++J_std)
         {
            double sixj = modelspace->GetSixJ(ja,jb,J_cc,jc,jd,J_std);
            if (std::abs(sixj) < 1e-8) continue;
            double tbme = TwoBody.GetTBME_J(J_std,a,d,c,b);
            Xbar -= (2*J_std+1) * sixj * tbme ;
         }
         if (orientation=="normal")
         {
           TwoBody_CC_ph(ibra,iket_cc) = Xbar;
         }
         else // "transpose"
         {
           TwoBody_CC_ph(iket_cc,ibra) = herm * Xbar * na_nb_factor;
         }

         // Exchange (a <-> b) to account for the (n_a - n_b) term
         jmin = max(std::abs(jb-jd),std::abs(jc-ja));
         jmax = min(jb+jd,jc+ja);
         Xbar = 0;
         for (int J_std=jmin; J_std<=jmax; ++J_std)
         {
            double sixj = modelspace->GetSixJ(jb,ja,J_cc,jc,jd,J_std);
            if (std::abs(sixj) < 1e-8) continue;
            double tbme = TwoBody.GetTBME_J(J_std,b,d,c,a);
            Xbar -= (2*J_std+1) * sixj * tbme ;
         }
         if (orientation=="normal")
         {
           TwoBody_CC_ph(ibra+nph_kets,iket_cc) = Xbar;
         }
         else  // "transpose"
         {
           TwoBody_CC_ph(iket_cc,ibra+nph_kets) = herm * Xbar * -na_nb_factor;
         }

      }
   }
}


void Operator::DoPandyaTransformation(deque<arma::mat>& TwoBody_CC_ph, string orientation="normal") const
{
   // loop over cross-coupled channels
   int n_nonzero = modelspace->SortedTwoBodyChannels_CC.size();
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->scalar_transform_first_pass)
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = modelspace->SortedTwoBodyChannels_CC[ich];
      DoPandyaTransformation_SingleChannel( TwoBody_CC_ph[ch_cc] , ch_cc, orientation);
   }
}



void Operator::AddInversePandyaTransformation_SingleChannel(arma::mat& Zbar, int ch_cc)
{
//   double t_start = omp_get_wtime();
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   TwoBodyChannel_CC& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
   int Jcc = tbc_cc.J;
   int Jhat2 = 2*Jcc+1;
   int parity_cc = tbc_cc.parity;
   int Tz_cc = tbc_cc.Tz;
   int nkets_cc = tbc_cc.GetNumberKets();
   int n_nonzeroChannels = modelspace->SortedTwoBodyChannels.size();
//   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->scalar_transform_first_pass)
   for (int ich = 0; ich < n_nonzeroChannels; ++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels[ich];
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int nKets = tbc.GetNumberKets();
      auto& Zmat = TwoBody.GetMatrix(ch,ch);

      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = 0.5*oi.j2;
         double jj = 0.5*oj.j2;
         int ketmin = IsHermitian() ? ibra : ibra+1;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = 0.5*ok.j2;
            double jl = 0.5*ol.j2;

            double commij = 0;
            double commji = 0;

            int jmin = max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int jmax = min(int(ji+jl),int(jk+jj));
            if ( ((oi.l+ol.l)%2==parity_cc)  and  (std::abs(oi.tz2+ol.tz2)==Tz_cc*2) and (Jcc>=jmin) and (Jcc<=jmax) )
            {
               double sixj = modelspace->GetSixJ(ji,jj,J,jk,jl,Jcc);
               int indx_il = tbc_cc.GetLocalIndex(i,l) ;
               int indx_kj = tbc_cc.GetLocalIndex(min(j,k),max(j,k)) +(k>j?nkets_cc:0);
               commij += Jhat2 * sixj * Zbar(indx_il,indx_kj) ;

            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              jmin = max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              jmax = min(int(jj+jl),int(jk+ji));
              if ( (oi.l+ok.l)%2==parity_cc  and  std::abs(oi.tz2+ok.tz2)==Tz_cc*2 and Jcc>=jmin and Jcc<=jmax)
              {
                 double sixj = modelspace->GetSixJ(jj,ji,J,jk,jl,Jcc);
                 int indx_ik = tbc_cc.GetLocalIndex(i,k) ;
                 int indx_lj = tbc_cc.GetLocalIndex(min(l,j),max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 commji += Jhat2 *  sixj *  Zbar(indx_ik, indx_lj) ;


              }
            }
            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            #pragma omp atomic
            Zmat(ibra,iket) += (commij - modelspace->phase(jk+jl-J ) * commji) / norm;
         }
      }
   }

//   profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;
}





void Operator::AddInversePandyaTransformation(const deque<arma::mat>& Zbar)
{
    // Do the inverse Pandya transform
    // Only go parallel if we've previously calculated the SixJ's. Otherwise, it's not thread safe.
   int n_nonzeroChannels = modelspace->SortedTwoBodyChannels.size();
//   int hZ = IsHermitian() ? 1 : -1;
//   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->SixJ_is_empty())

   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->scalar_transform_first_pass)
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
            int Tz_cc = std::abs(oi.tz2+ol.tz2)/2;
            int jmin = max(std::abs(int(ji-jl)),std::abs(int(jk-jj)));
            int jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(ji,jj,J,jk,jl,Jprime);
               if (std::abs(sixj)<1e-8) continue;
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               TwoBodyChannel_CC& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
               int nkets_cc = tbc_cc.GetNumberKets();
               int indx_il = tbc_cc.GetLocalIndex(min(i,l),max(i,l)) +(i>l?nkets_cc:0);
               int indx_kj = tbc_cc.GetLocalIndex(min(j,k),max(j,k)) +(k>j?nkets_cc:0);
               double me1 = Zbar.at(ch_cc)(indx_il,indx_kj);
               commij += (2*Jprime+1) * sixj * me1;
//               if (J==1 and i==0 and j==1 and k==0 and l==9)
//               if ( ch==0 and ibra+iket==1)
//               if ( ch==0 and ibra==2 and iket==4)
//               {
//                  cout << "commij: ch_cc = " << ch_cc << "  adding   " << 2*Jprime+1 << " * " << scientific << sixj << " * " << me1 << "  ->  " << commij << endl;
//                  cout << "ijkl = " << i << " " << j << " " << k << " " << l << "   ibra,iket =   " << ibra << " " << iket << endl;
//               }
            }

            if (k==l)
            {
              commji = commij;
            }
            else if (i==j)
            {
              commji = modelspace->phase(ji+jj+jk+jl) * commij;
            }
            else
            {
              // now loop over the cross coupled TBME's
              parity_cc = (oi.l+ok.l)%2;
              Tz_cc = std::abs(oi.tz2+ok.tz2)/2;
              jmin = max(std::abs(int(jj-jl)),std::abs(int(jk-ji)));
              jmax = min(int(jj+jl),int(jk+ji));
              for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
              {
                 double sixj = modelspace->GetSixJ(jj,ji,J,jk,jl,Jprime);
                 if (std::abs(sixj)<1e-8) continue;
                 int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
                 TwoBodyChannel_CC& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
                 int nkets_cc = tbc_cc.GetNumberKets();
                 int indx_ik = tbc_cc.GetLocalIndex(min(i,k),max(i,k)) +(i>k?nkets_cc:0);
                 int indx_lj = tbc_cc.GetLocalIndex(min(l,j),max(l,j)) +(l>j?nkets_cc:0);
                 // we always have i<=k so we should always flip Z_jlki = (-1)^{i+j+k+l} Z_iklj
                 double me1 = Zbar.at(ch_cc)(indx_ik, indx_lj) ;//* modelspace->phase(ji+jj+jk+jl);
                 commji += (2*Jprime+1) *  sixj * me1;

//                 if (J==1 and i==0 and j==1 and k==0 and l==9)
//               if ( ch==0 and ibra+iket==1)
//               if ( ch==0 and ibra==2 and iket==4)
//                 {
//                    cout <<  "commji: ch_cc = " << ch_cc << "  adding   " << 2*Jprime+1 << " x " << scientific << sixj << " x Z(" << indx_ik << "," << indx_lj << ") = " << me1 << "  ->  " << commji << endl;
//                 }
              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
//            TwoBody.GetMatrix(ch,ch)(ibra,iket) += (commij - modelspace->phase(ji+jj-J)*commji) / norm;
//            if (J==1 and i==0 and j==1 and k==0 and l==9)
//               if ( ch==0 and ibra+iket==1)
//            {
//              cout << "TBME was " << TwoBody.GetMatrix(ch,ch)(ibra,iket) << endl;
//            }
            TwoBody.GetMatrix(ch,ch)(ibra,iket) += (commij - modelspace->phase(jk+jl-J ) * commji) / norm;
//            if (J==1 and i==0 and j==1 and k==0 and l==9)
//               if ( ch==0 and ibra+iket==1)
//            {
//              cout << "added  (" << commij << " - " << modelspace->phase(jk+jl-J ) << " * " << commji << " ) / " << norm << " = " <<  (commij - modelspace->phase(jk+jl-J ) * commji) / norm << endl;
//              cout << "set TBME to " << TwoBody.GetMatrix(ch,ch)(ibra,iket) << endl;
//            }
         }
      }
   }
}

///*************************************
/// convenience function
/// called by comm222_phss
///*************************************
deque<arma::mat> Operator::InitializePandya(size_t nch, string orientation="normal")
{
   deque<arma::mat> X(nch);
   int n_nonzero = modelspace->SortedTwoBodyChannels_CC.size();
   for (int ich=0; ich<n_nonzero; ++ich)
   {
      int ch_cc = modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();
      if (orientation=="normal")
         X[ch_cc] = arma::mat(2*nph_kets,   nKets_cc, arma::fill::zeros);
      else if (orientation=="transpose")
         X[ch_cc] = arma::mat(nKets_cc, 2*nph_kets,   arma::fill::zeros);
   }
   return X;
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
void Operator::comm222_phss( const Operator& X, const Operator& Y ) 
{

//   cout << "start comm222_phss ******************************************************" << endl;
   Operator& Z = *this;
   int hy = Y.IsHermitian() ? 1 : -1;
//   Operator Z_debug(Z);
   // Create Pandya-transformed hp and ph matrix elements
   double t_start = omp_get_wtime();
//   deque<arma::mat> Y_bar_ph_all (InitializePandya( nChannels, "normal"));
//   deque<arma::mat> Xt_bar_ph_all (InitializePandya( nChannels, "transpose"));

//   Y.DoPandyaTransformation(Y_bar_ph_all, "normal" );
//   X.DoPandyaTransformation(Xt_bar_ph_all ,"transpose");
//   profiler.timer["DoPandyaTransformation"] += omp_get_wtime() - t_start;

   // Construct the intermediate matrix Z_bar
   const auto& pandya_lookup = modelspace->GetPandyaLookup(rank_J, rank_T, parity);
   int nch = modelspace->SortedTwoBodyChannels_CC.size();
   t_start = omp_get_wtime();
   deque<arma::mat> Z_bar (nChannels );
   vector<bool> lookup_empty(nChannels,true);
   for (int ich=0;ich<nch;++ich)
   {
      int ch = modelspace->SortedTwoBodyChannels_CC[ich];
//      if ( pandya_lookup.at({ch,ch})[0].size()<1 ) continue;
      index_t nKets_cc = modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
      Z_bar[ch].zeros( nKets_cc, 2*nKets_cc );
      if ( pandya_lookup.at({ch,ch})[0].size()>0 ) lookup_empty[ich] = false;
//      lookup_empty[ich] = false;
   }


   #ifndef OPENBLAS_NOUSEOMP
//   #pragma omp parallel for schedule(dynamic,1)
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->scalar_transform_first_pass)
   #endif
   for (int ich=0; ich<nch; ++ich )
   {
      if (lookup_empty.at(ich)) continue;
      int ch = modelspace->SortedTwoBodyChannels_CC.at(ich);
//      if ( pandya_lookup.find({ch,ch}) == pandya_lookup.end()) continue;
//      if ( pandya_lookup.at({ch,ch})[0].size()<1 ) continue;
      const TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch);
      index_t nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.GetKetIndex_hh().size() + tbc_cc.GetKetIndex_ph().size();

//      arma::mat Y_bar_ph(2*nph_kets,   nKets_cc, arma::fill::zeros);
//      arma::mat Xt_bar_ph(nKets_cc, 2*nph_kets,   arma::fill::zeros);
      arma::mat Y_bar_ph;
      arma::mat Xt_bar_ph;

      Y.DoPandyaTransformation_SingleChannel(Y_bar_ph,ch,"normal");
      X.DoPandyaTransformation_SingleChannel(Xt_bar_ph,ch,"transpose");
      auto& Zbar_ch = Z_bar.at(ch);

//      auto& Xt_bar_ph = Xt_bar_ph_all[ch];
//      auto& Y_bar_ph = Y_bar_ph_all[ch];

      if (Y_bar_ph.size()<1 or Xt_bar_ph.size()<1)
      {
//        Z_bar[ch] = arma::zeros( Xt_bar_ph.n_rows, Y_bar_ph.n_cols*2);
        Zbar_ch = arma::zeros( Xt_bar_ph.n_rows, Y_bar_ph.n_cols*2);
//        Zbar_ch.zeros();
        continue;
      }

      // get the phases for taking the transpose
      arma::mat PhaseMat(nKets_cc, nKets_cc, arma::fill::ones );
      for (index_t iket=0;iket<nKets_cc;iket++)
      {
         const Ket& ket = tbc_cc.GetKet(iket);
         if ( modelspace->phase( (ket.op->j2 + ket.oq->j2)/2 ) > 0) continue;
         PhaseMat.col( iket ) *= -1;
         PhaseMat.row( iket ) *= -1;
      }
      arma::uvec phkets = arma::join_cols(tbc_cc.GetKetIndex_hh(), tbc_cc.GetKetIndex_ph() );
      auto PhaseMatY = PhaseMat.rows(phkets) * hy;

//      Z_bar[ch] =  (Xt_bar_ph[ch] * Y_bar_ph[ch]);
//      arma::mat Z_bar =  (Xt_bar_ph[ch] * Y_bar_ph[ch]);

//                                           [      |     ]
//     create full Y matrix from the half:   [  Yhp | Y'ph]   where the prime indicates multiplication by (-1)^(i+j+k+l) h_y
//                                           [      |     ]   Flipping hp <-> ph and multiplying by the phase is equivalent to
//                                           [  Yph | Y'hp]   having kets |kj> with k>j.
//      int halfnry = Y_bar_ph.n_rows/2;
//      arma::mat Z_bar =  Xt_bar_ph * join_horiz(Y_bar_ph, join_vert(Y_bar_ph.tail_rows(halfnry)%PhaseMatY,
//                                                                    Y_bar_ph.head_rows(halfnry)%PhaseMatY) );
//      Z_bar[ch] =  Xt_bar_ph * join_horiz(Y_bar_ph, join_vert( Y_bar_ph.tail_rows(nph_kets)%PhaseMatY,
      Zbar_ch =  Xt_bar_ph * join_horiz(Y_bar_ph, join_vert(   Y_bar_ph.tail_rows(nph_kets)%PhaseMatY,
                                                               Y_bar_ph.head_rows(nph_kets)%PhaseMatY) );



      // If Z is hermitian, then XY is anti-hermitian, and so XY - YX = XY + (XY)^T
      if ( Z.IsHermitian() )
      {
//         Z_bar.cols(0,nKets_cc-1) += Z_bar.cols(0,nKets_cc-1).t();
//         Z_bar[ch].head_cols(nKets_cc) += Z_bar[ch].head_cols(nKets_cc).t();
         Zbar_ch.head_cols(nKets_cc) += Zbar_ch.head_cols(nKets_cc).t();
      }
      else
      {
//         Z_bar.cols(0,nKets_cc-1) -= Z_bar.cols(0,nKets_cc-1).t();
//         Z_bar[ch].head_cols(nKets_cc) -= Z_bar[ch].head_cols(nKets_cc).t();
         Zbar_ch.head_cols(nKets_cc) -= Zbar_ch.head_cols(nKets_cc).t();
      }
//      Z_bar.cols(nKets_cc,2*nKets_cc-1) += Z_bar.cols(nKets_cc,2*nKets_cc-1).t()%PhaseMat;
//      Z_bar[ch].tail_cols(nKets_cc) += Z_bar[ch].tail_cols(nKets_cc).t()%PhaseMat;
      Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t()%PhaseMat;

//     cout << "ch = " << ch << " --> nKets_cc = " << nKets_cc << "  size of Z_bar = " << Z_bar[ch].n_rows << " x " << Z_bar[ch].n_cols << endl;
//     Z_debug.AddInversePandyaTransformation_SingleChannel(Z_bar[ch],ch);
//     Z.AddInversePandyaTransformation_SingleChannel(Z_bar[ch],ch);
//     Z.AddInversePandyaTransformation_SingleChannel(Z_bar,ch);

   }

   profiler.timer["Build Z_bar"] += omp_get_wtime() - t_start;

   // Perform inverse Pandya transform on Z_bar to get Z
   t_start = omp_get_wtime();
   Z.AddInversePandyaTransformation(Z_bar);
//   for (auto& itch : Z_debug.TwoBody.MatEl)
//   {
//     cout << "ch = " << itch.first[0] << endl;
//     cout << itch.second << endl << endl << Z.TwoBody.GetMatrix(itch.first) << endl << endl << endl << endl;
//   }
//   cout << "end of commutator: " << endl << Z.TwoBody.GetMatrix(0).submat(0,0,1,1) << endl << endl << Z_debug.TwoBody.GetMatrix(0).submat(0,0,1,1);
//   cout << "done with ph commutator" << endl;
   modelspace->scalar_transform_first_pass = false;
   profiler.timer["InversePandyaTransformation"] += omp_get_wtime() - t_start;

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
   double tstart = omp_get_wtime();
   comm111ss(X,Y);
   profiler.timer["comm111st"] += omp_get_wtime() - tstart;
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
// X is scalar one-body, Y is tensor two-body
// There must be a better way to do this looping. 
//
//void Operator::comm121st( Operator& Y, Operator& Z) 
void Operator::comm121st( const Operator& X, const Operator& Y) 
{

   double tstart = omp_get_wtime();
   Operator& Z = *this;
   int norbits = modelspace->GetNumberOrbits();
   int Lambda = Z.GetJRank();
   modelspace->PreCalculateSixJ();
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(Z.rank_J))
   for (int i=0;i<norbits;++i)
   {
      Orbit &oi = modelspace->GetOrbit(i);
      double ji = 0.5*oi.j2;
      for (int j : Z.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) ) 
      {
          Orbit &oj = modelspace->GetOrbit(j);
          double jj = 0.5*oj.j2;
          if (j<i) continue; // only calculate upper triangle
          double& Zij = Z.OneBody(i,j);
          for (auto& a : modelspace->holes)  // C++11 syntax
          {
             Orbit &oa = modelspace->GetOrbit(a);
             double ja = 0.5*oa.j2;
             for (auto& b : X.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {
                Orbit &ob = modelspace->GetOrbit(b);
                double nanb = oa.occ * (1-ob.occ);
                  int J1min = std::abs(ji-ja);
                  int J1max = ji + ja;
                  for (int J1=J1min; J1<=J1max; ++J1)
                  {
                    int phasefactor = modelspace->phase(jj+ja+J1+Lambda);
                    int J2min = max(std::abs(Lambda - J1),std::abs(int(ja-jj)));
                    int J2max = min(Lambda + J1,int(ja+jj));
                    for (int J2=J2min; J2<=J2max; ++J2)
                    {
                      if ( ! ( J2>=std::abs(ja-jj) and J2<=ja+jj )) continue;
                      double prefactor = nanb*phasefactor * sqrt((2*J1+1)*(2*J2+1)) * modelspace->GetSixJ(J1,J2,Lambda,jj,ji,ja);
                      Zij +=  prefactor * ( X.OneBody(a,b) * Y.TwoBody.GetTBME_J(J1,J2,b,i,a,j) - X.OneBody(b,a) * Y.TwoBody.GetTBME_J(J1,J2,a,i,b,j ));
                    }
                  }
             }
             // Now, X is scalar two-body and Y is tensor one-body
             for (auto& b : Y.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
             {

                Orbit &ob = modelspace->GetOrbit(b);
                double jb = 0.5*ob.j2;
                if (std::abs(ob.occ-1) < OCC_CUT) continue;
                double nanb = oa.occ * (1-ob.occ);
                int J1min = max(std::abs(ji-jb),std::abs(jj-ja));
                int J1max = min(ji+jb,jj+ja);
                double zij = 0;
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij -= modelspace->phase(ji+jb+J1) * (2*J1+1) * modelspace->GetSixJ(ja,jb,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,b,i,a,j);
                }

                J1min = max(std::abs(ji-ja),std::abs(jj-jb));
                J1max = min(ji+ja,jj+jb);
                for (int J1=J1min; J1<=J1max; ++J1)
                {
                  zij += modelspace->phase(ji+jb+J1) * (2*J1+1) * modelspace->GetSixJ(jb,ja,Lambda,ji,jj,J1) * X.TwoBody.GetTBME_J(J1,J1,a,i,b,j) ;
                }

                Zij += nanb * Y.OneBody(a,b) * zij;
             }

             
          }
      }
   }
   
   profiler.timer["comm121st"] += omp_get_wtime() - tstart;
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
   double tstart = omp_get_wtime();
   Operator& Z = *this;
   int Lambda = Z.rank_J;

    vector< int > bra_channels;
    vector< int > ket_channels;
//    vector< array<int,2> > channels;
//    for ( auto& itmat : Z.TwoBody.MatEl ) channels.push_back( itmat.first );
    for ( auto& itmat : Z.TwoBody.MatEl )
    {
      bra_channels.push_back( itmat.first[0] );
      ket_channels.push_back( itmat.first[1] );
    }
    int nmat = bra_channels.size();
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(Z.rank_J))
    for (int ii=0; ii<nmat; ++ii)
    {
     int ch_bra = bra_channels[ii];
     int ch_ket = ket_channels[ii];

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
            double c1 = 0;
            double c2 = 0;
            double c3 = 0;
            double c4 = 0;

            for ( int a : X.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
              c1 += X.OneBody(i,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,a,j,k,l);
            }
            if (i==j)
            {
              c2 = c1; // there should be a phase here, but if the ket exists, it'd better be +1.
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 c2 += X.OneBody(j,a) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,a,k,l);
              }
            }
            for ( int a : X.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               c3 += X.OneBody(a,k) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,a,l);
            }
            if (k==l)
            {
              c4 = c3 ;
            }
            else
            {
              for ( int a : X.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 c4 += X.OneBody(a,l) * Y.TwoBody.GetTBME(ch_bra,ch_ket,i,j,k,a);
              }
            }

            cijkl = c1 + c2 - c3 - c4;


            c1=0;
            c2=0;
            c3=0;
            c4=0;
            int phase1 = modelspace->phase(ji+jj+J2+Lambda);
            int phase2 = modelspace->phase(J1-J2+Lambda);
            int phase3 = modelspace->phase(J1-J2+Lambda);
            int phase4 = modelspace->phase(jk+jl-J1+Lambda);


            for ( int a : Y.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2*0.5;
               c1 -=   modelspace->GetSixJ(J2,J1,Lambda,ji,ja,jj) * Y.OneBody(i,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,j,k,l) ;
            }
            if (i==j)
            {
              c2 = -c1;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({oj.l,oj.j2,oj.tz2}) )
              {
                 double ja = modelspace->GetOrbit(a).j2*0.5;
                 c2 +=  modelspace->GetSixJ(J2,J1,Lambda,jj,ja,ji) * Y.OneBody(j,a) * X.TwoBody.GetTBME(ch_ket,ch_ket,a,i,k,l);
              }
            }
            for ( int a : Y.OneBodyChannels.at({ok.l,ok.j2,ok.tz2}) )
            {
               double ja = modelspace->GetOrbit(a).j2*0.5;
               c3 -=  modelspace->GetSixJ(J1,J2,Lambda,jk,ja,jl) * Y.OneBody(a,k) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,l,a) ;
            }
            if (k==l)
            {
              c4 = -c3;
            }
            else
            {
              for ( int a : Y.OneBodyChannels.at({ol.l,ol.j2,ol.tz2}) )
              {
                 double ja = modelspace->GetOrbit(a).j2*0.5;
                 c4 +=  modelspace->GetSixJ(J1,J2,Lambda,jl,ja,jk) * Y.OneBody(a,l) * X.TwoBody.GetTBME(ch_bra,ch_bra,i,j,k,a) ;
              }
            }
            cijkl += hatfactor*(phase1*c1+phase2*c2+phase3*c3+phase4*c4);


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Z2(ibra,iket) += cijkl /norm;
         }
      }
   }
   profiler.timer["comm122st"] += omp_get_wtime() - tstart;
}





// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
// X is a scalar, Y is a tensor
//void Operator::comm222_pp_hh_221st( Operator& Y, Operator& Z )  
void Operator::comm222_pp_hh_221st( const Operator& X, const Operator& Y )  
{

   double tstart = omp_get_wtime();
   Operator& Z = *this;
   int Lambda = Z.GetJRank();

   TwoBodyME Mpp = Y.TwoBody;
   TwoBodyME Mhh = Y.TwoBody;
   TwoBodyME Mff = Y.TwoBody;

   vector<int> vch_bra;
   vector<int> vch_ket;
   vector<const arma::mat*> vmtx;
   for ( auto& itmat : Y.TwoBody.MatEl )
   {
     vch_bra.push_back(itmat.first[0]);
     vch_ket.push_back(itmat.first[1]);
     vmtx.push_back(&(itmat.second));
   }
   size_t nchan = vch_bra.size();
//   for ( auto& itmat : Y.TwoBody.MatEl )
   #pragma omp parallel for schedule(dynamic,1)
   for (size_t i=0;i<nchan; ++i)
   {
    int ch_bra = vch_bra[i];
    int ch_ket = vch_ket[i];
//    auto& Z2 = Z.TwoBody.GetMatrix(ch_bra,ch_ket);

    TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);

    auto& LHS1 = X.TwoBody.GetMatrix(ch_bra,ch_bra);
    auto& LHS2 = X.TwoBody.GetMatrix(ch_ket,ch_ket);

    auto& RHS  =  *vmtx[i];

    arma::mat& Matrixpp =  Mpp.GetMatrix(ch_bra,ch_ket);
    arma::mat& Matrixhh =  Mhh.GetMatrix(ch_bra,ch_ket);
   
    const arma::uvec& bras_pp = tbc_bra.GetKetIndex_pp();
    const arma::uvec& bras_hh = tbc_bra.GetKetIndex_hh();
    const arma::uvec& bras_ph = tbc_bra.GetKetIndex_ph();
    const arma::uvec& kets_pp = tbc_ket.GetKetIndex_pp();
    const arma::uvec& kets_hh = tbc_ket.GetKetIndex_hh();
    const arma::uvec& kets_ph = tbc_ket.GetKetIndex_ph();

    // the complicated-looking construct after the % signs just multiply the matrix elements by the proper occupation numbers (nanb, etc.)

    arma::mat MLeft  = join_horiz( LHS1.cols(bras_hh) , -RHS.cols(kets_hh) );
    arma::mat MRight = join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_occ_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros ) ),
                                 LHS2.rows(kets_hh)  % tbc_ket.Ket_occ_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ));

    Matrixhh = MLeft * MRight;


    MLeft  = join_horiz( LHS1.cols(join_vert(bras_pp,join_vert(bras_hh,bras_ph))), -RHS.cols(join_vert(kets_pp,join_vert(kets_hh,kets_ph))) );
    MRight = join_vert( join_vert(     RHS.rows(bras_pp), 
                          join_vert( RHS.rows(bras_hh)  % tbc_bra.Ket_unocc_hh.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) )  ,
                                     RHS.rows(bras_ph)  % tbc_bra.Ket_unocc_ph.cols( arma::uvec(RHS.n_cols,arma::fill::zeros) ) )),
                      join_vert(     LHS2.rows(kets_pp),
                          join_vert( LHS2.rows(kets_hh) % tbc_ket.Ket_unocc_hh.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) ),
                                     LHS2.rows(kets_ph) % tbc_ket.Ket_unocc_ph.cols( arma::uvec(LHS2.n_cols,arma::fill::zeros) )))
                     );

    Matrixpp = MLeft * MRight;
                                

    if (Z.GetParticleRank()>1)
    {
      Z.TwoBody.GetMatrix(ch_bra,ch_ket) += Matrixpp - Matrixhh;
    }
//    Z2 += Matrixpp - Matrixhh;

   }// for itmat

      // The one body part takes some additional work

   int norbits = modelspace->GetNumberOrbits();
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(Z.rank_J))
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
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = max(int(std::abs(jc-jj)),std::abs(Lambda-J1));
               int j2max = min(int(jc+jj),J1+Lambda); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * modelspace->phase(jj + jc + J1 + Lambda)
                        * ( oc.occ * Mpp.GetTBME_J(J1,J2,c,i,c,j) + (1-oc.occ) * Mhh.GetTBME_J(J1,J2,c,i,c,j));
//                cijJ += hatfactor * sixj * modelspace->phase(jj + jc + J1 + Lambda) * (1-oc.occ) * Mhh.GetTBME_J(J1,J2,c,i,c,j);  // This is probably right???
               }
              }
           // Sum c over particles and include the n_a * n_b terms
           }
           for (auto& c : modelspace->particles)
           {
              Orbit &oc = modelspace->GetOrbit(c);
              double jc = oc.j2/2.0;
              int j1min = std::abs(jc-ji);
              int j1max = jc+ji;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
               int j2min = max(int(std::abs(jc-jj)),std::abs(Lambda-J1));
               int j2max = min(int(jc+jj),J1+Lambda); 
               for (int J2=j2min; J2<=j2max; ++J2)
               {
                double hatfactor = sqrt( (2*J1+1)*(2*J2+1) );
                double sixj = modelspace->GetSixJ(J1, J2, Lambda, jj, ji, jc);
                cijJ += hatfactor * sixj * modelspace->phase(jj + jc + J1 + Lambda) * Mhh.GetTBME_J(J1,J2,c,i,c,j);
               }
              }
           }
//         #pragma omp critical
         Z.OneBody(i,j) += cijJ ;
      } // for j
    } // for i
    profiler.timer["comm222_pp_hh_221st"] += omp_get_wtime() - tstart;
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
//void Operator::DoTensorPandyaTransformation(map<array<int,2>,arma::mat>& TwoBody_CC_hp, map<array<int,2>,arma::mat>& TwoBody_CC_ph) const
void Operator::DoTensorPandyaTransformation( map<array<index_t,2>,arma::mat>& TwoBody_CC_ph) const
{
   int Lambda = rank_J;
   // loop over cross-coupled channels
   index_t nch = modelspace->SortedTwoBodyChannels_CC.size();

   // Allocate map for matrices -- this needs to be serial.
   for ( index_t ch_bra_cc : modelspace->SortedTwoBodyChannels_CC )
   {
      TwoBodyChannel& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
      index_t nph_bras = bras_ph.n_rows;
      for ( index_t ch_ket_cc : modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

         TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}] =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);
      }
   }

   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(Lambda))
   for (index_t ich=0;ich<nch;++ich)
   {
      index_t ch_bra_cc = modelspace->SortedTwoBodyChannels_CC[ich];
      TwoBodyChannel& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      int Jbra_cc = tbc_bra_cc.J;
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
//      arma::uvec& bras_ph = tbc_bra_cc.GetKetIndex_ph();
      index_t nph_bras = bras_ph.size();

      for ( index_t ch_ket_cc : modelspace->SortedTwoBodyChannels_CC )
      {
        TwoBodyChannel& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
        int Jket_cc = tbc_ket_cc.J;
        if ( (Jbra_cc+Jket_cc < rank_J) or std::abs(Jbra_cc-Jket_cc)>rank_J ) continue;
        if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + parity)%2>0 ) continue;

        index_t nKets_cc = tbc_ket_cc.GetNumberKets();

//        arma::mat& MatCC_hp = TwoBody_CC_hp[{ch_bra_cc,ch_ket_cc}];
        arma::mat& MatCC_ph = TwoBody_CC_ph[{ch_bra_cc,ch_ket_cc}];
        // loop over ph bras <ad| in this channel
        for (index_t ibra=0; ibra<nph_bras; ++ibra)
        {
           Ket & bra_cc = tbc_bra_cc.GetKet( bras_ph[ibra] );
           index_t a = bra_cc.p;
           index_t b = bra_cc.q;
           Orbit & oa = modelspace->GetOrbit(a);
           Orbit & ob = modelspace->GetOrbit(b);
           double ja = oa.j2*0.5;
           double jb = ob.j2*0.5;

           // loop over kets |bc> in this channel
           index_t iket_max =  nKets_cc ;
           for (index_t iket_cc=0; iket_cc<iket_max; ++iket_cc)
           {
              Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
              index_t c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
              index_t d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
              Orbit & oc = modelspace->GetOrbit(c);
              Orbit & od = modelspace->GetOrbit(d);
              double jc = oc.j2*0.5;
              double jd = od.j2*0.5;


              int j1min = std::abs(ja-jd);
              int j1max = ja+jd;
              double sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = max(int(std::abs(jc-jb)),std::abs(J1-Lambda));
                int j2max = min(int(jc+jb),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = TwoBody.GetTBME_J(J1,J2,a,d,c,b);
                  sm -= hatfactor * modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
//              MatCC_hp(ibra,iket_cc) = sm;
              MatCC_ph(ibra,iket_cc) = sm;

              // Exchange (a <-> b) to account for the (n_a - n_b) term
              // Get Tz,parity and range of J for <bd || ca > coupling
              j1min = std::abs(jb-jd);
              j1max = jb+jd;
              sm = 0;
              for (int J1=j1min; J1<=j1max; ++J1)
              {
                int j2min = max(int(std::abs(jc-ja)),std::abs(J1-Lambda));
                int j2max = min(int(jc+ja),J1+Lambda);
                for (int J2=j2min; J2<=j2max; ++J2)
                {
                  double ninej = modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
                  double tbme = TwoBody.GetTBME_J(J1,J2,b,d,c,a);
                  sm -= hatfactor * modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
                }
              }
//              MatCC_ph(ibra,iket_cc) = sm;
              MatCC_ph(ibra+nph_bras,iket_cc) = sm;

           }
        }
    }
   }
}



// This happens inside an OMP loop, and so everything here needs to be thread safe
//
//void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& TwoBody_CC_ph, int ch_bra_cc, int ch_ket_cc) const
void Operator::DoTensorPandyaTransformation_SingleChannel( arma::mat& MatCC_ph, int ch_bra_cc, int ch_ket_cc) const
{
   int Lambda = rank_J;

   TwoBodyChannel& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
   arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
   int nph_bras = bras_ph.n_rows;

   TwoBodyChannel& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
   int nKets_cc = tbc_ket_cc.GetNumberKets();

   // The Pandya-transformed (formerly cross-coupled) particle-hole type matrix elements
   // (this is the output of this method)
   MatCC_ph =  arma::mat(2*nph_bras,   nKets_cc, arma::fill::zeros);

   int Jbra_cc = tbc_bra_cc.J;
   int Jket_cc = tbc_ket_cc.J;
   if ( (Jbra_cc+Jket_cc < rank_J) or std::abs(Jbra_cc-Jket_cc)>rank_J ) return;
   if ( (tbc_bra_cc.parity + tbc_ket_cc.parity + parity)%2>0 ) return;


//   arma::mat& MatCC_ph = TwoBody_CC_ph;
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
      int iket_max =  nKets_cc ;
      for (int iket_cc=0; iket_cc<iket_max; ++iket_cc)
      {
         Ket & ket_cc = tbc_ket_cc.GetKet(iket_cc%nKets_cc);
         int c = iket_cc < nKets_cc ? ket_cc.p : ket_cc.q;
         int d = iket_cc < nKets_cc ? ket_cc.q : ket_cc.p;
         Orbit & oc = modelspace->GetOrbit(c);
         Orbit & od = modelspace->GetOrbit(d);
         double jc = oc.j2*0.5;
         double jd = od.j2*0.5;


         int j1min = std::abs(ja-jd);
         int j1max = ja+jd;
         double sm = 0;
         for (int J1=j1min; J1<=j1max; ++J1)
         {
           int j2min = max(int(std::abs(jc-jb)),std::abs(J1-Lambda));
           int j2max = min(int(jc+jb),J1+Lambda);
           for (int J2=j2min; J2<=j2max; ++J2)
           {
             double ninej = modelspace->GetNineJ(ja,jd,J1,jb,jc,J2,Jbra_cc,Jket_cc,Lambda);
             if (std::abs(ninej) < 1e-8) continue;
             double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
             double tbme = TwoBody.GetTBME_J(J1,J2,a,d,c,b);
             sm -= hatfactor * modelspace->phase(jb+jd+Jket_cc+J2) * ninej * tbme ;
           }
         }
         MatCC_ph(ibra,iket_cc) = sm;

         // Exchange (a <-> b) to account for the (n_a - n_b) term
         if (a==b)
         {
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
         else
         {

           // Get Tz,parity and range of J for <bd || ca > coupling
           j1min = std::abs(jb-jd);
           j1max = jb+jd;
           sm = 0;
           for (int J1=j1min; J1<=j1max; ++J1)
           {
             int j2min = max(int(std::abs(jc-ja)),std::abs(J1-Lambda));
             int j2max = min(int(jc+ja),J1+Lambda);
             for (int J2=j2min; J2<=j2max; ++J2)
             {
               double ninej = modelspace->GetNineJ(jb,jd,J1,ja,jc,J2,Jbra_cc,Jket_cc,Lambda);
               if (std::abs(ninej) < 1e-8) continue;
               double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*Jbra_cc+1)*(2*Jket_cc+1) );
               double tbme = TwoBody.GetTBME_J(J1,J2,b,d,c,a);
               sm -= hatfactor * modelspace->phase(ja+jd+Jket_cc+J2) * ninej * tbme ;
             }
           }
           MatCC_ph(ibra+nph_bras,iket_cc) = sm;
         }
      }
   }
}





// Take Pandya-transformed matrix Zbar for a single channel, invert the Pandya transformation and add the result to the current operator.
void Operator::AddInverseTensorPandyaTransformation_SingleChannel(arma::mat& Zbar, int ch_bra_cc, int ch_ket_cc)
{
    // Do the inverse Pandya transform
   if (ch_bra_cc > ch_ket_cc)  // hopefully this won't happen
   {
      cout << "WARNING: Called Operator::AddInverseTensorPandyaTransformation_SingleChannel with ch_bra_cc > ch_ket_cc : " << ch_bra_cc << " > " << ch_ket_cc << endl;
      cout << " Skipping this channel." << endl;
      return;
   }
   int n_channels_kept = 0;
//   cout << " [ " << ch_bra_cc << " , " << ch_ket_cc << " ] " << endl;
   Operator& Z = *this;
   const auto& pandya_lookup = modelspace->GetPandyaLookup(rank_J, rank_T, parity)[{ch_bra_cc,ch_ket_cc}];
   int Lambda = Z.rank_J;
   int hZ = Z.IsHermitian() ? 1 : -1;
   TwoBodyChannel_CC& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
   TwoBodyChannel_CC& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
   int J3 = tbc_bra_cc.J;
   int J4 = tbc_ket_cc.J;
   int Tz_bra_cc = tbc_bra_cc.Tz;
   int Tz_ket_cc = tbc_ket_cc.Tz;
   int parity_bra_cc = tbc_bra_cc.parity;
   int parity_ket_cc = tbc_ket_cc.parity;
//   int nbras_cc = tbc_bra_cc.GetNumberKets();
   int nkets_cc = tbc_ket_cc.GetNumberKets();

//   for ( auto& iter : Z.TwoBody.MatEl )
   int nchannels = pandya_lookup[0].size();
//   for ( auto& itchan : pandya_lookup )
   for (int ich=0;ich<nchannels;++ich)
   {
//      int ch_bra = iter.first[0];
//      int ch_ket = iter.first[1];
//      int ch_bra = itchan[0];
//      int ch_ket = itchan[1];
      int ch_bra = pandya_lookup[0][ich];
      int ch_ket = pandya_lookup[1][ich];
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      int nBras = tbc_bra.GetNumberKets();
      int nKets = tbc_ket.GetNumberKets();
//      arma::mat& Zijkl = iter.second;
//      arma::mat& Zijkl = Z.TwoBody.GetMatrix(itchan);
      arma::mat& Zijkl = Z.TwoBody.GetMatrix(ch_bra,ch_ket);
      bool inner_loop = false;

      double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );

      for (int ibra=0; ibra<nBras; ++ibra)
      {
         Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         double ji = 0.5*oi.j2;
         double jj = 0.5*oj.j2;
         int ketmin = ch_bra==ch_ket ? ibra : 0;
         for (int iket=ketmin; iket<nKets; ++iket)
         {
            Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            Orbit & ok = modelspace->GetOrbit(k);
            Orbit & ol = modelspace->GetOrbit(l);
            double jk = 0.5*ok.j2;
            double jl = 0.5*ol.j2;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int j3min = std::abs(int(ji-jl));
            int j3max = ji+jl;
            int j4min = std::abs(int(jk-jj));
            int j4max = int(jk+jj);

            if (   (oi.l+ol.l)%2==parity_bra_cc             and (ok.l+oj.l)%2==parity_ket_cc
                      and std::abs(oi.tz2+ol.tz2)==2*Tz_bra_cc   and std::abs(ok.tz2+oj.tz2)==2*Tz_ket_cc
                      and j3min<=J3 and J3<=j3max           and j4min<=J4 and J4<=j4max )
            {
               int indx_il = tbc_bra_cc.GetLocalIndex(min(i,l),max(i,l));
               int indx_kj = tbc_ket_cc.GetLocalIndex(min(j,k),max(j,k));

               double ninej = modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
               double tbme = 0;

               if (i<=l) tbme = Zbar( indx_il ,indx_kj+(k>j?nkets_cc:0) );
               else      tbme = Zbar( indx_il ,indx_kj+(k>j?0:nkets_cc) ) * hZ * modelspace->phase(J3-J4 + ji+jj+jk+jl);
               commij += hatfactor * modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;
               inner_loop = true;
            }

            if (  ch_bra_cc != ch_ket_cc  and (oi.l+ol.l)%2==parity_ket_cc        and (ok.l+oj.l)%2==parity_bra_cc
                                            and std::abs(oi.tz2+ol.tz2)==2*Tz_ket_cc   and std::abs(ok.tz2+oj.tz2)==2*Tz_bra_cc
                                            and j3min<=J4 and J4<=j3max           and j4min<=J3 and J3<=j4max )
              {
                 int indx_kj = tbc_bra_cc.GetLocalIndex(min(j,k),max(j,k));
                 int indx_il = tbc_ket_cc.GetLocalIndex(min(i,l),max(i,l));

                 double ninej = modelspace->GetNineJ(ji,jl,J4,jj,jk,J3,J1,J2,Lambda);
                 double tbme = 0;

                 if(k<=j) tbme = Zbar(indx_kj, indx_il+(i>l?nkets_cc:0) ) * hZ * modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                 else     tbme = Zbar(indx_kj, indx_il+(i>l?0:nkets_cc) ) * modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)

                 commij += hatfactor * modelspace->phase(jj+jl+J2+J3) * ninej * tbme ;
               inner_loop = true;
              }


            if (i==j)
            {
              commji = commij;
            }
            else
            {
              j3min = std::abs(int(jj-jl));
              j3max = jj+jl;
              j4min = std::abs(int(jk-ji));
              j4max = int(jk+ji);
              if (   (oj.l+ol.l)%2==parity_bra_cc and (ok.l+oi.l)%2==parity_ket_cc
                      and std::abs(oj.tz2+ol.tz2)==2*Tz_bra_cc and std::abs(ok.tz2+oi.tz2)==2*Tz_ket_cc
                      and J3>=j3min and J3<=j3max and J4>=j4min and J4<=j4max)
              {
                // Transform Z_jlki
                int indx_jl = tbc_bra_cc.GetLocalIndex(min(j,l),max(j,l));
                int indx_ki = tbc_ket_cc.GetLocalIndex(min(k,i),max(k,i));
                double ninej = modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                double tbme = 0;
                if (j<=l) tbme = Zbar( indx_jl ,indx_ki+(k>i?nkets_cc:0) );
                else      tbme = Zbar( indx_jl ,indx_ki+(k>i?0:nkets_cc) ) * hZ * modelspace->phase(J3-J4 + ji+jj+jk+jl);
                commji += hatfactor * modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
               inner_loop = true;
              }
              if (ch_bra_cc!=ch_ket_cc and (oj.l+ol.l)%2==parity_ket_cc and (ok.l+oi.l)%2==parity_bra_cc
                                       and std::abs(oj.tz2+ol.tz2)==2*Tz_ket_cc and std::abs(ok.tz2+oi.tz2)==2*Tz_bra_cc
                                       and J4>=j3min and J4<=j3max and J3>=j4min and J3<=j4max)
              {
                int indx_jl = tbc_ket_cc.GetLocalIndex(min(j,l),max(j,l));
                int indx_ki = tbc_bra_cc.GetLocalIndex(min(k,i),max(k,i));
                double ninej = modelspace->GetNineJ(jj,jl,J4,ji,jk,J3,J1,J2,Lambda);
                double tbme = 0;

                if(k<=i) tbme = Zbar(indx_ki, indx_jl+(j>l?nkets_cc:0) ) * hZ * modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                else     tbme = Zbar(indx_ki, indx_jl+(j>l?0:nkets_cc) ) * modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                commji += hatfactor * modelspace->phase(ji+jl+J2+J3) * ninej * tbme ;
               inner_loop = true;
              }
            }


            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            #pragma omp atomic
            Zijkl(ibra,iket) +=  (commij - modelspace->phase(ji+jj-J1)*commji) / norm;
            if (ch_bra==ch_ket) 
            {
               #pragma omp atomic write
               Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
            }
         }
      }
            if (inner_loop)
            {
//              cout << "    ( " << ch_bra << " , " << ch_ket  << " ) " << endl;
              n_channels_kept++;
            }
   }
//   cout << "Kept " << n_channels_kept << " out of " << Z.TwoBody.MatEl.size() << " channels.  Size of Zbar is " << Zbar.n_rows << "x" << Zbar.n_cols << endl;
}





void Operator::AddInverseTensorPandyaTransformation( const map<array<index_t,2>,arma::mat>&  Zbar )
{
    // Do the inverse Pandya transform
   Operator& Z = *this;
   int Lambda = Z.rank_J;
   vector<map<array<int,2>,arma::mat>::iterator> iteratorlist;
   for (map<array<int,2>,arma::mat>::iterator iter= Z.TwoBody.MatEl.begin(); iter!= Z.TwoBody.MatEl.end(); ++iter) iteratorlist.push_back(iter);
   int niter = iteratorlist.size();
   int hZ = Z.IsHermitian() ? 1 : -1;
    // Only go parallel if we've previously calculated the SixJs/NineJs. Otherwise, it's not thread safe.
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(Lambda))
   for (int i=0; i<niter; ++i)
   {
      const auto iter = iteratorlist[i];
      int ch_bra = iter->first[0];
      int ch_ket = iter->first[1];
      arma::mat& Zijkl = iter->second;
      const TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      const TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J1 = tbc_bra.J;
      int J2 = tbc_ket.J;
      index_t nBras = tbc_bra.GetNumberKets();
      index_t nKets = tbc_ket.GetNumberKets();

      for (index_t ibra=0; ibra<nBras; ++ibra)
      {
         const Ket & bra = tbc_bra.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         const Orbit & oi = modelspace->GetOrbit(i);
         const Orbit & oj = modelspace->GetOrbit(j);
         double ji = oi.j2/2.;
         double jj = oj.j2/2.;
         index_t ketmin = ch_bra==ch_ket ? ibra : 0;
         for (index_t iket=ketmin; iket<nKets; ++iket)
         {
            const Ket & ket = tbc_ket.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            const Orbit & ok = modelspace->GetOrbit(k);
            const Orbit & ol = modelspace->GetOrbit(l);
            double jk = ok.j2/2.;
            double jl = ol.j2/2.;

            double commij = 0;
            double commji = 0;

            // Transform Z_ilkj
            int parity_bra_cc = (oi.l+ol.l)%2;
            int parity_ket_cc = (ok.l+oj.l)%2;
            int Tz_bra_cc = std::abs(oi.tz2+ol.tz2)/2;
            int Tz_ket_cc = std::abs(ok.tz2+oj.tz2)/2;
            int j3min = std::abs(int(ji-jl));
            int j3max = ji+jl;
            for (int J3=j3min; J3<=j3max; ++J3)
            {
              index_t ch_bra_cc = modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
              const TwoBodyChannel_CC& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
              index_t nbras = tbc_bra_cc.GetNumberKets();
              index_t indx_il = tbc_bra_cc.GetLocalIndex(min(i,l),max(i,l));
              int j4min = max(std::abs(int(jk-jj)),std::abs(J3-Lambda));
              int j4max = min(int(jk+jj),J3+Lambda);
              for (int J4=j4min; J4<=j4max; ++J4)
              {
                 index_t ch_ket_cc = modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                 const TwoBodyChannel_CC& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                 index_t nkets = tbc_ket_cc.GetNumberKets();
                 index_t indx_kj = tbc_ket_cc.GetLocalIndex(min(j,k),max(j,k));

                  double ninej = modelspace->GetNineJ(ji,jl,J3,jj,jk,J4,J1,J2,Lambda);
                  if (std::abs(ninej) < 1e-8) continue;
                  double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
                  double tbme = 0;
                  index_t ch_lo = min(ch_bra_cc,ch_ket_cc);
                  index_t ch_hi = max(ch_bra_cc,ch_ket_cc);
                  auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                  if (zbar_iter == Zbar.end()) continue;
                  const auto& Zmat = zbar_iter->second;

                  if (ch_bra_cc <= ch_ket_cc)
                  {
                    if (i<=l) tbme = Zmat( indx_il ,indx_kj+(k>j?nkets:0) );
                    else      tbme = Zmat( indx_il ,indx_kj+(k>j?0:nkets) ) * hZ * modelspace->phase(J3-J4 + ji+jj+jk+jl);
                  }
                  else
                  {
                      if(k<=j) tbme = Zmat(indx_kj, indx_il+(i>l?nbras:0) ) * hZ * modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                      else     tbme = Zmat(indx_kj, indx_il+(i>l?0:nbras) ) * modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                  }
                  commij += hatfactor * modelspace->phase(jj+jl+J2+J4) * ninej * tbme ;

//                  if (J1==0 and J2==0 and i==0 and j==18 and k==1 and l==5 and std::abs(tbme)>1e-7 and std::abs(ninej)>1e-7)
//                  {
//                    cout << "    " << J3 << " " << J4 << "  " << hatfactor << " "  << modelspace->phase(jj+jl+J2+J4) << " " << ninej << " " << tbme << " " << commij << endl;
//                  }
              }
            }

            if (i==j)
            {
              commji = commij;
            }
            else
            {
              // Transform Z_jlki
              parity_bra_cc = (oj.l+ol.l)%2;
              parity_ket_cc = (ok.l+oi.l)%2;
              Tz_bra_cc = std::abs(oj.tz2+ol.tz2)/2;
              Tz_ket_cc = std::abs(ok.tz2+oi.tz2)/2;
              j3min = std::abs(int(jj-jl));
              j3max = jj+jl;
  
              for (int J3=j3min; J3<=j3max; ++J3)
              {
                int ch_bra_cc = modelspace->GetTwoBodyChannelIndex(J3,parity_bra_cc,Tz_bra_cc);
                const TwoBodyChannel_CC& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
                int nbras = tbc_bra_cc.GetNumberKets();
                int indx_jl = tbc_bra_cc.GetLocalIndex(min(j,l),max(j,l));
                int j4min = max(std::abs(int(jk-ji)),std::abs(J3-Lambda));
                int j4max = min(int(jk+ji),J3+Lambda);
                for (int J4=j4min; J4<=j4max; ++J4)
                {
                   int ch_ket_cc = modelspace->GetTwoBodyChannelIndex(J4,parity_ket_cc,Tz_ket_cc);
                   const TwoBodyChannel_CC& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
                   int nkets = tbc_ket_cc.GetNumberKets();
                   int indx_ki = tbc_ket_cc.GetLocalIndex(min(k,i),max(k,i));
                    double ninej = modelspace->GetNineJ(jj,jl,J3,ji,jk,J4,J1,J2,Lambda);
                    if (std::abs(ninej) < 1e-8) continue;
                    double hatfactor = sqrt( (2*J1+1)*(2*J2+1)*(2*J3+1)*(2*J4+1) );
  
                    index_t ch_lo = min(ch_bra_cc,ch_ket_cc);
                    index_t ch_hi = max(ch_bra_cc,ch_ket_cc);
                    auto zbar_iter = Zbar.find({ch_lo,ch_hi});
                    if (zbar_iter == Zbar.end()) continue;
                    const auto& Zmat = zbar_iter->second;
                    double tbme = 0;
                    if (ch_bra_cc <= ch_ket_cc)
                    {
                      if (j<=l) tbme = Zmat( indx_jl ,indx_ki+(k>i?nkets:0) );
                      else      tbme = Zmat( indx_jl ,indx_ki+(k>i?0:nkets) ) * hZ * modelspace->phase(J3-J4 + ji+jj+jk+jl);
                    }
                    else
                    {
                        if(k<=i) tbme = Zmat(indx_ki, indx_jl+(j>l?nbras:0) ) * hZ * modelspace->phase(J3-J4); // Z_ilkj = Z_kjil * (phase)
                        else     tbme = Zmat(indx_ki, indx_jl+(j>l?0:nbras) ) * modelspace->phase( ji+jj+jk+jl) ; // Z_ilkj = Z_kjil * (phase)
                    }
  
  
                      commji += hatfactor * modelspace->phase(ji+jl+J2+J4) * ninej * tbme ;
                }
              }
            }

            double norm = bra.delta_pq()==ket.delta_pq() ? 1+bra.delta_pq() : SQRT2;
            Zijkl(ibra,iket) +=  (commij - modelspace->phase(ji+jj-J1)*commji) / norm;
            if (ch_bra==ch_ket) Zijkl(iket,ibra) = hZ * Zijkl(ibra,iket);
//            if (J1==0 and J2==0 and i==0 and j==18 and k==1 and l==5 )
//            {
//              cout << "debug: adding term with commij = " << commij << "  and commji = " << commji << " :   (" << ch_bra << ", " << ch_ket << ") " << ibra << " " << iket << "  norm = " << norm << " -> " << Zijkl(ibra,iket) << endl;
//            }
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
//   Operator Z_debug1(Z);
//   Z_debug1.TwoBody.Scale(0.0);
//   Operator Z_debug2(Z_debug1);
   int hX = X.IsHermitian() ? 1 : -1;
   int hY = Y.IsHermitian() ? 1 : -1;
//   int hZ = Z.IsHermitian() ? 1 : -1;
   // Create Pandya-transformed hp and ph matrix elements
//   deque<arma::mat> X_bar_hp = InitializePandya( nChannels, "transpose");

   double t_start = omp_get_wtime();
   // We reuse Xt_bar multiple times, so it makes sense to calculate them once and store them in a deque.
   deque<arma::mat> Xt_bar_ph = InitializePandya( nChannels, "transpose"); // We re-use the scalar part multiple times, so there's a significant speed gain for saving it
   map<array<index_t,2>,arma::mat> Y_bar_ph;
   X.DoPandyaTransformation(Xt_bar_ph, "transpose" );
//   Y.DoTensorPandyaTransformation(Y_bar_ph );
   profiler.timer["DoTensorPandyaTransformation"] += omp_get_wtime() - t_start;


   t_start = omp_get_wtime();
   // Construct the intermediate matrix Z_bar
   // First, we initialize the map Z_bar with empty matrices
   // to avoid problems in the parallel loop -- (do we even want a parallel loop here?)
   map<array<index_t,2>,arma::mat> Z_bar;
//   vector<int> ybras(Y_bar_ph.size());
//   vector<int> ykets(Y_bar_ph.size());
//   int counter = 0;
//   for (auto& iter : Y_bar_ph )
//   {
//      ybras[counter] = iter.first[0];
//      ykets[counter] = iter.first[1];
//      if (iter.first[0] > iter.first[1]) continue;
//      Z_bar[iter.first].zeros( Xt_bar_ph[iter.first[0]].n_rows, 2*Xt_bar_ph[iter.first[1]].n_rows); // important to initialize this for parallelization
//      counter++;
//   }

   t_start = omp_get_wtime();
//   auto& pandya_lookup = modelspace->GetPandyaLookup(rank_J, rank_T, parity);
   const auto& pandya_lookup = modelspace->GetPandyaLookup(rank_J, rank_T, parity);
   profiler.timer["PandyaLookup"] += omp_get_wtime() - t_start;

   vector<index_t> ybras;
   vector<index_t> ykets;
//   cout << "start loop over TBC_CC" << endl;
   for (auto ich_bra : modelspace->SortedTwoBodyChannels_CC)
   {
     int n_rows = modelspace->GetTwoBodyChannel_CC(ich_bra).GetNumberKets();
     for (auto ich_ket : modelspace->SortedTwoBodyChannels_CC)
     {
       if (ich_bra>ich_ket) continue;
       if (pandya_lookup.at({(int)ich_bra,(int)ich_ket})[0].size()<1) continue;
         int n_cols = 2*modelspace->GetTwoBodyChannel_CC(ich_ket).GetNumberKets();
         ybras.push_back(ich_bra);
         ykets.push_back(ich_ket);
         Z_bar[{ich_bra,ich_ket}] = arma::mat(n_rows,n_cols);
     }
   }
   int counter = ybras.size();
//  cout << "Done allocating" << endl;


   profiler.timer["Allocate Z_bar_tensor"] += omp_get_wtime() - t_start;

   t_start = omp_get_wtime();

   #ifndef OPENBLAS_NOUSEOMP
   #pragma omp parallel for schedule(dynamic,1) if (not modelspace->tensor_transform_first_pass.at(rank_J))
   #endif
   for(int i=0;i<counter;++i)
   {
//      double t_start2 = omp_get_wtime();
      index_t ch_bra_cc = ybras[i];
      index_t ch_ket_cc = ykets[i];
//      cout << i << " " << ch_bra_cc << " " << ch_ket_cc << endl;
      const auto plookup = pandya_lookup.find({(int)ch_bra_cc,(int)ch_ket_cc});
      if ( plookup == pandya_lookup.end() or plookup->second[0].size()<1 )
      {
//       profiler.timer["BuildZbarTensor_setup"] += omp_get_wtime() - t_start2;
       continue;
      }

      const auto& tbc_bra_cc = modelspace->GetTwoBodyChannel_CC(ch_bra_cc);
      const auto& tbc_ket_cc = modelspace->GetTwoBodyChannel_CC(ch_ket_cc);
      int Jbra = tbc_bra_cc.J;
      int Jket = tbc_ket_cc.J;

//      arma::mat XJ1;
//      arma::mat XJ2;
      arma::mat YJ1J2;
      arma::mat YJ2J1;
      const auto& XJ1 = Xt_bar_ph[ch_bra_cc];
      const auto& XJ2 = Xt_bar_ph[ch_ket_cc];

      arma::uvec kets_ph = arma::join_cols( tbc_ket_cc.GetKetIndex_hh(), tbc_ket_cc.GetKetIndex_ph() );
      arma::uvec bras_ph = arma::join_cols( tbc_bra_cc.GetKetIndex_hh(), tbc_bra_cc.GetKetIndex_ph() );
//      profiler.timer["BuildZbarTensor_setup"] += omp_get_wtime() - t_start2;

//      t_start2 = omp_get_wtime();
      Y.DoTensorPandyaTransformation_SingleChannel(YJ1J2,ch_bra_cc,ch_ket_cc);
      if (ch_bra_cc==ch_ket_cc)
      {
//         XJ2 = XJ1;
         YJ2J1 = YJ1J2;
      }
      else
      {
         Y.DoTensorPandyaTransformation_SingleChannel(YJ2J1,ch_ket_cc,ch_bra_cc);
      }
//      profiler.timer["DoTensorPandyaTransformation_SingleChannel"] += omp_get_wtime() - t_start2;

//      t_start2 = omp_get_wtime();
      int flipphaseY = hY * modelspace->phase( Jbra - Jket ) ;
      // construct a matrix of phases (-1)^{k+j+p+h} used below to generate X_phkj for k>j
      arma::mat PhaseMatXJ2( tbc_ket_cc.GetNumberKets(), kets_ph.size(), arma::fill::ones) ;
      arma::mat PhaseMatYJ1J2( bras_ph.size(), tbc_ket_cc.GetNumberKets(), arma::fill::ones) ;
      for ( index_t iket=0;iket<(index_t)tbc_ket_cc.GetNumberKets();iket++)
      {
        const Ket& ket = tbc_ket_cc.GetKet(iket);
        if ( modelspace->phase((ket.op->j2+ket.oq->j2)/2)<0)
        {
           PhaseMatXJ2.row(iket) *=-1;
           PhaseMatYJ1J2.col(iket) *=-1;
        }
      }
      for (index_t iph=0;iph<kets_ph.size();iph++)
      {
        const Ket& ket_ph = tbc_ket_cc.GetKet( kets_ph[iph] );
        if ( modelspace->phase((ket_ph.op->j2+ket_ph.oq->j2)/2)<0)   PhaseMatXJ2.col(iph) *=-1;
      }
      for (index_t iph=0;iph<bras_ph.size();iph++)
      {
        const Ket& bra_ph = tbc_bra_cc.GetKet( bras_ph[iph] );
        if ( modelspace->phase((bra_ph.op->j2+bra_ph.oq->j2)/2)<0)   PhaseMatYJ1J2.row(iph) *=-1;
      }
      PhaseMatYJ1J2 *= flipphaseY;

//      profiler.timer["MakePhaseMat"] += omp_get_wtime() - t_start2;



//                J2                       J1         J2                       J2          J2
//             k<=j     k>=j                hp  -ph    hp   ph                 k<=j       k<=j
//      J1   [       |       ]       J1   [           |          ]         [hp        |ph        ]
//     i<=j  [  Zbar | Zbar  ]  =   i<=j  [   Xbar    | -Ybar    ]   * J1  [   Ybar   |   Ybar'  ]
//           [       |       ]            [           |          ]         [ph        |hp        ]      where Ybar'_phkj = Ybar_hpkj * (-1)^{p+h+k+j}*(-1)^{J1-J2}*hY
//                                                                         [----------|----------]       and
//                                                                     J2  [hp        |ph        ]            Xbar'_phkj = Xbar_hpkj * (-1)^{p+h+k+j}*hX
//                                                                         [   Xbar   |   Xbar'  ]
//                                                                         [-ph       |-hp       ]
//    
//
//      t_start2 = omp_get_wtime();
      int halfncx2 = XJ2.n_cols/2;
      int halfnry12 = YJ1J2.n_rows/2;
//      auto& Zmat = Z_bar[{ch_bra_cc,ch_ket_cc}];
      auto& Zmat = Z_bar.at({ch_bra_cc,ch_ket_cc});
//      arma::mat Zmat ;

      arma::mat Mleft = join_horiz( XJ1,  -flipphaseY * YJ2J1.t() );
      arma::mat Mright = join_vert( join_horiz( YJ1J2 ,  join_vert( YJ1J2.tail_rows(halfnry12)%PhaseMatYJ1J2 ,
                                                                    YJ1J2.head_rows(halfnry12)%PhaseMatYJ1J2  )   ), 
                                  hX*join_vert( XJ2,    join_horiz(   XJ2.tail_cols(halfncx2)%PhaseMatXJ2 ,
                                                                      XJ2.head_cols(halfncx2)%PhaseMatXJ2     )   ).t() );

      Zmat = Mleft * Mright;
//      profiler.timer["Multiply_Mleft_Mright"] += omp_get_wtime() - t_start2;

//      Z.AddInverseTensorPandyaTransformation_SingleChannel(Zmat,ch_bra_cc,ch_ket_cc); 
//      cout << "PANDYA LOOKUP size = " << plookup.size() << endl;
//      Z_debug1.AddInverseTensorPandyaTransformation_SingleChannel(Zmat,ch_bra_cc,ch_ket_cc); 

   }
   profiler.timer["Build Z_bar_tensor"] += omp_get_wtime() - t_start;

//   cout << "Done with parallel loop" << endl;

   t_start = omp_get_wtime();
//   cout << "Adding INversePanyaTransformation" << endl;
   Z.AddInverseTensorPandyaTransformation(Z_bar); // TODO: Do this one channel at a time <-- done did it, and it's sloooowwww.
//   cout << "done." << endl;

   profiler.timer["InverseTensorPandyaTransformation"] += omp_get_wtime() - t_start;

   modelspace->tensor_transform_first_pass.at( rank_J ) = false;

}










