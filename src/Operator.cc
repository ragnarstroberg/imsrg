
#include "Operator.hh"
#include "AngMom.hh"
#include "IMSRGProfiler.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <deque>
#include <array>
#include <gsl/gsl_math.h>
#include <math.h>
#include "omp.h"

#ifndef SQRT2
  #define SQRT2 1.4142135623730950488L
#endif

//using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

//double  Operator::bch_transform_threshold = 1e-6;
//double  Operator::bch_transform_threshold = 1e-9;
//double  Operator::bch_product_threshold = 1e-4;
//bool Operator::use_brueckner_bch = false;
//bool Operator::use_goose_tank_correction = false;
//bool Operator::use_goose_tank_correction_titus = false;

//IMSRGProfiler Operator::IMSRGProfiler::

//Operator& Operator::TempOp(size_t n)
//{
//  static deque<Operator> TempArray;
//  if (n >= TempArray.size()) TempArray.resize(n+1,*this);
//  return TempArray[n];
//}

//vector<arma::mat>& Operator::TempMatVec(size_t n)
//{
//  static deque<vector<arma::mat>> TempMatVecArray;
//  if (n>= TempMatVecArray.size()) TempMatVecArray.resize(max(n,(size_t)5));
//  return TempMatVecArray[n];
//}

//////////////////// DESTRUCTOR //////////////////////////////////////////
Operator::~Operator()
{
  IMSRGProfiler::counter["N_Operators"] --;
}

/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator()
 :   modelspace(NULL), 
    rank_J(0), rank_T(0), parity(0), particle_rank(2), legs(4), 
    hermitian(true), antihermitian(false), nChannels(0)
{
  IMSRGProfiler::counter["N_Operators"] ++;
}


// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank) : 
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms,Jrank,Trank,p),  ThreeBody(&ms),
    rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank), legs(2*part_rank),
    E3max(ms.GetE3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels()) , Q_space_orbit(-1)
{
  SetUpOneBodyChannels();
  if (particle_rank >=3) ThreeBody.Allocate();
  IMSRGProfiler::counter["N_Operators"] ++;
}

Operator::Operator(ModelSpace& ms) :
    modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(),arma::fill::zeros),
    TwoBody(&ms),  ThreeBody(&ms),
    rank_J(0), rank_T(0), parity(0), particle_rank(2), legs(4),
    E3max(ms.GetE3max()),
    hermitian(true), antihermitian(false),  
    nChannels(ms.GetNumberTwoBodyChannels()), Q_space_orbit(-1)
{
  SetUpOneBodyChannels();
  IMSRGProfiler::counter["N_Operators"] ++;
}

Operator::Operator(const Operator& op)
: modelspace(op.modelspace),  ZeroBody(op.ZeroBody),
  OneBody(op.OneBody), TwoBody(op.TwoBody) ,ThreeBody(op.ThreeBody),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank), legs(op.legs),
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels), Q_space_orbit(op.Q_space_orbit)
{
  IMSRGProfiler::counter["N_Operators"] ++;
}

Operator::Operator(Operator&& op)
: modelspace(op.modelspace), ZeroBody(op.ZeroBody),
  OneBody(std::move(op.OneBody)), TwoBody(std::move(op.TwoBody)) , ThreeBody(std::move(op.ThreeBody)),
  rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank), legs(op.legs), 
  E2max(op.E2max), E3max(op.E3max), 
  hermitian(op.hermitian), antihermitian(op.antihermitian),
  nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels), Q_space_orbit(op.Q_space_orbit)
{
  IMSRGProfiler::counter["N_Operators"] ++;
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
  for ( size_t i=0; i<modelspace->GetNumberOrbits(); ++i )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    // The +-1 comes from the spin [LxS](J)
    int lmin = std::max( oi.l - rank_J-1, 0);
    int lmax = std::min( oi.l + rank_J+1, modelspace->GetEmax() );
    for (int l=lmin; l<=lmax; l+=1)
    {
      if ((l + oi.l + parity)%2>0) continue;
      int j2min = std::max(std::max(oi.j2 - 2*rank_J, 2*l-1),1);
      int j2max = std::min(oi.j2 + 2*rank_J, 2*l+1);
      for (int j2=j2min; j2<=j2max; j2+=2)
      {
        int tz2min = std::max( oi.tz2 - 2*rank_T, -1);
        int tz2max = std::min( oi.tz2 + 2*rank_T, 1);
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

double Operator::GetTwoBody(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket)
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


void Operator::WriteBinary(std::ofstream& ofs)
{
  double tstart = omp_get_wtime();
  ofs.write((char*)&rank_J,sizeof(rank_J));
  ofs.write((char*)&rank_T,sizeof(rank_T));
  ofs.write((char*)&parity,sizeof(parity));
  ofs.write((char*)&particle_rank,sizeof(particle_rank));
  ofs.write((char*)&legs,sizeof(legs));
  ofs.write((char*)&E2max,sizeof(E2max));
  ofs.write((char*)&E3max,sizeof(E3max));
  ofs.write((char*)&hermitian,sizeof(hermitian));
  ofs.write((char*)&antihermitian,sizeof(antihermitian));
  ofs.write((char*)&nChannels,sizeof(nChannels));
  ofs.write((char*)&ZeroBody,sizeof(ZeroBody));
  ofs.write((char*)OneBody.memptr(),OneBody.size()*sizeof(double));
//  if (particle_rank > 1)
  if (legs > 3)
    TwoBody.WriteBinary(ofs);
//  if (particle_rank > 2)
  if (legs > 5)
    ThreeBody.WriteBinary(ofs);
  IMSRGProfiler::timer["Write Binary Op"] += omp_get_wtime() - tstart;
}


void Operator::ReadBinary(std::ifstream& ifs)
{
  double tstart = omp_get_wtime();
  ifs.read((char*)&rank_J,sizeof(rank_J));
  ifs.read((char*)&rank_T,sizeof(rank_T));
  ifs.read((char*)&parity,sizeof(parity));
  ifs.read((char*)&particle_rank,sizeof(particle_rank));
  ifs.read((char*)&legs,sizeof(legs));
  ifs.read((char*)&E2max,sizeof(E2max));
  ifs.read((char*)&E3max,sizeof(E3max));
  ifs.read((char*)&hermitian,sizeof(hermitian));
  ifs.read((char*)&antihermitian,sizeof(antihermitian));
  ifs.read((char*)&nChannels,sizeof(nChannels));
  SetUpOneBodyChannels();
  ifs.read((char*)&ZeroBody,sizeof(ZeroBody));
  ifs.read((char*)OneBody.memptr(),OneBody.size()*sizeof(double));
//  if (particle_rank > 1)
  if (legs > 3)
    TwoBody.ReadBinary(ifs);
//  if (particle_rank > 2)
  if (legs > 5)
    ThreeBody.ReadBinary(ifs);
  IMSRGProfiler::timer["Read Binary Op"] += omp_get_wtime() - tstart;
}





////////////////// MAIN INTERFACE METHODS //////////////////////////

Operator Operator::DoNormalOrdering() const
{
   if (legs%2>0)
      return DoNormalOrderingDagger(+1);
   if (legs>5)
      return DoNormalOrdering3(+1);
   else
      return DoNormalOrdering2(+1);
}

//*************************************************************
///  Normal ordering of a 2body operator
///  set up for scalar or tensor operators, but
///  the tensor part hasn't been tested
//*************************************************************
Operator Operator::DoNormalOrdering2(int sign) const
{
   Operator opNO(*this);
   bool scalar = (opNO.rank_J==0 and opNO.rank_T==0 and opNO.parity==0);
   if (scalar)
   {
     for (auto& k : modelspace->holes) // loop over hole orbits
     {
        Orbit& ok = modelspace->GetOrbit(k);
        opNO.ZeroBody += (ok.j2+1) * sign*ok.occ * OneBody(k,k);
     }
   }
//   std::cout << "OneBody contribution: " << opNO.ZeroBody << std::endl;

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
          // We have two occupations (na*nb), so if we're undoing the normal ordering the signs cancel out and we get no minus sign here.
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
                   opNO.OneBody(a,b) += hatfactor /(2*ja+1.0) * sign*oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
                }
                else
                {
                   double jh = oh.j2*0.5;
                   if ((ja+jh < J_bra) or (abs(ja-jh)>J_bra) or (jb+jh < J_ket) or (abs(jb-jh)>J_ket) ) continue;
                   if ((oa.l + oh.l + tbc_bra.parity)%2 >0) continue;
                   if ((ob.l + oh.l + tbc_ket.parity)%2 >0) continue;
                   if ((oa.tz2 + oh.tz2) != tbc_bra.Tz*2) continue;
                   if ((ob.tz2 + oh.tz2) != tbc_ket.Tz*2) continue;
                   double ME = hatfactor  * sign*oh.occ *modelspace->phase(ja+jh-J_ket-opNO.rank_J)
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
//     std::cout << "------------------------------------------" << std::endl;
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
//Operator Operator::DoNormalOrdering3()
Operator Operator::DoNormalOrdering3(int sign) const
{
   Operator opNO3 = Operator(*modelspace);
//   #pragma omp parallel for
   for ( auto& itmat : opNO3.TwoBody.MatEl )
   {
      int ch = itmat.first[0]; // assume ch_bra = ch_ket for 3body...
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& Gamma = (arma::mat&) itmat.second;
      for (size_t ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         Orbit & oi = modelspace->GetOrbit(i);
         Orbit & oj = modelspace->GetOrbit(j);
         for (size_t iket=ibra; iket<tbc.GetNumberKets(); ++iket)
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
                  Gamma(ibra,iket) += (K2+1) * sign*oa.occ * ThreeBody.GetME_pn(tbc.J,tbc.J,K2,i,j,a,k,l,a); // This is unnormalized, but it should be normalized!!!!
               }
            }
            Gamma(ibra,iket) /= (2*tbc.J+1)* sqrt((1+bra.delta_pq())*(1+ket.delta_pq()));
         }
      }
   }
   opNO3.Symmetrize();
   Operator opNO2 = opNO3.DoNormalOrdering2(sign);
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);

   // Also normal order the 1 and 2 body pieces
   opNO2 += DoNormalOrdering2(sign);
   return opNO2;

}


///  The normal ordering is slightly different if the operator is a
///  dagger operator.
///
Operator Operator::DoNormalOrderingDagger( int sign) const
{
  Operator opNO(*this);
 
  index_t Q = opNO.GetQSpaceOrbit();
  Orbit& oQ = modelspace->GetOrbit(Q);
//  double jQ = oQ.j2*0.5;

  for ( auto& itmat : TwoBody.MatEl )
  {
     int ch_bra = itmat.first[0];
     int ch_ket = itmat.first[1];
     
     TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch_bra);
     int J = tbc.J;
     double hatfactor = 2*J+1.0;

     // One body part

     for ( auto a : OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
     {
        Orbit &oa = modelspace->GetOrbit(a);
        double ja = oa.j2*0.5;
        for (auto& h : modelspace->holes)  // C++11 syntax
        {
          Orbit& oh = modelspace->GetOrbit(h);

          if (opNO.rank_J==0)
          {
             opNO.OneBody(a,Q) -= hatfactor/(2*ja+1) * sign*oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,h,a,h,Q);
          }
        }
     }
  } // loop over channels

  return opNO;

}




Operator Operator::UndoNormalOrdering() const
{
  if (legs%2>0)
    return UndoNormalOrderingDagger();
  else if (legs < 5)
    return UndoNormalOrdering2();
  else
  {
    return UndoNormalOrdering3();
//    std::cout << "WARNING: calling Operator::UndoNormalOrdering on a 3-body operator. Not yet implemented." << std::endl;
//    return UndoNormalOrdering2();
  }
}

/// Convert to a basis normal ordered wrt the vacuum.
/// This doesn't handle 3-body terms. In that case,
/// the 2-body piece is unchanged.
//Operator Operator::UndoNormalOrdering2() const
//{
//   Operator opNO = *this;
////   std::cout << "Undoing Normal ordering. Initial ZeroBody = " << opNO.ZeroBody << std::endl;
//
//   if (opNO.GetJRank()==0 and opNO.GetTRank()==0 and opNO.GetParity()==0)
//   {
//     for (auto& k : modelspace->holes) // loop over hole orbits
//     {
//        Orbit& ok = modelspace->GetOrbit(k);
//        opNO.ZeroBody -= (ok.j2+1) * ok.occ * OneBody(k,k);
//     }
//   }
//
//   index_t norbits = modelspace->GetNumberOrbits();
//
//   for ( auto& itmat : TwoBody.MatEl )
//   {
//      int ch_bra = itmat.first[0];
//      int ch_ket = itmat.first[1];
//      auto& matrix = itmat.second;
//      
//      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
//      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
//      int J_bra = tbc_bra.J;
//      int J_ket = tbc_ket.J;
//      double hatfactor = sqrt((2*J_bra+1.0)*(2*J_ket+1.0));
//
//      // Zero body part
//      if (opNO.GetJRank()==0 and opNO.GetTRank()==0 and opNO.GetParity()==0)
//      {
//        arma::vec diagonals = matrix.diag();
//        auto hh = tbc_ket.GetKetIndex_hh();
//        auto hocc = tbc_ket.Ket_occ_hh;
//        opNO.ZeroBody +=  arma::sum( hocc % diagonals.elem(hh) ) *hatfactor;
//      }
//
//
//      // One body part
//      for (index_t a=0;a<norbits;++a)
//      {
//         Orbit &oa = modelspace->GetOrbit(a);
//         double ja = oa.j2*0.5;
////         index_t bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
//           index_t bstart = (IsNonHermitian() or ch_bra!=ch_ket )? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
//         for ( auto& b : opNO.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) ) 
//         {
//            if (b < bstart) continue;
//            Orbit &ob = modelspace->GetOrbit(b);
//            double jb = ob.j2*0.5;
//            for (auto& h : modelspace->holes)  // C++11 syntax
//            {
//              Orbit& oh = modelspace->GetOrbit(h);
//
//              if (opNO.rank_J==0)
//              {
//                 opNO.OneBody(a,b) -= hatfactor/(2*ja+1) * oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
//              }
//              else
//              {
//                 double jh = oh.j2*0.5;
//                 if ((ja+jh < J_bra) or (abs(ja-jh)>J_bra) or (jb+jh < J_ket) or (abs(jb-jh)>J_ket) ) continue;
//
//                 if ((oa.l + oh.l + tbc_bra.parity)%2 >0) continue;
//                 if ((ob.l + oh.l + tbc_ket.parity)%2 >0) continue;
//                 if ((oa.tz2 + oh.tz2) != tbc_bra.Tz*2) continue;
//                 if ((ob.tz2 + oh.tz2) != tbc_ket.Tz*2) continue;
//                 double ME = hatfactor  * oh.occ *modelspace->phase(ja+jh-J_ket-opNO.rank_J)
//                                             * modelspace->GetSixJ(J_bra,J_ket,opNO.rank_J,jb,ja,jh) * TwoBody.GetTBME(ch_bra,ch_ket,a,h,b,h);
//                 if (a>b)
//                 {
//                   int herm = IsHermitian() ? 1 : -1;
//                   opNO.OneBody(b,a) -= herm * modelspace->phase(ja-jb) * ME;
//                 }
//                 else
//                 {
//                   opNO.OneBody(a,b) -= ME;
//                 }
//
//              }
//            }
//
//         }
//      }
//   } // loop over channels
//
//   if (hermitian) opNO.Symmetrize();
//   if (antihermitian) opNO.AntiSymmetrize();
//
//   return opNO;
//
//}


//Operator Operator::UndoNormalOrderingDagger() const
//{
//   Operator opNO(*this);
//
////   index_t norbits = modelspace->GetNumberOrbits();
//
//   index_t Q = opNO.GetQSpaceOrbit();
//   Orbit &oQ = modelspace->GetOrbit(Q);
////   double jQ = oQ.j2*0.5;
//
//   for ( auto& itmat : TwoBody.MatEl )
//   {
//      int ch_bra = itmat.first[0];
//      int ch_ket = itmat.first[1];
////      auto& matrix = itmat.second;
//      
//      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch_bra);
////      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
////      int J_bra = tbc_bra.J;
//      int J = tbc.J;
//      double hatfactor = 2*J+1.0;
//
//      // One body part
//
////      for (index_t a=0;a<norbits;++a)
//      for ( auto a : OneBodyChannels.at({oQ.l,oQ.j2,oQ.tz2}) )
//      {
//         Orbit &oa = modelspace->GetOrbit(a);
//         double ja = oa.j2*0.5;
//            for (auto& h : modelspace->holes)  // C++11 syntax
//            {
//              Orbit& oh = modelspace->GetOrbit(h);
//
//              if (opNO.rank_J==0)
//              {
//                 opNO.OneBody(a,Q) += hatfactor/(2*ja+1) * oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,h,a,h,Q);
//              }
//            }
//      }
//   } // loop over channels
//
//
//  return opNO;
//
//}



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
    std::cout << "Error: Cannot truncate an operator with emax = " << modelspace->GetEmax() << " to one with emax = " << new_emax << std::endl;
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



ModelSpace* Operator::GetModelSpace() const
{
   return modelspace;
}


void Operator::Erase()
{
  EraseZeroBody();
  EraseOneBody();
  TwoBody.Erase();
//  if (particle_rank >=3)
  if (legs >=6)
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
    std::cout << "Trying to reduce an operator with J rank = " << rank_J << ". Not good!!!" << std::endl;
    return;
  }
  for ( size_t a=0;a<modelspace->GetNumberOrbits();++a )
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for (size_t b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
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
    std::cout << "Trying to un-reduce an operator with J rank = " << rank_J << ". Not good!!!" << std::endl;
    return;
  }
  for ( size_t a=0;a<modelspace->GetNumberOrbits();++a )
  {
    Orbit& oa = modelspace->GetOrbit(a);
    for ( size_t b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
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
//void Operator::Eye()
//{
//   ZeroBody = 1;
//   OneBody.eye();
//   TwoBody.Eye();
//}


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
       if (abs(OneBody(i,a))>1e-6)
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
           int Jmin = std::max(std::abs(oi.j2-oj.j2),std::abs(oa.j2-ob.j2))/2;
           int Jmax = std::min(oi.j2+oj.j2,oa.j2+ob.j2)/2;
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
   IMSRGProfiler::timer["GetMP2_Energy"] += omp_get_wtime() - t_start;
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
//double Operator::GetMP3_Energy()
std::array<double,3> Operator::GetMP3_Energy()
{
   // So far, the pp and hh parts seem to work. No such luck for the ph.
   double t_start = omp_get_wtime();
//   double Emp3 = 0;
   double Epp = 0;
   double Ehh = 0;
   double Eph = 0;
   // This can certainly be optimized, but I'll wait until this is the bottleneck.
   int nch = modelspace->GetNumberTwoBodyChannels();

//   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Emp3)
   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Epp,Ehh)
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
//           Emp3 += (2*J+1)*Mat(iket_ab,iket_ij) * Mat(iket_ij,iket_cd) * Mat(iket_cd,iket_ab) / (Delta_abij * Delta_cdij);
           Ehh += (2*J+1)*Mat(iket_ab,iket_ij) * Mat(iket_ij,iket_cd) * Mat(iket_cd,iket_ab) / (Delta_abij * Delta_cdij);
         }

       // pp term
         for (auto iket_kl : tbc.GetKetIndex_pp() )
         {
           Ket& ket_kl = tbc.GetKet(iket_kl);
           index_t k = ket_kl.p;
           index_t l = ket_kl.q;
           double Delta_abkl = OneBody(a,a) + OneBody(b,b) - OneBody(k,k) - OneBody(l,l);
//           Emp3 += (2*J+1)*Mat(iket_ab,iket_ij)*Mat(iket_ij,iket_kl)*Mat(iket_kl,iket_ab) / (Delta_abij * Delta_abkl);
           Epp += (2*J+1)*Mat(iket_ab,iket_ij)*Mat(iket_ij,iket_kl)*Mat(iket_kl,iket_ab) / (Delta_abij * Delta_abkl);
         }

       } // for ij
     } // for ab
   } // for ich
//   cout << "done with pp and hh. E(3) = " << Emp3 << endl;




   index_t nparticles = modelspace->particles.size();
   modelspace->PreCalculateSixJ();
//   #pragma omp parallel for schedule(dynamic,1)  reduction(+:Emp3)
   #pragma omp parallel for schedule(dynamic,1)  reduction(+:Eph)
   for (index_t ii=0;ii<nparticles;ii++)
   {
     auto i = modelspace->particles[ii];
     double ji = 0.5*modelspace->GetOrbit(i).j2;
     for (auto a : modelspace->holes)
     {
      double ja = 0.5*modelspace->GetOrbit(a).j2;
      int J_min = abs(ja-ji);
      int J_max = ja+ji;
      for (int J_tot=J_min;J_tot<=J_max;++J_tot)
      {
       double Jfactor = (2*J_tot + 1)*(2*J_tot + 1) ; // I don't yet understand why it's (2J+1)**2, but this is what came from Johannes.
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
//             Emp3 +=  Jfactor * tbme_abij * tbme_cjkb * tbme_ikac / (Delta_abij * Delta_acik);
             Eph +=  Jfactor * tbme_abij * tbme_cjkb * tbme_ikac / (Delta_abij * Delta_acik);
            } // for k
          } // for c
        } // for j
       } // for b
      } // for J_tot
     } // for a
   } // for i


   IMSRGProfiler::timer["GetMP3_Energy"] += omp_get_wtime() - t_start;
//   return Emp3;
   return {Epp,Ehh,Eph};
}




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
   for ( size_t p=0; p<modelspace->GetNumberOrbits(); ++p)
   {
     Orbit& op = modelspace->GetOrbit(p);
     for ( auto q : OneBodyChannels.at({op.l, op.j2, op.tz2}) )
     {
       Orbit& oq = modelspace->GetOrbit(q);
       int degeneracy_factor = (op.j2+1) * ( (std::min(oq.j2,op.j2+rank_J)-std::max(-oq.j2,op.j2-rank_J))/2+1 );
       nrm += OneBody(p,q)*OneBody(p,q) * degeneracy_factor * degeneracy_factor;
     }
   }
   return sqrt(nrm);
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

  for ( size_t ch=0; ch<modelspace->GetNumberTwoBodyChannels(); ++ch )
  {
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    double trVijij = arma::trace( OpVac.TwoBody.GetMatrix(ch,ch) );
    switch ( tbc.Tz )
    {
      case ( -1 ) :  trace += (tbc.J*2+1)* trVijij * norm_pp; break;
      case (  0 ) :  trace += (tbc.J*2+1)* trVijij * norm_pn; break;
      case (  1 ) :  trace += (tbc.J*2+1)* trVijij * norm_nn; break;
      default: std::cout << "AAAHHH blew the switch statement. tbc.Tz = " << tbc.Tz << std::endl;
    }
  }
  IMSRGProfiler::timer["Operator::Trace"] += omp_get_wtime() - t_start;
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
    for (size_t ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
    {
      Ket& bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      double ei = OneBody(i,i);
      double ej = OneBody(j,j);
      for (size_t iket=0; iket<tbc.GetNumberKets(); ++iket)
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





