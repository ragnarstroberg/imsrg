
#include "Operator.hh"
#include "AngMom.hh"
#include "IMSRGProfiler.hh"
#include "PhysicalConstants.hh" // for SQRT2
#include <cmath>
#include <iostream>
#include <iomanip>
#include <deque>
#include <array>
#include <gsl/gsl_math.h>
#include <math.h>
#include "omp.h"

// #ifndef SQRT2
//   #define SQRT2 1.4142135623730950488L
// #endif

// using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

// double  Operator::bch_transform_threshold = 1e-6;
// double  Operator::bch_transform_threshold = 1e-9;
// double  Operator::bch_product_threshold = 1e-4;
// bool Operator::use_brueckner_bch = false;
// bool Operator::use_goose_tank_correction = false;
// bool Operator::use_goose_tank_correction_titus = false;

// IMSRGProfiler Operator::IMSRGProfiler::

// Operator& Operator::TempOp(size_t n)
//{
//   static deque<Operator> TempArray;
//   if (n >= TempArray.size()) TempArray.resize(n+1,*this);
//   return TempArray[n];
// }

// vector<arma::mat>& Operator::TempMatVec(size_t n)
//{
//   static deque<vector<arma::mat>> TempMatVecArray;
//   if (n>= TempMatVecArray.size()) TempMatVecArray.resize(max(n,(size_t)5));
//   return TempMatVecArray[n];
// }

//////////////////// DESTRUCTOR //////////////////////////////////////////
Operator::~Operator()
{
  IMSRGProfiler::counter["N_Operators"]--;
}

/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator::Operator()
    : modelspace(NULL),
      rank_J(0), rank_T(0), parity(0), particle_rank(2), legs(4),
      hermitian(true), antihermitian(false), nChannels(0), is_reduced(false)
{
  IMSRGProfiler::counter["N_Operators"]++;
}

// Create a zero-valued operator in a given model space
Operator::Operator(ModelSpace &ms, int Jrank, int Trank, int p, int part_rank) : modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(), arma::fill::zeros),
                                                                                 //    TwoBody(&ms,Jrank,Trank,p),  ThreeBody(&ms,Jrank,Trank,p), ThreeLeg(&ms), ThreeBodyNO2B(),
                                                                                 //    TwoBody(&ms,Jrank,Trank,p),  ThreeBody(&ms,Jrank,Trank,p), ThreeLeg(&ms),
                                                                                 TwoBody(), ThreeBody(&ms, Jrank, Trank, p), ThreeLeg(&ms),
                                                                                 rank_J(Jrank), rank_T(Trank), parity(p), particle_rank(part_rank), legs(2 * part_rank),
                                                                                 E3max(ms.GetE3max()),
                                                                                 hermitian(true), antihermitian(false),
                                                                                 nChannels(ms.GetNumberTwoBodyChannels()), Q_space_orbit(-1),
                                                                                 is_reduced(Jrank > 0 or Trank > 0 or p > 0) // by default, Hamiltonian-like operators are not reduced, all others are reduced.
{
  if (part_rank >= 2)
  {
    TwoBody = TwoBodyME(&ms, Jrank, Trank, p);
  }
  SetUpOneBodyChannels();
  //  if (particle_rank >=3) ThreeBody.Allocate();   // Don't allocate automatically. Wait until we're sure we want it.
  IMSRGProfiler::counter["N_Operators"]++;
}

Operator::Operator(ModelSpace &ms) : modelspace(&ms), ZeroBody(0), OneBody(ms.GetNumberOrbits(), ms.GetNumberOrbits(), arma::fill::zeros),
                                     TwoBody(&ms), ThreeBody(&ms), ThreeLeg(&ms),
                                     //    TwoBody(&ms),  ThreeBody(&ms), ThreeLeg(&ms), ThreeBodyNO2B(),
                                     rank_J(0), rank_T(0), parity(0), particle_rank(2), legs(4),
                                     E3max(ms.GetE3max()),
                                     hermitian(true), antihermitian(false),
                                     nChannels(ms.GetNumberTwoBodyChannels()), Q_space_orbit(-1),
                                     is_reduced(false)
{
  SetUpOneBodyChannels();
  IMSRGProfiler::counter["N_Operators"]++;
}

Operator::Operator(const Operator &op)
    : modelspace(op.modelspace), ZeroBody(op.ZeroBody),
      OneBody(op.OneBody), TwoBody(op.TwoBody), ThreeBody(op.ThreeBody), ThreeLeg(op.ThreeLeg),
      //  OneBody(op.OneBody), TwoBody(op.TwoBody) ,ThreeBody(op.ThreeBody), ThreeLeg(op.ThreeLeg), ThreeBodyNO2B(op.ThreeBodyNO2B),
      rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank), legs(op.legs),
      E2max(op.E2max), E3max(op.E3max),
      hermitian(op.hermitian), antihermitian(op.antihermitian),
      nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels), Q_space_orbit(op.Q_space_orbit),
      OneBodyChannels_vec(op.OneBodyChannels_vec),
      is_reduced(op.is_reduced)
{
  IMSRGProfiler::counter["N_Operators"]++;
}

Operator::Operator(Operator &&op)
    : modelspace(op.modelspace), ZeroBody(op.ZeroBody),
      OneBody(std::move(op.OneBody)), TwoBody(std::move(op.TwoBody)), ThreeBody(std::move(op.ThreeBody)), ThreeLeg(std::move(op.ThreeLeg)),
      //  ThreeBodyNO2B(std::move(op.ThreeBodyNO2B)),
      rank_J(op.rank_J), rank_T(op.rank_T), parity(op.parity), particle_rank(op.particle_rank), legs(op.legs),
      E2max(op.E2max), E3max(op.E3max),
      hermitian(op.hermitian), antihermitian(op.antihermitian),
      nChannels(op.nChannels), OneBodyChannels(op.OneBodyChannels), Q_space_orbit(op.Q_space_orbit),
      OneBodyChannels_vec(op.OneBodyChannels_vec),
      is_reduced(op.is_reduced)
{
  IMSRGProfiler::counter["N_Operators"]++;
}

/////////////// OVERLOADED OPERATORS =,+,-,*,etc ////////////////////

Operator &Operator::operator=(const Operator &rhs) = default;
Operator &Operator::operator=(Operator &&rhs) = default;

// multiply operator by a scalar
Operator &Operator::operator*=(const double rhs)
{
  ZeroBody *= rhs;
  OneBody *= rhs;
  TwoBody *= rhs;
  ThreeLeg *= rhs;
  if (particle_rank > 2)
    ThreeBody *= rhs;
  //   if (particle_rank > 2)  ThreeBodyNO2B *= rhs;
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
Operator operator*(const double lhs, const Operator &rhs)
{
  return rhs * lhs;
}
Operator operator*(const double lhs, const Operator &&rhs)
{
  return rhs * lhs;
}

// divide operator by a scalar
Operator &Operator::operator/=(const double rhs)
{
  return *this *= (1.0 / rhs);
}

Operator Operator::operator/(const double rhs) const
{
  Operator opout = Operator(*this);
  opout *= (1.0 / rhs);
  return opout;
}

// Add operators
Operator &Operator::operator+=(const Operator &rhs)
{
  ZeroBody += rhs.ZeroBody;
  OneBody += rhs.OneBody;
  if (rhs.GetParticleRank() > 1)
    TwoBody += rhs.TwoBody;
  if (rhs.GetParticleRank() > 2)
    ThreeBody += rhs.ThreeBody;
  //   if (rhs.GetParticleRank() >2 )
  //     ThreeBodyNO2B += rhs.ThreeBodyNO2B;
  if (rhs.GetNumberLegs() % 2 == 1)
    ThreeLeg += rhs.ThreeLeg;
  return *this;
}

Operator Operator::operator+(const Operator &rhs) const
{
  if (GetParticleRank() >= rhs.GetParticleRank())
    return (Operator(*this) += rhs);
  else
    return (Operator(rhs) += *this);
}

Operator &Operator::operator+=(const double &rhs)
{
  ZeroBody += rhs;
  return *this;
}

Operator Operator::operator+(const double &rhs) const
{
  return (Operator(*this) += rhs);
}

// Subtract operators
Operator &Operator::operator-=(const Operator &rhs)
{
  ZeroBody -= rhs.ZeroBody;
  OneBody -= rhs.OneBody;
  if (rhs.GetParticleRank() > 1)
    TwoBody -= rhs.TwoBody;
  if (rhs.GetParticleRank() > 2)
    ThreeBody -= rhs.ThreeBody;
  //   if (rhs.GetParticleRank() > 2)
  //     ThreeBodyNO2B -= rhs.ThreeBodyNO2B;
  if (rhs.GetNumberLegs() % 2 == 1)
    ThreeLeg -= rhs.ThreeLeg;
  return *this;
}

Operator Operator::operator-(const Operator &rhs) const
{
  return (Operator(*this) -= rhs);
}

Operator &Operator::operator-=(const double &rhs)
{
  ZeroBody -= rhs;
  return *this;
}

Operator Operator::operator-(const double &rhs) const
{
  return (Operator(*this) -= rhs);
}

// Negation operator
Operator Operator::operator-() const
{
  return (*this) * -1.0;
}

void Operator::SetUpOneBodyChannels()
{
  for ( auto i : modelspace->all_orbits )
  {
    Orbit& oi = modelspace->GetOrbit(i);
    if ( OneBodyChannels.find( {oi.l,oi.j2,oi.tz2} ) == OneBodyChannels.end() ) OneBodyChannels[{oi.l,oi.j2,oi.tz2}] = {};
    // The +-1 comes from the spin [LxS](J)
    int lmin = std::max( oi.l - rank_J-1, 0);
    int lmax = std::min( oi.l + rank_J+1, modelspace->GetEmax() );
    for (int l=lmin; l<=lmax; l++)
    {
      if ((l + oi.l + parity)%2>0) continue;

      int j2min = std::max( std::abs(oi.j2 - 2*rank_J), 2*l-1);
      int j2max = std::min(          oi.j2 + 2*rank_J,  2*l+1);

      for (int j2=j2min; j2<=j2max; j2+=2)
      {
        for ( int tz2=-1; tz2<=1; tz2+=2)
        {
          if (std::abs(oi.tz2-tz2) == 2*rank_T)
          {
              OneBodyChannels[ {l, j2, tz2} ].insert(i);
              //std::cout << oi.n << "  " << oi.l << "  " << oi.j2 << "  " << oi.tz2 << "    ||    " << l << "  " << j2 << "  " << tz2 << std::endl;
          }
        }
      }
    }
  }


  OneBodyChannels_vec.resize(4 * (modelspace->Emax + 1) + 1, {});
  for (auto &it : OneBodyChannels)
  {
    int l = it.first[0];
    int twoj = it.first[1];
    int twotz = it.first[2];
    size_t indx = l * 4 + (twoj + 1 - 2 * l) + (twotz + 1) / 2;
    OneBodyChannels_vec[indx] = it.second;
  }
}

// l runs from 0 to emax, tz is -1,1, so (tz+1)/2 runs 0 to 1
// twoj is 2l-1, 2l+1,  so we can use (twoj+1-2*l) 0,2
// twoj runs from 1 to 2emax+1, so (twoj-1)/2 runs from 0 to emax
// std::set<index_t>& Operator::GetOneBodyChannel( int l, int twoj, int twotz )
const std::set<index_t> &Operator::GetOneBodyChannel(int l, int twoj, int twotz) const
{
  size_t indx = l * 4 + (twoj + 1 - 2 * l) + (twotz + 1) / 2;
  //  if (indx >= OneBodyChannels_vec.size() )
  //  {
  //     std::cout << "AHH ljt = " << l << " " << twoj << " " << twotz << "   -> index = " << indx << "  but size of OneBodyChannels_vec = " << OneBodyChannels_vec.size() << std::endl;
  //  }

  return OneBodyChannels_vec[indx];
}

size_t Operator::Size()
{
  return sizeof(ZeroBody) + OneBody.size() * sizeof(double) + TwoBody.size() + ThreeBody.size();
}

void Operator::SetOneBody(int i, int j, double val)
{
  OneBody(i, j) = val;
  if (IsNonHermitian())
    return;
  int flip_phase = IsHermitian() ? 1 : -1;
  if (rank_J > 0)
    flip_phase *= modelspace->phase((modelspace->GetOrbit(i).j2 - modelspace->GetOrbit(j).j2) / 2);
  OneBody(j, i) = flip_phase * val;
}

void Operator::SetTwoBody(int J1, int p1, int T1, int J2, int p2, int T2, int i, int j, int k, int l, double v)
{
  TwoBody.SetTBME(J1, p1, T1, J2, p2, T2, i, j, k, l, v);
}

double Operator::GetTwoBody(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket)
{
  if (ch_bra <= ch_ket or IsNonHermitian())
  {
    return TwoBody.GetMatrix(ch_bra, ch_ket)(ibra, iket);
  }
  else
  {
    TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int flipphase = modelspace->phase(tbc_bra.J - tbc_ket.J) * (IsHermitian() ? 1 : -1);
    return flipphase * TwoBody.GetMatrix(ch_ket, ch_bra)(iket, ibra);
  }
}

void Operator::WriteBinary(std::ofstream &ofs)
{
  double tstart = omp_get_wtime();
  ofs.write((char *)&rank_J, sizeof(rank_J));
  ofs.write((char *)&rank_T, sizeof(rank_T));
  ofs.write((char *)&parity, sizeof(parity));
  ofs.write((char *)&particle_rank, sizeof(particle_rank));
  ofs.write((char *)&legs, sizeof(legs));
  ofs.write((char *)&E2max, sizeof(E2max));
  ofs.write((char *)&E3max, sizeof(E3max));
  ofs.write((char *)&hermitian, sizeof(hermitian));
  ofs.write((char *)&antihermitian, sizeof(antihermitian));
  ofs.write((char *)&nChannels, sizeof(nChannels));
  ofs.write((char *)&ZeroBody, sizeof(ZeroBody));
  ofs.write((char *)OneBody.memptr(), OneBody.size() * sizeof(double));
  //  if (particle_rank > 1)
  if (legs > 3)
    TwoBody.WriteBinary(ofs);
  //  if (particle_rank > 2)
  if (legs > 5)
    ThreeBody.WriteBinary(ofs);
  IMSRGProfiler::timer["Write Binary Op"] += omp_get_wtime() - tstart;
}

void Operator::ReadBinary(std::ifstream &ifs)
{
  double tstart = omp_get_wtime();
  ifs.read((char *)&rank_J, sizeof(rank_J));
  ifs.read((char *)&rank_T, sizeof(rank_T));
  ifs.read((char *)&parity, sizeof(parity));
  ifs.read((char *)&particle_rank, sizeof(particle_rank));
  ifs.read((char *)&legs, sizeof(legs));
  ifs.read((char *)&E2max, sizeof(E2max));
  ifs.read((char *)&E3max, sizeof(E3max));
  ifs.read((char *)&hermitian, sizeof(hermitian));
  ifs.read((char *)&antihermitian, sizeof(antihermitian));
  ifs.read((char *)&nChannels, sizeof(nChannels));
  SetUpOneBodyChannels();
  ifs.read((char *)&ZeroBody, sizeof(ZeroBody));
  ifs.read((char *)OneBody.memptr(), OneBody.size() * sizeof(double));
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
  if (legs % 2 > 0)
    return DoNormalOrderingDagger(+1, modelspace->holes);
  if (legs > 5)
    return DoNormalOrdering3(+1, modelspace->holes);
  else
    return DoNormalOrdering2(+1, modelspace->holes);
}

Operator Operator::UndoNormalOrdering() const
{
  std::cout << " IN " << __func__ << "   legs = " << legs << std::endl;
  if (legs % 2 > 0)
    return DoNormalOrderingDagger(-1, modelspace->holes);
  //    return UndoNormalOrderingDagger();
  else if (legs < 5)
    return DoNormalOrdering2(-1, modelspace->holes);
  //    return UndoNormalOrdering2();
  else
  {
    return DoNormalOrdering3(-1, modelspace->holes);
    //    return UndoNormalOrdering3();
  }
}

// Operator Operator::UndoNormalOrdering2() const
//{
//    return this->DoNormalOrdering2(-1, modelspace->holes);
// }

// Operator Operator::UndoNormalOrdering3() const
//{
//    return this->DoNormalOrdering3(-1, modelspace->holes);
// }
// Operator Operator::UndoNormalOrderingDagger() const
//{
//    return this->DoNormalOrderingDagger(-1, modelspace->holes);
// }

Operator Operator::DoNormalOrderingCore() const
{
  std::cout << " IN " << __func__ << "   legs = " << legs << std::endl;
  if (legs % 2 > 0)
    return DoNormalOrderingDagger(+1, modelspace->core);
  if (legs > 5)
    return DoNormalOrdering3(+1, modelspace->core);
  else
    return DoNormalOrdering2(+1, modelspace->core);
}


Operator Operator::DoNormalOrderingFilledValence() const
{
  std::cout << " IN " << __func__ << "   legs = " << legs << std::endl;
   if (legs%2>0)
      return DoNormalOrderingDagger(+1, modelspace->valence);
   if (legs>5)
      return DoNormalOrdering3(+1, modelspace->valence);
   else
      return DoNormalOrdering2(+1, modelspace->valence);
}

//*************************************************************
///  Normal ordering of a 2body operator
///  set up for scalar or tensor operators, but
///  the tensor part hasn't been tested
//*************************************************************
// Operator Operator::DoNormalOrdering2(int sign) const
Operator Operator::DoNormalOrdering2(int sign, std::set<index_t> occupied) const
{
  //   for ( auto o : occupied ) std::cout << o << " ";
  //   std::cout << std::endl;

  Operator opNO(*this);
  bool scalar = (opNO.rank_J == 0 and opNO.rank_T == 0 and opNO.parity == 0);
  if (scalar)
  {
    //     for (auto& k : modelspace->holes) // loop over hole orbits
    for (auto &k : occupied) // loop over hole orbits
    {
      Orbit &ok = modelspace->GetOrbit(k);
      opNO.ZeroBody += (ok.j2 + 1) * sign * ok.occ * OneBody(k, k);

      if (particle_rank > 1)
      {
        for (auto &l : occupied)
        {
          if (l < k)
            continue;
          Orbit &ol = modelspace->GetOrbit(l);
          int Jmin = std::abs(ok.j2 - ol.j2) / 2;
          int Jmax = (ok.j2 + ol.j2) / 2;
          for (int J = Jmin; J <= Jmax; J++)
          {
            opNO.ZeroBody += (2 * J + 1) * ok.occ * ol.occ * TwoBody.GetTBME_J_norm(J, J, k, l, k, l);
          }
        }
      }
    }
  }
  //   std::cout << "OneBody contribution: " << opNO.ZeroBody << std::endl;

  index_t norbits = modelspace->GetNumberOrbits();
  if (TwoBody.Norm() > 1e-7)
  {
    for (auto &itmat : TwoBody.MatEl)
    {
      int ch_bra = itmat.first[0];
      int ch_ket = itmat.first[1];
      //        auto& matrix = itmat.second;

      TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
      int J_bra = tbc_bra.J;
      int J_ket = tbc_ket.J;
      double hatfactor = sqrt((2 * J_bra + 1.0) * (2 * J_ket + 1.0));

      // One body part
      for (index_t a = 0; a < norbits; ++a)
      {
        Orbit &oa = modelspace->GetOrbit(a);
        double ja = oa.j2 / 2.0;
        index_t bstart = (IsNonHermitian() or ch_bra != ch_ket) ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
        for (auto &b : opNO.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
        {
          if (b < bstart)
            continue;
          Orbit &ob = modelspace->GetOrbit(b);
          double jb = ob.j2 / 2.0;
          //              for (auto& h : modelspace->holes)  // C++11 syntax
          for (auto &h : occupied) // C++11 syntax
          {
            Orbit &oh = modelspace->GetOrbit(h);
            if (opNO.rank_J == 0)
            {
              opNO.OneBody(a, b) += hatfactor / (2 * ja + 1.0) * sign * oh.occ * TwoBody.GetTBME(ch_bra, ch_ket, a, h, b, h);
            }
            else
            {
              double jh = oh.j2 * 0.5;
              if ((ja + jh < J_bra) or (abs(ja - jh) > J_bra) or (jb + jh < J_ket) or (abs(jb - jh) > J_ket))
                continue;
              if ((oa.l + oh.l + tbc_bra.parity) % 2 > 0)
                continue;
              if ((ob.l + oh.l + tbc_ket.parity) % 2 > 0)
                continue;
              if ((oa.tz2 + oh.tz2) != tbc_bra.Tz * 2)
                continue;
              if ((ob.tz2 + oh.tz2) != tbc_ket.Tz * 2)
                continue;
              double ME = hatfactor * sign * oh.occ * modelspace->phase(ja + jh - J_ket - opNO.rank_J) * modelspace->GetSixJ(J_bra, J_ket, opNO.rank_J, jb, ja, jh) * TwoBody.GetTBME(ch_bra, ch_ket, a, h, b, h);
              if (a > b)
              {
                int herm = IsHermitian() ? 1 : -1;
                opNO.OneBody(b, a) += herm * modelspace->phase(ja - jb) * ME;
              }
              else
              {
                opNO.OneBody(a, b) += ME;
              }
            }
          }
        }
      }
    } // loop over channels
    //     std::cout << "------------------------------------------" << std::endl;
  }

  if (hermitian)
    opNO.Symmetrize();
  if (antihermitian)
    opNO.AntiSymmetrize();

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
// Operator Operator::DoNormalOrdering3()
Operator Operator::DoNormalOrdering3(int sign, std::set<index_t> occupied) const
{
  double t_start = omp_get_wtime();
  std::cout << "begin " << __func__ << "   norm of 3b is " << ThreeBodyNorm() << std::endl;
  //   Operator opNO3 = Operator(*modelspace);
  if (rank_J > 0)
  {
    std::cout << " Uh oh. Trying to call " << __func__ << "  on an operator with rank_J = " << rank_J << "   you should probably implement that first..." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  //    double vread = ThreeBody.GetME_pn(0,0,3,10,10,3,11,11,3);
  //    std::cout << " IN " << __func__ << "   vread =  " << vread << std::endl;
  Operator opNO3 = Operator(*modelspace, rank_J, rank_T, parity, 2);
  std::vector<int> ch_bra_list, ch_ket_list;
  //   std::vector<arma::mat *> mat_ptr_list;
  for (auto &itmat : opNO3.TwoBody.MatEl)
  {
    ch_bra_list.push_back(itmat.first[0]);
    ch_ket_list.push_back(itmat.first[1]);
  }
  int niter = ch_bra_list.size();
  //   for ( auto& itmat : opNO3.TwoBody.MatEl )
#pragma omp parallel for schedule(dynamic, 1)
  for (int iter = 0; iter < niter; iter++)
  {
    int ch_bra = ch_bra_list[iter];
    int ch_ket = ch_ket_list[iter];
    //      auto& Gamma = *(mat_ptr_list[iter]);
    auto &Gamma = opNO3.TwoBody.GetMatrix(ch_bra, ch_ket);
    //      int ch_bra = itmat.first[0]; // assume ch_bra = ch_ket for 3body...
    //      int ch_ket = itmat.first[1]; // assume ch_bra = ch_ket for 3body...
    TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    //      auto& Gamma =  itmat->second;
    for (size_t ibra = 0; ibra < tbc_bra.GetNumberKets(); ++ibra)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      Orbit &oi = modelspace->GetOrbit(i);
      Orbit &oj = modelspace->GetOrbit(j);
      size_t iket_min = ch_bra == ch_ket ? ibra : 0;
      for (size_t iket = iket_min; iket < tbc_ket.GetNumberKets(); ++iket)
      {
        Ket &ket = tbc_ket.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        Orbit &ok = modelspace->GetOrbit(k);
        Orbit &ol = modelspace->GetOrbit(l);
        //            for (auto& a : modelspace->holes)
        for (auto &a : occupied)
        {
          Orbit &oa = modelspace->GetOrbit(a);
          if ((2 * (oi.n + oj.n + oa.n) + oi.l + oj.l + oa.l) > E3max)
            continue;
          if ((2 * (ok.n + ol.n + oa.n) + ok.l + ol.l + oa.l) > E3max)
            continue;

          //               int kmin2 = abs(2*tbc_bra.J-oa.j2);
          //               int kmax2 = 2*tbc_bra.J+oa.j2;
          //               for (int K2=kmin2; K2<=kmax2; K2+=2)
          //               {
          //                  Gamma(ibra,iket) += (K2+1) * sign*oa.occ * ThreeBody.GetME_pn(tbc_bra.J,tbc_ket.J,K2,i,j,a,k,l,a); // This is unnormalized.
          //                   std::cout << " accessing 3bme   "<< tbc_bra.J << " " << tbc_ket.J << " " << K2 << "    " << i << " " << j << " " << a << "  " << k << " "  << l << " " << a << "       " << ThreeBody.GetME_pn(tbc_bra.J,tbc_ket.J,K2,i,j,a,k,l,a) << "  ->  " << Gamma(ibra,iket) << std::endl;
          //                                                   }

          Gamma(ibra, iket) += sign * oa.occ * ThreeBody.GetME_pn_no2b(i, j, a, k, l, a, tbc_bra.J);
        }
        Gamma(ibra, iket) /= (2 * tbc_bra.J + 1) * sqrt((1 + bra.delta_pq()) * (1 + ket.delta_pq()));
        if (opNO3.GetTRank() != 0 or opNO3.GetParity() != 0)
        {
          Gamma(ibra, iket) *= sqrt(2 * tbc_bra.J + 1); // reduced matrix element
        }
      }
    }
  }
  opNO3.Symmetrize();
  Operator opNO2 = opNO3.DoNormalOrdering2(sign, occupied);
  opNO2.ScaleZeroBody(1. / 3.);
  opNO2.ScaleOneBody(1. / 2.);
  //   std::cout << "IN " << __func__ << "  line " << __LINE__ << "   norms of NO 3b pieces are " << opNO2.ZeroBody << "   " << opNO2.OneBodyNorm() << "   " << opNO2.TwoBodyNorm() << "  and thie original 3b norm was  " << ThreeBody.Norm() << "  which produced a no2b with norm " << opNO3.TwoBodyNorm() << std::endl;
  //   std::cout << " opNO2 has storage mode " << opNO2.ThreeBody.GetStorageMode() << "  and this has storage mode " << ThreeBody.GetStorageMode() << "  and opNO3 has " << opNO3.ThreeBody.GetStorageMode() << std::endl;
  //   std::cout << "Are they allocated? " << opNO2.ThreeBody.IsAllocated() << "  " << ThreeBody.IsAllocated() << "  " << opNO3.ThreeBody.IsAllocated() << std::endl;
  std::cout << __func__ << "  contributed " << opNO2.ZeroBody << "  to the zero body part" << std::endl;
  // Also normal order the 1 and 2 body pieces
  opNO2 += DoNormalOrdering2(sign, occupied);
  opNO2.ThreeBody.SetMode("pn");

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return opNO2;
}

///  The normal ordering is slightly different if the operator is a
///  dagger operator.
///
Operator Operator::DoNormalOrderingDagger(int sign, std::set<index_t> occupied) const
{
  Operator opNO(*this);

  index_t Q = opNO.GetQSpaceOrbit();
  Orbit &oQ = modelspace->GetOrbit(Q);

  for (auto &itmat : ThreeLeg.MatEl)
  {
    int ch_bra = itmat.first;
    //     int ch_ket = itmat.first[1];

    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch_bra);
    int J = tbc.J;
    double hatfactor = 2 * J + 1.0;

    // a dagger part

    for (auto a : OneBodyChannels.at({oQ.l, oQ.j2, oQ.tz2}))
    {
      Orbit &oa = modelspace->GetOrbit(a);
      double ja = oa.j2 * 0.5;
      //        for (auto& h : modelspace->holes)  // C++11 syntax
      for (auto &h : occupied) // C++11 syntax
      {
        Orbit &oh = modelspace->GetOrbit(h);
        // TODO: Confirm that this should be a minus sign rather than plus
        opNO.OneBody(a, 0) -= hatfactor / (2 * ja + 1) * sign * oh.occ * ThreeLeg.GetME(ch_bra, h, a, h); // The GetTBME returns an unnormalized matrix element.
                                                                                                          //          opNO.OneBody(a,Q) -= hatfactor/(2*ja+1) * sign*oh.occ * ThreeLeg.GetME(ch_bra,h,a,h);  // The GetTBME returns an unnormalized matrix element.
                                                                                                          //          opNO.OneBody(a,Q) -= hatfactor/(2*ja+1) * sign*oh.occ * TwoBody.GetTBME(ch_bra,ch_ket,h,a,h,Q);  // The GetTBME returns an unnormalized matrix element.
      }
    }
  } // loop over channels

  return opNO;
}

//********************************************
/// Truncate an operator to a smaller emax (and E3max)
/// A corresponding ModelSpace object must be
/// created at the appropriate scope.
//********************************************
Operator Operator::Truncate(ModelSpace &ms_new)
{
  Operator OpNew(ms_new, rank_J, rank_T, parity, particle_rank);

  int new_emax = ms_new.GetEmax();
  if (new_emax > modelspace->GetEmax())
  {
    std::cout << "Error: Cannot truncate an operator with emax = " << modelspace->GetEmax() << " to one with emax = " << new_emax << std::endl;
    return OpNew;
  }

  OpNew.ZeroBody = ZeroBody;
  OpNew.hermitian = hermitian;
  OpNew.antihermitian = antihermitian;
  size_t norb = ms_new.GetNumberOrbits();
  arma::uvec old_orbs(norb);
  for (size_t i = 0; i < norb; i++)
  {
    Orbit &oi = ms_new.GetOrbit(i);
    size_t iold = modelspace->GetOrbitIndex(oi.n, oi.l, oi.j2, oi.tz2);
    old_orbs(i) = iold;
  }
  //  OpNew.OneBody = OneBody.submat(0,0,norb-1,norb-1);
  OpNew.OneBody = OneBody.submat(old_orbs, old_orbs);
  //  std::cout << "Done truncating one body " << std::endl << OneBody << std::endl << std::endl << OpNew.OneBody << std::endl;

  for (auto &itmat : OpNew.TwoBody.MatEl)
  {

    int ch_bra_new = itmat.first[0];
    int ch_ket_new = itmat.first[1];
    TwoBodyChannel &tbc_bra_new = ms_new.GetTwoBodyChannel(ch_bra_new);
    TwoBodyChannel &tbc_ket_new = ms_new.GetTwoBodyChannel(ch_ket_new);
    auto &Mat_new = itmat.second;

    int ch_bra_old = modelspace->GetTwoBodyChannelIndex(tbc_bra_new.J, tbc_bra_new.parity, tbc_bra_new.Tz);
    int ch_ket_old = modelspace->GetTwoBodyChannelIndex(tbc_ket_new.J, tbc_ket_new.parity, tbc_ket_new.Tz);
    TwoBodyChannel &tbc_bra_old = modelspace->GetTwoBodyChannel(ch_bra_old);
    TwoBodyChannel &tbc_ket_old = modelspace->GetTwoBodyChannel(ch_ket_old);
    auto &Mat_old = TwoBody.GetMatrix(ch_bra_old, ch_ket_old);

    int nbras_new = tbc_bra_new.GetNumberKets();
    int nkets_new = tbc_ket_new.GetNumberKets();
    arma::uvec ibra_old(nbras_new);
    arma::uvec iket_old(nkets_new);

    for (int ibra = 0; ibra < nbras_new; ++ibra)
    {
      auto &bra_new = tbc_bra_new.GetKet(ibra);
      ibra_old(ibra) = tbc_bra_old.GetLocalIndex(old_orbs(bra_new.p), old_orbs(bra_new.q));
    }

    for (int iket = 0; iket < nkets_new; ++iket)
    {
      auto ket_new = tbc_ket_new.GetKet(iket);
      iket_old(iket) = tbc_ket_old.GetLocalIndex(old_orbs(ket_new.p), old_orbs(ket_new.q));
    }

    Mat_new = Mat_old.submat(ibra_old, iket_old);
  }

  // We may also want to truncate an operator with a 3N part
  if (particle_rank >= 3)
  {
    auto storage_mode = ThreeBody.GetStorageMode();
    OpNew.ThreeBody.SetMode(storage_mode);
    size_t nch = ms_new.GetNumberThreeBodyChannels();
    for (size_t ch3 = 0; ch3 < nch; ch3++)
    {
      ThreeBodyChannel Tbc_new = ms_new.GetThreeBodyChannel(ch3);
      size_t ch3_old = modelspace->GetThreeBodyChannelIndex(Tbc_new.twoJ, Tbc_new.parity, Tbc_new.twoTz);
      ThreeBodyChannel Tbc_old = modelspace->GetThreeBodyChannel(ch3_old);
      size_t nkets = Tbc_new.GetNumberKets();
      for (size_t ibra = 0; ibra < nkets; ibra++)
      {
        auto &bra = Tbc_new.GetKet(ibra);
        size_t ibra_old = Tbc_old.GetLocalIndex(bra.p, bra.q, bra.r, bra.Jpq);
        // we store bra <= ket
        for (size_t iket = ibra; iket < nkets; iket++)
        {
          auto &ket = Tbc_new.GetKet(iket);
          size_t iket_old = Tbc_old.GetLocalIndex(ket.p, ket.q, ket.r, ket.Jpq);
          if (storage_mode == "pn")
          {
            auto matel = ThreeBody.GetME_pn_ch(ch3_old, ch3_old, ibra_old, iket_old);
            OpNew.ThreeBody.SetME_pn_ch(ch3, ch3, ibra, iket, matel);
          }
          else if (storage_mode == "iso")
          {
            for (int tab = 0; tab <= 1; tab++)
            {
              for (int tde = 0; tde <= 1; tde++)
              {
                for (int twoT = std::abs(Tbc_new.twoTz); twoT <= 3; twoT += 2)
                {
                  auto matel = ThreeBody.GetME_iso(bra.Jpq, ket.Jpq, Tbc_new.twoJ, tab, tde, twoT, bra.p, bra.q, bra.r, ket.p, ket.q, ket.r);
                  OpNew.ThreeBody.SetME_iso(bra.Jpq, ket.Jpq, Tbc_new.twoJ, tab, tde, twoT, bra.p, bra.q, bra.r, ket.p, ket.q, ket.r, matel);
                }
              }
            }
          }

        } // for iket
      } // for ibra
    } // for ch3

  } // if particle_rank>=3

  return OpNew;
}


Operator Operator::DoIsospinAveraging() const
{
   Operator OpIso = 0*(*this);
   OpIso.ZeroBody = this->ZeroBody;

   for (auto p : modelspace->proton_orbits )
   {
      Orbit& op = modelspace->GetOrbit(p);
      int n = modelspace->GetOrbitIndex( op.n, op.l, op.j2, -op.tz2 );
      for ( auto pp : modelspace->proton_orbits )
      {
         Orbit& opp = modelspace->GetOrbit(pp);
         int nn = modelspace->GetOrbitIndex( opp.n, opp.l, opp.j2, -opp.tz2 );
         double vavg = (this->OneBody(p,pp) + this->OneBody(n,nn) )/2;
         OpIso.OneBody(p,pp) = vavg;
         OpIso.OneBody(n,nn) = vavg;
      }
   }


   // We loop over only the Tz=0 (proton-neutron) channel
   // and then work out the pp and nn terms
   int nch = modelspace->GetNumberTwoBodyChannels();
   for ( int ch=0; ch<nch; ch++)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      if ( tbc.Tz != 0 ) continue;
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets; ibra++)
      {
        Ket& bra = tbc.GetKet( ibra );
        int ap = bra.p;
        int bp = bra.q;
        int an = modelspace->GetOrbitIndex( bra.op->n, bra.op->l, bra.op->j2, -bra.op->tz2 );
        int bn = modelspace->GetOrbitIndex( bra.oq->n, bra.oq->l, bra.oq->j2, -bra.oq->tz2 );
        if ( bra.op->tz2==1 )  std::swap(ap,an);
        if ( bra.oq->tz2==1 )  std::swap(bp,bn);
//        for (int iket=0;iket<nkets; iket++)
        for (int iket=ibra;iket<nkets; iket++)
        {
           Ket& ket = tbc.GetKet( iket );
           int cp = ket.p;
           int dp = ket.q;
           int cn = modelspace->GetOrbitIndex( ket.op->n, ket.op->l, ket.op->j2, -ket.op->tz2 );
           int dn = modelspace->GetOrbitIndex( ket.oq->n, ket.oq->l, ket.oq->j2, -ket.oq->tz2 );
           if ( ket.op->tz2==1 )  std::swap(cp,cn);
           if ( ket.oq->tz2==1 )  std::swap(dp,dn);
//           double Vpppp = this->TwoBody.GetTBME(ch,ch,ibra,iket);
           double Vpppp = this->TwoBody.GetTBME_J(tbc.J,tbc.J, ap,bp,cp,dp);
           double Vnnnn = this->TwoBody.GetTBME_J(tbc.J,tbc.J, an,bn,cn,dn);
           double Vpnpn = this->TwoBody.GetTBME_J(tbc.J,tbc.J, ap,bn,cp,dn);
           double Vpnnp = this->TwoBody.GetTBME_J(tbc.J,tbc.J, ap,bn,cn,dp);
           double Vnpnp = this->TwoBody.GetTBME_J(tbc.J,tbc.J, an,bp,cn,dp);
           double Vnppn = this->TwoBody.GetTBME_J(tbc.J,tbc.J, an,bp,cp,dn);
           double VT1 = ( Vpppp + Vnnnn + 0.5*(Vpnpn + Vnpnp + Vpnnp + Vnppn) ) / 3;
           double VT0 = 0.5*(Vpnpn + Vnpnp - Vpnnp - Vnppn);
         
           if ( not ( (ap==bp or cp==dp) and tbc.J%2==1 )  ) // Only set these if T=1 channel exists
           {
             double Norm_pp = 1.0;
             if ( ap==bp) Norm_pp /= PhysConst::SQRT2;
             if ( cp==dp) Norm_pp /= PhysConst::SQRT2;
             OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  ap,bp,cp,dp,  VT1*Norm_pp); // pppp
             OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  an,bn,cn,dn,  VT1*Norm_pp); // nnnn
           }
           OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  ap,bn,cp,dn,  0.5*(VT1+VT0) ); // pnpn
           OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  an,bp,cn,dp,  0.5*(VT1+VT0) ); // npnp
           OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  an,bp,cp,dn,  0.5*(VT1-VT0) ); // nppn
           OpIso.TwoBody.SetTBME_J( tbc.J,tbc.J,  ap,bn,cn,dp,  0.5*(VT1-VT0) ); // pnnp
           
        }
      }
   }
   return OpIso;
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
  if (legs >= 6)
    ThreeBody.Erase();
  if ((legs % 2) > 0)
    ThreeLeg.Erase();
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
  ThreeBody.Erase();
  //  ThreeBody = ThreeBodyME();
}

void Operator::EraseThreeLeg()
{
  ThreeLeg.Erase();
}

void Operator::SetHermitian()
{
  hermitian = true;
  antihermitian = false;
  TwoBody.SetHermitian();
  ThreeBody.SetHermitian();
}

void Operator::SetAntiHermitian()
{
  hermitian = false;
  antihermitian = true;
  TwoBody.SetAntiHermitian();
  ThreeBody.SetAntiHermitian();
}

void Operator::SetNonHermitian()
{
  hermitian = false;
  antihermitian = false;
  TwoBody.SetNonHermitian();
}

void Operator::SetNumberLegs(int l)
{
  int old_legs = legs;
  legs = l;
  if (l == old_legs)
    return;
  if (legs % 2 == 0)
  {
    if (old_legs < 4)
      TwoBody = TwoBodyME(modelspace, rank_J, rank_T, parity);
    if (TwoBody.MatEl.size() < 1)
      TwoBody.Allocate();
    ThreeLeg.Deallocate();
    //    if (legs>5 and (not ThreeBody.is_allocated)) ThreeBody.Allocate();
    if (legs > 5 and (not ThreeBody.IsAllocated()))
      ThreeBody.Allocate();
    particle_rank = l / 2;
  }
  else
  {
    TwoBody.Deallocate();
    ThreeBody.Deallocate();
    OneBody.zeros(modelspace->GetNumberOrbits(), 1); // reduce it to a single column
    ThreeLeg.Allocate();
    OneBody.zeros(modelspace->GetNumberOrbits(), 1);
  }
}

void Operator::MakeReduced()
{
  if (is_reduced)
  {
    std::cout << "Calling MakeReduced(), but this operator is already reduced." << std::endl;
  }
  if (rank_J > 0)
  {
    std::cout << "Trying to reduce an operator with J rank = " << rank_J << ". Not good!!!" << std::endl;
    return;
  }
  for (size_t a = 0; a < modelspace->GetNumberOrbits(); ++a)
  {
    Orbit &oa = modelspace->GetOrbit(a);
    for (size_t b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
    {
      //      if (b<a) continue;
      OneBody(a, b) *= sqrt(oa.j2 + 1);
      //      OneBody(b,a) *= sqrt(oa.j2+1);
    }
  }
  for (auto &itmat : TwoBody.MatEl)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(itmat.first[0]);
    itmat.second *= sqrt(2 * tbc.J + 1);
  }

  if (particle_rank > 2)
  {
    for (auto &it : ThreeBody.Get_ch_start())
    {
      size_t ThCH_bra = it.first.ch_bra;
      size_t ThCH_ket = it.first.ch_ket;
      ThreeBodyChannel &Tbc_bra = modelspace->GetThreeBodyChannel(ThCH_bra);
      ThreeBodyChannel &Tbc_ket = modelspace->GetThreeBodyChannel(ThCH_ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        size_t nket3 = Tbc_bra.GetNumberKets();
        for (size_t iket = ibra; iket < nket3; iket++)
        {
          double ME3b = ThreeBody.GetME_pn_ch(ThCH_bra, ThCH_ket, ibra, iket);
          ThreeBody.SetME_pn_ch(ThCH_bra, ThCH_ket, ibra, iket, ME3b * sqrt(Tbc_bra.twoJ + 1) );
        }
      }
    }
  }

  is_reduced = true;
}

void Operator::MakeNotReduced()
{
  if (not is_reduced)
  {
    std::cout << "Calling MakeNotReduced(), but this operator is already not reduced." << std::endl;
  }
  if (rank_J > 0)
  {
    std::cout << "Trying to un-reduce an operator with J rank = " << rank_J << ". Not good!!!" << std::endl;
    return;
  }
  for (size_t a = 0; a < modelspace->GetNumberOrbits(); ++a)
  {
    Orbit &oa = modelspace->GetOrbit(a);
    for (size_t b : OneBodyChannels.at({oa.l, oa.j2, oa.tz2}))
    {
      //      if (b<a) continue;
      OneBody(a, b) /= sqrt(oa.j2 + 1);
      //      OneBody(b,a) = OneBody(a,b);
    }
  }
  for (auto &itmat : TwoBody.MatEl)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(itmat.first[0]);
    itmat.second /= sqrt(2 * tbc.J + 1);
  }

  if (particle_rank > 2)
  {
    for (auto &it : ThreeBody.Get_ch_start())
    {
      size_t ThCH_bra = it.first.ch_bra;
      size_t ThCH_ket = it.first.ch_ket;
      ThreeBodyChannel &Tbc_bra = modelspace->GetThreeBodyChannel(ThCH_bra);
      ThreeBodyChannel &Tbc_ket = modelspace->GetThreeBodyChannel(ThCH_ket);
      size_t nbras3 = Tbc_bra.GetNumberKets();
      for (size_t ibra = 0; ibra < nbras3; ibra++)
      {
        size_t nket3 = Tbc_bra.GetNumberKets();
        for (size_t iket = ibra; iket < nket3; iket++)
        {
          double ME3b = ThreeBody.GetME_pn_ch(ThCH_bra, ThCH_ket, ibra, iket);
          ThreeBody.SetME_pn_ch(ThCH_bra, ThCH_ket, ibra, iket, ME3b / sqrt(Tbc_bra.twoJ + 1) );
        }
      }
    }
  }

  is_reduced = false;
}

//// this routine then multiplies the TBME <ab|Op|cd> by coeff if a==b, and again if c==d
// void Operator::ChangeNormalization( double coeff )
//{
//   for (auto& it_mat : TwoBody.MatEl )
//   {
//     int ch_bra = it_mat.first[0];
//     int ch_ket = it_mat.first[1];
//     TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
//     TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
//     int nbras = tbc_bra.GetNumberKets();
//     int nkets = tbc_ket.GetNumberKets();
//     for (int ibra=0; ibra<nbras; ++ibra)
//     {
//       Ket& bra = tbc_bra.GetKet(ibra);
//       if ( bra.p == bra.q ) it_mat.second.row(ibra) *= coeff;
//     }
//     for (int iket=0; iket<nkets; ++iket)
//     {
//       Ket& ket = tbc_ket.GetKet(iket);
//       if ( ket.p == ket.q ) it_mat.second.col(iket) *= coeff;
//     }
//   }
//
// }

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
// void Operator::Eye()
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
  //   std::cout << "  a    b    i    j    J    na     nb    tbme    denom     dE" << std::endl;
  double t_start = omp_get_wtime();
  double Emp2 = 0;
  int nparticles = modelspace->particles.size();
  std::vector<index_t> particles_vec(modelspace->particles.begin(), modelspace->particles.end()); // convert set to vector for OMP looping
                                                                                                  //   for ( auto& i : modelspace->particles)
                                                                                                  //   #pragma omp parallel for reduction(+:Emp2)
  for (int ii = 0; ii < nparticles; ++ii)
  {
    //     std::cout << " i = " << i << std::endl;
    //     index_t i = modelspace->particles[ii];
    index_t i = particles_vec[ii];
    double ei = OneBody(i, i);
    Orbit &oi = modelspace->GetOrbit(i);
    for (auto &a : modelspace->holes)
    {
      Orbit &oa = modelspace->GetOrbit(a);
      double ea = OneBody(a, a);
      if (abs(OneBody(i, a)) > 1e-9)
        Emp2 += (oa.j2 + 1) * oa.occ * OneBody(i, a) * OneBody(i, a) / (OneBody(a, a) - OneBody(i, i));
      for (index_t j : modelspace->particles)
      {
        if (j < i)
          continue;
        double ej = OneBody(j, j);
        Orbit &oj = modelspace->GetOrbit(j);
        for (auto &b : modelspace->holes)
        {
          if (b < a)
            continue;
          Orbit &ob = modelspace->GetOrbit(b);
          if ((oi.l + oj.l + oa.l + ob.l) % 2 > 0)
            continue;
          if ((oi.tz2 + oj.tz2) != (oa.tz2 + ob.tz2))
            continue;
          double eb = OneBody(b, b);
          double denom = ea + eb - ei - ej;
          int Jmin = std::max(std::abs(oi.j2 - oj.j2), std::abs(oa.j2 - ob.j2)) / 2;
          int Jmax = std::min(oi.j2 + oj.j2, oa.j2 + ob.j2) / 2;
          int dJ = 1;
          if (a == b or i == j)
          {
            Jmin += Jmin % 2;
            dJ = 2;
          }
          for (int J = Jmin; J <= Jmax; J += dJ)
          {
            double tbme = TwoBody.GetTBME_J_norm(J, a, b, i, j);
            if (std::abs(tbme) > 1e-9)
            {
              Emp2 += (2 * J + 1) * oa.occ * ob.occ * tbme * tbme / denom; // no factor 1/4 because of the restricted sum
                                                                           //              std::cout << "MBPT2 " << a << " " << b << " " << i << " " << j << "    " << J << "  " << oa.occ << " " << ob.occ << "  "
                                                                           //                        << std::setw(12) << std::setprecision(6) << tbme << " "
                                                                           //                        << std::setw(12) << std::setprecision(6) << denom << "   "
                                                                           //                        << std::setw(12) << std::setprecision(6) << (2*J+1) * oa.occ*ob.occ*tbme*tbme/denom  << std::endl;
            }
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
// double Operator::GetMP3_Energy()
std::array<double, 3> Operator::GetMP3_Energy()
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
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Epp, Ehh)
  for (int ich = 0; ich < nch; ++ich)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    auto &Mat = TwoBody.GetMatrix(ich, ich);

    size_t n_hh = tbc.GetKetIndex_hh().size();
    size_t n_pp = tbc.GetKetIndex_pp().size();
    arma::mat M_hhpp(n_hh, n_pp, arma::fill::zeros);
    arma::mat M_hhhh(n_hh, n_hh, arma::fill::zeros);
    arma::mat M_pppp(n_pp, n_pp, arma::fill::zeros);

    size_t I_hh = 0;
    for (auto iket_ij : tbc.GetKetIndex_hh())
    {
      Ket &ket_ij = tbc.GetKet(iket_ij);
      index_t i = ket_ij.p;
      index_t j = ket_ij.q;

      size_t II_pp = 0;
      for (auto iket_ab : tbc.GetKetIndex_pp())
      {
        Ket &ket_ab = tbc.GetKet(iket_ab);
        index_t a = ket_ab.p;
        index_t b = ket_ab.q;
        double Delta_ijab = OneBody(i, i) + OneBody(j, j) - OneBody(a, a) - OneBody(b, b);
        M_hhpp(I_hh, II_pp) = Mat(iket_ij, iket_ab) / Delta_ijab;
        II_pp++;
      }
      size_t II_hh = 0;
      for (auto iket_kl : tbc.GetKetIndex_hh())
      {
        //         Ket& ket_kl = tbc.GetKet(iket_kl);
        //         index_t k = ket_kl.p;
        //         index_t l = ket_kl.q;
        M_hhhh(I_hh, II_hh) = Mat(iket_ij, iket_kl);
        II_hh++;
      }
      I_hh++;
    }

    size_t I_pp = 0;
    for (auto iket_ab : tbc.GetKetIndex_pp())
    {
      //       Ket& ket_ab = tbc.GetKet(iket_ab);
      //       index_t a = ket_ab.p;
      //       index_t b = ket_ab.q;
      size_t II_pp = 0;
      for (auto iket_cd : tbc.GetKetIndex_pp())
      {
        //         Ket& ket_cd = tbc.GetKet(iket_cd);
        //         index_t c = ket_cd.p;
        //         index_t d = ket_cd.q;
        M_pppp(I_pp, II_pp) = Mat(iket_ab, iket_cd);
        II_pp++;
      }
      I_pp++;
    }

    int J = tbc.J;
    Ehh += (2 * J + 1) * arma::trace(M_hhpp.t() * M_hhhh * M_hhpp);
    Epp += (2 * J + 1) * arma::trace(M_hhpp * M_pppp * M_hhpp.t());

  } // for ich
  ////   cout << "done with pp and hh. E(3) = " << Emp3 << endl;

  //   index_t nparticles = modelspace->particles.size();
  modelspace->PreCalculateSixJ();

  int nch_CC = modelspace->GetNumberTwoBodyChannels_CC();

  //   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Emp3)
  //   #pragma omp parallel for  schedule(dynamic,1) reduction(+:Epp,Ehh)
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Eph)
  for (int ich_CC = 0; ich_CC < nch_CC; ++ich_CC)
  {
    TwoBodyChannel_CC &tbc_CC = modelspace->GetTwoBodyChannel_CC(ich_CC);
    size_t nkets_ph = tbc_CC.GetKetIndex_ph().size();
    arma::mat Vbar_iabj(nkets_ph, nkets_ph, arma::fill::zeros);
    arma::mat Vbar_bjck(nkets_ph, nkets_ph, arma::fill::zeros);
    int Jph = tbc_CC.J;

    size_t I_ph = 0;
    for (auto iket_ai : tbc_CC.GetKetIndex_ph())
    {
      Ket &ket_ai = tbc_CC.GetKet(iket_ai);
      index_t a = ket_ai.p;
      index_t i = ket_ai.q;
      double ja = 0.5 * modelspace->GetOrbit(a).j2;
      double ji = 0.5 * modelspace->GetOrbit(i).j2;

      int phase_ai = 1;
      int phase_ia = -AngMom::phase(ja + ji - Jph);
      if (ket_ai.op->occ < ket_ai.oq->occ)
      {
        std::swap(a, i);
        std::swap(ja, ji);
        std::swap(phase_ai, phase_ia);
      }

      size_t II_ph = 0;
      for (auto iket_bj : tbc_CC.GetKetIndex_ph())
      {
        Ket &ket_bj = tbc_CC.GetKet(iket_bj);
        index_t b = ket_bj.p;
        index_t j = ket_bj.q;

        double jb = 0.5 * modelspace->GetOrbit(b).j2;
        double jj = 0.5 * modelspace->GetOrbit(j).j2;

        int phase_bj = 1;
        int phase_jb = -AngMom::phase(jb + jj - Jph);
        if (ket_bj.op->occ < ket_bj.oq->occ)
        {
          std::swap(b, j);
          std::swap(jb, jj);
          std::swap(phase_bj, phase_jb);
        }

        double Delta_ijab = OneBody(i, i) + OneBody(j, j) - OneBody(a, a) - OneBody(b, b);
        int J1min = std::max(std::abs(ja - jb), std::abs(ji - jj));
        int J1max = std::min(ja + jb, ji + jj);
        double tbme_iabj = 0;
        double tbme_bjck = 0;
        //         double tbme_ckia = 0;

        if (AngMom::Triangle(jj, jb, Jph) and AngMom::Triangle(ji, ja, Jph))
        {
          for (int J1 = J1min; J1 <= J1max; ++J1) // Pandya 1: <ai`| V |jb`>_Jtot
          {
            tbme_iabj -= modelspace->GetSixJ(ja, ji, Jph, jj, jb, J1) * (2 * J1 + 1) * TwoBody.GetTBME_J(J1, i, j, b, a);
          }
        }

        J1min = std::max(std::abs(ji - jb), std::abs(ja - jj));
        J1max = std::min(ji + jb, ja + jj);

        if (AngMom::Triangle(jj, jb, Jph) and AngMom::Triangle(ji, ja, Jph))
        {
          for (int J1 = J1min; J1 <= J1max; ++J1) // Pandya 1: <ai`| V |jb`>_Jtot
          {
            tbme_bjck -= modelspace->GetSixJ(jb, jj, Jph, ja, ji, J1) * (2 * J1 + 1) * TwoBody.GetTBME_J(J1, b, i, a, j);
          }
        }

        Vbar_iabj(I_ph, II_ph) = tbme_iabj * phase_ia * phase_bj / Delta_ijab;
        Vbar_bjck(II_ph, I_ph) = tbme_bjck * phase_bj * phase_ai;
        II_ph++;
      }
      I_ph++;
    }
    auto Vbar_ckia = Vbar_iabj.t();
    Eph += (2 * Jph + 1) * arma::trace(Vbar_iabj * Vbar_bjck * Vbar_ckia);

  } // for ich_CC

  IMSRGProfiler::timer["GetMP3_Energy"] += omp_get_wtime() - t_start;
  //   return Emp3;
  return {Epp, Ehh, Eph};
}

//*************************************************************
/// The second order MBPT correction due to the 3B interaction
/// using Moller-Plesset energy denominators
//*************************************************************
double Operator::GetMP2_3BEnergy()
{
  double t_start = omp_get_wtime();
  double Emp2 = 0;
  if (legs < 6)
    return 0;
  //   if ( not ThreeBody.is_allocated ) return 0;
  if (not ThreeBody.IsAllocated())
    return 0;
  size_t nch3 = modelspace->GetNumberThreeBodyChannels();
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : Emp2)
  for (size_t ch3 = 0; ch3 < nch3; ch3++)
  {
    ThreeBodyChannel &Tbc = modelspace->GetThreeBodyChannel(ch3);
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumberKets();
    for (size_t ibra = 0; ibra < nkets; ibra++)
    {
      Ket3 &bra = Tbc.GetKet(ibra);
      double occ_bra = (bra.op->occ) * (bra.oq->occ) * (bra.oR->occ);
      if (std::abs(occ_bra) < 1e-9)
        continue;
      size_t i = bra.p;
      size_t j = bra.q;
      size_t k = bra.r;
      double symm_ijk = 6;
      if (i == j and i == k)
        symm_ijk = 1;
      else if (i == j or i == k)
        symm_ijk = 3;
      double Eijk = OneBody(i, i) + OneBody(j, j) + OneBody(k, k);
      for (size_t iket = 0; iket < nkets; iket++)
      {
        Ket3 &ket = Tbc.GetKet(iket);
        double unocc_ket = (1 - ket.op->occ) * (1 - ket.oq->occ) * (1 - ket.oR->occ);
        if (std::abs(unocc_ket) < 1e-9)
          continue;
        size_t a = ket.p;
        size_t b = ket.q;
        size_t c = ket.r;
        double symm_abc = 6;
        if (a == b and a == c)
          symm_abc = 1;
        else if (a == b or a == c)
          symm_abc = 3;
        double Eabc = OneBody(a, a) + OneBody(b, b) + OneBody(c, c);
        //         double V = ThreeBody.GetME_pn_PN_ch(ch3,ch3,ibra,iket);
        double V = ThreeBody.GetME_pn_ch(ch3, ch3, ibra, iket);
        Emp2 += 1. / 36 * symm_ijk * symm_abc * (twoJ + 1) * occ_bra * unocc_ket * V * V / (Eijk - Eabc);
      } // for iket
    } // for ibra
  } // for ch3

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;
  return Emp2;
}

//////////////////////////////////////////////////////////
///  All-order resummation of pp ladders and hh ladders
////////////////////////////////////////////////////////////
std::array<double, 2> Operator::GetPPHH_Ladders()
{
  // So far, the pp and hh parts seem to work. No such luck for the ph.
  double t_start = omp_get_wtime();
  //   double Emp3 = 0;
  double Epp = 0;
  double Ehh = 0;
  //   double Eph = 0;
  // This can certainly be optimized, but I'll wait until this is the bottleneck.
  int nch = modelspace->GetNumberTwoBodyChannels();

  for (int ich = 0; ich < nch; ++ich)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
    auto &Mat = TwoBody.GetMatrix(ich, ich);

    size_t n_hh = tbc.GetKetIndex_hh().size();
    size_t n_pp = tbc.GetKetIndex_pp().size();

    for (auto iket_ij : tbc.GetKetIndex_hh())
    {
      Ket &ket_ij = tbc.GetKet(iket_ij);
      index_t i = ket_ij.p;
      index_t j = ket_ij.q;

      arma::mat ONE(n_pp, n_pp, arma::fill::eye);
      arma::mat M_pppp(n_pp, n_pp, arma::fill::zeros);
      arma::mat M_ijpp(1, n_pp, arma::fill::zeros);
      arma::mat M_ppij(n_pp, 1, arma::fill::zeros);

      size_t I_pp = 0;
      for (auto iket_ab : tbc.GetKetIndex_pp())
      {
        Ket &ket_ab = tbc.GetKet(iket_ab);
        index_t a = ket_ab.p;
        index_t b = ket_ab.q;
        double Delta_ijab = OneBody(i, i) + OneBody(j, j) - OneBody(a, a) - OneBody(b, b);
        M_ijpp(0, I_pp) = Mat(iket_ij, iket_ab);
        M_ppij(I_pp, 0) = Mat(iket_ab, iket_ij) / Delta_ijab;
        size_t II_pp = 0;
        for (auto iket_cd : tbc.GetKetIndex_pp())
        {
          //           Ket& ket_cd = tbc.GetKet(iket_cd);
          //           index_t c = ket_cd.p;
          //           index_t d = ket_cd.q;
          M_pppp(I_pp, II_pp) += Mat(iket_ab, iket_cd) / Delta_ijab;
          II_pp++;
        }
        I_pp++;
      }
      int J = tbc.J;
      Epp += (2 * J + 1) * arma::trace(M_ijpp * (arma::inv(ONE - M_pppp) - ONE) * M_ppij);

    } // for iket_ij

    for (auto iket_ab : tbc.GetKetIndex_pp())
    {
      Ket &ket_ab = tbc.GetKet(iket_ab);
      index_t a = ket_ab.p;
      index_t b = ket_ab.q;

      arma::mat ONE(n_hh, n_hh, arma::fill::eye);
      arma::mat M_hhhh(n_hh, n_hh, arma::fill::zeros);
      arma::mat M_abhh(1, n_hh, arma::fill::zeros);
      arma::mat M_hhab(n_hh, 1, arma::fill::zeros);

      size_t I_hh = 0;
      for (auto iket_ij : tbc.GetKetIndex_hh())
      {
        Ket &ket_ij = tbc.GetKet(iket_ij);
        index_t i = ket_ij.p;
        index_t j = ket_ij.q;
        double Delta_ijab = OneBody(i, i) + OneBody(j, j) - OneBody(a, a) - OneBody(b, b);
        M_abhh(0, I_hh) = Mat(iket_ab, iket_ij);
        M_hhab(I_hh, 0) = Mat(iket_ij, iket_ab) / Delta_ijab;
        size_t II_hh = 0;
        for (auto iket_kl : tbc.GetKetIndex_hh())
        {
          //           Ket& ket_kl = tbc.GetKet(iket_kl);
          //           index_t k = ket_kl.p;
          //           index_t l = ket_kl.q;
          M_hhhh(I_hh, II_hh) += Mat(iket_ij, iket_kl) / Delta_ijab;
          II_hh++;
        }
        I_hh++;
      }
      int J = tbc.J;
      Ehh += (2 * J + 1) * arma::trace(M_abhh * (arma::inv(ONE - M_hhhh) - ONE) * M_hhab);
      //       Ehh += (2*J+1) * arma::trace( M_abhh * M_hhhh  * M_hhab );

    } // for iket_ij

  } // for ich

  IMSRGProfiler::timer[__func__] += omp_get_wtime() - t_start;

  return {Epp, Ehh};
}

//*************************************************************
/// Evaluate first order perturbative correction to the operator's
/// ground-state expectation value. A HF basis is assumed.
/// \f[
///  \mathcal{O}^{(1)} = 2\sum_{abij} \frac{H_{abij}\mathcal{O}_{ijab}}{\Delta_{abij}}
/// \f]
///
//*************************************************************
double Operator::MP1_Eval(Operator &H)
{
  auto &Op = *this;
  double opval = 0;
  if (Op.GetJRank() == 0 and Op.GetTRank() == 0 and Op.GetParity() == 0)
  {
    int nch = modelspace->GetNumberTwoBodyChannels();
    for (int ich = 0; ich < nch; ++ich)
    {
      //      std::cout << "     ich = " << ich << std::endl;
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ich);
      int J = tbc.J;
      auto &Hmat = H.TwoBody.GetMatrix(ich, ich);
      auto &Opmat = Op.TwoBody.GetMatrix(ich, ich);
      for (auto ibra : tbc.GetKetIndex_hh())
      {
        Ket &bra = tbc.GetKet(ibra);
        for (auto iket : tbc.GetKetIndex_pp())
        {
          Ket &ket = tbc.GetKet(iket);
          double Delta_abij = H.OneBody(bra.p, bra.p) + H.OneBody(bra.q, bra.q) - H.OneBody(ket.p, ket.p) - H.OneBody(ket.q, ket.q);
          opval += 2 * (2 * J + 1) * Hmat(ibra, iket) * Opmat(iket, ibra) / Delta_abij;
        }
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
  if (legs % 2 == 0)
  {
    double n1 = OneBodyNorm();
    double n2 = TwoBody.Norm();
    double n3 = 0.;
    if (legs > 5)
      n3 = ThreeBody.Norm();
    //      return sqrt(n1*n1+n2*n2);
    return sqrt(n1 * n1 + n2 * n2 + n3 * n3);
  }
  else
  {
    double n1 = OneLegNorm();
    double n2 = ThreeLeg.Norm();
    return sqrt(n1 * n1 + n2 * n2);
  }
}

double Operator::OneBodyNorm() const
{
  double nrm = 0;
  //   for ( size_t p=0; p<modelspace->GetNumberOrbits(); ++p)
  for (auto p : modelspace->all_orbits)
  {
    Orbit &op = modelspace->GetOrbit(p);
    for (auto q : OneBodyChannels.at({op.l, op.j2, op.tz2}))
    {
      Orbit &oq = modelspace->GetOrbit(q);
      int degeneracy_factor = (op.j2 + 1) * ((std::min(oq.j2, op.j2 + rank_J) - std::max(-oq.j2, op.j2 - rank_J)) / 2 + 1);
      //       nrm += OneBody(p,q)*OneBody(p,q) * degeneracy_factor * degeneracy_factor;
      nrm += OneBody(p, q) * OneBody(p, q) * degeneracy_factor;
    }
  }
  return sqrt(nrm);
}

double Operator::TwoBodyNorm() const
{
  return TwoBody.Norm();
}

double Operator::ThreeBodyNorm() const
{
  return ThreeBody.Norm();
}

double Operator::OneLegNorm() const
{
  double nrm = 0;
  //   for ( size_t p=0; p<modelspace->GetNumberOrbits(); ++p)
  for (auto p : modelspace->all_orbits)
  {
    Orbit &op = modelspace->GetOrbit(p);
    int degeneracy_factor = (op.j2 + 1);
    nrm += OneBody(p, 0) * OneBody(p, 0) * degeneracy_factor * degeneracy_factor;
  }
  return sqrt(nrm);
}

double Operator::ThreeLegNorm() const
{
  return ThreeLeg.Norm();
}

// void Operator::MakeNormalized(){ ChangeNormalization( 1./SQRT2)  ;}
// void Operator::MakeUnNormalized(){ ChangeNormalization( SQRT2)  ;}
void Operator::MakeNormalized() { ChangeNormalization(PhysConst::INVSQRT2); }
void Operator::MakeUnNormalized() { ChangeNormalization(PhysConst::SQRT2); }
void Operator::ChangeNormalization(double factor)
{
  for (auto &itmat : TwoBody.MatEl)
  {
    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];
    auto &TBME = itmat.second;
    TwoBodyChannel &tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel &tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = tbc_bra.GetNumberKets();
    int nkets = tbc_ket.GetNumberKets();
    for (int ibra = 0; ibra < nbras; ibra++)
    {
      Ket &bra = tbc_bra.GetKet(ibra);
      if (bra.p == bra.q)
      {
        TBME.col(ibra) *= factor;
      }
    }
    for (int iket = 0; iket < nkets; iket++)
    {
      Ket &ket = tbc_ket.GetKet(iket);
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
  double M = (emax + 1) * (emax + 2) * (emax + 3) / 3; // number of m-scheme orbits for the proton or neutron partitions
  double norm_p = Ztrace / M;
  double norm_n = Ntrace / M;
  double norm_pp = Ztrace * (Ztrace - 1) / (M * (M - 1));
  double norm_nn = Ntrace * (Ntrace - 1) / (M * (M - 1));
  double norm_pn = norm_p * norm_n;

  //  for (int i=0; i<modelspace->GetNumberOrbits(); ++i)
  for (auto i : modelspace->proton_orbits)
  {
    Orbit &oi = modelspace->GetOrbit(i);
    trace += OpVac.OneBody(i, i) * (oi.j2 + 1) * norm_p;
  }
  for (auto i : modelspace->neutron_orbits)
  {
    Orbit &oi = modelspace->GetOrbit(i);
    trace += OpVac.OneBody(i, i) * (oi.j2 + 1) * norm_n;
  }

  for (size_t ch = 0; ch < modelspace->GetNumberTwoBodyChannels(); ++ch)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
    double trVijij = arma::trace(OpVac.TwoBody.GetMatrix(ch, ch));
    switch (tbc.Tz)
    {
    case (-1):
      trace += (tbc.J * 2 + 1) * trVijij * norm_pp;
      break;
    case (0):
      trace += (tbc.J * 2 + 1) * trVijij * norm_pn;
      break;
    case (1):
      trace += (tbc.J * 2 + 1) * trVijij * norm_nn;
      break;
    default:
      std::cout << "AAAHHH blew the switch statement. tbc.Tz = " << tbc.Tz << std::endl;
    }
  }
  IMSRGProfiler::timer["Operator::Trace"] += omp_get_wtime() - t_start;
  return trace;
}

void Operator::ScaleFermiDirac(Operator &H, double T, double Efermi)
{
  //  int norb = modelspace->GetNumberOrbits();

  //  for (int i=0; i<norb; ++i)
  for (auto i : modelspace->all_orbits)
  {
    Orbit &oi = modelspace->GetOrbit(i);
    double ei = OneBody(i, i);
    for (auto j : OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
    {
      double ej = OneBody(j, j);
      OneBody(i, j) *= 1. / (1 + exp((ei + ej - 2 * Efermi) / (2 * T)));
    }
  }
  for (auto itmat : TwoBody.MatEl)
  {
    TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(itmat.first[0]);
    for (size_t ibra = 0; ibra < tbc.GetNumberKets(); ++ibra)
    {
      Ket &bra = tbc.GetKet(ibra);
      int i = bra.p;
      int j = bra.q;
      double ei = OneBody(i, i);
      double ej = OneBody(j, j);
      for (size_t iket = 0; iket < tbc.GetNumberKets(); ++iket)
      {
        Ket &ket = tbc.GetKet(iket);
        int k = ket.p;
        int l = ket.q;
        double ek = OneBody(k, k);
        double el = OneBody(l, l);
        itmat.second.row(ibra) *= 1. / (1 + exp((ei + ej + ek + el - 4 * Efermi) / (2 * T)));
      }
    }
  }
}

void Operator::Symmetrize()
{
  if (rank_J == 0)
    OneBody = arma::symmatu(OneBody);
  else
  {
    //     int norb = modelspace->GetNumberOrbits();
    //     for (int i=0;i<norb; ++i)
    for (auto i : modelspace->all_orbits)
    {
      Orbit &oi = modelspace->GetOrbit(i);
      for (size_t j : OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
      {
        if (j <= i)
          continue;
        Orbit &oj = modelspace->GetOrbit(j);
        OneBody(j, i) = modelspace->phase((oi.j2 - oj.j2) / 2) * OneBody(i, j);
      }
    }
  }
  TwoBody.Symmetrize();
}

void Operator::AntiSymmetrize()
{
  if (rank_J == 0)
  {
    OneBody = arma::trimatu(OneBody) - arma::trimatu(OneBody).t();
  }
  else
  {
    //     int norb = modelspace->GetNumberOrbits();
    //     for (int i=0;i<norb; ++i)
    for (auto i : modelspace->all_orbits)
    {
      Orbit &oi = modelspace->GetOrbit(i);
      for (size_t j : OneBodyChannels.at({oi.l, oi.j2, oi.tz2}))
      {
        if (j <= i)
          continue;
        if (rank_J == 0)
          OneBody(j, i) = -OneBody(i, j);
        else
        {
          Orbit &oj = modelspace->GetOrbit(j);
          OneBody(j, i) = -modelspace->phase((oi.j2 - oj.j2) / 2) * OneBody(i, j);
        }
      }
    }
  }
  TwoBody.AntiSymmetrize();
}

/*
// Modified version of GetMP2_Energy. Determines each orbital's impact on the total MP2 energy
// (i.e., by how much would EMP2 change if this single orbital were removed)
arma::vec Operator::GetMP2_Impacts() const
{
//   std::cout << "  a    b    i    j    J    na     nb    tbme    denom     dE" << std::endl;
   double t_start = omp_get_wtime();
   double de = 0;

   int nparticles = modelspace->particles.size();
   arma::vec orbit_impacts(modelspace->all_orbits.size(), arma::fill::zeros);

   std::vector<index_t> particles_vec(modelspace->particles.begin(),modelspace->particles.end()); // convert set to vector for OMP looping
//   for ( auto& i : modelspace->particles)
//   #pragma omp parallel for reduction(+:Emp2)
   for ( int ii=0;ii<nparticles;++ii)
   {
//     index_t i = modelspace->particles[ii];
     index_t i = particles_vec[ii];
     //  std::cout << " i = " << i << std::endl;

     double ei = OneBody(i,i);
     Orbit& oi = modelspace->GetOrbit(i);
     for (auto& a : modelspace->holes)
     {
       Orbit& oa = modelspace->GetOrbit(a);
       double ea = OneBody(a,a);
       if (abs(OneBody(i,a))>1e-6)
        de = (oa.j2+1) * oa.occ * OneBody(i,a)*OneBody(i,a)/(OneBody(a,a)-OneBody(i,i));
        orbit_impacts(a) += de;
        orbit_impacts(i) += de;

       for (index_t j : modelspace->particles)
       {
         if (j<i) continue;
         double ej = OneBody(j,j);
         Orbit& oj = modelspace->GetOrbit(j);
         for ( auto& b: modelspace->holes)
         {
           if (b<a) continue;
           Orbit& ob = modelspace->GetOrbit(b);
           if ( (oi.l+oj.l+oa.l+ob.l)%2 >0) continue;
           if ( (oi.tz2 + oj.tz2) != (oa.tz2 +ob.tz2) ) continue;
           double eb = OneBody(b,b);
           double denom = ea+eb-ei-ej;
           int Jmin = std::max(std::abs(oi.j2-oj.j2),std::abs(oa.j2-ob.j2))/2;
           int Jmax = std::min(oi.j2+oj.j2,oa.j2+ob.j2)/2;
           int dJ = 1;
           if (a==b or i==j)
           {
             Jmin += Jmin%2;
             dJ=2;
           }
           for (int J=Jmin; J<=Jmax; J+=dJ)
           {
             double tbme = TwoBody.GetTBME_J_norm(J,a,b,i,j);
             if (std::abs(tbme)>1e-6)
             {
              de = (2*J+1)* oa.occ * ob.occ * tbme*tbme/denom; // no factor 1/4 because of the restricted sum

              orbit_impacts(a) += de;
              if (a!=b) orbit_impacts(b) += de;

              orbit_impacts(i) += de;
              if (i!=j) orbit_impacts(j) += de;

              }
           }
         }
       }
     }
   }
   IMSRGProfiler::timer["GetMP2_Impacts"] += omp_get_wtime() - t_start;
   return orbit_impacts;
}
*/
