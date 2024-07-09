
#include "IMSRGSolver.hh"
#include "Commutator.hh"
#include "BCH.hh"
#include "Operator.hh"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

IMSRGSolver::~IMSRGSolver()
{
  CleanupScratch();
}

IMSRGSolver::IMSRGSolver()
    : s(0), ds(0.1), ds_max(0.5),
      norm_domega(0.1), omega_norm_max(2.0), eta_criterion(1e-6), method("magnus_euler"),
      flowfile(""), n_omega_written(0), max_omega_written(500), magnus_adaptive(true), hunter_gatherer(false), perturbative_triples(false),
      /*pert_triples_this_omega(0),pert_triples_sum(0),*/ ode_monitor(*this), ode_mode("H"), ode_e_abs(1e-6), ode_e_rel(1e-6)
{
}

// Constructor
IMSRGSolver::IMSRGSolver(Operator &H_in)
    : modelspace(H_in.GetModelSpace()), H_0(&H_in), FlowingOps(1, H_in), Eta(H_in),
      istep(0), s(0), ds(0.1), ds_max(0.5),
      smax(2.0), norm_domega(0.1), omega_norm_max(2.0), eta_criterion(1e-6), method("magnus_euler"),
      flowfile(""), n_omega_written(0), max_omega_written(500), magnus_adaptive(true), hunter_gatherer(false), perturbative_triples(false),
      /*pert_triples_this_omega(0),pert_triples_sum(0),*/ ode_monitor(*this), ode_mode("H"), ode_e_abs(1e-6), ode_e_rel(1e-6)
{
  Eta.Erase();
  Eta.SetAntiHermitian();
//  Eta.ThreeBody.SetMode("pn");  // DONT DO THIS, it allocates a 3N structure even if you don't want it.
  Omega.emplace_back(Eta);
}

void IMSRGSolver::NewOmega()
{
  H_saved = FlowingOps[0];
  std::cout << "pushing back another Omega. Omega.size = " << Omega.size()
            << " , operator size = " << Omega.front().Size() / 1024. / 1024. << " MB"
            << ",  memory usage = " << profiler.CheckMem()["RSS"] / 1024. / 1024. << " GB";
  //  if ( perturbative_triples )
  //  {
  //       pert_triples_this_omega = GetPerturbativeTriples();
  //       pert_triples_sum += pert_triples_this_omega;
  //       std::cout << "  pert. triples = " << pert_triples_this_omega << "   sum = " << pert_triples_sum;
  //  }
  std::cout << std::endl;
  if (scratchdir != "")
  {

    if (scratchdir.find("/dev/null") == std::string::npos)
    {
      std::ostringstream filename;
      filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << n_omega_written;
      std::ofstream ofs(filename.str(), std::ios::binary);
      Omega.back().WriteBinary(ofs);
      std::cout << "Omega written to file " << filename.str() << "  written " << n_omega_written << " so far." << std::endl;
      if (n_omega_written > max_omega_written)
      {
        std::cout << "n_omega_written > max_omega_written.  (" << n_omega_written << " > " << max_omega_written
                  << " ) deleting OMEGA files and calling terminate." << std::endl;
        CleanupScratch();
        std::terminate();
      }
    }
    if (Omega.back().GetModelSpace() != Eta.GetModelSpace())
      Omega.back() = Eta;
    n_omega_written++;
  }
  else
  {
    Omega.emplace_back(Eta);
  }
  Omega.back().Erase();
}

// Use a hunter-gatherer approach to finding Omega.
// We only store two Omega operators, the "hunter" and the "gatherer".
// The hunter is updated by the generator eta with the BCH formula
// exp[Omega_h(s+ds)] = exp[eta(s)]*exp[Omega_h(s)]
// When the norm of the hunter gets to a threshold set by omega_norm_max,
// we gather it, also using the BCH formula
// exp[Omega_g] = exp[Omega_h] * exp[Omega_g].
// We then clear out the hunter, and update H(s) according to the gathered Omega.
// And we continue hunting, using H(s) = exp[Omega_g] H(0) exp[-Omega_g] as our starting point.
// This is essentially a compromise between the "split" and "no split" approaches, combining
// the advantages of both. As we hunt, we have a relatively small Omega, so we don't need to
// evaluate as many nested commutators, but in the end we have just one Omega so that if
// we want to transform some operator we don't need to do a bunch of transformations.
void IMSRGSolver::GatherOmega()
{
  std::cout << "gathering Omega. " << std::endl;
  if (Omega.size() < 2 ) 
  {
    auto &last = Omega.back();
    Omega.emplace_back(last);
  }
  // the last omega in the list is the hunter. the one just preceeding it is the gatherer.
  auto &hunter = Omega.back();
  auto &gatherer = Omega[Omega.size() - 2];
  if (hunter.Norm() > 1e-6)
  {
    gatherer = BCH::BCH_Product(hunter, gatherer);
  }
  hunter.Erase();
  H_saved = *H_0;
  for (size_t i = 0; i < Omega.size() - 1; i++)
  {
    H_saved = BCH::BCH_Transform(H_saved, Omega[i]);
  }
}

void IMSRGSolver::SetHin(Operator &H_in)
{
  modelspace = H_in.GetModelSpace();
  H_0 = &H_in;
  FlowingOps[0] = H_in;
  Eta = Operator(H_in);
  Eta.Erase();
  Eta.SetAntiHermitian();
  if (Omega.back().Norm() > 1e-6)
  {
    NewOmega();
  }
  else
  {
    Omega.back() = Eta;
  }
}

void IMSRGSolver::SetOmega(size_t i, Operator &om)
{
  if ((i + 1) > Omega.size())
  {
    Omega.resize(i + 1);
  }
  Omega[i] = om;
}

void IMSRGSolver::Reset()
{
  s = 0;
  Eta.Erase();
  Omega.resize(0);
  NewOmega();
}

void IMSRGSolver::SetGenerator(std::string gen)
{
  generator.SetType(gen);
  if (Omega.back().Norm() > 1e-6)
  {
    Eta.Erase();
    // NewOmega();   B.C. He
    if (hunter_gatherer)
    {
      GatherOmega();
    }
    else
    {
      NewOmega();
    }
  }
  if (magnus_adaptive)
  {
    ds = ds_0;
    in_soft_landing_phase_ = false;
  }
}

void IMSRGSolver::SetDenominatorPartitioning(std::string dp)
{
  generator.SetDenominatorPartitioning(dp);
  if (Omega.back().Norm() > 1e-6)
  {
    Eta.Erase();
    NewOmega();
  }
}

void IMSRGSolver::SetFlowFile(std::string str)
{
  flowfile = str;
  std::ofstream flowf;
  if (flowfile != "")
  {
    flowf.open(flowfile, std::ofstream::out);
    flowf.close();
  }
}

void IMSRGSolver::Solve()
{

  if (s < 1e-4)
    WriteFlowStatusHeader(std::cout);

  if (method == "magnus_euler" or method == "magnus")
    Solve_magnus_euler();
  else if (method == "magnus_backoff")
    Solve_magnus_backoff();
  else if (method == "magnus_modified_euler")
    Solve_magnus_modified_euler();
  else if (method == "flow_adaptive" or method == "flow")
    Solve_ode_adaptive();
  else if (method == "magnus_adaptive")
    Solve_ode_magnus();
  else if (method == "flow_euler")
    Solve_ode();
  else if (method == "flow_RK4")
    Solve_flow_RK4();
  else if (method == "restore_4th_order")
  {
    FlowingOps.emplace_back(Operator(*(FlowingOps[0].GetModelSpace()), 0, 0, 0, 1));
    Solve_ode_adaptive();
  }
  else
    std::cout << "IMSRGSolver: I don't know method " << method << std::endl;
}

void IMSRGSolver::UpdateEta()
{
  generator.Update(FlowingOps[0], Eta);
}

// This is the default solver
void IMSRGSolver::Solve_magnus_euler()
{
  istep = 0;

  generator.Update(FlowingOps[0], Eta);

  // SRS noticed this on June 12 2024. If these two parameters are equal, and especially if we're using the hunter-gatherer mode, then we become sensitive to
  // numerical precision when deciding if we should split omega, leading to machine-dependent behavior.
  if ( std::abs( omega_norm_max - norm_domega)<1e-6 )
  {
     norm_domega += 1e-4;
     std::cout << __func__ << ":  adjusting norm_domega to " << norm_domega << "  to avoid numerical trouble, since omega_norm_max = " << omega_norm_max << std::endl;
  }

  Elast = H_0->ZeroBody;
  cumulative_error = 0;
  // Write details of the flow
  WriteFlowStatus(flowfile);
  WriteFlowStatus(std::cout);

  for (istep = 1; s < smax; ++istep)
  {

    double norm_eta = Eta.Norm();
    if (norm_eta < eta_criterion)
    {
      break;
    }
    if (norm_eta > 1e12 or std::abs(Elast) > 1e9) // This is obviously going nowhere...
    {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "!!!!!!!!!!!  Norm of eta is " << norm_eta << " E0 = " << Elast << "  things are clearly broken. Giving up." << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      FlowingOps[0] *= 1.0 / 0.0;
      break;
    }
    double norm_omega = Omega.back().Norm();
    if (norm_omega > omega_norm_max)
    {
      if (hunter_gatherer)
      {
        GatherOmega();
      }
      else
      {
        NewOmega();
      }
      norm_omega = 0;
    }
    // ds should never be more than 1, as this is over-rotating
    // Also, since we check if ||Omega|| < omega_norm_max, if we choose ds so that ||Omega|| = omega_norm_max, then we become sensitive
    // to numerical precision details when evaluating the inequality and behavior becomes machine dependent. So we add 1e-5 to the omega_norm_max
    // option to ensure that we're definitely on one side of the inequality.
    if (magnus_adaptive)
      ds = std::min( { norm_domega / norm_eta,   norm_domega / norm_eta / (norm_omega + 1.0e-9),    (omega_norm_max+1e-5) / norm_eta, ds_max   });
    ds = std::min(ds, smax - s);

    s += ds;
    Eta *= ds; // Here's the Euler step.

    // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
    Omega.back() = BCH::BCH_Product(Eta, Omega.back());

    // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
    if ((Omega.size() + n_omega_written) < 2)
    {
      FlowingOps[0] = BCH::BCH_Transform(*H_0, Omega.back());
    }
    else
    {
      FlowingOps[0] = BCH::BCH_Transform(H_saved, Omega.back());
    }

    if (norm_eta < 1.0 and generator.GetType() == "shell-model-atan")
    {
      generator.SetDenominatorCutoff(1e-6);
    }

    generator.Update(FlowingOps[0], Eta);

    // Write details of the flow
    WriteFlowStatus(flowfile);
    WriteFlowStatus(std::cout);
    Elast = FlowingOps[0].ZeroBody;
  }
}

/// Modification added by Matthias
void IMSRGSolver::Solve_magnus_backoff()
{
  istep = 0;

  generator.Update(FlowingOps[0], Eta);

  Elast = H_0->ZeroBody;
  cumulative_error = 0;
  // Write details of the flow
  WriteFlowStatus(flowfile);
  WriteFlowStatus(std::cout);

  for (istep = 1; s < smax; ++istep)
  {
    double norm_eta = Eta.Norm();
    // Factor 20.0 is somewhat arbitrary, but this should be relative to eta_criterion
    if ((!in_soft_landing_phase_) && (std::abs(norm_eta) < eta_criterion * 20.0))
    {
      in_soft_landing_phase_ = true;
      std::cout << "Entering soft landing phase.\n";
    }

    if (norm_eta < eta_criterion)
    {
      break;
    }
    while (magnus_adaptive && (ds * norm_eta > M_PI_4))
    {
      double new_ds = ds * ds_backoff_factor_;
      std::cout << "Backing off ds because ds * norm_eta > pi / 4: ds = " << ds << " -> " << new_ds << "\n";
      ds = new_ds;
    }
    double norm_omega = Omega.back().Norm();
    if (norm_omega > omega_norm_max)
    {
      //               if ( perturbative_triples )
      //               {
      ////                 GetPerturbativeTriples();
      //                 pert_triples_this_omega = GetPerturbativeTriples();
      //                 pert_triples_sum += pert_triples_this_omega;
      //               }
      if (hunter_gatherer)
      {
        GatherOmega();
      }
      else
      {
        NewOmega();
      }
      norm_omega = 0;
    }
    // ds should never be more than 1, as this is over-rotating
    // if (magnus_adaptive)
    //    ds = std::min( std::min( std::min(norm_domega/norm_eta, norm_domega /
    //    norm_eta / (norm_omega+1.0e-9)), omega_norm_max/norm_eta), ds_max);
    ds = std::min(ds, smax - s);

    s += ds;
    Eta *= ds; // Here's the Euler step.

    // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) *
    // exp(Omega_last)
    Omega.back() = BCH::BCH_Product(Eta, Omega.back());

    // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
    if ((Omega.size() + n_omega_written) < 2)
    {
      FlowingOps[0] = BCH::BCH_Transform(*H_0, Omega.back());
    }
    else
    {
      FlowingOps[0] = BCH::BCH_Transform(H_saved, Omega.back());
    }

    //     if (norm_eta < 1.0 and generator.GetType() == "shell-model-atan") {
    //       generator.SetDenominatorCutoff(1e-6);
    //     }

    //      if ( generator.GetType() == "rspace" ) {
    //      generator.SetRegulatorLength(s); };
    //      generator.Update(&FlowingOps[0],&Eta);
    generator.Update(FlowingOps[0], Eta);
    //      cumulative_error += EstimateStepError();

    if (magnus_adaptive)
    {
      if (Eta.Norm() > norm_eta)
      {
        double ds_new = ds * ds_backoff_factor_;
        std::cout << "Backing off ds because of Eta norm growth, ds = " << ds << " -> " << ds_new
                  << "\n";
        ds = ds_new;
      }
      else
      {
        if (!in_soft_landing_phase_)
        {
          double ds_new = std::min(ds_max, ds * ds_max_growth_factor_);
          if (ds_new > ds)
          {
            std::cout << "Growing ds, ds = " << ds << " -> " << ds_new << "\n";
            ds = ds_new;
          }
        }
      }
    }

    // Write details of the flow
    WriteFlowStatus(flowfile);
    WriteFlowStatus(std::cout);
    Elast = FlowingOps[0].ZeroBody;

    //      Operator DHDS = Commutator::Commutator( Eta, FlowingOps[0] );
    //      size_t ch =
    //      FlowingOps[0].modelspace->GetTwoBodyChannelIndex(0,0,1); auto MAT
    //      = FlowingOps[0].TwoBody.GetMatrix(ch,ch); auto dMAT  =
    //      DHDS.TwoBody.GetMatrix(ch,ch); std::cout << " ======>  " <<
    //      MAT(0,0) << "     " << MAT(1,0) << "   " << MAT(2,0) << "  " <<
    //      MAT(1,1) << " " << MAT(2,1) << "  " << MAT(2,2) << std::endl;
    //      std::cout << " ______>  " << dMAT(0,0) << "     " << dMAT(1,0) << "
    //      " << dMAT(2,0) << "  " << dMAT(1,1) << " " << dMAT(2,1) << "  " <<
    //      dMAT(2,2) << std::endl;
  }
}

void IMSRGSolver::Solve_magnus_modified_euler()
{
  istep = 0;
  //   generator.Update(&FlowingOps[0],&Eta);
  generator.Update(FlowingOps[0], Eta);

  Operator H_temp;
  // Write details of the flow
  WriteFlowStatus(flowfile);
  WriteFlowStatus(std::cout);

  for (istep = 1; s < smax; ++istep)
  {
    double norm_eta = Eta.Norm();
    double norm_omega = Omega.back().Norm();
    if (norm_omega > omega_norm_max)
    {
      NewOmega();
      norm_omega = 0;
    }
    // ds should never be more than 1, as this is over-rotating
    ds = std::min(std::min(norm_domega / norm_eta, norm_domega / norm_eta / (norm_omega + 1.0e-9)), ds_max);
    if (s + ds > smax)
      ds = smax - s;
    s += ds;

    H_temp = FlowingOps[0] + ds * Commutator::Commutator(Eta, FlowingOps[0]);
    //      generator.AddToEta(&H_temp,&Eta);
    generator.AddToEta(H_temp, Eta);

    Eta *= ds * 0.5; // Here's the modified Euler step.

    // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
    //      Omega.back() = Eta.BCH_Product( Omega.back() );
    Omega.back() = BCH::BCH_Product(Eta, Omega.back());

    if ((Omega.size() + n_omega_written) < 2)
    {
      //        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
      FlowingOps[0] = BCH::BCH_Transform(*H_0, Omega.back());
    }
    else
    {
      //        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
      FlowingOps[0] = BCH::BCH_Transform(H_saved, Omega.back());
    }

    //      generator.Update(&FlowingOps[0],&Eta);
    generator.Update(FlowingOps[0], Eta);

    // Write details of the flow
    WriteFlowStatus(flowfile);
    WriteFlowStatus(std::cout);
  }
}

// Solve with fixed-step 4th-order Runge-Kutta
void IMSRGSolver::Solve_flow_RK4()
{
  istep = 0;

  //   if ( generator.GetType() == "rspace" ) { generator.modelspace = (Eta.modelspace); generator.SetRegulatorLength(800005.0); };

  //   generator.Update(&FlowingOps[0],&Eta);
  generator.Update(FlowingOps[0], Eta);

  if (generator.GetType() == "shell-model-atan")
  {
    generator.SetDenominatorCutoff(1.0); // do we need this?
  }

  Elast = H_0->ZeroBody;
  cumulative_error = 0;
  // Write details of the flow
  WriteFlowStatus(flowfile);
  WriteFlowStatus(std::cout);

  Operator goosetank_chi(*modelspace, 0, 0, 0, 1);  // for use if we do IMSRG2* to mock up the goose tanks
  Operator goosetank_dchi(*modelspace, 0, 0, 0, 1); // for use if we do IMSRG2* to mock up the goose tanks

  for (istep = 1; s < smax; ++istep)
  {

    double norm_eta = Eta.Norm();
    if (norm_eta < eta_criterion)
    {
      break;
    }

    ds = std::min(ds_max, smax - s);
    s += ds;

    int nops = FlowingOps.size();
    std::vector<Operator> K1(nops);
    std::vector<Operator> K2(nops);
    std::vector<Operator> K3(nops);
    std::vector<Operator> K4(nops);
    std::vector<Operator> Ktmp(nops);

    //      Operator& Hs = FlowingOps[0];   // this is not used explicitly
    for (int i = 0; i < nops; i++)
    {
      if (i == 0)
        K1[i] = Commutator::Commutator(Eta, FlowingOps[i] + goosetank_chi);
      else
        K1[i] = Commutator::Commutator(Eta, FlowingOps[i]);
      Ktmp[i] = FlowingOps[i] + 0.5 * ds * K1[i];
    }
    //      Operator K1 = Commutator::Commutator( Eta, Hs );
    //      Operator Htmp = Hs + 0.5*ds*K1[0];
    //      generator.Update(&Htmp,&Eta);
    //      generator.Update(&Ktmp[0],&Eta);
    generator.Update(Ktmp[0], Eta);
    for (int i = 0; i < nops; i++)
    {
      //        K2[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
      if (i == 0)
        K2[i] = Commutator::Commutator(Eta, Ktmp[i] + goosetank_chi);
      else
        K2[i] = Commutator::Commutator(Eta, Ktmp[i]);
      Ktmp[i] = FlowingOps[i] + 0.5 * ds * K2[i];
    }
    //      Operator K2 = Commutator::Commutator( Eta, Hs+Htmp );
    //      Htmp = Hs + 0.5*ds*K2;
    //      generator.Update(&Htmp,&Eta);
    //      generator.Update(&Ktmp[0],&Eta);
    generator.Update(Ktmp[0], Eta);
    for (int i = 0; i < nops; i++)
    {
      //        K3[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
      if (i == 0)
        K3[i] = Commutator::Commutator(Eta, Ktmp[i] + goosetank_chi);
      else
        K3[i] = Commutator::Commutator(Eta, Ktmp[i]);
      Ktmp[i] = FlowingOps[i] + 1.0 * ds * K3[i];
    }
    //      Operator K3 = Commutator::Commutator( Eta, Hs+Htmp );
    //      Htmp = Hs + 1.0*ds*K3;
    //      generator.Update(&Htmp,&Eta);
    //      generator.Update(&Ktmp[0],&Eta);
    generator.Update(Ktmp[0], Eta);
    for (int i = 0; i < nops; i++)
    {
      //        K4[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
      if (i == 0)
        K4[i] = Commutator::Commutator(Eta, Ktmp[i] + goosetank_chi);
      else
        K4[i] = Commutator::Commutator(Eta, Ktmp[i]);
      //        Ktmp[i] = FlowingOps[i] + 1.0*ds*K2[i];
      FlowingOps[i] += ds / 6.0 * (K1[i] + 2 * K2[i] + 2 * K3[i] + K4[i]);
    }
    //      Operator K4 = Commutator::Commutator( Eta, Hs+Htmp );
    //      Hs += ds/6.0 * ( K1 + 2*K2 + 2*K3 + K4);

    if (norm_eta < 1.0 and generator.GetType() == "shell-model-atan")
    {
      generator.SetDenominatorCutoff(1e-6);
    }

    if (BCH::use_goose_tank_correction)
    {
      goosetank_dchi.EraseOneBody();
      Commutator::comm221ss(Eta, FlowingOps[0], goosetank_dchi); // update chi.
      for (auto i : modelspace->all_orbits)                      // enforce n_in_j + nbar_i nbar_j
      {
        Orbit &oi = modelspace->GetOrbit(i);
        for (auto j : modelspace->all_orbits)
        {
          Orbit &oj = modelspace->GetOrbit(j);
          goosetank_dchi.OneBody(i, j) *= oi.occ * oj.occ + (1.0 - oi.occ) * (1.0 - oj.occ);
        }
      }
      goosetank_chi += goosetank_dchi * ds;
      //        std::cout << " " << __FILE__ << "  line " << __LINE__ << s << "  " << goosetank_chi.OneBody(1,1) << std::endl;;
    }

    //      if ( generator.GetType() == "rspace" ) { generator.SetRegulatorLength(s); };
    //      generator.Update(&FlowingOps[0],&Eta);
    generator.Update(FlowingOps[0], Eta);
    //      cumulative_error += EstimateStepError();

    // Write details of the flow
    WriteFlowStatus(flowfile);
    WriteFlowStatus(std::cout);
    //      profiler.PrintMemory();
    Elast = FlowingOps[0].ZeroBody;
  }
}

#ifndef NO_ODE

// Implement element-wise division and abs and reduce for Operators.
// This is required for adaptive steppers
// vector<Operator> operator/ (const vector<Operator>& num, const vector<Operator>& denom)
std::deque<Operator> operator/(const std::deque<Operator> &num, const std::deque<Operator> &denom)
{
  //   vector<Operator> quotient = num;
  std::deque<Operator> quotient = num;
  for (size_t i = 0; i < num.size(); ++i)
  {
    quotient[i].ZeroBody /= denom[i].ZeroBody;
    quotient[i].OneBody /= denom[i].OneBody;
    for (auto &itmat : quotient[i].TwoBody.MatEl)
      itmat.second /= denom[i].TwoBody.GetMatrix(itmat.first[0], itmat.first[1]);
  }
  return quotient;
}

// vector<Operator> operator* (const double a, const vector<Operator>& X)
std::deque<Operator> operator*(const double a, const std::deque<Operator> &X)
{
  //  vector<Operator> Y = X;
  std::deque<Operator> Y = X;
  for (auto &y : Y)
    y *= a;
  return Y;
}

// vector<Operator> operator+ ( const vector<Operator>& X, const vector<Operator>& Y)
std::deque<Operator> operator+(const std::deque<Operator> &X, const std::deque<Operator> &Y)
{
  //  vector<Operator> Z = X;
  std::deque<Operator> Z = X;
  for (size_t i = 0; i < Z.size(); ++i)
    Z[i] += Y[i];
  return Z;
}

// Also need the dubious operation of adding a double to an operator.
// vector<Operator> operator+ (const double a, const vector<Operator>& X)
std::deque<Operator> operator+(const double a, const std::deque<Operator> &X)
{
  //   vector<Operator> Y = X;
  std::deque<Operator> Y = X;
  for (auto &y : Y)
  {
    y.ZeroBody += a;
    y.OneBody += a;
    //     for( auto& v : y.OneBody ) v += a;
    for (auto &itmat : y.TwoBody.MatEl)
      itmat.second += a;
    //       for ( auto& v : itmat.second ) v += a;
  }
  return Y;
}

// Return the element-wise absolute value of an operator
// this is needed for ODE adaptive solver
// vector<Operator> abs(const vector<Operator>& OpIn)
std::deque<Operator> abs(const std::deque<Operator> &OpIn)
{
  //   vector<Operator> OpOut = OpIn;
  std::deque<Operator> OpOut = OpIn;
  for (auto &opout : OpOut)
  {
    opout.ZeroBody = std::abs(opout.ZeroBody);
    opout.OneBody = arma::abs(opout.OneBody);
    for (auto &itmat : opout.TwoBody.MatEl)
      itmat.second = arma::abs(itmat.second);
  }
  return OpOut;
}

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION < 1.56
#ifdef OLD_BOOST
namespace boost
{
  namespace numeric
  {
    namespace odeint
    {
      template <>
      struct vector_space_reduce<std::deque<Operator>>
      {
        template <class Op>
        double operator()(const std::deque<Operator> &X, Op op, double init)
        {
          for (auto &x : X)
          {
            init = op(init, x.ZeroBody);
            for (auto &v : x.OneBody)
              init = op(init, v);
            for (auto &itmat : x.TwoBody.MatEl)
            {
              for (auto &v : itmat.second)
                init = op(init, v);
            }
          }
          return init;
        }
      };
    }
  }
}
#endif

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION >= 1.56
// struct vector_space_norm_inf< vector<Operator> >
#ifndef OLD_BOOST
namespace boost
{
  namespace numeric
  {
    namespace odeint
    {
      template <>
      struct vector_space_norm_inf<std::deque<Operator>>
      {
        typedef double result_type;
        //   double operator()(const vector<Operator>& X)
        double operator()(const std::deque<Operator> &X)
        {
          double norm = 0;
          for (auto &x : X)
            norm += x.Norm();
          return norm;
        }
      };
    }
  }
}
#endif

void IMSRGSolver::Solve_ode()
{

  ode_mode = "H";
  WriteFlowStatusHeader(std::cout);
  WriteFlowStatus(flowfile);
  //   using namespace boost::numeric::odeint;
  namespace odeint = boost::numeric::odeint;
  //   runge_kutta4< vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
  odeint::runge_kutta4<std::deque<Operator>, double, std::deque<Operator>, double, odeint::vector_space_algebra> stepper;
  auto system = *this;
  auto monitor = ode_monitor;
  //   size_t steps = integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
  odeint::integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
  monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
  ode_mode = "H";
  if (method == "restore_4th_order")
    ode_mode = "Restored";
  WriteFlowStatusHeader(std::cout);
  WriteFlowStatus(flowfile);
  std::cout << "done writing header and status" << std::endl;
  //   using namespace boost::numeric::odeint;
  namespace odeint = boost::numeric::odeint;
  auto system = *this;
  //   typedef runge_kutta_dopri5< vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
  typedef odeint::runge_kutta_dopri5<std::deque<Operator>, double, std::deque<Operator>, double, odeint::vector_space_algebra> stepper;
  //   typedef adams_bashforth_moulton< 4, vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
  auto monitor = ode_monitor;
  //   size_t steps = integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
  odeint::integrate_adaptive(odeint::make_controlled<stepper>(ode_e_abs, ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
  monitor.report();
}

// Evaluate dx/dt for boost ode
// void IMSRGSolver::operator()( const Operator& x, Operator& dxdt, const double t)
// void IMSRGSolver::operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t)
void IMSRGSolver::operator()(const std::deque<Operator> &x, std::deque<Operator> &dxdt, const double t)
{
  s = t;
  if (ode_mode == "H")
  {
    FlowingOps[0] = x[0];
    if (dxdt.size() < x.size())
      dxdt.resize(x.size());
    auto &H_s = FlowingOps[0];
    //     generator.Update(&H_s,&Eta);
    generator.Update(H_s, Eta);
    if (Eta.Norm() < eta_criterion)
    {
      for (size_t i = 0; i < x.size(); ++i)
      {
        dxdt[i] = 0 * x[i];
      }
    }
    else
    {
      for (size_t i = 0; i < x.size(); ++i)
      {
        dxdt[i] = Commutator::Commutator(Eta, x[i]);
      }
    }
  }
  else if (ode_mode == "Omega")
  {

    double norm_omega = Omega.back().Norm();
    if (norm_omega > omega_norm_max) // This doesn't seem to works so well just yet...
    {
      NewOmega();
      norm_omega = 0;
    }
    Omega.back() = x.back();
    auto &Omega_s = x.back();
    Operator &H_s = FlowingOps[0];
    if ((Omega.size() + n_omega_written) > 1)
      H_s = BCH::BCH_Transform(H_saved, Omega_s);
    //       H_s = H_saved.BCH_Transform(Omega_s);
    else
      H_s = BCH::BCH_Transform(*H_0, Omega_s);
    //       H_s = H_0->BCH_Transform(Omega_s);
    //     generator.Update(&H_s,&Eta);
    generator.Update(H_s, Eta);
    if (dxdt.size() < x.size())
      dxdt.resize(x.size());
    dxdt.back() = Eta - 0.5 * Commutator::Commutator(Omega_s, Eta);
  }
  else if (ode_mode == "Restored")
  {
    FlowingOps[0] = x[0];
    FlowingOps[1] = x[1];
    if (dxdt.size() < x.size())
      dxdt.resize(x.size());
    dxdt[1] = Operator(x[1]);
    auto &H_s = FlowingOps[0];
    //     generator.Update(&H_s,&Eta);
    generator.Update(H_s, Eta);
    if (Eta.Norm() < eta_criterion)
    {
      for (size_t i = 0; i < x.size(); ++i)
      {
        dxdt[i] = 0 * x[i];
      }
    }
    else
    {
      dxdt[0] = Commutator::Commutator(Eta, x[0] + x[1]);
      dxdt[1].Erase();
      //       dxdt[1].comm221ss(Eta,x[0]);
      Commutator::comm221ss(Eta, x[0], dxdt[1]);
      // keep only pp and hh parts of d chi/ ds
      for (auto &a : modelspace->holes)
      {
        for (auto &i : modelspace->particles)
        {
          dxdt[1].OneBody(a, i) = 0;
          dxdt[1].OneBody(i, a) = 0;
        }
      }
      for (size_t i = 2; i < x.size(); ++i)
      {
        dxdt[i] = Commutator::Commutator(Eta, x[i]);
      }
    }
  }
  WriteFlowStatus(flowfile);
  WriteFlowStatus(std::cout);
}

void IMSRGSolver::Solve_ode_magnus()
{
  ode_mode = "Omega";
  WriteFlowStatus(std::cout);
  WriteFlowStatus(flowfile);
  //   using namespace boost::numeric::odeint;
  namespace odeint = boost::numeric::odeint;
  namespace pl = std::placeholders;
  //   runge_kutta4<vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
  odeint::runge_kutta4<std::deque<Operator>, double, std::deque<Operator>, double, odeint::vector_space_algebra> stepper;
  auto system = *this;
  auto monitor = ode_monitor;
  //   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
  odeint::integrate_const(stepper, system, Omega, s, smax, ds, monitor);
  monitor.report();
}

#endif

/// Returns \f$ e^{Omega} \mathcal{O} e^{-Omega} \f$
Operator IMSRGSolver::Transform(Operator &OpIn)
{
  return Transform_Partial(OpIn, 0);
}

Operator IMSRGSolver::Transform(Operator &&OpIn)
{
  return Transform_Partial(OpIn, 0);
}

/// Returns \f$ e^{-Omega} \mathcal{O} e^{Omega} \f$
Operator IMSRGSolver::InverseTransform(Operator &OpIn)
{
  //  if (OpIn.GetJRank()+OpIn.GetTRank()+OpIn.GetParity()>0)
  //  {
  //    OpIn.ResetTensorTransformFirstPass();
  //  }
  Operator OpOut = OpIn;
  for (auto omega = Omega.rbegin(); omega != Omega.rend(); ++omega)
  {
    Operator negomega = -(*omega);
    //    OpOut = OpOut.BCH_Transform( negomega );
    OpOut = BCH::BCH_Transform(OpOut, negomega);
  }
  return OpOut;
}

/// Returns \f$ e^{\Omega} \mathcal{O} e^{-\Omega} \f$
/// for the \f$\Omega_i\f$s with index greater than or equal to n.
Operator IMSRGSolver::Transform_Partial(Operator &OpIn, int n)
{
  Operator OpOut = OpIn;
  if (OpOut.GetParticleRank() == 1)
    OpOut.SetParticleRank(2);

  //  if ((rw != NULL) and rw->GetScratchDir() != "")
  if (scratchdir != "")
  {
    //    char tmp[512];
    for (int i = n; i < n_omega_written; i++)
    {
      //    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
      //     std::string fname(tmp);
      //    Operator omega(OpIn);
      Operator omega(Eta);
      std::ostringstream filename;
      //     filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
      filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
      std::cout << "Transforming using " << filename.str() << std::endl;
      std::ifstream ifs(filename.str(), std::ios::binary);
      omega.ReadBinary(ifs);
      //     if (OpIn.GetJRank()>0) cout << "step " << i << endl;
      //     OpOut = OpOut.BCH_Transform( omega );
      OpOut = BCH::BCH_Transform(OpOut, omega);
      std::cout << "norm of omega = " << omega.Norm() << std::endl;
      std::cout << " op zero body = " << OpOut.ZeroBody << std::endl;
      //     if (OpIn.GetJRank()>0)cout << "done" << endl;
    }
  }

  for (size_t i = std::max(n - n_omega_written, 0); i < Omega.size(); ++i)
  {
    //     if (OpIn.GetJRank()>0) std::cout << "step " << i << endl;
    //     std::cout << "  line " << __LINE__ << "   i= " << i << std::endl;
    //     std::cout << " Omega.size() = " << Omega.size() << std::endl;
    //     std::cout << " norm of omega[i] = " << Omega[i].Norm() << std::endl;
    //     std::cout << " norm of op = " << OpOut.Norm() << std::endl;
    //     std::cout << " op zero body = " << OpOut.ZeroBody << std::endl;
    //    OpOut = OpOut.BCH_Transform( Omega[i] );
    OpOut = BCH::BCH_Transform(OpOut, Omega[i]);
    //     if (OpIn.GetJRank()>0)cout << "done" << endl;
  }

  return OpOut;
}

Operator IMSRGSolver::Transform_Partial(Operator &&OpIn, int n)
{
  //  cout << "Calling r-value version of Transform_Partial, n = " << n << endl;
  Operator OpOut = OpIn;
  //  if ((rw != NULL) and rw->GetScratchDir() != "")
  if (scratchdir != "")
  {
    //    char tmp[512];
    for (int i = n; i < n_omega_written; i++)
    {
      //     sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
      //     std::string fname(tmp);
      //     Operator omega(OpIn);
      Operator omega(Eta);
      std::ostringstream filename;
      //     filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
      filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
      std::ifstream ifs(filename.str(), std::ios::binary);
      omega.ReadBinary(ifs);
      //     OpOut = OpOut.BCH_Transform( omega );
      OpOut = BCH::BCH_Transform(OpOut, omega);
    }
  }

  for (size_t i = std::max(n - n_omega_written, 0); i < Omega.size(); ++i)
  {
    //    OpOut = OpOut.BCH_Transform( Omega[i] );
    OpOut = BCH::BCH_Transform(OpOut, Omega[i]);
  }
  return OpOut;
}

// count number of equations to be solved
int IMSRGSolver::GetSystemDimension()
{
  int dim = 1; // zero-body part

  int N = H_0->OneBody.n_cols;
  dim += N * (N + 1) / 2;
  dim += H_0->TwoBody.Dimension();
  return dim;
}

void IMSRGSolver::FlushOmegaToScratch()
{
  if (scratchdir.find("/dev/null") != std::string::npos)
    return; // if we're writing to /dev/null, we'd better not try to erase things
  if (scratchdir == "")
    return;
  for (size_t i = 0; i < Omega.size(); i++)
  {
    std::stringstream filename;
    filename << scratchdir << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i + n_omega_written;
    std::ofstream ofs(filename.str(), std::ios::binary);
    Omega[i].WriteBinary(ofs);
    n_omega_written++;
    std::cout << "Omega written to file " << filename.str() << "  written " << i + n_omega_written << " so far." << std::endl;
  }
}

//// It's important that we don't refer to rw in this routine, since we call on destruction
//// and it may be that rw gets destroyed before this instance, making the pointer to rw invalid.
void IMSRGSolver::CleanupScratch()
{
  if (n_omega_written <= 0)
    return;
  //  if ( rw->GetScratchDir().find("/dev/null") != std::string::npos ) return; // if we're writing to /dev/null, we'd better not try to erase things
  if (scratchdir.find("/dev/null") != std::string::npos)
    return; // if we're writing to /dev/null, we'd better not try to erase things
            //  std::cout << "Cleaning up files written to scratch space" << std::endl;
            //  char tmp[512];
  for (int i = 0; i < n_omega_written; i++)
  {
    std::ostringstream filename;
    //    filename << rw->GetScratchDir() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
    filename << scratchdir << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
    //    filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
    //    std::cout << "  scratchdir:  " << rw->GetScratchDir() << "  ,   " << rw->GetScratchDir().c_str() << std::endl;
    //    std::cout << "  scratchdir:  " << scratchdir << "  ,   " << scratchdir.c_str() << std::endl;
    //    std::cout << "  pid :  " << std::setw(6) << std::setfill('0') << getpid() << std::endl;
    //    std::cout << "    i :  " << std::setw(3) << std::setfill('0') << i << std::endl;
    //    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
    //    std::string fname(tmp);
    //    if ( remove(tmp) !=0 )
    if (remove(filename.str().c_str()) != 0)
    {
      std::cout << "Error when attempting to delete " << filename.str() << std::endl;
    }
  }
}

// This doesn't really work all that well. Probably shouldn't use it.
double IMSRGSolver::EstimateStepError()
{
  double err = 0;
  double Ecurrent = FlowingOps[0].ZeroBody;
  std::string gtype = generator.generator_type;
  //  std::cout << "gtype = " << gtype << std::endl;
  if (std::abs(s) > 1e-9 and gtype.find("qtransfer") != std::string::npos)
  {
    int n;
    std::istringstream(gtype.substr(gtype.find("_") + 1)) >> n;
    double dEds = (Ecurrent - Elast);
    double kF = 1.33;
    double R = pow(s, 1.0 / n);
    //    double suppression = std::min( (modelspace->GetTargetMass()-2.)/3.0, pow( (kF*R), 3) ); // Not sure if this is actually legit
    double suppression = pow((kF * R), 3);
    err = std::abs(dEds) * suppression;
    //    std::cout << std::endl << "Made it inside, suppression = " << suppression << "  dEds = " << dEds << " ->  "
    //              << Elast << " " << Ecurrent << "  err = " << err << std::endl;
  }
  //  Elast = Ecurrent;
  return err;
}

double IMSRGSolver::EstimateBCHError()
{
  double err = 0;
  int counter = 0;
  for (auto &omega : Omega)
  {
    err += BCH::EstimateBCHError(omega, *H_0);
    std::cout << counter++ << "  " << err << std::endl;
  }
  return err;
}

/// Compute the perturbative triples correction to the energy
/// \f[ \Delta E_3 = \sum |W_{ijkabc}|^2/\Delta_{ijkabc} \f]
/// where
/// \f[ W=[\Omega,H]_3b \f]
double IMSRGSolver::CalculatePerturbativeTriples()
{

  Operator Wbar((*modelspace), 0, 0, 0, 2);
  // Wbar.ThreeBody.SetMode("pn");  // Dont do this. It automatically allocates and we don't want that.

  // If we've split the Omegas up, we combine them here, implicitly summing the 3N generated by each.
  Operator omega = Omega[0];
//  for (size_t n=1; n<Omega.size(); n++) omega += Omega[n];
  for (size_t n=1; n<Omega.size(); n++) omega = BCH::BCH_Product(omega,Omega[n]);

  Operator &Hs = FlowingOps[0];

  BCH::SetBCHSkipiEq1(true);
  Operator Htilde = Transform(*H_0);
  BCH::SetBCHSkipiEq1(false);

  // We double the first commutator account for [O,[O,H]_3]_3 diagrams which produce identical Wod and which we would otherwise miss
  // Leave this off for now.
//  Htilde.TwoBody  += 0.5*Commutator::Commutator(omega,*H_0).TwoBody ;

  // Need to put the one-body part of H into Wbar so we can get the denominators. I'm not sure this is the best way to do that...
  Wbar.OneBody = Hs.OneBody;
  Wbar.TwoBody = Hs.TwoBody;

  Commutator::perturbative_triples = true;
  Commutator::comm223ss(omega, Htilde, Wbar);
  Commutator::perturbative_triples = false; // turn it back off in case we want to do any more transformations

  double pert_triples = Wbar.ZeroBody;

  return pert_triples;
}

double IMSRGSolver::CalculatePerturbativeTriples(Operator &Op_0)
{
  Operator Wbar((*modelspace), 0, 0, 0, 2);
  // Wbar.ThreeBody.SetMode("pn");  // Dont do this. It automatically allocates and we don't want that.

  // If we've split the Omegas up, we combine them here, implicitly summing the 3N generated by each.
  Operator omega = Omega[0];
//  for (size_t n=1; n<Omega.size(); n++) omega += Omega[n];
  for (size_t n=1; n<Omega.size(); n++) omega = BCH::BCH_Product(omega,Omega[n]);

  Operator &Hs = FlowingOps[0];

  BCH::SetBCHSkipiEq1(true);
  Operator Otilde = Transform(Op_0);
  BCH::SetBCHSkipiEq1(false);

  // Need to put the one-body part of H into Wbar so we can get the denominators. I'm not sure this is the best way to do that...
  Wbar.OneBody = Hs.OneBody;
  Wbar.TwoBody = Hs.TwoBody;

  Commutator::perturbative_triples = true;
  Commutator::comm223ss(omega, Otilde, Wbar);
  Commutator::perturbative_triples = false; // turn it back off in case we want to do any more transformations

  double pert_triples = Wbar.ZeroBody;

  return pert_triples;
}

void IMSRGSolver::WriteFlowStatus(std::string fname)
{
  if (fname != "")
  {
    std::ofstream ff(fname, std::ios::app);
    WriteFlowStatus(ff);
  }
}
void IMSRGSolver::WriteFlowStatus(std::ostream &f)
{
  if (f.good())
  {
    int fwidth = 16;
    int fprecision = 9;
    auto &H_s = FlowingOps[0];
    f.setf(std::ios::fixed);
    f << std::fixed << std::setw(5) << istep
      << std::setw(12) << std::setprecision(5) << s
      << std::setw(fwidth) << std::setprecision(fprecision) << H_s.ZeroBody
      << std::setw(fwidth) << std::setprecision(fprecision) << H_s.Norm()
      << std::setw(fwidth) << std::setprecision(fprecision) << cumulative_error
      << std::setw(fwidth) << std::setprecision(fprecision) << Omega.back().OneBodyNorm()
      << std::setw(fwidth) << std::setprecision(fprecision) << Omega.back().TwoBodyNorm()
      << std::setw(fwidth) << std::setprecision(fprecision) << Omega.back().ThreeBodyNorm()
      << std::setw(fwidth) << std::setprecision(fprecision) << Eta.OneBodyNorm()
      << std::setw(fwidth) << std::setprecision(fprecision) << Eta.TwoBodyNorm()
      << std::setw(fwidth) << std::setprecision(fprecision) << Eta.ThreeBodyNorm()
      << std::setw(7) << std::setprecision(0) << profiler.counter["N_ScalarCommutators"] + profiler.counter["N_TensorCommutators"]
      << std::setw(fwidth) << std::setprecision(fprecision) << H_s.GetMP2_Energy()
      << std::setw(7) << std::setprecision(0) << profiler.counter["N_Operators"]
      << std::setprecision(fprecision)
      << std::setw(12) << std::setprecision(3) << profiler.GetTimes()["real"]
      << std::setw(12) << std::setprecision(3) << profiler.CheckMem()["RSS"] / 1024. << " / " << std::skipws << profiler.MaxMemUsage() / 1024. << std::fixed
      << std::endl;
  }
}

void IMSRGSolver::WriteFlowStatusHeader(std::string fname)
{
  std::ofstream ff;
  if (fname != "")
    ff.open(fname, std::ios::app);
  WriteFlowStatusHeader(ff);
}
void IMSRGSolver::WriteFlowStatusHeader(std::ostream &f)
{
  if (f.good())
  {
    int fwidth = 16;
    int fprecision = 9;
    f.setf(std::ios::fixed);
    f << std::fixed << std::setw(5) << "i"
      << std::setw(12) << std::setprecision(5) << "s"
      << std::setw(fwidth) << std::setprecision(fprecision) << "E0"
      //        << std::setw(fwidth) << std::setprecision(fprecision) << "||H_1||"
      //        << std::setw(fwidth) << std::setprecision(fprecision) << "||H_2||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||H||"
      //        << std::setw(fwidth) << std::setprecision(fprecision) << "Tr(H)/Tr(1)"
      << std::setw(fwidth) << std::setprecision(fprecision) << " estim. err"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Omega_1||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Omega_2||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Omega_3||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_1||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_2||"
      << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_3||"
      << std::setw(7) << std::setprecision(fprecision) << "Ncomm"
      << std::setw(16) << std::setprecision(fprecision) << "E(MP2)"
      << std::setw(7) << std::setprecision(fprecision) << "N_Ops"
      << std::setw(16) << std::setprecision(fprecision) << "Walltime (s)"
      << std::setw(19) << std::setprecision(fprecision) << "Memory (MB)"
      << std::endl;
    for (int x = 0; x < 175; x++)
      f << "-";
    f << std::endl;
  }
}
