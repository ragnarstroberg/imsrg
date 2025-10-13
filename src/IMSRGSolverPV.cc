#include "IMSRGSolverPV.hh"
#include "Commutator.hh"
#include "BCH.hh"
// #include "Operator.hh"
// #include <algorithm>
// #include <cmath>
// #include <iomanip>
// #include <iostream>
// #include <string>
// #include <sstream>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif



IMSRGSolverPV::IMSRGSolverPV(Operator &H_in, Operator &VPT_in)
    : IMSRGSolver(H_in), VPT_0(&VPT_in), FlowingOpsPV(1, VPT_in), Etapv(VPT_in)                                                                                                                                                                                                                   
{   
    Etapv.Erase();
    //Note the Eta is already added to Omega deque in the IMSRGSolver Constructor
    if (VPT_in.IsHermitian())
    {
        Etapv.SetAntiHermitian();
    }
    else
    {
        Etapv.SetHermitian();
    }
    OmegaPV.emplace_back(Etapv);
}


void IMSRGSolverPV::SetFlowFilePV(std::string str)
{
      flowfile = str;
      std::ofstream flowf;
      if (flowfile != "")
      {
            flowf.open(flowfile, std::ofstream::out);
            flowf.close();
      }
}


void IMSRGSolverPV::Solve_flow_RK4_PV()
{
      istep = 0;
      
      generatorPV.Update(FlowingOps[0], FlowingOpsPV[0], Eta, Etapv);

      Elast = H_0->ZeroBody;
      cumulative_error = 0;
      // Write details of the flow
      WriteFlowStatusHeaderPV(std::cout); // commented beatriz 03/05/2025
      WriteFlowStatusPV(std::cout);
      for (istep = 1; s < smax; ++istep)
      {
            double norm_eta = Eta.Norm();
            double norm_etapv = Etapv.Norm();
            if (sqrt(norm_eta*norm_eta+norm_etapv*norm_etapv) < eta_criterion)
            {       
                  break;
            }

            ds = std::min(ds_max, smax - s);
            s += ds;
            int nops = FlowingOps.size();
            int nops_pv = FlowingOpsPV.size();
            if (nops != nops_pv)
            {
                std::cout<<"Problem, number of PC and PV operators isn't the same..."<<std::endl;
                exit(0);
            }
            std::vector<Operator> K1H(nops);
            std::vector<Operator> K2H(nops);
            std::vector<Operator> K3H(nops);
            std::vector<Operator> K4H(nops);
            std::vector<Operator> KtmpH(nops);
            std::vector<Operator> K1V(nops);
            std::vector<Operator> K2V(nops);
            std::vector<Operator> K3V(nops);
            std::vector<Operator> K4V(nops);
            std::vector<Operator> KtmpV(nops);
            
            for (int i = 0; i < nops; i++)
            {
                K1V[i] = Commutator::Commutator(Eta, FlowingOpsPV[i]) + Commutator::Commutator(Etapv, FlowingOps[i]);
                K1H[i] = Commutator::Commutator(Eta, FlowingOps[i]) + Commutator::Commutator(Etapv, FlowingOpsPV[i]);
                KtmpV[i] = FlowingOpsPV[i] + 0.5 * ds * K1V[i];
                KtmpH[i] = FlowingOps[i] + 0.5 * ds * K1H[i];
            }
            generatorPV.Update(KtmpH[0], KtmpV[0], Eta, Etapv);
            for (int i = 0; i < nops; i++)
            {
                K2V[i] = Commutator::Commutator(Eta, KtmpV[i]) + Commutator::Commutator(Etapv, KtmpH[i]);
                K2H[i] = Commutator::Commutator(Eta, KtmpH[i]) + Commutator::Commutator(Etapv, KtmpV[i]);
                KtmpV[i] = FlowingOpsPV[i] + 0.5 * ds * K2V[i];
                KtmpH[i] = FlowingOps[i] + 0.5 * ds * K2H[i];
            }
            generatorPV.Update(KtmpH[0], KtmpV[0], Eta, Etapv);
            for (int i = 0; i < nops; i++)
            {
                K3V[i] = Commutator::Commutator(Eta, KtmpV[i]) + Commutator::Commutator(Etapv, KtmpH[i]);
                K3H[i] = Commutator::Commutator(Eta, KtmpH[i]) + Commutator::Commutator(Etapv, KtmpV[i]);
                KtmpV[i] = FlowingOpsPV[i] + 1.0 * ds * K3V[i];
                KtmpH[i] = FlowingOps[i] + 1.0 * ds * K3H[i];
     
            }
            generatorPV.Update(KtmpH[0], KtmpV[0], Eta, Etapv);
            for (int i = 0; i < nops; i++)
            {
                K4V[i] = Commutator::Commutator(Eta, KtmpV[i]) + Commutator::Commutator(Etapv, KtmpH[i]);
                K4H[i] = Commutator::Commutator(Eta, KtmpH[i]) + Commutator::Commutator(Etapv, KtmpV[i]);
                FlowingOps[i] += ds / 6.0 * (K1H[i] + 2 * K2H[i] + 2 * K3H[i] + K4H[i]);
                FlowingOpsPV[i] += ds / 6.0 * (K1V[i] + 2 * K2V[i] + 2 * K3V[i] + K4V[i]);
            }
            // FlowingOps[0].PrintOneBody();
            generatorPV.Update(FlowingOps[0], FlowingOpsPV[0], Eta, Etapv);
            // Write details of the flow
            WriteFlowStatusPV(flowfile);
            WriteFlowStatusPV(std::cout);
            Elast = FlowingOps[0].ZeroBody;
      }
    //   FlowingOps[0].PrintOneBody();
}


void IMSRGSolverPV::NewOmega_PV()
{
    H_saved = FlowingOps[0];
    VPT_saved = FlowingOpsPV[0];
    std::cout << "pushing back another Omega. Omega.size = " << Omega.size()
                << " , operator size = " << Omega.front().Size() / 1024. / 1024. << " MB"
                << ",  memory usage = " << profiler.CheckMem()["RSS"] / 1024. / 1024. << " GB";
    std::cout << std::endl;
    std::cout << "pushing back another OmegaPV. OmegaPV.size = " << OmegaPV.size()
              << " , operator size = " << OmegaPV.front().Size() / 1024. / 1024. << " MB"
              << ",  memory usage = " << profiler.CheckMem()["RSS"] / 1024. / 1024. << " GB";
    std::cout << std::endl;
    if (scratchdir != "")
    {
        if (scratchdir.find("/dev/null") == std::string::npos)
        {
            std::ostringstream filename;
            filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << n_omega_written;
            std::ofstream ofs(filename.str(), std::ios::binary);
            Omega.back().WriteBinary(ofs);
            std::ostringstream filenamePV;
            filenamePV << scratchdir.c_str() << "/OMEGAPV_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << n_omega_written;
            std::ofstream ofsPV(filenamePV.str(), std::ios::binary);
            OmegaPV.back().WriteBinary(ofsPV);
            std::cout << "Omega written to file " << filenamePV.str() << "  written " << n_omega_written << " so far." << std::endl;
            if (n_omega_written > max_omega_written)
            {
                std::cout << "n_omega_written > max_omega_written.  (" << n_omega_written << " > " << max_omega_written
                        << " ) deleting OMEGA files and calling terminate." << std::endl;
                CleanupScratch();
                std::terminate();
            }
        }
        if (Omega.back().GetModelSpace() != Eta.GetModelSpace() or OmegaPV.back().GetModelSpace() != Etapv.GetModelSpace())
        {
            Omega.back() = Eta;
            OmegaPV.back() = Etapv;
        }
        n_omega_written++;
    }
    else
    {
        Omega.emplace_back(Eta);
        OmegaPV.emplace_back(Etapv);
    }
    Omega.back().Erase();
    OmegaPV.back().Erase();
}


void IMSRGSolverPV::Solve_magnus_euler_PV()
{
  istep = 0;
  generatorPV.Update(FlowingOps[0], FlowingOpsPV[0], Eta, Etapv);
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
  //   WriteFlowStatusPV(flowfile);
  WriteFlowStatusHeaderPV(std::cout);
  WriteFlowStatusPV(std::cout);

  for (istep = 1; s < smax; ++istep)
  {
    double norm_eta = Eta.Norm();
    double norm_etaPV = Etapv.Norm();
    if (sqrt(norm_eta*norm_eta+norm_etaPV*norm_etaPV)< eta_criterion)
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
    double norm_omegaPV = OmegaPV.back().Norm();
    if (sqrt(norm_omega*norm_omega+norm_omegaPV*norm_omegaPV) > omega_norm_max)
    {
      if (hunter_gatherer)
      {
        std::cout<<"Hunter-Gatherer not implemented for PV solver."<<std::endl;
        exit(0);
      }
      else
      {
        NewOmega_PV();
        std::cout<<"new omega"<<std::endl;
      }
      norm_omega = 0;
      norm_omegaPV = 0;
    }
    // ds should never be more than 1, as this is over-rotating
    // Also, since we check if ||Omega|| < omega_norm_max, if we choose ds so that ||Omega|| = omega_norm_max, then we become sensitive
    // to numerical precision details when evaluating the inequality and behavior becomes machine dependent. So we add 1e-5 to the omega_norm_max
    // option to ensure that we're definitely on one side of the inequality.
    if (magnus_adaptive)
        ds = std::min({norm_domega / sqrt(norm_eta * norm_eta + norm_etaPV * norm_etaPV), norm_domega / sqrt(norm_eta * norm_eta + norm_etaPV * norm_etaPV) / (sqrt(norm_omega * norm_omega + norm_omegaPV * norm_omegaPV) + 1.0e-9), (omega_norm_max + 1e-5) / sqrt(norm_eta * norm_eta + norm_etaPV * norm_etaPV), ds_max});
    ds = std::min(ds, smax - s);

    s += ds;
    Eta *= ds; // Here's the Euler step.
    Etapv *= ds;
    // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
    std::tuple<Operator, Operator> Omega_tmp = BCH::BCH_ProductPV(Eta, Etapv, Omega.back(), OmegaPV.back());
    Omega.back() = std::get<0>(Omega_tmp);
    OmegaPV.back() = std::get<1>(Omega_tmp);
    // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
    std::tuple<Operator, Operator> H_tmp; 
    if ((Omega.size() + n_omega_written) < 2)
    {
        H_tmp = BCH::BCH_TransformPV(*H_0, *VPT_0, Omega.back(), OmegaPV.back());
        FlowingOps[0] = std::get<0>(H_tmp);
        FlowingOpsPV[0] = std::get<1>(H_tmp);
    }
    else
    {
        H_tmp = BCH::BCH_TransformPV(H_saved, VPT_saved, Omega.back(), OmegaPV.back());
        FlowingOps[0] = std::get<0>(H_tmp);
        FlowingOpsPV[0] = std::get<1>(H_tmp);
    }

    // if (norm_eta < 1.0 and generator.GetType() == "shell-model-atan")
    // {
    //   generatorPV.SetDenominatorCutoff(1e-6);
    // }

    generatorPV.Update(FlowingOps[0], FlowingOpsPV[0], Eta, Etapv);
    // Etapv.PrintOneBody();
    // Etapv.PrintTwoBody();

    // Write details of the flow
    // WriteFlowStatusPV(flowfile);
    WriteFlowStatusPV(std::cout);
    Elast = FlowingOps[0].ZeroBody;
  }
  //   WriteFlowStatusPV(std::cout);
}


void IMSRGSolverPV::SetGeneratorPV(std::string gen)
{
      generatorPV.SetType(gen);
      if (Omega.back().Norm() > 1e-6 or OmegaPV.back().Norm() >1e-6)
      {
        Eta.Erase();
        Etapv.Erase();
        NewOmega_PV();
      }
      if (magnus_adaptive)
      {
          ds = ds_0;
          in_soft_landing_phase_ = false;
      }
}

void IMSRGSolverPV::WriteFlowStatusPV(std::string fname)
{
      if (fname != "")
      {
            std::ofstream ff(fname, std::ios::app);
            WriteFlowStatusPV(ff);
      }
}
void IMSRGSolverPV::WriteFlowStatusPV(std::ostream &f)
{
      if (f.good())
      {
            int fwidth = 18;
            int fprecision = 8;
            auto &H_s = FlowingOps[0];
            auto &V_s = FlowingOpsPV[0];


            f.setf(std::ios::fixed);
            f << std::setw(5) << std::setw(5) << istep
              << std::setw(10) << std::setprecision(3) << s
              << std::setw(fwidth) << std::setprecision(fprecision) << H_s.ZeroBody
              << std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << H_s.TwoBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << H_s.Norm()
              << std::setw(fwidth) << std::setprecision(fprecision) << V_s.OneBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << V_s.TwoBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << V_s.Norm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Eta.Norm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Eta.OneBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Eta.TwoBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Etapv.Norm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Etapv.OneBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << Etapv.TwoBodyNorm()
              << std::setw(fwidth) << std::setprecision(fprecision) << sqrt(Etapv.Norm()*Etapv.Norm()+Eta.Norm()*Eta.Norm())
              << std::setw(fwidth) << std::setprecision(fprecision) << H_s.GetMP2_Energy()
              << std::setprecision(fprecision)
              << std::setw(7) << std::setprecision(1) << profiler.GetTimes()["real"]
              << std::endl;
      }
}

void IMSRGSolverPV::WriteFlowStatusHeaderPV(std::string fname)
{
      std::ofstream ff;
      if (fname != "")
            ff.open(fname, std::ios::app);
      WriteFlowStatusHeaderPV(ff);
}
void IMSRGSolverPV::WriteFlowStatusHeaderPV(std::ostream &f)
{
      if (f.good())
      {
            for (int x = 0; x < 220; x++)
                  f << "-";
            f << std::endl;
            int fwidth = 18;
            int fprecision = 8;
            f.setf(std::ios::fixed);
            f << std::fixed << std::setw(5) << "i"
              << std::setw(8) << std::setprecision(3) << "s"
              << std::setw(fwidth) << std::setprecision(fprecision) << "E0"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||H1||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||H2||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||H||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Vpt1||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Vpt2||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Vpt||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_1||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_2||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_pv||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_pv_1||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_pv_2||"
              << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_total||"
              << std::setw(15) << std::setprecision(fprecision) << "E(MP2)"
              << std::setw(10) << std::setprecision(fprecision) << "t(s)"
              << std::endl;
            for (int x = 0; x < 220; x++)
                  f << "-";
            f << std::endl;
      }
}

std::tuple<Operator, Operator> IMSRGSolverPV::Transform(Operator &OpIn, Operator &OpInPV)
{
    return Transform_Partial(OpIn, OpInPV,  0);
}

std::tuple<Operator, Operator> IMSRGSolverPV::Transform(Operator &&OpIn, Operator &&OpInPV)
{
    return Transform_Partial(OpIn, OpInPV, 0);
}

std::tuple<Operator, Operator> IMSRGSolverPV::Transform_Partial(Operator &OpIn, Operator &OpInPV, int n)
{
    Operator OpOut = OpIn;
    Operator OpOutPV = OpInPV;
    if (OpOut.GetParticleRank() == 1)
    {
        OpOut.SetParticleRank(2);
    }
       
    if (OpOutPV.GetParticleRank() == 1)
    {
        OpOutPV.SetParticleRank(2);
    }
        
    if (scratchdir != "")
    {
        for (int i = n; i < n_omega_written; i++)
        {
            Operator omega(Eta);
            Operator omegaPV(Etapv);
            // Read Omega
            std::ostringstream filename;
            filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
            std::cout << "Transforming using " << filename.str() << std::endl;
            std::ifstream ifs(filename.str(), std::ios::binary);
            omega.ReadBinary(ifs);
            // Read OmegaPV
            std::ostringstream filenamePV;
            filenamePV << scratchdir.c_str() << "/OMEGAPV_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
            std::cout << "Transforming using " << filenamePV.str() << std::endl;
            std::ifstream ifsPV(filenamePV.str(), std::ios::binary);
            omegaPV.ReadBinary(ifsPV);
            std::tuple<Operator, Operator> Op_tmp = BCH::BCH_TransformPV(OpOut, OpOutPV, omega, omegaPV);
            OpOut = std::get<0>(Op_tmp);
            OpOutPV = std::get<1>(Op_tmp);
            std::cout << "norm of omega = " << omega.Norm() << "norm of omegaPV = "<<omegaPV.Norm()<< std::endl;
            std::cout << " op zero body = " << OpOut.ZeroBody << " opPV zero body = " << OpOutPV.ZeroBody << std::endl;
        }
    }

    for (size_t i = std::max(n - n_omega_written, 0); i < Omega.size(); ++i)
    {
        std::tuple<Operator, Operator> Op_tmp = BCH::BCH_TransformPV(OpOut, OpOutPV, Omega[i], OmegaPV[i]);
        OpOut = std::get<0>(Op_tmp);
        OpOutPV = std::get<1>(Op_tmp);
    }

    return {OpOut, OpOutPV};
}

std::tuple<Operator, Operator> IMSRGSolverPV::Transform_Partial(Operator &&OpIn, Operator &&OpInPV, int n)
{
    
    Operator OpOut = OpIn;
    Operator OpOutPV = OpInPV;
    Operator Optmp;
    Operator OptmpPV;
    if (scratchdir != "")
    {
        for (int i = n; i < n_omega_written; i++)
        {
            Operator omega(Eta);
            Operator omegaPV(Etapv);
            // Read Omega
            std::ostringstream filename;
            filename << scratchdir.c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
            std::ifstream ifs(filename.str(), std::ios::binary);
            omega.ReadBinary(ifs);
            // Read OmegaPV
            std::ostringstream filenamePV;
            filenamePV << scratchdir.c_str() << "/OMEGAPV_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
            std::ifstream ifsPV(filenamePV.str(), std::ios::binary);
            omegaPV.ReadBinary(ifsPV);
            std::tuple<Operator, Operator> Op_tmp = BCH::BCH_TransformPV(OpOut, OpOutPV, omega, omegaPV);
            OpOut = std::get<0>(Op_tmp);
            OpOutPV = std::get<1>(Op_tmp);
            std::cout << "norm of omega = " << omega.Norm() << "norm of omegaPV = " << omegaPV.Norm() << std::endl;
            std::cout << " op zero body = " << OpOut.ZeroBody << " opPV zero body = " << OpOutPV.ZeroBody << std::endl;
        }
    }

    for (size_t i = std::max(n - n_omega_written, 0); i < Omega.size(); ++i)
    {
        std::tuple<Operator, Operator> Op_tmp = BCH::BCH_TransformPV(OpOut, OpOutPV, Omega[i], OmegaPV[i]);
        OpOut = std::get<0>(Op_tmp);
        OpOutPV = std::get<1>(Op_tmp);
    }
    return {OpOut, OpOutPV};
}