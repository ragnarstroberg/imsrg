#include "IMSRGSolverPV.hh"
#include "Commutator.hh"
//#include "Operator.hh"
//#include <algorithm>
//#include <cmath>
//#include <iomanip>
//#include <iostream>
//#include <string>
//#include <sstream>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

/*
IMSRGSolverPV::~IMSRGSolverPV()
{
  CleanupScratch();
}

IMSRGSolverPV::IMSRGSolverPV()
    : s(0),ds(0.1),ds_max(0.5),
     norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
     flowfile(""), n_omega_written(0),max_omega_written(500),magnus_adaptive(true),hunter_gatherer(false),perturbative_triples(false),
     pert_triples_this_omega(0),pert_triples_sum(0),ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{}
*/

IMSRGSolverPV::IMSRGSolverPV( Operator &H_in, Operator &VPT_in, Operator &Schiff_in, Operator &Schiffpp_in)
   : IMSRGSolver(H_in), VPT_0(&VPT_in), Schiff_0(&Schiff_in), Schiffpp_0(&Schiffpp_in), FlowingOpsH(1,H_in), FlowingOpsV(1,VPT_in), FlowingOpsSchiff(1,Schiff_in), FlowingOpsSchiffpp(1,Schiffpp_in), Etapv(VPT_in)//modelspace(H_in.GetModelSpace()), H_0(&H_in), VPT_0(&VPT_in), FlowingOpsH(1,H_in),FlowingOpsV(1,VPT_in), Eta(H_in), Etapv(VPT_in), istep(0), s(0),ds(0.1),ds_max(0.5), smax(2.0) , norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
    //flowfile(""), n_omega_written(0),max_omega_written(500),magnus_adaptive(true),hunter_gatherer(false),perturbative_triples(false),
    //pert_triples_this_omega(0),pert_triples_sum(0),ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{
   Eta.Erase();
   Etapv.Erase();
   Eta.SetAntiHermitian();
   Etapv.SetAntiHermitian();
 //  Omega.emplace_back( Eta);
}

//void IMSRGSolverPV::UpdateEta()
//{
//   generatorPV.Update(FlowingOpsH[0],FlowingOpsV[0],Eta,Etapv);
//}

void IMSRGSolverPV::Solve_flow_RK4_PV()
{
   istep = 0;

//   generator.Update(&FlowingOps[0],&Eta);
   generatorPV.Update(FlowingOpsH[0],FlowingOpsV[0],Eta,Etapv);

 //  if (generator.GetType() == "shell-model-atan")
//   {
//     generator.SetDenominatorCutoff(1.0); // do we need this?
//   }

   Elast = H_0->ZeroBody;
   cumulative_error = 0;
    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(std::cout);

//   Operator goosetank_chi( *modelspace, 0,0,0,1);// for use if we do IMSRG2* to mock up the goose tanks
//   Operator goosetank_dchi( *modelspace, 0,0,0,1);// for use if we do IMSRG2* to mock up the goose tanks

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
      double norm_etapv = Etapv.Norm();
      if (norm_eta < eta_criterion and norm_etapv < eta_criterion )
      {
        break;
      }

      ds = std::min(ds_max,smax-s);
      s += ds;

      int nops = FlowingOps.size();
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
      std::vector<Operator> K1S(nops); //Schiff
      std::vector<Operator> K2S(nops);
      std::vector<Operator> K3S(nops);
      std::vector<Operator> K4S(nops);
      std::vector<Operator> KtmpS(nops);
      std::vector<Operator> K1Spp(nops); //Schiffpp
      std::vector<Operator> K2Spp(nops);
      std::vector<Operator> K3Spp(nops);
      std::vector<Operator> K4Spp(nops);
      std::vector<Operator> KtmpSpp(nops);


//      Operator& Hs = FlowingOps[0];   // this is not used explicitly
      for ( int i=0;i<nops; i++ )
      {
//        if (i==0)
//           K1[i] =  Commutator::Commutator( Eta, FlowingOps[i] + goosetank_chi ) ;
//        else
           K1V[i] =  Commutator::Commutator( Eta, FlowingOpsV[i] ) + Commutator::Commutator( Etapv, FlowingOpsH[i] );
           K1H[i] =  Commutator::Commutator( Eta, FlowingOpsH[i] );
           K1Spp[i] = Commutator::Commutator(Eta, FlowingOpsSchiffpp[i]) + Commutator::Commutator(Etapv, FlowingOpsSchiff[i]);
           K1S[i] = Commutator::Commutator(Eta, FlowingOpsSchiff[i]);
           KtmpV[i] = FlowingOpsV[i] + 0.5*ds*K1V[i];
           KtmpH[i] = FlowingOpsH[i] + 0.5*ds*K1H[i];
           KtmpS[i] = FlowingOpsSchiff[i] + 0.5*ds*K1S[i];
           KtmpSpp[i] = FlowingOpsSchiffpp[i] + 0.5*ds*K1Spp[i];
      }
//      Operator K1 = Commutator::Commutator( Eta, Hs );
//      Operator Htmp = Hs + 0.5*ds*K1[0];
//      generator.Update(&Htmp,&Eta);
//      generator.Update(&Ktmp[0],&Eta);
      generatorPV.Update(KtmpH[0],KtmpV[0],Eta,Etapv);
      for (int i=0; i<nops; i++ )
      {
//        K2[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
//        if (i==0)
//            K2[i] = Commutator::Commutator( Eta, Ktmp[i] + goosetank_chi );
//        else
           K2V[i] =  Commutator::Commutator( Eta, KtmpV[i] ) + Commutator::Commutator( Etapv, KtmpH[i] );
           K2H[i] =  Commutator::Commutator( Eta, KtmpH[i] );
           K2Spp[i] = Commutator::Commutator(Eta, FlowingOpsSchiffpp[i]) + Commutator::Commutator(Etapv, FlowingOpsSchiff[i]);
           K2S[i] = Commutator::Commutator(Eta, FlowingOpsSchiff[i]);
           KtmpV[i] = FlowingOpsV[i] + 0.5*ds*K2V[i];
           KtmpH[i] = FlowingOpsH[i] + 0.5*ds*K2H[i];  
           KtmpS[i] = FlowingOpsSchiff[i] + 0.5*ds*K2S[i];
           KtmpSpp[i] = FlowingOpsSchiffpp[i] + 0.5*ds*K2Spp[i]; 
   }

//      Operator K2 = Commutator::Commutator( Eta, Hs+Htmp );
//      Htmp = Hs + 0.5*ds*K2;
//      generator.Update(&Htmp,&Eta);
//      generator.Update(&Ktmp[0],&Eta);
      generatorPV.Update(KtmpH[0],KtmpV[0],Eta,Etapv);
      for (int i=0; i<nops; i++ )
      {
//        K3[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
//        if (i==0)
//           K3[i] = Commutator::Commutator( Eta, Ktmp[i] + goosetank_chi);
//        else
           K3V[i] =  Commutator::Commutator( Eta, KtmpV[i] ) + Commutator::Commutator( Etapv, KtmpH[i] );
           K3H[i] =  Commutator::Commutator( Eta, KtmpH[i] );
           K3Spp[i] = Commutator::Commutator(Eta, FlowingOpsSchiffpp[i]) + Commutator::Commutator(Etapv, FlowingOpsSchiff[i]);
           K3S[i] = Commutator::Commutator(Eta, FlowingOpsSchiff[i]);
           KtmpV[i] = FlowingOpsV[i] + 1.0*ds*K3V[i];
           KtmpH[i] = FlowingOpsH[i] + 1.0*ds*K3H[i];    
           KtmpS[i] = FlowingOpsSchiff[i] + 1.0*ds*K3S[i];
           KtmpSpp[i] = FlowingOpsSchiffpp[i] + 1.0*ds*K3Spp[i];
  }
//      Operator K3 = Commutator::Commutator( Eta, Hs+Htmp );
//      Htmp = Hs + 1.0*ds*K3;
//      generator.Update(&Htmp,&Eta);
//      generator.Update(&Ktmp[0],&Eta);
     generatorPV.Update(KtmpH[0],KtmpV[0],Eta,Etapv);
      for (int i=0; i<nops; i++ )
      {
//        K4[i] = Commutator::Commutator( Eta, FlowingOps[i]+Ktmp[i]);
//        if (i==0)
//           K4[i] = Commutator::Commutator( Eta, Ktmp[i] + goosetank_chi);
//        else
           K4V[i] =  Commutator::Commutator( Eta, KtmpV[i] ) + Commutator::Commutator( Etapv, KtmpH[i] );
           K4H[i] =  Commutator::Commutator( Eta, KtmpH[i] );
           K4Spp[i] = Commutator::Commutator(Eta, FlowingOpsSchiffpp[i]) + Commutator::Commutator(Etapv, FlowingOpsSchiff[i]);
           K4S[i] = Commutator::Commutator(Eta, FlowingOpsSchiff[i]);
     FlowingOpsH[i] += ds/6.0 * ( K1H[i] + 2*K2H[i] + 2*K3H[i] + K4H[i] );
     FlowingOpsV[i] += ds/6.0 * ( K1V[i] + 2*K2V[i] + 2*K3V[i] + K4V[i] );
     FlowingOpsSchiff[i] += ds/6.0 * ( K1S[i] + 2*K2S[i] + 2*K3S[i] + K4S[i] );
     FlowingOpsSchiffpp[i] += ds/6.0 * ( K1Spp[i] + 2*K2Spp[i] + 2*K3Spp[i] + K4Spp[i] );     
 }
//      Operator K4 = Commutator::Commutator( Eta, Hs+Htmp );
//      Hs += ds/6.0 * ( K1 + 2*K2 + 2*K3 + K4);

/*
      if (norm_eta<1.0 and generator.GetType() == "shell-model-atan")
      {
        generator.SetDenominatorCutoff(1e-6);
      }

      if ( Commutator::use_goose_tank_correction )
      {
         goosetank_dchi.EraseOneBody();
         Commutator::comm221ss( Eta, FlowingOps[0], goosetank_dchi );  // update chi.
         for (auto i : modelspace->all_orbits )  // enforce n_in_j + nbar_i nbar_j
         {
           Orbit &oi = modelspace->GetOrbit(i);
           for ( auto j : modelspace->all_orbits )
           {
            Orbit &oj = modelspace->GetOrbit(j);
            goosetank_dchi.OneBody(i,j) *=  oi.occ*oj.occ + (1.0-oi.occ)*(1.0-oj.occ) ;
           }
         }
         goosetank_chi += goosetank_dchi * ds;
//        std::cout << " " << __FILE__ << "  line " << __LINE__ << s << "  " << goosetank_chi.OneBody(1,1) << std::endl;;
      }
*/
//      if ( generator.GetType() == "rspace" ) { generator.SetRegulatorLength(s); };
//      generator.Update(&FlowingOps[0],&Eta);
      generatorPV.Update(FlowingOpsH[0],FlowingOpsV[0],Eta,Etapv);
//      cumulative_error += EstimateStepError();

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(std::cout);
//      profiler.PrintMemory();
      Elast = FlowingOps[0].ZeroBody;

   }
}

void IMSRGSolverPV::SetGeneratorPV(std::string gen)
{     
  generatorPV.SetType(gen);

}



