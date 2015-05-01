
#include "IMSRGSolver.hh"
#include <iomanip>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

IMSRGSolver::~IMSRGSolver()
{
//   cout << "In IMSRGSolver destructor." << endl;
}

IMSRGSolver::IMSRGSolver()
#ifndef NO_ODE
    :ode_monitor(*this)
#endif
{
   method = "BCH";
   generator = "white";
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   flowfile = "";
}

// Constructor
//IMSRGSolver::IMSRGSolver(const Operator &H_in)
IMSRGSolver::IMSRGSolver( Operator &H_in)
//   : H_0(H_in), H_s(H_in), Eta(H_in), Omega(H_in)// ,dOmega(H_in)
   : H_0(&H_in), H_s(H_in), Eta(H_in), Omega(H_in)// ,dOmega(H_in)
#ifndef NO_ODE
    ,ode_monitor(*this)
#endif
{
   method = "BCH";
   generator = "white";
   istep = 0;
   s = 0;
   ds = 0.1;
   ds_max = 0.5;
   smax  = 2.0;
   norm_domega = 0.1;
   flowfile = "";
   modelspace = H_in.GetModelSpace();
   Eta.Erase();
   Eta.SetAntiHermitian();
   Omega.Erase();
   Omega.SetAntiHermitian();
}

void IMSRGSolver::SetHin( Operator & H_in)
{
   modelspace = H_in.GetModelSpace();
   H_0 = &H_in;
   H_s = H_in;
   Eta = H_in;
   Eta.Erase();
   Eta.SetAntiHermitian();
   Omega = Eta;
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
   Omega.Erase();
}

void IMSRGSolver::Solve()
{
   // If we have a flow output file, open it up and write to it here.
   ofstream flowf;
   if (flowfile != "")
   {
      flowf.open(flowfile,ofstream::out);
      flowf.close();
   }
   WriteFlowStatusHeader(cout);

   istep = 0;
   UpdateEta();

    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(cout);

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
      double norm_omega = Omega.Norm();
      // ds should never be more than 1, as this is over-rotating
      //ds = min(norm_domega / norm_eta / (norm_omega+1.0e-9), ds_max); 
      ds = min(min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max); 
      if (ds == ds_max) norm_domega /=2;
      if (s+ds > smax) ds = smax-s;
      s += ds;
//      dOmega = Eta * ds; // Here's the Euler step.
      Eta *= ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
//      Omega = dOmega.BCH_Product( Omega ); 
      Omega = Eta.BCH_Product( Omega ); 

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
//      H_s = H_0.BCH_Transform( Omega );
      H_s = H_0->BCH_Transform( Omega );
        
      UpdateEta();

      // Write details of the flow
//      WriteFlowStatus(flowf);
//      WriteFlowStatus(cout);
      WriteFlowStatus(flowfile);
      WriteFlowStatus(cout);

   }

//   if (flowfile != "")
//      flowf.close();

}


#ifndef NO_ODE
// Implement element-wise division and abs and reduce for Operators.
// This is required for adaptive steppers


/*
Operator operator/ (const Operator& num, const Operator& denom)
{
   Operator quotient = num;
   quotient.ZeroBody /= denom.ZeroBody;
   quotient.OneBody /= denom.OneBody;
   for (int ch=0;ch<quotient.nChannels;++ch)
   {
      for (auto &twobody : quotient.TwoBody.MatEl[ch] )
      {
        int chket = twobody.first;
        twobody.second /= denom.TwoBody.GetMatrix(ch,chket);
      }
   }

}

Operator abs(const Operator& opin)
{
   Operator opout = opin;
   opout.ZeroBody = abs(opout.ZeroBody);
   opout.OneBody = arma::abs(opout.OneBody);
   for (int ch=0;ch<opout.nChannels;++ch)
   {
      for (auto &twobody : opout.TwoBody.MatEl[ch] )
      {
        twobody.second = arma::abs(twobody.second);
      }
   }

}
// Black magic...
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_reduce<Operator&>
{
   template<class Value, class Op>
   Value operator()(const Operator& X, Op op, Value init)
   {
      init = op(init, X.ZeroBody());
      for ( double& v : X.OneBody )
      {
         init = op(init, v) ;
      }
      for (int ch=0;ch<X.nChannels;++ch)
      {
        for (auto &twobody : X.TwoBody[ch] )
        {
          for (double& v : twobody )
          {
             init = op(init, v) ;
          }
        }
      }
      return init;
   }
};
}}}
*/

void IMSRGSolver::Solve_ode()
{

   WriteFlowStatusHeader(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemH, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   size_t steps = integrate_const(stepper, system, H_s, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
/*
   WriteFlowStatusHeader(cout);
     WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
//   auto stepper = make_controlled<runge_kutta_dopri5<Operator, double, Operator, double, vector_space_algebra> >(1e-8,1e-8);
   auto system = std::bind( &IMSRGSolver::ODE_systemH, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   auto stepper = make_controlled(1e-8,1e-8, runge_kutta_dopri5<Operator, double, Operator, double, vector_space_algebra>() );
   size_t steps = integrate_adaptive(stepper, system, H_s, s, smax, ds, monitor);
   monitor.report();
*/
}


void IMSRGSolver::ODE_systemH( Operator& x, Operator& dxdt, const double t)
{
   H_s = x;
   s = t;
   UpdateEta();
   dxdt = Eta.Commutator(x);
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
}


void IMSRGSolver::Solve_ode_magnus()
{
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
   runge_kutta4<Operator, double, Operator, double, vector_space_algebra> stepper;
   auto system = std::bind( &IMSRGSolver::ODE_systemOmega, *this, pl::_1, pl::_2, pl::_3);
   auto monitor = ode_monitor;
   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::ODE_systemOmega( Operator& x, Operator& dxdt, const double t)
{
   s = t;
   Omega = x;
//   H_s = H_0.BCH_Transform(Omega);
   H_s = H_0->BCH_Transform(Omega);
   UpdateEta();
   dxdt = Eta - 0.5*Omega.Commutator(Eta);
   WriteFlowStatus(cout);
   WriteFlowStatus(flowfile);

}

#endif


/// Returns \f$ e^{Omega} \mathcal{O} e^{-Omega} \f$
Operator IMSRGSolver::Transform(Operator& OpIn)
{
   return OpIn.BCH_Transform( Omega );
}

/// Returns \f$ e^{-Omega} \mathcal{O} e^{Omega} \f$
Operator IMSRGSolver::InverseTransform(Operator& OpIn)
{
   Operator negomega = -Omega;
   return OpIn.BCH_Transform( negomega );
}

void IMSRGSolver::UpdateEta()
{
   Eta.EraseOneBody();
   Eta.EraseTwoBody();
   if (generator == "wegner")
   {
     ConstructGenerator_Wegner();
   } 
   else if (generator == "white")
   {
     ConstructGenerator_White();
   } 
   else if (generator == "atan")
   {
     ConstructGenerator_Atan();
   } 
   else if (generator == "shell-model")
   {
     ConstructGenerator_ShellModel();
   }
   else if (generator == "shell-model-atan")
   {
     ConstructGenerator_ShellModel_Atan();
   }
   else if (generator == "hartree-fock")
   {
     ConstructGenerator_HartreeFock();
   }
   else
   {
      cout << "Error. Unkown generator: " << generator << endl;
   }

}



// Epstein-Nesbet energy denominators for White-type generators
// i=particle, j=hole
double IMSRGSolver::Get1bDenominator_ph(int i, int j) 
{
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j) - H_s.TwoBody.GetTBMEmonopole(i,j,i,j);
   return denominator;
}

// i=p j=p'
double IMSRGSolver::Get1bDenominator_pp(int i, int j)
{
   double denominator = H_s.OneBody(i,i) - H_s.OneBody(j,j);
   if (abs(denominator) < 1e-6) return 1e0;
   return denominator;
}



// This could likely be sped up by constructing and storing the monopole matrix
// bra=pp'  ket = hh'
double IMSRGSolver::Get2bDenominator_pphh(int ch, int ibra, int iket) 
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int i = bra.p;
   int j = bra.q;
   int a = ket.p;
   int b = ket.q;
   double denominator = H_s.OneBody(i,i)+ H_s.OneBody(j,j) - H_s.OneBody(a,a) - H_s.OneBody(b,b);

//   denominator       += H_s.GetTBME(ch,ch,i,j,i,j); // pp'pp'
//   denominator       += H_s.GetTBME(ch,ch,a,b,a,b); // hh'hh'
   denominator       += H_s.TwoBody.GetTBMEmonopole(i,j,i,j); // pp'pp'
   denominator       += H_s.TwoBody.GetTBMEmonopole(a,b,a,b); // hh'hh'
   denominator       -= H_s.TwoBody.GetTBMEmonopole(i,a,i,a); // phph
   denominator       -= H_s.TwoBody.GetTBMEmonopole(i,b,i,b); // ph'ph'
   denominator       -= H_s.TwoBody.GetTBMEmonopole(j,a,j,a); // p'hp'h
   denominator       -= H_s.TwoBody.GetTBMEmonopole(j,b,j,b); // p'h'p'h'

   return denominator;
}


// bra = pp'  ket=hp
double IMSRGSolver::Get2bDenominator_pphp(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int p1 = bra.p;
   int p2 = bra.q;
   int h = ket.p;
   int p3 = ket.q;
   double denominator = H_s.OneBody(p1,p1)+ H_s.OneBody(p2,p2) - H_s.OneBody(p3,p3) - H_s.OneBody(h,h);

//   denominator       += H_s.GetTBME(ch,ch,p1,p2,p1,p2); // phph
   denominator       += H_s.TwoBody.GetTBMEmonopole(p1,p2,p1,p2); // phph
   denominator       -= H_s.TwoBody.GetTBMEmonopole(p1,h,p1,h); // pp'pp'
   denominator       -= H_s.TwoBody.GetTBMEmonopole(p2,h,p2,h); // hh'hh'

   return denominator;
}

// bra = pq  ket=vv'
double IMSRGSolver::Get2bDenominator_pppp(int ch, int ibra, int iket)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   Ket & bra = tbc.GetKet(ibra);
   Ket & ket = tbc.GetKet(iket);
   int p1 = bra.p;
   int p2 = bra.q;
   int p3 = ket.p;
   int p4 = ket.q;
   double denominator = H_s.OneBody(p1,p1)+ H_s.OneBody(p2,p2) - H_s.OneBody(p3,p3) - H_s.OneBody(p4,p4);

//   denominator       += H_s.GetTBME(ch,ch,p1,p2,p1,p2); // pp'pp'
//   denominator       -= H_s.GetTBME(ch,ch,p3,p4,p3,p4); // hh'hh'
   denominator       += H_s.TwoBody.GetTBMEmonopole(p1,p2,p1,p2); // pp'pp'
   denominator       -= H_s.TwoBody.GetTBMEmonopole(p3,p4,p3,p4); // hh'hh'

   return denominator;
}



// I haven't used this, so I don't know if it's right.
void IMSRGSolver::ConstructGenerator_Wegner()
{
   Operator H_diag = H_s;
   H_diag.ZeroBody = 0;
   for (auto& a : modelspace->holes)
   {
      for (auto& b : modelspace->valence)
      {
         H_diag.OneBody(a,b) =0;
         H_diag.OneBody(b,a) =0;
      }
   }

   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {  // Note, should also decouple the v and q spaces
      // This is wrong. The projection operator should be different.
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_ph() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_ph() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_pp() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_ph(), tbc.GetKetIndex_hh() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_pp(), tbc.GetKetIndex_hh() ).zeros();
      H_diag.TwoBody.GetMatrix(ch).submat(tbc.GetKetIndex_hh(), tbc.GetKetIndex_pp() ).zeros();
   }
//   H_diag.Symmetrize();

   Eta = H_diag.Commutator(H_s);
}



void IMSRGSolver::ConstructGenerator_White()
{
   // One body piece -- eliminate ph bits
   for ( auto& i : modelspace->particles)
   {
      for ( auto& a : modelspace->holes)
      {
         //double denominator = GetEpsteinNesbet1bDenominator(i,a);
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }

   // Two body piece -- eliminate pp'hh' bits. Note that the hh'hp pieces are accounted
   // for in the normal ordered one-body part.
   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
//      auto& ETA2 = Eta.TwoBody[ch].at(ch);
//      auto& H2 = H_s.TwoBody[ch].at(ch);
      auto& ETA2 = Eta.TwoBody.GetMatrix(ch);
      auto& H2 = H_s.TwoBody.GetMatrix(ch);
//      for ( auto& ibra : tbc.KetIndex_pp)
      for ( auto& ibra : tbc.GetKetIndex_pp() )
      {
//         for ( auto& iket : tbc.KetIndex_hh)
         for ( auto& iket : tbc.GetKetIndex_hh() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
//            double denominator = GetEpsteinNesbet2bDenominator(ch,ibra,iket);
            ETA2(ibra,iket) =  H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
}


void IMSRGSolver::ConstructGenerator_Atan()
{
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->hole_qspace)
   {
      for ( auto& i : modelspace->particles)
      {
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = 0.5*atan(2*H_s.OneBody(i,a)/denominator);
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
      // now decouple valence holes from core holes
      if (modelspace->holes.size() != modelspace->hole_qspace.size() )
      {
        for ( auto& b : modelspace->valence)
        {
           Orbit & ob = modelspace->GetOrbit(b);
           if (ob.ph == 0) continue;
           double denominator = Get1bDenominator_pp(b,a);
           Eta.OneBody(b,a) = 0.5*atan(2*H_s.OneBody(b,a)/denominator);
           Eta.OneBody(a,b) = - Eta.OneBody(b,a);
        }
      }
   }

   // Two body piece -- eliminate pp'hh' bits
   for (int ch=0;ch<Eta.nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta.TwoBody.GetMatrix(ch);
      arma::mat& H2 = H_s.TwoBody.GetMatrix(ch);
      for ( auto& iket : tbc.GetKetIndex_holeq_holeq() )
      {
         for ( auto& ibra : tbc.GetKetIndex_pp() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);

            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // If there are hole orbits in the valence space, decouple those too
         if (modelspace->holes.size() != modelspace->hole_qspace.size() )
         {
           for ( auto& ibra : tbc.GetKetIndex_vv() )
           {
              Ket& bra = modelspace->GetKet(ibra);
              int php = modelspace->GetOrbit(bra.p).ph;
              int phq = modelspace->GetOrbit(bra.q).ph;
              if (php==0 and phq==0) continue;
              double denominator;
              if (php+phq==2) // hhhh
                 denominator = Get2bDenominator_pppp(ch,ibra,iket);
              else
                 denominator = Get2bDenominator_pppp(ch,ibra,iket);
//                 denominator = Get2bDenominator_ppph(ch,ibra,iket);
  
              ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
              ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
           }
         }
      }
    }
}






void IMSRGSolver::ConstructGenerator_ShellModel()
{
   // One body piece -- make sure the valence one-body part is diagonal
   // no excitations out of the core

   Eta.Erase();

   for ( auto& i : modelspace->particles )
   {
      for ( auto& j : modelspace->holes )
      {
         double denominator = Get1bDenominator_ph(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
      for ( auto& j : modelspace->particles )
      {
         if (i==j) continue;
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   }
   for ( auto& i : modelspace->holes )
   {
      for ( auto& j : modelspace->holes )
      {
         if (i==j) continue;
         double denominator = -Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   }


   // Two body piece 
   // we want to drive to zero any terms that take 
   //   |hh> to |pp>
   //   |vh> to |pp>
   //   |vv> to |qq> or |vq>

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);

//      auto& ETA2 = Eta.TwoBody[ch].at(ch);
//      auto& H2 = H_s.TwoBody[ch].at(ch);
      auto& ETA2 = Eta.TwoBody.GetMatrix(ch);
      auto& H2 = H_s.TwoBody.GetMatrix(ch);


      // Decouple vv states from pq states
//      for ( auto& iket : tbc.KetIndex_vv)
      for ( auto& iket : tbc.GetKetIndex_vv() )
      {
         // < qq' | vv' >
//         for ( auto& ibra : tbc.KetIndex_particleq_particleq) 
         for ( auto& ibra : tbc.GetKetIndex_particleq_particleq() ) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vq | vv' > 
//         for ( auto& ibra : tbc.KetIndex_v_particleq) 
         for ( auto& ibra : tbc.GetKetIndex_v_particleq() ) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }


      // Decouple hh states

//      for ( auto& iket : tbc.KetIndex_holeq_holeq)
      for ( auto& iket : tbc.GetKetIndex_holeq_holeq() )
      {

         // < qq' | hh' >
//         for ( auto& ibra : tbc.KetIndex_particleq_particleq)
         for ( auto& ibra : tbc.GetKetIndex_particleq_particleq() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
         // < vq | hh' >
//         for ( auto& ibra : tbc.KetIndex_v_particleq)
         for ( auto& ibra : tbc.GetKetIndex_v_particleq() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vv | hh' >
//         for ( auto& ibra : tbc.KetIndex_vv)
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }


      // Decouple vh states

//      for ( auto& iket : tbc.KetIndex_v_holeq)
      for ( auto& iket : tbc.GetKetIndex_v_holeq() )
      {
         // < qq | vh >
//         for ( auto& ibra : tbc.KetIndex_particleq_particleq)
         for ( auto& ibra : tbc.GetKetIndex_particleq_particleq() )
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vq | vh >
//         for ( auto& ibra : tbc.KetIndex_v_particleq)
         for ( auto& ibra : tbc.GetKetIndex_v_particleq() )
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         // < vv | vh >
//         for ( auto& ibra : tbc.KetIndex_vv)
         for ( auto& ibra : tbc.GetKetIndex_vv() )
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = H2(ibra,iket) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

      }


    }
}



void IMSRGSolver::ConstructGenerator_ShellModel_Atan()
{
   // One body piece -- make sure the valence one-body part is diagonal
   for ( auto& i : modelspace->valence)
   {
      for ( auto& j : modelspace->holes)
      {
         double denominator = Get1bDenominator_ph(i,j);
         Eta.OneBody(i,j) = 0.5*atan(2*H_s.OneBody(i,j)/denominator);
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
      for ( auto& j : modelspace->particles)
      {
         if (i==j) continue;
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = 0.5*atan(2*H_s.OneBody(i,j)/denominator);
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   
   }
   // Two body piece -- eliminate ppvh and pqvv  

   int nchan = modelspace->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat& ETA2 =  Eta.TwoBody.GetMatrix(ch);
      arma::mat& H2 =  H_s.TwoBody.GetMatrix(ch);
      // Decouple vv from qq and qv

      for ( auto& ibra : tbc.GetKetIndex_vv() )
      {
         for ( auto& iket : tbc.GetKetIndex_particleq_particleq() ) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& iket : tbc.GetKetIndex_holeq_holeq() ) 
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& iket : tbc.GetKetIndex_v_particleq() ) 
         {
            double denominator = Get2bDenominator_pppp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& iket : tbc.GetKetIndex_v_holeq() ) 
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }


      // Decouple hh states

      for ( auto& iket : tbc.GetKetIndex_holeq_holeq() )
      {
         for ( auto& ibra : tbc.GetKetIndex_particleq_particleq() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& ibra : tbc.GetKetIndex_v_particleq() )
         {
            double denominator = Get2bDenominator_pphh(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

      // Decouple vh states

      for ( auto& iket : tbc.GetKetIndex_v_holeq() )
      {
         for ( auto& ibra : tbc.GetKetIndex_particleq_particleq() )
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }

         for ( auto& ibra : tbc.GetKetIndex_v_particleq() )
         {
            double denominator = Get2bDenominator_pphp(ch,ibra,iket);
            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket)) / denominator;
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }

    }
}

void IMSRGSolver::ConstructGenerator_HartreeFock()
{
   Eta.SetParticleRank(1);
   // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->holes)
   {
      for ( auto& b : modelspace->holes)
      {
         // Note that for the hole-hole case, the denominator
         // has a minus sign relative to the pp case.
         if (a==b) continue;
         double denominator = -Get1bDenominator_pp(a,b);
         Eta.OneBody(a,b) = H_s.OneBody(a,b)/denominator;
         Eta.OneBody(b,a) = - Eta.OneBody(a,b);
      }
   }
   for ( auto& i : modelspace->particles)
   {
      for ( auto& a : modelspace->holes)
      {
         double denominator = Get1bDenominator_ph(i,a);
         Eta.OneBody(i,a) = H_s.OneBody(i,a)/denominator;
         Eta.OneBody(a,i) = - Eta.OneBody(i,a);
      }
   }  
   for ( auto& i : modelspace->particles)
   {
      for ( auto& j : modelspace->particles)
      {
         if (i==j) continue;
         double denominator = Get1bDenominator_pp(i,j);
         Eta.OneBody(i,j) = H_s.OneBody(i,j)/denominator;
         Eta.OneBody(j,i) = - Eta.OneBody(i,j);
      }
   }  
}



// count number of equations to be solved
int IMSRGSolver::GetSystemDimension()
{
   int dim = 1; // zero-body part

   // one-body part
   int N = H_0->OneBody.n_cols;
   dim += N*(N+1)/2;
   dim += H_0->TwoBody.Dimension();
/*
   for( auto channels : H_0.TwoBody.MatEl )
   {
     for( auto mtx : channels )
     {
        N = mtx.second.n_cols;
        dim += N*(N+1)/2;
     }
   }
*/
   return dim;
}



void IMSRGSolver::WriteFlowStatus(string fname)
{
   ofstream ff;
   if (fname !="") ff.open(fname,ios::app);
   WriteFlowStatus(ff);
}
void IMSRGSolver::WriteFlowStatus(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(ios::fixed);
      f << setw(5) << istep
        << setw(fwidth) << setprecision(3) << s
        << setw(fwidth) << setprecision(fprecision) << H_s.ZeroBody 
        << setw(fwidth) << setprecision(fprecision) << H_s.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << H_s.TwoBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Omega.Norm()
        << setw(fwidth) << setprecision(fprecision) << Eta.OneBodyNorm()
        << setw(fwidth) << setprecision(fprecision) << Eta.TwoBodyNorm()
        << endl;
   }

}

void IMSRGSolver::WriteFlowStatusHeader(string fname)
{
   ofstream ff;
   if (fname !="") ff.open(fname,ios::app);
   WriteFlowStatusHeader(ff);
}
void IMSRGSolver::WriteFlowStatusHeader(ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(ios::fixed);
      f << setw(5) << "i"
        << setw(fwidth) << setprecision(3) << "s"
        << setw(fwidth) << setprecision(fprecision) << "E0"
        << setw(fwidth) << setprecision(fprecision) << "||H_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||H_2||" 
        << setw(fwidth) << setprecision(fprecision) << "||Omega||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_1||" 
        << setw(fwidth) << setprecision(fprecision) << "||Eta_2||" 
        << endl;
      f << "-----------------------------------------------------------------------------------------------------------------------" << endl;
   }

}
