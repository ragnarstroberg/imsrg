///////////////////////////////////////////////////////////////////////////////////
//    IMSRGSolver.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef IMSRGSolver_h
#define IMSRGSolver_h 1

#include <fstream>
#include <string>
#include <deque>
#include "Operator.hh"
#include "Generator.hh"
#include "IMSRGProfiler.hh"
#include "ReadWrite.hh"

using namespace std;


class IMSRGSolver
{

  public:

//  private:
  ModelSpace* modelspace;
  ReadWrite* rw;
  Operator* H_0; 
  deque<Operator> FlowingOps;
  Operator H_saved;
  Operator Eta;
  deque<Operator> Omega;
  Generator generator;
  int istep;
  double s;
  double ds;
  double ds_max;
  double smax;
  double norm_domega;
  double omega_norm_max;
  double eta_criterion;
  string method;
  string flowfile;
  IMSRGProfiler profiler;
  int n_omega_written;
  int max_omega_written;
  bool magnus_adaptive;



  ~IMSRGSolver();
  IMSRGSolver();
  IMSRGSolver( Operator& H_in);
  void NewOmega();
  void SetHin( Operator& H_in);
  void SetReadWrite( ReadWrite& r){rw = &r;};
  void Reset();
  void AddOperator(Operator& Op){FlowingOps.push_back(Op);};
  void UpdateEta(); // Force eta to be calculated. For debugging.

  void SetMethod(string m){method=m;};
  void Solve();
  void Solve_magnus_euler();
  void Solve_magnus_modified_euler();

  Operator Transform(Operator& OpIn);
  Operator Transform(Operator&& OpIn);
  Operator InverseTransform(Operator& OpIn);
  Operator GetOmega(int i){return Omega[i];};
  void SetOmega(size_t i, Operator& om);
  size_t GetOmegaSize(){return Omega.size();};
  int GetNOmegaWritten(){return n_omega_written;};
  Operator Transform_Partial(Operator& OpIn, int n);
  Operator Transform_Partial(Operator&& OpIn, int n);

  void SetFlowFile(string s);
  void SetDs(double d){ds = d;};
  void SetDsmax(double d){ds_max = d;};
  void SetdOmega(double d){norm_domega = d;};
  void SetSmax(double d){smax = d;};
  void SetGenerator(string g);
  void SetOmegaNormMax(double x){omega_norm_max = x;};
  void SetODETolerance(float x){ode_e_abs=x;ode_e_rel=x;};
  void SetEtaCriterion(float x){eta_criterion = x;};
  void SetMagnusAdaptive(bool b){magnus_adaptive = b;};

  int GetSystemDimension();
  Operator& GetH_s(){return FlowingOps[0];};
  Operator& GetEta(){return Eta;};
  Generator& GetGenerator(){return generator;};

  void UpdateOmega();
  void UpdateH();

  void WriteFlowStatus(ostream&);
  void WriteFlowStatusHeader(ostream&);
  void WriteFlowStatus(string);
  void WriteFlowStatusHeader(string);

  void SetDenominatorCutoff(double c){generator.SetDenominatorCutoff(c);};
  void SetDenominatorDelta(double d){generator.SetDenominatorDelta(d);};
  void SetDenominatorDeltaIndex(int i){generator.SetDenominatorDeltaIndex(i);};
  void SetDenominatorDeltaOrbit(string o){generator.SetDenominatorDeltaOrbit(o);};

  void CleanupScratch();


  // This is used to get flow info from odeint
  class ODE_Monitor
  {
    public:
     ODE_Monitor(IMSRGSolver& solver)
           : imsrgsolver(solver), times(solver.times), E0(solver.E0),
              eta1(solver.eta1),eta2(solver.eta2) {};
     IMSRGSolver& imsrgsolver;
     vector<double>& times;
     vector<double>& E0;
     vector<double>& eta1;
     vector<double>& eta2;
//     void operator() (const vector<Operator>& x, double t)
     void operator() (const deque<Operator>& x, double t)
     {
        times.push_back(t);
        E0.push_back(x.front().ZeroBody);
        eta1.push_back(imsrgsolver.Eta.OneBodyNorm());
        eta2.push_back(imsrgsolver.Eta.TwoBodyNorm());
     }
     void report()
     {
        for (size_t i=0; i<times.size(); ++i)
        {
           cout << times[i] << "  " << E0[i] << "  "
                << eta1[i]  << "  " << eta2[i] << endl;
        }
     }
  };


  ODE_Monitor ode_monitor;

  vector<double> times;
  vector<double> E0;
  vector<double> eta1;
  vector<double> eta2;
  string ode_mode;
  float ode_e_abs;
  float ode_e_rel;


//  void operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t);
  void operator()( const deque<Operator>& x, deque<Operator>& dxdt, const double t);
  void Solve_ode();
  void Solve_ode_adaptive();
  void Solve_ode_magnus();


};





#endif

