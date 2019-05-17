/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///        
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///       
///                                               \___\/                                          ///        
///                                                                                               ///        
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//    imsrg++.cc, part of  imsrg++
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


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
#ifdef BUILDVERSION
  cout << "######  imsrg++ build version: " << BUILDVERSION << endl;
#endif
  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  string inputtbme = parameters.s("2bme");
  string input3bme = parameters.s("3bme");
  string reference = parameters.s("reference");
  string valence_space = parameters.s("valence_space");
  string custom_valence_space = parameters.s("custom_valence_space");
  string basis = parameters.s("basis");
  string method = parameters.s("method");
  string flowfile = parameters.s("flowfile");
  string intfile = parameters.s("intfile");
  string core_generator = parameters.s("core_generator");
  string valence_generator = parameters.s("valence_generator");
  string fmt2 = parameters.s("fmt2");
  string fmt3 = parameters.s("fmt3");
  string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  string LECs = parameters.s("LECs");
  string scratch = parameters.s("scratch");
  string use_brueckner_bch = parameters.s("use_brueckner_bch");
  string valence_file_format = parameters.s("valence_file_format");
  string occ_file = parameters.s("occ_file");
  string goose_tank = parameters.s("goose_tank");
  string write_omega = parameters.s("write_omega");
  string nucleon_mass_correction = parameters.s("nucleon_mass_correction");

  int eMax = parameters.i("emax");
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");

  double hw = parameters.d("hw");
  double smax = parameters.d("smax");
  double ode_tolerance = parameters.d("ode_tolerance");
  double dsmax = parameters.d("dsmax");
  double ds_0 = parameters.d("ds_0");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max"); 
  double denominator_delta = parameters.d("denominator_delta");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double eta_criterion = parameters.d("eta_criterion");

  vector<string> opnames = parameters.v("Operators");
  vector<string> opsfromfile = parameters.v("OperatorsFromFile");

  vector<Operator> ops;
  vector<string> spwf = parameters.v("SPWF");


  ifstream test;
  // test 2bme file
  if (inputtbme != "none")
  {
    test.open(inputtbme);
//    if( not test.good() and fmt2!="oakridge_binary")
    if( not test.good() and  fmt2.find("oakridge")== string::npos)
    {
      cout << "trouble reading " << inputtbme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }
  // test 3bme file
  if (input3bme != "none")
  {
    test.open(input3bme);
    if( not test.good() )
    {
      cout << "trouble reading " << input3bme << " exiting. " << endl;
      return 1;
    }
    test.close();
  }



  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );

//  ModelSpace modelspace;

  if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
  {
    if (valence_space=="") // if no name is given, then just name it "custom"
    {
      parameters.string_par["valence_space"] = "custom";
      flowfile = parameters.DefaultFlowFile();
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  modelspace.SetE3max(E3max);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);
  
  cout << "Making the operator..." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


  if ( goose_tank == "true" or goose_tank == "True")
  {
    Hbare.SetUseGooseTank(true);
  }

  cout << "Reading interactions..." << endl;


  if (inputtbme != "none")
  {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2.find("oakridge") != string::npos )
    { // input format should be: singleparticle.dat,vnn.dat
      size_t comma_pos = inputtbme.find_first_of(",");
      if ( fmt2.find("bin") != string::npos )
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "binary");
      else
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "ascii");
    }
    else if (fmt2 == "takayuki" )
      rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
    else if (fmt2 == "nushellx" )
      rw.ReadNuShellX_int( Hbare, inputtbme );
  
    cout << "done reading 2N" << endl;
  }

  if (fmt2 != "nushellx")  // Don't need to add kinetic energy if we read a shell model interaction
  {
    Hbare += Trel_Op(modelspace);
  }

  if ( nucleon_mass_correction == "true" or nucleon_mass_correction == "True" )  
  {  // correction to kinetic energy because M_proton != M_neutron
    Hbare += Trel_Masscorrection_Op(modelspace);
  }
  
  if (Hbare.particle_rank >=3)
  {
    rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
    cout << "done reading 3N" << endl;
  }  


  // Add a Lawson term. If hwBetaCM is specified, use that frequency
  if (std::abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }

  cout << "Creating HF" << endl;
  HartreeFock hf(Hbare);
  cout << "Solving" << endl;
  hf.Solve();
//  cout << "EHF = " << hf.EHF << endl;
  
//  Operator HNO;
  Operator& HNO = Hbare;
  if (basis == "HF" and method !="HF")
    HNO = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    HNO = Hbare.DoNormalOrdering();


  int n_radial_points = 40;
  double Rmax = 10.0;
  vector<index_t> spwf_indices = modelspace.String2Index(spwf);
  vector<double> R(n_radial_points);
  vector<double> PSI(n_radial_points);
  for ( index_t i=0; i< spwf.size(); ++i)
  {
    for (int rstep=0;rstep<n_radial_points;++rstep) R[rstep] = Rmax/n_radial_points * rstep;
    hf.GetRadialWF(spwf_indices[i], R, PSI);
    ofstream wf_file (intfile + "_spwf_" + spwf[i] + ".dat");
    for ( index_t rstep=0; rstep<R.size(); ++rstep)  wf_file << fixed << setw(10) << setprecision(7) << R[rstep] << "   " << setw(10) << setprecision(7) << PSI[rstep] << endl;
    cout << "About to close wf file" << endl;
//    wf_file.close();
  }
  if (spwf.size() > 0)   cout << "Done with SPWF" << endl;

  HNO -= BetaCM * 1.5*hwBetaCM;
  cout << "Hbare 0b = " << HNO.ZeroBody << endl;

  if (method != "HF")
  {
    cout << "Perturbative estimates of gs energy:" << endl;
    double EMP2 = HNO.GetMP2_Energy();
    cout << "EMP2 = " << EMP2 << endl; 
    double EMP3 = HNO.GetMP3_Energy();
    cout << "EMP3 = " << EMP3 << endl; 
    cout << "To 3rd order, E = " << HNO.ZeroBody+EMP2+EMP3 << endl;
  }



  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
    ops.emplace_back( imsrg_util::OperatorFromString(modelspace,opname) );
  }


  // the format should look like OpName^j_t_p_r^/path/to/file
  for (auto& tag : opsfromfile)
  {
    istringstream ss(tag);
    string opname,qnumbers,fname;
    vector<int> qn(4);
    
    getline(ss,opname,'^');
    getline(ss,qnumbers,'^');
    getline(ss,fname,'^');
    ss.str(qnumbers);
    ss.clear();
//    cout << " ss.str = " << ss.str() << endl;
    for (int i=0;i<4;i++)
    {
      string tmp;
      getline(ss,tmp,'_');
      istringstream(tmp) >> qn[i];
//      cout << i << " [" << tmp << "] " << qn[i] << endl;
    }
//    ss >> j; ss.ignore();
//    ss >> t; ss.ignore();
//    ss >> p; ss.ignore();
//    ss >> r;
    int j,t,p,r;
    j = qn[0];
    t = qn[1];
    p = qn[2];
    r = qn[3];
//    cout << "Parsed tag. opname = " << opname << "  qnumbers = " << qnumbers << "  " << j << " " << t << " " << p << " " << r << "   file = " << fname << endl;
    Operator op(modelspace,j,t,p,r);
    rw.Read2bCurrent_Navratil( fname, op );
    ops.push_back( op );
    opnames.push_back( opname );
  }


//  // This is only for testing and should be deleted
//  if (opsfromfile.size()>1)
//  {
//    Operator ogt = imsrg_util::OperatorFromString( modelspace, "GamowTeller");
//    imsrg_util::Embed1BodyIn2Body( ogt, modelspace.GetTargetMass());
//    ogt.EraseOneBody();
//    ops.push_back(ogt);
//    opnames.push_back("GTembed");
//  }


  for (auto& op : ops)
  {
     if (basis == "HF") op = hf.TransformToHFBasis(op);
     op = op.DoNormalOrdering();
     if (method == "MP3")
     {
       double dop = op.MP1_Eval( HNO );
       cout << "Operator 1st order correction  " << dop << "  ->  " << op.ZeroBody + dop << endl;
     }
//     cout << endl << op.OneBody << endl;
  }
  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
    cout << " HF charge radius = " << ( std::abs(Rp2.ZeroBody)<1e-6 ? 0.0 : sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) ) << endl; 
  }
  for (index_t i=0;i<ops.size();++i)
  {
    cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
  }


  cout << "HF Single particle energies:" << endl;
//  hf.PrintSPE();
  hf.PrintSPEandWF();
  cout << endl;
  
  if ( method == "HF" or method == "MP3")
  {
    HNO.PrintTimes();
    return 0;
  }


  if (method == "FCI")
  {
    HNO = HNO.UndoNormalOrdering();
    rw.WriteNuShellX_int(HNO,intfile+".int");
    rw.WriteNuShellX_sps(HNO,intfile+".sp");

    for (index_t i=0;i<ops.size();++i)
    {
      ops[i] = ops[i].UndoNormalOrdering();
      if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
      {
        rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
        rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
      }
    }
    HNO.PrintTimes();
    return 0;
  }

//  Operator HlowT = HNO;
//  double Temp = hw;
//  double Efermi = 0;
//  Operator Eye = HNO;
//  Eye.Eye();
//  HlowT.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  Eye.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  cout << "Initial low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) <<"  with normalization  " << Eye.Trace( modelspace.GetAref(),modelspace.GetZref() ) << endl; 

  IMSRGSolver imsrgsolver(HNO);
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  bool brueckner_restart = false;
  
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=500;
    method = "magnus";
  }
  if (method.find("brueckner") != string::npos)
  {
    if (method=="brueckner2") brueckner_restart=true;
    if (method=="brueckner1step")
    { 
       nsteps = 1;
       core_generator = valence_generator;
    }
    use_brueckner_bch = "true";
    omega_norm_max=500;
    method = "magnus";
  }

  if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
  {
//    Hbare.SetUseBruecknerBCH(true);
    HNO.SetUseBruecknerBCH(true);
    cout << "Using Brueckner flavor of BCH" << endl;
  }

  imsrgsolver.SetMethod(method);
//  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetHin(HNO);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=string::npos)
  {
   if (ds_0>1e-2)
   {
     ds_0 = 1e-4;
     dsmax = 1e-2;
     imsrgsolver.SetDs(ds_0);
     imsrgsolver.SetDsmax(dsmax);
   }
  }
  imsrgsolver.Solve();

//  HlowT = imsrgsolver.Transform(HlowT);
//  cout << "After Solve, low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) << endl; 

  if (method == "magnus")
  {
//    for (size_t i=0;i<ops.size();++i)
//    {
//      Operator tmp = imsrgsolver.Transform(ops[i]);
////      rw.WriteOperatorHuman(tmp,intfile+opnames[i]+"_step1.op");
//    }
//    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  if (brueckner_restart)
  {
     arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = hf.GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;

    imsrgsolver.SetGenerator(valence_generator);
    modelspace.ResetFirstPass();
    if (valence_generator.find("imaginary")!=string::npos)
    {
     if (ds_0>1e-2)
     {
       ds_0 = 1e-4;
       dsmax = 1e-2;
       imsrgsolver.SetDs(ds_0);
       imsrgsolver.SetDsmax(dsmax);
     }
    }
    imsrgsolver.SetSmax(smax);
    imsrgsolver.Solve();
  }



  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) cout << "transforming operators" << endl;
    for (size_t i=0;i<ops.size();++i)
    {
      cout << opnames[i] << " " << endl;
      ops[i] = imsrgsolver.Transform(ops[i]);
      cout << " (" << ops[i].ZeroBody << " ) " << endl; 
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
    }
    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted/ensemble normal ordering 
  // we now re-normal order wrt to the core
  // and do any remaining flow.
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
  {
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (not renormal_order)
    {
      for (auto c : modelspace.core)
      {
         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {

    HNO = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    HNO = HNO.UndoNormalOrdering();

    ms2.SetReference(ms2.core); // change the reference
    HNO.SetModelSpace(ms2);

    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << endl;
    HNO = HNO.DoNormalOrdering();

    imsrgsolver.SetHin(HNO);
    imsrgsolver.SetEtaCriterion(1e-4);
    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    cout << "Final transformation on the operators..." << endl;
    int iop = 0;
    for (auto& op : ops)
    {
      cout << opnames[iop++] << endl;
      op = op.UndoNormalOrdering();
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
    }
    cout << "Writing files: " << intfile << endl;
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");

    if (method == "magnus")
    {
       for (index_t i=0;i<ops.size();++i)
       {
          if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
          {
            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
          }
          else
          {
            rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
            rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    cout << "Core Energy = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;
    for (index_t i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << endl; 
         cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
      }
      if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
      {
        cout << "Writing operator to " << intfile+opnames[i]+".op" << endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }
  }


//  cout << "Made it here and write_omega is " << write_omega << endl;
  if (write_omega == "true" or write_omega == "True")
  {
    cout << "writing Omega to " << intfile << "_omega.op" << endl;
    rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");
  }


  Hbare.PrintTimes();
 
  return 0;
}

