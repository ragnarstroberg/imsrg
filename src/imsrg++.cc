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
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"
#include "PhysicalConstants.hh"


int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
#ifdef BUILDVERSION
  std::cout << "######  imsrg++ build version: " << BUILDVERSION << std::endl;
#endif

  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string input3bme_type = parameters.s("3bme_type");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string basis = parameters.s("basis");
  std::string method = parameters.s("method");
  std::string flowfile = parameters.s("flowfile");
  std::string intfile = parameters.s("intfile");
  std::string core_generator = parameters.s("core_generator");
  std::string valence_generator = parameters.s("valence_generator");
  std::string fmt2 = parameters.s("fmt2");
  std::string fmt3 = parameters.s("fmt3");
  std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  std::string LECs = parameters.s("LECs");
  std::string scratch = parameters.s("scratch");
  std::string valence_file_format = parameters.s("valence_file_format");
  std::string occ_file = parameters.s("occ_file");
  std::string physical_system = parameters.s("physical_system");
  bool use_brueckner_bch = parameters.s("use_brueckner_bch") == "true";
  bool nucleon_mass_correction = parameters.s("nucleon_mass_correction") == "true";
  bool relativistic_correction = parameters.s("relativistic_correction") == "true";
  bool IMSRG3 = parameters.s("IMSRG3") == "true";
  bool write_omega = parameters.s("write_omega") == "true";
  bool freeze_occupations = parameters.s("freeze_occupations")=="true";
  bool hunter_gatherer = parameters.s("hunter_gatherer") == "true";
  bool goose_tank = parameters.s("goose_tank") == "true";
  bool discard_residual_input3N = parameters.s("discard_residual_input3N")=="true";
  bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
  bool store_3bme_pn = (parameters.s("store_3bme_pn")=="true");

  int eMax = parameters.i("emax");
  int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
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
  int atomicZ = parameters.i("atomicZ");
  int emax_unocc = parameters.i("emax_unocc");
  int dE3max = parameters.i("dE3max");

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
  double hw_trap = parameters.d("hw_trap");

  std::vector<std::string> opnames = parameters.v("Operators");
  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
  std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
  std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
  std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

  std::vector<Operator> ops;
  std::vector<std::string> spwf = parameters.v("SPWF");


  std::ifstream test;
  // test 2bme file
  if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
  {
    test.open(inputtbme);
////    if( not test.good() and fmt2!="oakridge_binary")
//    if( not test.good() and  fmt2.find("oakridge")== std::string::npos)
    if( not test.good() )
    {
      std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
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
      std::cout << "trouble reading " << input3bme << " exiting. " << std::endl;
      return 1;
    }
    test.close();
  }



  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );

  // Test whether the scratch directory exists and we can write to it.
  // This is necessary because otherwise you get garbage for transformed operators and it's
  // not obvious what went wrong.
  if ( method=="magnus" and  scratch != "" and scratch!= "/dev/null" and scratch != "/dev/null/")
  {
    std::string testfilename = scratch + "/_this_is_a_test_delete_me";
    std::ofstream testout(testfilename);
    testout << "PASSED" << std::endl;
    testout.close();
    std::remove( testfilename.c_str() );
    if ( not testout.good() )
    {
      std::cout << "WARNING in " << __FILE__ <<  " failed test write to scratch directory " << scratch;
      if (opnames.size()>0 )
      {
      std::cout << "   dying now. " << std::endl;
      exit(EXIT_FAILURE);
      }
      else std::cout << std::endl;
    }
  }
  if ( (method=="magnus") and (scratch=="/dev/null" or scratch=="/dev/null/") )
  {
    if ( opnames.size() > 0 )
    {
      std::cout << "WARNING!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators: ";
      for (auto opn : opnames ) std::cout << opn << " ";
      std::cout << "   dying now." << std::endl;
      exit(EXIT_FAILURE);
    }
  }


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

  modelspace.SetE3max(E3max);
  modelspace.SetLmax(lmax);
//  if (lmax!= 99999)
//  {
//    modelspace.ClearVectors();
//    modelspace.Init(eMax, reference,valence_space);  
//  }

  if (emax_unocc>0)
  {
    modelspace.SetEmaxUnocc(emax_unocc);
  }

  if (physical_system == "atomic")
  {
    modelspace.InitSingleSpecies(eMax, reference, valence_space);
  }

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }


  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);



// For both dagger operators and single particle wave functions, it's convenient to
// just get every orbit in the valence space. So if SPWF="valence" ,  we append all valence orbits
  if ( std::find( spwf.begin(), spwf.end(), "valence" ) != spwf.end() )
  {
    // this erase/remove idiom is needed because remove just shuffles things around rather than actually removing it.
    spwf.erase( std::remove( spwf.begin(), spwf.end(), "valence" ), std::end(spwf) );
    for ( auto v : modelspace.valence )
    {
      spwf.push_back( modelspace.Index2String(v) );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerHF_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerHF_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerHF_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerHF_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerAlln_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerAlln_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerAlln_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerAlln_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }


//  std::cout << "Making the Hamiltonian..." << std::endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


  if ( goose_tank )
  {
//    Hbare.SetUseGooseTank(true);
    Commutator::SetUseGooseTank(true);
  }

  std::cout << "Reading interactions..." << std::endl;


  if (inputtbme != "none")
  {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2.find("oakridge") != std::string::npos )
    { // input format should be: singleparticle.dat,vnn.dat
      size_t comma_pos = inputtbme.find_first_of(",");
      if ( fmt2.find("bin") != std::string::npos )
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "binary");
      else
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "ascii");
    }
    else if (fmt2 == "takayuki" )
      rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
    else if (fmt2 == "nushellx" )
      rw.ReadNuShellX_int( Hbare, inputtbme );
    else if (fmt2 == "schematic" )
    {
      std::cout << "using schematic potential " << inputtbme << std::endl;
      if ( inputtbme == "Minnesota") Hbare += imsrg_util::MinnesotaPotential( modelspace );
    }

    std::cout << "done reading 2N" << std::endl;
  }

  if (inputtbme == "none" and physical_system == "atomic")
  {
//    std::cout << "||||| Calling InitSingleSpecies( " <<eMax << ", " << reference << ", " << valence_space << "  ||||||| " << std::endl;
//    std::cout << "*****************************************************************************************************" << std::endl;
//    modelspace.SetLmax(lmax);
//    modelspace.InitSingleSpecies(eMax, reference, valence_space);
//    Hbare = Operator(modelspace,0,0,0,particle_rank);
//    Hbare.SetHermitian();
    using PhysConst::M_ELECTRON;
    using PhysConst::M_NUCLEON;
    int Z = (atomicZ>=0) ?  atomicZ : modelspace.GetTargetZ() ;
//    const double HARTREE = 27.21138602; // 1 Hartree in eV
    Hbare -= Z*imsrg_util::VCentralCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;
//    std::cout << "After conversion, central coulomb 1-body looks like " << std::endl << Hbare.OneBody << std::endl << std::endl;
    Hbare += imsrg_util::VCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;  // convert oscillator length from fm with nucleon mass to nm with electon mass (in eV).
//    std::cout << "done with VCoulomb_Op" << std::endl;
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // Don't need to rescale this, because it's related to the oscillator frequency, which we input.
//    std::cout << "done with KineticEnergy" << std::endl;
    Hbare /= PhysConst::HARTREE; // Convert to Hartree
  }

  if (fmt2 != "nushellx" and physical_system != "atomic" and hw_trap < 0)  // Don't need to add kinetic energy if we read a shell model interaction
  {
    Hbare += imsrg_util::Trel_Op(modelspace);
    if (Hbare.OneBody.has_nan())
    {
       std::cout << "  Looks like the Trel op is hosed from the get go." << std::endl;
    }
  }

  if ( hw_trap > 0 )
  {
    Hbare += 0.5 * (PhysConst::M_NUCLEON * hw_trap * hw_trap)/(PhysConst::HBARC*PhysConst::HBARC) * imsrg_util::RSquaredOp(modelspace); // add lab-frame harmonic trap
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // use lab-frame kinetic energy
  }


  if ( nucleon_mass_correction)
  {  // correction to kinetic energy because M_proton != M_neutron
    Hbare += imsrg_util::Trel_Masscorrection_Op(modelspace);
  }
  if ( relativistic_correction)
  {
    Hbare += imsrg_util::KineticEnergy_RelativisticCorr(modelspace);
  }

  if (Hbare.particle_rank >=3)
  {
    if(input3bme_type == "full"){
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
      std::cout << "done reading 3N" << std::endl;
    }
    if(input3bme_type == "no2b"){
      double t_start = omp_get_wtime();
      Hbare.ThreeBodyNO2B.Allocate(modelspace, file3e1max, file3e2max, file3e3max, file3e1max, input3bme);
      Hbare.profiler.timer["ThreeBodyNO2B::Allocate"] += omp_get_wtime() - t_start;
      t_start = omp_get_wtime();
      Hbare.ThreeBodyNO2B.ReadFile();
      Hbare.profiler.timer["ThreeBodyNO2B::ReadFile"] += omp_get_wtime() - t_start;
      std::cout << "done reading 3N" << std::endl;
    }
  }

  if (store_3bme_pn)
  {
    Hbare.ThreeBody.TransformToPN();
  }


  // Add a Lawson term. If hwBetaCM is specified, use that frequency
  if (std::abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    std::ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }


  // If we're doing FCI, we want things normal ordered wrt the vacuum
  // and we want things diagonal when normal ordered wrt the vacuum.
//  if ( method == "FCI" )
//  {
//    modelspace.SetReference("vacuum");
//  }


  std::cout << "Creating HF" << std::endl;
//  HartreeFock hf(Hbare);
  HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so no harm done here.

  if (not freeze_occupations )  hf.UnFreezeOccupations();
  std::cout << "Solving" << std::endl;
//  if (basis=="HF")
  hf.Solve();

  int hno_particle_rank = 2;
  if ((IMSRG3) and (Hbare.ThreeBodyNorm() > 1e-5))  hno_particle_rank = 3;
  if (discard_residual_input3N) hno_particle_rank = 2;
//  Operator HNO;
  Operator& HNO = Hbare;
  if (basis == "HF" and method !="HF")
  {
//    ThreeBodyME hf3b;
    HNO = hf.GetNormalOrderedH( hno_particle_rank );
  }
  else if (basis == "NAT") // we want to use the natural orbital basis
  {
    hf.UseNATOccupations( use_NAT_occupations );

// This calls GetDensityMatrix(), which computes the 1b density matrix up to MBPT2
// using the NO2B Hamiltonian in the HF basis, obtained with GetNormalOrderedH().
// Then it calls DiagonalizeRho() which diagonalizes the density matrix, yielding the natural orbital basis.
    hf.GetNaturalOrbitals();
    HNO = hf.GetNormalOrderedHNAT( hno_particle_rank );
//    HNO = hf.GetNormalOrderedHNAT();

    // For now, even if we use the NAT occupations, we switch back to naive occupations after the normal ordering
    // This should be investigated in more detail.
    if (use_NAT_occupations)
    {
      hf.FillLowestOrbits();
      std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
      HNO = HNO.UndoNormalOrdering();
      hf.UpdateReference();
      modelspace.SetReference(modelspace.core); // change the reference
      std::cout << "Doing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
      HNO = HNO.DoNormalOrdering();
    }

  }
  else if (basis == "oscillator")
  {
    HNO = Hbare.DoNormalOrdering();
  }


  if (IMSRG3)
  {
    modelspace.SetdE3max(dE3max);
    std::cout << "You have chosen IMSRG3. good luck..." << std::endl;

    if (hno_particle_rank<3 )
    {
      Operator H3(modelspace,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      if ( Hbare.ThreeBodyNorm() <1e-6 )
      {
        HNO.ThreeBody.TransformToPN();
      }
    }
  }



  if ( spwf.size() > 0 )
  {
    imsrg_util::WriteSPWaveFunctions( spwf, hf, intfile);
  }



  HNO -= BetaCM * 1.5*hwBetaCM;
  std::cout << "Hbare 0b = " << HNO.ZeroBody << std::endl;

  if (method != "HF")
  {
    std::cout << "Perturbative estimates of gs energy:" << std::endl;
    double EMP2 = HNO.GetMP2_Energy();
    double EMP2_3B = HNO.GetMP2_3BEnergy();
    std::cout << "EMP2 = " << EMP2 << std::endl;
    std::cout << "EMP2_3B = " << EMP2_3B << std::endl;
//    double EMP3 = HNO.GetMP3_Energy();
    std::array<double,3> Emp_3 = HNO.GetMP3_Energy();
    double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
    std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
//    cout << "EMP3 = " << EMP3 << endl;
    std::cout << "To 3rd order, E = " << HNO.ZeroBody + EMP2 + EMP3 + EMP2_3B << std::endl;
  }



  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
    ops.emplace_back( imsrg_util::OperatorFromString(modelspace,opname) );
  }
  // Calculate first order perturbative correction to some operators, if that's what we asked for.
  // Strictly speaking, it doesn't make much sense to do this and then proceed with the IMSRG calculation,
  // but I'm not here to tell people what to do...
  for (auto& opnamept1 : opnamesPT1 )
  {
    ops.emplace_back( imsrg_util::FirstOrderCorr_1b( imsrg_util::OperatorFromString(modelspace,opnamept1)   , HNO ) );
    opnames.push_back( opnamept1+"PT1" );
  }
  for (auto& opnametda : opnamesTDA )
  {  // passing the argument "TDA" just sets the phhp and hpph blocks to zero in the RPA calculation
    ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnametda)   , HNO, "TDA" ) );
    opnames.push_back( opnametda+"TDA" );
  }
  for (auto& opnamerpa : opnamesRPA )
  {
    ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnamerpa)   , HNO, "RPA" ) );
    opnames.push_back( opnamerpa+"RPA" );
  }


  // the format should look like OpName^j_t_p_r^/path/to/file
  for (auto& tag : opsfromfile)
  {
    std::istringstream ss(tag);
    std::string opname,qnumbers,fname;
    std::vector<int> qn(4);

    getline(ss,opname,'^');
    getline(ss,qnumbers,'^');
    getline(ss,fname,'^');
    ss.str(qnumbers);
    ss.clear();
    for (int i=0;i<4;i++)
    {
      std::string tmp;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> qn[i];
    }

    int j,t,p,r;
    j = qn[0];
    t = qn[1];
    p = qn[2];
    r = qn[3];
//    std::cout << "Parsed tag. opname = " << opname << "  qnumbers = " << qnumbers << "  " << j << " " << t << " " << p << " " << r << "   file = " << fname << std::endl;
    Operator op(modelspace,j,t,p,r);
    rw.Read2bCurrent_Navratil( fname, op );
    ops.push_back( op );
    opnames.push_back( opname );
  }



//  for (auto& op : ops)
  for (size_t i=0;i<ops.size();++i)
  {
     // We don't transform a DaggerHF, because we want the a^dagger to already refer to the HF basis.
    if ((basis == "HF") and (opnames[i].find("DaggerHF") == std::string::npos)  )
    {
      ops[i] = hf.TransformToHFBasis(ops[i]);
    }
    else if ((basis == "NAT") and (opnames[i].find("DaggerHF") == std::string::npos)  )
    {
      ops[i] = hf.TransformHOToNATBasis(ops[i]);
    }
    ops[i] = ops[i].DoNormalOrdering();
    if (method == "MP3")
    {
      double dop = ops[i].MP1_Eval( HNO );
      std::cout << "Operator 1st order correction  " << dop << "  ->  " << ops[i].ZeroBody + dop << std::endl;
    }
  }


  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    std::cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << std::endl;
    std::cout << " HF charge radius = " << ( abs(Rp2.ZeroBody)<1e-6 ? 0.0 : sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) ) << std::endl;
  }
  for (index_t i=0;i<ops.size();++i)
  {
    std::cout << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
  }


  if (basis=="HF" or basis=="NAT")
  {
    std::cout << basis << " Single particle energies and wave functions:" << std::endl;
    hf.PrintSPEandWF();
    std::cout << std::endl;
  }

  if ( method == "HF" or method == "MP3")
  {
    HNO.PrintTimes();
    return 0;
  }


  if (method == "FCI")
  {
   // we want the 1b piece to be diagonal in the vacuum NO representation
    HNO = HNO.UndoNormalOrdering();
    double previous_zero_body = HNO.ZeroBody;
    modelspace.SetReference("vacuum");
    HartreeFock hfvac(HNO);
    hfvac.Solve();

//    Operator Hvac = hfvac.GetNormalOrderedH();
    HNO = hfvac.GetNormalOrderedH();
    std::cout << "HNO had zero body = " << HNO.ZeroBody << "  and I add " << previous_zero_body << " to it. " << std::endl;
    HNO.ZeroBody += previous_zero_body;

    rw.WriteNuShellX_int(HNO,intfile+".int");
    rw.WriteNuShellX_sps(HNO,intfile+".sp");

    std::cout << "NO wrt vacuum. One Body term is hopfully still diagonal?" << std::endl << HNO.OneBody << std::endl;

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


  IMSRGSolver imsrgsolver(HNO);
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  imsrgsolver.max_omega_written = 500;
  bool brueckner_restart = false;
  if (hunter_gatherer) imsrgsolver.SetHunterGatherer( true);

  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=50000;
    method = "magnus";
  }
  if (method.find("brueckner") != std::string::npos)
  {
    if (method=="brueckner2") brueckner_restart=true;
    if (method=="brueckner1step")
    {
       nsteps = 1;
       core_generator = valence_generator;
    }
    use_brueckner_bch = true;
    omega_norm_max=500;
    method = "magnus";
  }

  if (use_brueckner_bch)
  {
//    Hbare.SetUseBruecknerBCH(true);
//    HNO.SetUseBruecknerBCH(true);
    Commutator::SetUseBruecknerBCH(true);
    std::cout << "Using Brueckner flavor of BCH" << std::endl;
  }
  if (IMSRG3)
  {
    Commutator::SetUseIMSRG3(true);
    std::cout << "Using IMSRG(3) commutators. This will probably be slow..." << std::endl;
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

  if (method == "flow" or method == "flow_RK4" )
  {
    for (auto& op : ops )  imsrgsolver.AddOperator( op );
  }

  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=std::string::npos or core_generator.find("wegner")!=std::string::npos )
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

  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }

//  HlowT = imsrgsolver.Transform(HlowT);
//  std::cout << "After Solve, low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) << std::endl;

//  if (method == "magnus")
//  {
////    for (size_t i=0;i<ops.size();++i)
////    {
////      Operator tmp = imsrgsolver.Transform(ops[i]);
//////      rw.WriteOperatorHuman(tmp,intfile+opnames[i]+"_step1.op");
////    }
////    std::cout << std::endl;
//    // increase smax in case we need to do additional steps
//    smax *= 1.5;
//    imsrgsolver.SetSmax(smax);
//  }


  if (brueckner_restart)
  {
     arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = hf.GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1 and valence_space != reference) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;

    imsrgsolver.SetGenerator(valence_generator);
    std::cout << "Setting generator to " << valence_generator << std::endl;
    modelspace.ResetFirstPass();
    if (valence_generator.find("imaginary")!=std::string::npos or valence_generator.find("wegner")!=std::string::npos)
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
    if (ops.size()>0) std::cout << "transforming operators" << std::endl;
    for (size_t i=0;i<ops.size();++i)
    {
      std::cout << opnames[i] << " " << std::endl;
      ops[i] = imsrgsolver.Transform(ops[i]);
      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
    }
    std::cout << std::endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }
  if (method == "flow" or method == "flow_RK4" )
  {
    for (size_t i=0;i<ops.size();++i)
    {
      ops[i] = imsrgsolver.GetOperator(i+1);  // the zero-th operator is the Hamiltonian
    }
  }


  // If we're doing targeted/ensemble normal ordering
  // we now re-normal order wrt to the core
  // and do any remaining flow.
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
//  if (modelspace.valence.size() > 0 or basis=="NAT")
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
    std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
    std::cout << "Before doing so, the spes are " << std::endl;
    for ( auto i : modelspace.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    if (IMSRG3)
    {
      std::cout << "Re-normal-ordering wrt the core. For now, we just throw away the 3N at this step." << std::endl;
      HNO.SetNumberLegs(4);
      HNO.SetParticleRank(2);
    }
    HNO = HNO.UndoNormalOrdering();

    ms2.SetReference(ms2.core); // change the reference
    HNO.SetModelSpace(ms2);


    std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
    HNO = HNO.DoNormalOrdering();

// More flowing is unnecessary, since things should stay decoupled.
    imsrgsolver.SetHin(HNO);
//    imsrgsolver.SetEtaCriterion(1e-4);
//    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    std::cout << "Final transformation on the operators..." << std::endl;
    int iop = 0;
    for (auto& op : ops)
    {
      std::cout << opnames[iop++] << std::endl;
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
    std::cout << "Writing files: " << intfile << std::endl;
    if (valence_file_format == "tokyo")
    {
     rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt", "");
    }
    else
    {
      rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
    }

    if (method == "magnus" or method=="flow_RK4")
    {
       for (index_t i=0;i<ops.size();++i)
       {
          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
          {
            if (valence_file_format == "tokyo")
            {
              rw.WriteTokyo(ops[i],intfile+opnames[i]+".snt", "op");
            }
            else
            {
              rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
            }
          }
          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
          {
//            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
          }
          else
          {
            if (valence_file_format == "tokyo")
            {
              rw.WriteTensorTokyo(intfile+opnames[i]+"_2b.snt",ops[i]);
            }
            else
            {
              rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
              rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
            }
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody << std::endl;
    for (index_t i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      std::cout << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
         std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << std::endl;
      }
      if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
      {
        std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }
  }


//  std::cout << "Made it here and write_omega is " << write_omega << std::endl;
  if (write_omega)
  {
    std::cout << "writing Omega to " << intfile << "_omega.op" << std::endl;
    rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");
  }


  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }
  Hbare.PrintTimes();

  return 0;
}

