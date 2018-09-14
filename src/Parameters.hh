///////////////////////////////////////////////////////////////////////////////////
//    Parameters.hh, part of  imsrg++
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

#include <map>
#include <string>
#include <sstream>

double r2p = 0.770; // effective proton intrinsic charge radius squared
double r2n = -0.1149; // effective neutron intrisic charge radius squared
double DF = 0.033; // Darwin-Foldy correction to the charge radius


//////////////////
class Parameters
{
 public:
  static map<string,string> string_par;
  static map<string,double> double_par;
  static map<string,int> int_par;
  static map<string,vector<string>> vec_par;

  Parameters(){};
  Parameters(int, char**);
  void ParseCommandLineArgs(int, char**);
  void PrintOptions();
  string s(string);
  double d(string);
  int i(string);
  vector<string> v(string);
  string DefaultFlowFile();
  string DefaultIntFile();
  bool help_mode;
};

map<string,string> Parameters::string_par = {
  {"2bme",			"none"},
  {"3bme",			"none"},
  {"core_generator",		"atan"},	// generator used for core part of 2-step decoupling
  {"valence_generator",		"shell-model-atan"},	// generator used for valence decoupling and 1-step (also single-ref)
  {"flowfile",			"default"},	// name of output flow file
  {"intfile",			"default"},	// name of output interaction fille
  {"fmt2",			"me2j"},	// can also be navratil or Navratil to read Petr's TBME format
  {"fmt3",			"me3j"},	// can also be navratil or Navratil to read Petr's TBME format
  {"reference",			"default"},	// nucleus used for HF and normal ordering.
  {"valence_space",		""},		// either valence space or nucleus for single reference
  {"custom_valence_space",      ""},		// if the provided valence spaces just aren't good enough for you
  {"basis",			"HF"},		// use HF basis or oscillator basis. HF is better.
  {"method",			"magnus"},	// can be magnus or flow or a few other things
  {"denominator_delta_orbit",	"none"},	// pick specific orbit to apply the delta
  {"LECs",			"EM2.0_2.0"},	// low energy constants for the interaction, only used with Johannes' hdf5 file format
  {"scratch",			""},		// scratch directory for writing operators in binary format
  {"use_brueckner_bch",          "false"},	// switch to Brueckner version of BCH
  {"valence_file_format",       "nushellx"},	// file format for valence space interaction
  {"occ_file",			"none"},	// name of file containing orbit occupations
  {"goose_tank",		"false"},	// do goose_tank correction to commutators
  {"write_omega",		"false"},	// write omega to disk
  {"nucleon_mass_correction",	"false"},	// include effect of proton-neutron mass splitting
};


map<string,double> Parameters::double_par = {
  {"hw",		20.0},
  {"smax",		200.0},	// maximum s. If we reach this,	terminate even if we're not converged.
  {"dsmax",		0.5},	// maximum step size
  {"ds_0",		0.5},	// initial step size
  {"domega",		0.2},	// max for norm of eta * ds
  {"omega_norm_max",	0.25},  // norm of omega before we do the splitting
  {"ode_tolerance",	1e-6},	// error tolerance for the ode solver
  {"denominator_delta",	   0},	// offset added to the denominator in the generator
  {"BetaCM",               0},  // Prefactor for Lawson-Glockner term
  {"hwBetaCM",            -1},  // Oscillator frequency used in the Lawson-Glockner term. Negative value means use the frequency of the basis
  {"eta_criterion",     1e-6},  // Threshold on ||eta|| for convergence in the flow

};

map<string,int> Parameters::int_par = {
  {"A",	-1},	// Aeff for kinetic energy. -1 means take A of reference
  {"e3max",		12},	
  {"emax",		6},
  {"lmax3",		-1}, // lmax for the 3body interaction
  {"nsteps",		-1},	// do the decoupling in 1 step or core-then-valence. -1 means default
  {"file2e1max",	12},
  {"file2e2max",	24},
  {"file2lmax",		10},
  {"file3e1max",	12},
  {"file3e2max",	24},
  {"file3e3max",	12},
};

map<string,vector<string>> Parameters::vec_par = {
 {"Operators", {} },
 {"OperatorsFromFile", {} },  // These will mostly be MECs for operators
 {"SPWF",{} }, // single-particle wave functions in HF basis
};


Parameters::Parameters(int argc, char** argv)
{
  help_mode = false;
  ParseCommandLineArgs(argc, argv);
} 

void Parameters::ParseCommandLineArgs(int argc, char** argv)
{
  cout << "====================  Parameters (defaults set in Parameters.hh) ===================" << endl;
  for (int iarg=1; iarg<argc; ++iarg)
  {
    string arg = argv[iarg];
    if (arg=="help" or arg=="-help" or arg=="--help")
    {
       cout << "\nUsage:\n\timsrg++ option1=variable1 option2=variable2...\n" << endl;
       cout << "At minimum, the 2bme file is required.\n" << endl;
       PrintOptions();
       help_mode = true;
       return;
    }
    size_t pos = arg.find("=");
    string var = arg.substr(0,pos);
    string val = arg.substr(pos+1);
    if (string_par.find(var) != string_par.end() )
    {
      if (val.size() > 0)
        istringstream(val) >> string_par[var];
      cout << var << " => " << string_par[var] << endl;
    }
    else if (double_par.find(var) != double_par.end() )
    {
      if (val.size() > 0)
        istringstream(val) >> double_par[var];
      cout << var << " => " << double_par[var] << endl;
    }
    else if (int_par.find(var) != int_par.end() )
    {
      if (val.size() > 0)
        istringstream(val) >> int_par[var];
      cout << var << " => " << int_par[var] << endl;
    }
    else if (vec_par.find(var) != vec_par.end() )
    {
      if (val.size() > 0)
      {
        istringstream ss(val);
        string tmp;
        while( getline(ss,tmp,',')) vec_par[var].push_back(tmp);
      }
      cout << var << " => ";
      for (auto x : vec_par[var]) cout << x << ",";
      cout << endl;
    }
    else
    {
      cout << "Unkown parameter: " << var << " => " << val << endl;
    }
    
  }
  if (string_par["flowfile"]=="default") string_par["flowfile"] = DefaultFlowFile();
  if (string_par["intfile"]=="default") string_par["intfile"] = DefaultIntFile();
  cout << "====================================================================================" << endl;
}

string Parameters::s(string key)
{
  return string_par[key];
}
double Parameters::d(string key)
{
  return double_par[key];
}
int Parameters::i(string key)
{
  return int_par[key];
}
vector<string> Parameters::v(string key)
{
  return vec_par[key];
}

string Parameters::DefaultFlowFile()
{
  char strbuf[200];
  sprintf(strbuf, "output/BCH_%s_%s_%s_hw%.0f_e%d_A%d.dat",string_par["method"].c_str(),string_par["reference"].c_str(),string_par["valence_space"].c_str(),double_par["hw"],int_par["emax"],int_par["A"]);
  return string(strbuf);
}

string Parameters::DefaultIntFile()
{
  char strbuf[200];
  sprintf(strbuf, "output/%s_%s_%s_hw%.0f_e%d_A%d",string_par["method"].c_str(),string_par["reference"].c_str(),string_par["valence_space"].c_str(),double_par["hw"],int_par["emax"],int_par["A"]);
  return string(strbuf);
}

void Parameters::PrintOptions()
{
  cout << "Input parameters and default values: " << endl;
  for (auto& strpar : string_par)
  {
    cout << "\t" << left << setw(30) << strpar.first << ":  " << strpar.second << endl;
  }
  for (auto& doublepar : double_par)
  {
    cout <<  "\t" << left << setw(30) <<doublepar.first << ":  " << doublepar.second << endl;
  }
  for (auto& intpar : int_par)
  {
    cout <<  "\t" << left << setw(30) <<intpar.first << ":  " << intpar.second << endl;
  }
  for (auto& vecpar : vec_par)
  {
    cout <<  "\t" << left << setw(30) << vecpar.first << ":  ";
    for (auto& op : vecpar.second) cout << op << ",";
    cout << endl;
  }

}


