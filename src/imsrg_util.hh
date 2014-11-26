#ifndef imsrg_util_hh
#define imsrg_util_hh 1

#include "ModelSpace.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>

#define HBARC 197.3269718 // hc in MeV * fm
#define M_NUCLEON 938.9185 // average nucleon mass in MeV

namespace imsrg_util
{
 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2);
 Operator PSquaredOp(ModelSpace& modelspace);

 double T_1body(ModelSpace& modelspace, int a, int b);
 double Calculate_p1p2(ModelSpace& modelspace, Ket* bra, Ket* ket, int J);
 double HO_density(int n, int l, double hw, double r);
 vector<double> GetOccupations(HartreeFock& hf);
 vector<double> GetOccupations(HartreeFock& hf, IMSRGSolver& imsrgsolver);
 vector<double> GetDensity(vector<double>& occ, vector<double>& R, vector<int>& orbits, ModelSpace& modelspace);
}

#endif
