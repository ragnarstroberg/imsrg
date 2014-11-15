#include "imsrg_util.hh"


namespace imsrg_util
{

 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(&modelspace);
   int indx = modelspace.Index1(n,l,j2,tz2);
   NumOp.ZeroBody = 0;
   NumOp.EraseOneBody();
   NumOp.EraseTwoBody();
   NumOp.OneBody(indx,indx) = 1;
   return NumOp;
 }


 double HO_density(int n, int l, double hw, double r)
 {
    double v = M_NUCLEON * hw / (HBARC*HBARC);
    double Norm = pow(v/2.,1.5+l) * M_SQRT2/M_SQRTPI * pow(2,n+2*l+3) * gsl_sf_fact(n) / gsl_sf_doublefact(2*n + 2*l + 1);
    double L = gsl_sf_laguerre_n(n, l+0.5, v*r*r);
    double rho = Norm * pow(r,2*l) * exp(-v * r*r) * L * L;
    return rho;

 }

 // Just do the HF transformation
 vector<double> GetOccupations(HartreeFock& hf)
 {
    ModelSpace* modelspace = hf.Hbare.modelspace;
    int norb = modelspace->GetNumberOrbits();
    vector<double> occupation(norb);

    for (int i=0; i<norb; ++i)
    {
      Orbit * oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi->n,oi->l,oi->j2,oi->tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      occupation[i] = N_NO.ZeroBody;
      cout << oi->n << " " << oi->l << " " << oi->j2 << "/2 " << occupation[i] << endl;
    }
    return occupation;
 }

 // Do the full IMSRG transformation
 vector<double> GetOccupations(HartreeFock& hf, IMSRGSolver& imsrgsolver)
 {
    ModelSpace* modelspace = imsrgsolver.modelspace;
    int norb = modelspace->GetNumberOrbits();
    vector<double> occupation(norb,0);

    for (int i=0; i<norb; ++i)
    {
      Orbit * oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi->n,oi->l,oi->j2,oi->tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      // Transform to the imsrg evolved basis
      Operator N_final = imsrgsolver.Transform(N_NO);
      occupation[i] = N_final.ZeroBody;
    }
    return occupation;
 }

 vector<double> GetDensity( vector<double>& occupation, vector<double>& R, vector<int>& orbits, ModelSpace& modelspace )
 {
     int nr_steps = R.size();
     double hw = modelspace.GetHbarOmega();
     vector<double> dens(nr_steps,0);
     for (int& i : orbits)
     {
       Orbit * oi = modelspace.GetOrbit(i);
       for (int ir=0;ir<nr_steps;++ir)
       {
          dens[ir] += HO_density(oi->n,oi->l,hw,R[ir]) * occupation[i];
       }
     }
     return dens;
 }


 Operator PSquaredOp(ModelSpace& modelspace)
 {
   Operator PsqOp = Operator(&modelspace);
   int nchan = modelspace.GetNumberTwoBodyChannels();
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket * bra = tbc.GetKet(ibra);
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            double mat_el = PSquaredMatEl(modelspace,bra,ket,tbc.J);
            PsqOp.TwoBody[ch](ibra,iket) = mat_el;
            PsqOp.TwoBody[ch](iket,ibra) = mat_el;
         }
      }
   }

 }

 double PSquaredMatEl(ModelSpace& modelspace, Ket* bra, Ket* ket, int J)
 {
   Orbit * oa = modelspace.GetOrbit(bra->p);
   Orbit * ob = modelspace.GetOrbit(bra->q);
   Orbit * oc = modelspace.GetOrbit(ket->p);
   Orbit * od = modelspace.GetOrbit(ket->q);
   

   return 0;
 }


}


