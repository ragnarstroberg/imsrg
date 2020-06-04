#include "M0nu.hh"
#include "AngMom.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>

//makusing namespace AngMom;
using namespace PhysConst;


namespace M0nu
{


  // Converts (a,b,c,d) in base (maxa+1,maxb+1,maxc+1,maxd+1) to an ordered decimal integer
  // eg: for (maxa,maxb,maxc,maxd) = (1,1,1,1) decimalgen would convert (0,1,1,0) to 6, ie: a binary to decimal converter
  // NOTE: to make comparisons between such decimals, keep (maxa,maxb,maxc,maxd) consistent
  // PS: I tots made this thing up, did some testing and it seemed to work...
  // ... hopefully it's good, but could be a source of weird bugs if anything comes up with missing integrals wrt dq, see GetIntegral(...)
  int decimalgen(int a, int b, int c, int d, int maxb, int maxc, int maxd)
  {
    int coeff1 = maxd + 1;
    int coeff2 = (maxc + 1)*coeff1;
    int coeff3 = (maxb + 1)*coeff2;
    return coeff3*a + coeff2*b + coeff1*c + d; // eg: (0,1,1,0) -> (2*2*2)*0 + (2*2)*1 + (2)*1 + 0 = 6
  }

  // My personal GSL integration error handling
  // NOTE: must put "#pragma omp critical" calls around this function when calling it in the code
  // NOTE: in order to use this, one must set "gsl_set_error_handler_off();" beforehand
  // NOTE: I'm setting the integral to 0 when it doesn't converge because for the all the case I've encountered so far, the resul was smaller than 1e-5
  // this might be the cause for weird bug later on so make sure to check the log file afterwards
  void GSLerror(std::string errstr, int status, size_t limit, double epsabs, double epsrel, double abserr, double result, int n, int l, int np, int lp, double q)
  {
    if (status)
    {
      std::cout<<"WARNING "<<errstr;
      if (status == GSL_EMAXITER)
      {
        std::cout<<" the maximum number of subdivisions was exceeded."<<std::endl;
      }
      else if (status == GSL_EROUND)
      {
        std::cout<<" cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table."<<std::endl;
      }
      else if (status == GSL_ESING)
      {
        std::cout<<" a non-integrable singularity or other bad integrand behavior was found in the integration interval."<<std::endl;
      }
      else if (status == GSL_EDIVERGE)
      { std::cout<<"result = "<<std::setprecision(12)<<result<<std::endl;
        std::cout<<" the integral is divergent, or too slowly convergent to be integrated numerically."<<std::endl;
        std::cout<<"Assuming value is 0"<<std::endl;
        std::cout<<"q ="<<q<<std::endl;
        result = 0;
      }
      else
      {
        std::cout<<"apparently I'm missing a GSL error status handle..."<<std::endl;
      }
      std::cout<<"limit  = "<<limit<<std::endl;
      std::cout<<"epsabs = "<<epsabs<<std::endl;
      std::cout<<"epsrel = "<<epsrel<<std::endl;
      std::cout<<"abserr = "<<abserr<<std::endl;
      std::cout<<"result = "<<std::setprecision(12)<<result<<std::endl;
      std::cout<<"n = "<<n<<", l = "<<l<<", np = "<<np<<", lp = "<<lp<<std::endl;
      double errat = abs(abserr/result);
      if (errat <= epsrel) // this if-else is a last resort, but only seems to rarely trigger for Tensor with high emax for a few n,l,np,lp
      {
        std::cout<<"moving forward since abs(abserr/result) = "<<errat<<" <= epsrel..."<<std::endl;
      }
      else if (abserr <= epsabs)
      {
        std::cout<<"moving forward since abserr <= epsabs..."<<std::endl;
      }
      else
      {
        std::cout<<"ERROR 37728: my GSL integration failed! :("<<std::endl;
        std::cout<<"Assuming value is 0"<<std::endl;
        result = 0;
      }
    }
  }

  //Intgrand function to compute the RBMEs
  double frRBME(double r, void *params)
  {
    struct frparams * p = (struct frparams *)params;
    int n       = (p->n);
    int l       = (p->l);
    int np      = (p->np);
    int lp      = (p->lp);
    double hw   = (p->hw);
    int rho     = (p->rho);
    double q    = (p->q);
    q /= HBARC; // put it into [fm^-1] since r is in [fm]
    q *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    return r*r*imsrg_util::HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,q*r)*imsrg_util::HO_Radial_psi(np,lp,hw,r);
  }

  //Integrand to compute the RBMEs with AV18 SRC
  double frRBME_AV18(double r, void *params)
  {
    struct frparams * p = (struct frparams *)params;
    int n       = (p->n);
    int l       = (p->l);
    int np      = (p->np);
    int lp      = (p->lp);
    double hw   = (p->hw);
    int rho     = (p->rho);
    double q    = (p->q);
    double a = 1.49;
    double b = 1.45; 
    double c = 0.92;
    q /= HBARC; // put it into [fm^-1] since r is in [fm]
    q *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    a *= 2.0; // " " " " " " " " " " " "
    b *= 2.0; // " " " " " " " " " " " "
    double rsq = r*r;
    double JTF = -c*exp(-a*rsq)*(1 - b*rsq); // the Jastrow-type correlation function
    return rsq*imsrg_util::HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,q*r)*imsrg_util::HO_Radial_psi(np,lp,hw,r)*(1 + JTF)*(1 + JTF); // see Equation (6) of PRC.86.067304(2012), for example
  }

  //Integrand to compute the RBMEs with CDB SRC
  double frRBME_CDB(double r, void *params)
  {
    struct frparams * p = (struct frparams *)params;
    int n       = (p->n);
    int l       = (p->l);
    int np      = (p->np);
    int lp      = (p->lp);
    double hw   = (p->hw);
    int rho     = (p->rho);
    double q    = (p->q);
    double a = 1.52;
    double b = 1.88; 
    double c = 0.46;
    q /= HBARC; // put it into [fm^-1] since r is in [fm]
    q *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    a *= 2.0; // " " " " " " " " " " " "
    b *= 2.0; // " " " " " " " " " " " "
    double rsq = r*r;
    double JTF = -c*exp(-a*rsq)*(1 - b*rsq); // the Jastrow-type correlation function
    return rsq*imsrg_util::HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,q*r)*imsrg_util::HO_Radial_psi(np,lp,hw,r)*(1 + JTF)*(1 + JTF); // see Equation (6) of PRC.86.067304(2012), for example
  }

    //Integrand to compute the RBMEs with MS SRC
  double frRBME_MS(double r, void *params)
  {
    struct frparams * p = (struct frparams *)params;
    int n       = (p->n);
    int l       = (p->l);
    int np      = (p->np);
    int lp      = (p->lp);
    double hw   = (p->hw);
    int rho     = (p->rho);
    double q    = (p->q);
    double a = 1.10;
    double b = 0.68; 
    double c = 1;
    q /= HBARC; // put it into [fm^-1] since r is in [fm]
    q *= SQRT2; // this is to stay consistent with r_{rel} = (r_1 - r_2)/SQRT2 in Moshinksy brackets
    a *= 2.0; // " " " " " " " " " " " "
    b *= 2.0; // " " " " " " " " " " " "
    double rsq = r*r;
    double JTF = -c*exp(-a*rsq)*(1 - b*rsq); // the Jastrow-type correlation function
    return rsq*imsrg_util::HO_Radial_psi(n,l,hw,r)*gsl_sf_bessel_jl(rho,q*r)*imsrg_util::HO_Radial_psi(np,lp,hw,r)*(1 + JTF)*(1 + JTF); // see Equation (6) of PRC.86.067304(2012), for example
  }


  // Calculates the RBMEs by doing the integral
  //with respect to dr using gsl_integration_qagiu
  double rbme(int n, int l, int np, int lp, double hw, int rho, double q, std::string src)
  {
    gsl_set_error_handler_off(); // set standard GSL error handling off, so I can do my own
    int statusr; // this will hold the GSL error status
    std::string errstrr="7828877-1: qagiu dr:"; // a unique GSL error string
    size_t limitr = 100000; // the number of double precision intervals that can be held by the workspace below
    gsl_integration_workspace * workspacer = gsl_integration_workspace_alloc(limitr); // GSL: "One workspace may be used multiple times...", "...workspaces should be allocated on a per-thread basis."
    struct frparams frparams = {n,l,np,lp,hw,rho,q};
    std::map<std::string, double(*) (double,void *)> frRBMElist = {{"none",frRBME},{"AV18",frRBME_AV18},{"CD-Bonn",frRBME_CDB},{"Miller-Spencer",frRBME_MS}};
    gsl_function Fr;
    Fr.function = frRBMElist[src]; // picks up the function to integrate for the src choice
    Fr.params = &frparams;
    double r1 = 0; // integrate over the range [r1,Inf)
    double epsabsr = 1e-7; // absolute error tolerance for the integration
    double epsrelr = 1e-5; // relative error tolerance for the adaptive integration
    double rbme = 0; // qagiu integration result
    double abserrr; // qagiu integration error estimate
    statusr = gsl_integration_qagiu(&Fr,r1,epsabsr,epsrelr,limitr,workspacer,&rbme,&abserrr); // perform the qagiu integration (over r)
    #pragma omp critical
    {
      GSLerror(errstrr,statusr,limitr,epsabsr,epsrelr,abserrr,rbme,n,l,np,lp,q); // check if the GSL integration worked
    }
    gsl_integration_workspace_free(workspacer); // free the allocated memory for the integration workspace, since GSL is written in C
    return rbme;
  }


  //Form factors of the neutrino potential
  double GTFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
    double Mprosq  = M_PROTON*M_PROTON; // the proton mass squared [MeV^2]
    double Mpionsq = M_PION_CHARGED*M_PION_CHARGED; // the pion mass squared [MeV^2]
    double gV      = NUCLEON_VECTOR_G/pow((1.0 + (qsq/(CUTV*CUTV))),2); // from Equation (4.14) of my thesis
    double gA      = NUCLEON_AXIAL_G/pow((1.0 + (qsq/(CUTA*CUTA))),2); // " " " " " "
    double gP      = (2*M_PROTON*gA)/(qsq + Mpionsq); // " " " " " "
    double gM      = MAGMOM*gV; // " " " " " "
    double ff      = ((gA*gA) - ((gA*gP*qsq)/(3*M_PROTON)) + (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(6*Mprosq)))/(NUCLEON_AXIAL_G*NUCLEON_AXIAL_G); //Equation (4.12) of Charlie's thesis
    return ff;
  }

  double FermiFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
//    double Mprosq  = M_PROTON*M_PROTON; // the proton mass squared [MeV^2]      SRS commented out because not used.
//    double Mpionsq = M_PION_CHARGED*M_PION_CHARGED; // the pion mass squared [MeV^2]
    double gV      = NUCLEON_VECTOR_G/pow((1.0 + (qsq/(CUTV*CUTV))),2); // from Equation (4.14) of my thesis
//    double gA      = NUCLEON_AXIAL_G/pow((1.0 + (qsq/(CUTA*CUTA))),2); // " " " " " "
//    double gP      = (2*M_PROTON*gA)/(qsq + Mpionsq); // " " " " " "
//    double gM      = MAGMOM*gV; // " " " " " "
    double ff      = (gV*gV)/(NUCLEON_VECTOR_G*NUCLEON_VECTOR_G); //Equation (4.11) of Charlie's thesis
    return ff;
  }

  double TensorFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
    double Mprosq  = M_PROTON*M_PROTON; // the proton mass squared [MeV^2]
    double Mpionsq = M_PION_CHARGED*M_PION_CHARGED; // the pion mass squared [MeV^2]
    double gV      = NUCLEON_VECTOR_G/pow((1.0 + (qsq/(CUTV*CUTV))),2); // from Equation (4.14) of my thesis
    double gA      = NUCLEON_AXIAL_G/pow((1.0 + (qsq/(CUTA*CUTA))),2); // " " " " " "
    double gP      = (2*M_PROTON*gA)/(qsq + Mpionsq); // " " " " " "
    double gM      = MAGMOM*gV; // " " " " " "
    double ff      = (((gA*gP*qsq)/(3*M_PROTON)) - (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(12*Mprosq)))/(NUCLEON_AXIAL_G*NUCLEON_AXIAL_G); //Equation (4.13) of Charlie's thesis
    return ff;
  }


  //Integrand for the intgeral with respect to dq
  //Uses QAGIU for everything, works up to n,l,nr,lr = 30
  //Haven't  been tested for higher values
  double fq(double q, void *params)
  { 
    struct fqparams * p = (struct fqparams *)params;
    int n                  = (p->n);
    int l                  = (p->l);
    int np                 = (p->np);
    int lp                 = (p->lp);
    double hw              = (p->hw);
    std::string transition = (p->transition); 
    int rho                = (p->rho);
    double Eclosure        = (p->Eclosure);
    std::string src        = (p->src); 
    std::map<std::string, std::function<double(double)> > FFlist = {{"F",FermiFormFactor},{"GT",GTFormFactor},{"T",TensorFormFactor}}; // make a map from std::string to a function which takes in q and returns a form factor     double ff  = FFlist[transition](q);
    double ff = FFlist[transition](q);
    double RBME = rbme(n, l, np, lp, hw, rho, q, src);
    return (q/(q + Eclosure))*ff*RBME;
  }

  // //Performs the integrals over q 
  double integrate_dq(int n, int l, int np, int lp, double hw, std::string transition, int rho, double Eclosure, std::string src)
  {
    gsl_set_error_handler_off(); // set standard GSL error handling off, so I can do my own
    int statusq; // this will hold the GSL error status
    std::string errstrq="7828877-1: qagiu dq:"; // a unique GSL error string
    size_t limitq = 10000; // the number of double precision intervals that can be held by the workspace below
    gsl_integration_workspace * workspacer = gsl_integration_workspace_alloc(limitq); // GSL: "One workspace may be used multiple times...", "...workspaces should be allocated on a per-thread basis."
    //Define param for the integrand function
    struct fqparams fqparams = {n,l,np,lp,hw,transition,rho,Eclosure,src};
    gsl_function Fq;
    Fq.function = &fq; // fq is defined above
    Fq.params = &fqparams;
    double r1 = 0; // integrate over the range [r1,Inf)
    double epsabsq = 1e-7; // absolute error tolerance for the integration
    double epsrelq = 1e-4; // relative error tolerance for the adaptive integration
    double val = 0; // qagiu integration result
    double abserrq; // qagiu integration error estimate
    double q = 0;
    statusq = gsl_integration_qagiu(&Fq,r1,epsabsq,epsrelq,limitq,workspacer,&val,&abserrq); // perform the qagiu integration (over r)
    #pragma omp critical
    {
      GSLerror(errstrq,statusq,limitq,epsabsq,epsrelq,abserrq,val,n,l,np,lp,q); // check if the GSL integration worked
    }
    gsl_integration_workspace_free(workspacer); // free the allocated memory for the integration workspace, since GSL is written in C
    return val;
  }


  

  // Hashtag for the IntList cache
  uint64_t IntHash(uint64_t n, uint64_t l, uint64_t np, uint64_t lp)
  {
    return   (n  << 21)
           + (l  << 15)
           + (np <<  8)
           +  lp;
  }

  // Inverse Hashtag for the IntList cache
  void IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp)
  {
    n  = (key >> 21) & 0x7FL;
    l  = (key >> 15) & 0x3FL;
    np = (key >>  8) & 0x7FL;
    lp = (key      ) & 0xFFL;
  }


  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals(int e2max, double hw, std::string transition, int rho, double Eclosure, std::string src)
  {
    IMSRGProfiler profiler;
    double t_start_pci = omp_get_wtime(); // profiling (s)
    //Map declaration to cache required integrals
    std::unordered_map<uint64_t,double> IntList;
    std::cout<<"calculating integrals wrt dq and dr..."<<std::endl;
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    int maxlp = e2max;
    std::vector<uint64_t> KEYS;
    if (transition == "F" or transition == "GT")
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=0; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
          for (int np=tempminnp; np<=maxnp; np++)
          {
            //long int key = IntHash(n,l,np);
            uint64_t key = IntHash(n,l,np,l);
            KEYS.push_back(key);
            IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
          }
        }
      }
    }
    else if (transition == "T")
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=0; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
          for (int np=tempminnp; np<=maxnp; np++)
          {
            int tempminlp = (n==np ? l : 0); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
            for (int lp=tempminlp; lp<=maxlp; lp++)
            {
              //long int key = IntHash(n,l,np,lp);
              uint64_t key = IntHash(n,l,np,lp);
              KEYS.push_back(key);
              IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }
          }
        }
      }
    }
    std::cout << " In " << __func__ << "  " << transition << "   size of KEYS = " << KEYS.size() << std::endl;
    #pragma omp parallel for schedule(dynamic, 1)// this works as long as the gsl_function handle is within this for-loop
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp;
      IntUnHash(key, n,l,np,lp);
      IntList[key] = integrate_dq(n,l,np,lp,hw,transition,rho,Eclosure,src); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
      //Uncomment if you want to verify integrals values
      // std::stringstream intvalue;
      // intvalue<<"n = "<<n<<", l = "<<l<<", np = "<<np<<", lp = "<<lp<<"; I = "<<IntList[key];
      // std::cout<<intvalue.str()<<std::endl; 
    }
    std::cout<<"...done calculating the integrals"<<std::endl;
    std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
    return IntList;
  }


  /// Get an integral from the IntList cache or calculate it (parallelization dependent)
//  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp, double hw, std::string transition, int rho, double Eclosure, std::string src, std::unordered_map<uint64_t,double> IntList)
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp, double hw, std::string transition, int rho, double Eclosure, std::string src, std::unordered_map<uint64_t,double>& IntList)
  {
//    int maxn = e2max/2;  // SRS commented out. not used.
    int maxl = e2max;
    int maxnp = e2max/2;
    int maxlp = e2max;
    int order1 = decimalgen(n,l,np,lp,maxl,maxnp,maxlp);
    int order2 = decimalgen(np,lp,n,l,maxl,maxnp,maxlp); // notice I was careful here with the order of maxl,maxnp,maxlp to make proper comparison
    if (order1 > order2)
    {
      std::swap(n,np); // using symmetry IntHash(n,l,np,lp) = IntHash(np,lp,n,l)
      std::swap(l,lp); // " " " " "
    }
    //long int key = IntHash(n,l,np,lp); // if I ever get that version working...
    uint64_t key = IntHash(n,l,np,lp);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double integral;
      integral = integrate_dq(n, l, np, lp, hw, transition, rho, Eclosure, src);
      if (omp_get_num_threads() >= 2)
      {
        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl;
        std::cout << "   I shouldn't be here in GetIntegral(" << n << ", " << l << ", " << np << ", " << lp << "):   key = " << key << "   integral = " << integral << std::endl;
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d):   key =%llx   integral=%f\n",n,l,np,lp,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;
      return integral;
    }
  }

  /// Hashtag for the T6jList cache
  uint64_t T6jHash(uint64_t l1, uint64_t L1, uint64_t R, uint64_t L2, uint64_t l2)
  {
    return   (l1 << 40)
           + (L1 << 30)
           + (R  << 20)
           + (L2 << 10)
           +  l2;
  }


  /// Inverse Hashtag for the T6jList cache
  void T6jUnHash(uint64_t key, uint64_t &l1, uint64_t &L1, uint64_t &R, uint64_t &L2, uint64_t &l2)
  {
    l1 = (key >> 40) & 0x3FFL;
    L1 = (key >> 30) & 0x3FFL;
    R  = (key >> 20) & 0x3FFL;
    L2 = (key >> 10) & 0x3FFL;
    l2 = (key      ) & 0x3FFL;
  }


  /// NDBD member function which caches the relevant 6j-symbols for the Tensor M0nu component
  /// { l1 L1 R }
  /// { L2 l2 2 } NOTE: the "2" in the bottom-right entry
  std::unordered_map<uint64_t,double> PreCalculateM0nuT6j(int e2max)
  {
    IMSRGProfiler profiler;
    double t_start_pctsj = omp_get_wtime(); // profiling (s)
    std::unordered_map<uint64_t,double> T6jList;
    std::cout<<"calculating 6j-symbols for the Tensor component..."<<std::endl;
    std::vector<uint64_t> KEYS;
    for (int l1=0; l1<=e2max; l1++)
    {
      for (int l2=0; l2<=e2max; l2++)
      {
        for (int L1=0; L1<=e2max; L1++)
        {
          for (int L2=0; L2<=e2max; L2++)
          {
            for (int R=0; R<=e2max; R++)
            {
              uint64_t key = T6jHash(l1,L1,R,L2,l2);
              KEYS.push_back(key);
              T6jList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }
          }
        }
      }
    }
    #pragma omp parallel for schedule(dynamic,1)
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t l1,l2,L1,L2,R;
      uint64_t key = KEYS[i];
      T6jUnHash(key,l1,L1,R,L2,l2);
      T6jList[key] = AngMom::SixJ(l1,L1,R,L2,l2,2);
    }
    std::cout<<"...done calculating the 6j-symbols"<<std::endl;
    std::cout<<"T6jList has "<<T6jList.bucket_count()<<" buckets and a load factor "<<T6jList.load_factor()
      <<", estimated storage ~= "<<((T6jList.bucket_count() + T6jList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalcT6j"] += omp_get_wtime() - t_start_pctsj; // profiling (r)
    return T6jList;
  }


  /// Get a 6j from the T6jList cache or calculate it (parallelization dependent)
  /// { l1 L1 R }
  /// { L2 l2 2 } NOTE: the "2" in the bottom-right entry
//  double GetM0nuT6j(int l1, int L1, int R, int L2, int l2, std::unordered_map<uint64_t,double> T6jList)
  double GetM0nuT6j(int l1, int L1, int R, int L2, int l2, std::unordered_map<uint64_t,double>& T6jList)
  {
    uint64_t key = T6jHash(l1,L1,R,L2,l2);
    auto it = T6jList.find(key);
    if (it != T6jList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double sixj = AngMom::SixJ(l1,L1,R,L2,l2,2);
      if (omp_get_num_threads() >= 2)
      {
        std::cout << "DANGER!!!!!!!  Updating T6jList inside a parellel loop breaks thread safety!" << std::endl;
        std::cout << "   I shouldn't be here in GetT6j(" << l1 << ", " << L1 << ", " << R << ", " << L2 << ", " << l2 << "):   key = " << key << "   sixj = " << sixj << std::endl;
//        printf("DANGER!!!!!!!  Updating T6jList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetT6j(%d, %d, %d, %d, %d):   key =%lx   sixj=%f\n",l1,L1,R,L2,l2,key,sixj);
        exit(EXIT_FAILURE);
      }
      T6jList[key] = sixj;
      return sixj;
    }
  }



  Operator GamowTeller(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "GT";
    int rho = 0;
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    std::cout<<"     reduced            =  "<<reduced<<std::endl;
    M0nuGT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = (R0/HBARC)*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/PI; // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, rho, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuGT_TBME.profiler.timer["M0nuGT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuGT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
      if (reduced == false)
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "R")
      {
        Jhat = sqrt(2*J + 1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
//      #pragma omp parallel for schedule(dynamic,100) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state orbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
            int Seval;
            Seval = 2*S*(S + 1) - 3;
            for (int L = abs(la-lb); L <= la+lb; L++) // ...and sum over angular momentum coupled to l_a and l_b
            {
              if ((abs(lc-ld) <= L) and (L <= lc+ld))
              {
                double bulk = 0; // this is a bulk product which will only be worth calculating if the Moshinsky brackets are non-zero
                double bulkas = 0; // (anti-symmetric part)
                double sumMT = 0; // for the Moshinsky transformation
                double sumMTas = 0; // (anti-symmetric part)
                int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
                int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
                  double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                  double intnpr = npr - floor(npr); // to check if npr is an integer
                  if ((npr >= 0) and (intnpr == 0))
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                        double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                        double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                        double nNJab = normab*modelspace.GetNineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
                        double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                        double nNJcd = normcd*modelspace.GetNineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
                        double nNJdc = normcd*modelspace.GetNineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
                        bulk = Seval*nNJab*nNJcd; // bulk product of the above
                        bulkas = Seval*nNJab*nNJdc; // (anti-symmetric part)
                        double Df,Di,asDi,integral; // these *will* be set below
                        #pragma omp critical // SRS WHY SHOULD THIS BE CRITICAL? IT *SHOULD* BE THREAD SAFE...
                        {
                          Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                          Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                          asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                          integral = GetM0nuIntegral(e2max,nr,lr,npr,lr,hw,transition,rho,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                        }
                        sumMT += Df*Di*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*integral; // (anti-symmetric part)
                      } // end of for-loop over: lr
                    } // end of for-loop over: Ncom
                  } // end of if: npr \in \Nat_0
                } // end of for-loop over: nr
                sumLS += bulk*sumMT; // perform the LS-coupling sum
                sumLSas += bulkas*sumMTas; // (anti-symmetric part)
              } // end of if: |lc - ld| <= L <= lc + ld
            } // end of for-loop over: L
          } // end of for-loop over: S
          double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuGT_TBME.profiler.timer["M0nuGT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    std::cout<<"...done calculating M0nu TBMEs"<<std::endl;
    M0nuGT_TBME.profiler.timer["M0nuGT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuGT_TBME.profiler.timer["M0nuGT_Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuGT_TBME;
  }

  Operator Fermi(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "F";
    int rho = 0;
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuF_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    std::cout<<"     reduced            =  "<<reduced<<std::endl;
    M0nuF_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = (R0/HBARC)*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/PI; // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, rho, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuF_TBME.profiler.timer["M0nuF_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuF_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      std::cout<<chbra<<" "<<chket<<std::endl;
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      double J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
      if (reduced == false)
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "R")
      {
        Jhat = sqrt(2*J + 1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
//      #pragma omp parallel for schedule(dynamic,100) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
            for (int L = abs(la-lb); L <= la+lb; L++) // ...and sum over angular momentum coupled to l_a and l_b
            {
              if ((abs(lc-ld) <= L) and (L <= lc+ld))
              {
                double bulk = 0; // this is a bulk product which will only be worth calculating if the Moshinsky brackets are non-zero
                double bulkas = 0; // (anti-symmetric part)
                double sumMT = 0; // for the Moshinsky transformation
                double sumMTas = 0; // (anti-symmetric part)
                int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
                int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
                  double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                  double intnpr = npr - floor(npr); // to check if npr is an integer
                  if ((npr >= 0) and (intnpr == 0))
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                        double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                        double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                        double nNJab = normab*modelspace.GetNineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
                        double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                        double nNJcd = normcd*modelspace.GetNineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
                        double nNJdc = normcd*modelspace.GetNineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
                        bulk = nNJab*nNJcd; // bulk product of the above
                        bulkas = nNJab*nNJdc; // (anti-symmetric part)
                        double Df,Di,asDi,integral; // these *will* be set below
                        #pragma omp critical
                        {
                          Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                          Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                          asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                          integral = GetM0nuIntegral(e2max,nr,lr,npr,lr,hw,transition,rho,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                        }
                        sumMT += Df*Di*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*integral; // (anti-symmetric part)
                      } // end of for-loop over: lr
                    } // end of for-loop over: Ncom
                  } // end of if: npr \in \Nat_0
                } // end of for-loop over: nr
                sumLS += bulk*sumMT; // perform the LS-coupling sum
                sumLSas += bulkas*sumMTas; // (anti-symmetric part)
              } // end of if: |lc - ld| <= L <= lc + ld
            } // end of for-loop over: L
          } // end of for-loop over: S
          double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuF_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuF_TBME.profiler.timer["M0nuF_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    std::cout<<"...done calculating M0nu TBMEs"<<std::endl;
    M0nuF_TBME.profiler.timer["M0nuF_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuF_TBME.profiler.timer["M0nuF_adpt_Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuF_TBME;
  }


  Operator Tensor(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "T";
    int rho = 2;
    // run through the initial set-up routine
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    std::cout<<"     reduced            =  "<<reduced<<std::endl;
    M0nuT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = (R0/HBARC)*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = 2*(2*Rnuc)/PI; // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> T6jList = PreCalculateM0nuT6j(e2max); // pre-calculate the expected 6j-symbols, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, rho, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency 
    // now we'll precalculate and cache the factors involving the Clebsch-Goradan coefficients, see "factSH" within the TMT
    double SQRT6 = sqrt(6);
    double listSH[PairFN(e2max,e2max)];
    int maxlr = e2max;
    int maxlpr = e2max;
    #pragma omp parallel for schedule(dynamic,1)
    for (int lr=0; lr<=maxlr; lr++)
    {
      for (int lpr=0; lpr<=maxlpr; lpr++)
      {
        int tempPF = PairFN(lr,lpr);
        if ((abs(lr-2) <= lpr and lpr <= lr+2) and ((lr+lpr)%2 == 0))
        {
          listSH[tempPF] = SQRT6*AngMom::CG(lr,0,2,0,lpr,0)*sqrt(2*lr + 1);
        }
        else
        {
          listSH[tempPF] = 0; // just to be safe, amirite!?
        }
      }
    }
    M0nuT_TBME.profiler.timer["M0nuT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    double Seval = 2.0*sqrt(5.0); // eigenvalue of S for "T", where the only non-zero one is for Sf = Si = 1
//    for (auto& itmat : M0nuT_TBME.TwoBody.MatEl)
    std::vector<int> chbralist,chketlist;
    for (auto& itmat : M0nuT_TBME.TwoBody.MatEl)
    {
      chbralist.push_back( itmat.first[0] );
      chketlist.push_back( itmat.first[1] );
    }
    size_t nch = chbralist.size();
     #pragma omp parallel for schedule(dynamic,1) 
    for (size_t ich=0; ich<nch; ich++)
    {
      int chbra = chbralist[ich];
      int chket = chketlist[ich];
//      int chbra = itmat.first[0]; // grab the channel count from auto
//      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat =  reduced  ?  sqrt(2*J+1.0)  :  1.0;

      t_start_omp = omp_get_wtime(); // profiling (s)
//      #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky(), PreCalcT6j, and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra=0; ibra<nbras; ibra++)
      {
        Ket& bra = tbc_bra.GetKet(ibra); // get the final state = <ab|
        int ia = bra.p; // get the integer label a
        int ib = bra.q; // get the integer label b
        Orbit& oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit& ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        for (int iket=0; iket<nkets; iket++)
        {
          Ket& ket = tbc_ket.GetKet(iket); // get the initial state = |cd>
          int ic = ket.p; // get the integer label c
          int id = ket.q; // get the integer label d
          Orbit& oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit& od = modelspace.GetOrbit(id); // get the |d> state orbit
          int na = oa.n; // this is just...
          int nb = ob.n;
          int nc = oc.n;
          int nd = od.n;
          int la = oa.l;
          int lb = ob.l;
          int lc = oc.l;
          int ld = od.l;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience

          int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
          int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          int Lf_min = std::max( std::abs(la-lb), std::abs(J-1));
          int Lf_max = std::min( la+lb, J+1);
          for (int Lf = Lf_min; Lf<=Lf_max; Lf++) // sum over angular momentum coupled to l_a and l_b
          {
            double tempLf = (2*Lf + 1); // just for the lines below
            double normab = sqrt(tempLf*3*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
            double nNJab = normab*modelspace.GetNineJ(la,lb,Lf,0.5,0.5,1,ja,jb,J); // the normalized 9j-symbol out front
            int Li_min = std::max( std::abs(lc-ld), std::max( std::abs(J-1),std::abs(Lf-2)));
            int Li_max = std::min( lc+ld, std::min( J+1, Lf+2) );
            for (int Li = Li_min; Li <= Li_max; Li++) // sum over angular momentum coupled to l_c and l_d
            {
              double tempLi = (2*Li + 1); // " " " " "
              double factLS = modelspace.phase(1 + J + Li)*GetM0nuT6j(Lf,1,J,1,Li,T6jList); // factor from transforming from jj-coupling to ls-coupling
              double normcd = sqrt(tempLi*3*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
              double nNJcd = normcd*modelspace.GetNineJ(lc,ld,Li,0.5,0.5,1,jc,jd,J); // the second normalized 9j-symbol
              double nNJdc = normcd*modelspace.GetNineJ(ld,lc,Li,0.5,0.5,1,jd,jc,J); // (anti-symmetric part)
              double tempMT = sqrt(tempLf*tempLi); // the hat factor of Lf and Li
              double bulk = Seval*factLS*nNJab*nNJcd*tempMT; // bulk product of the above
              double bulkas = Seval*factLS*nNJab*nNJdc*tempMT; // (anti-symmetric part)
              if ( std::abs(bulk)<1e-9  and std::abs(bulkas)<1e-9) continue;

              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)

              for ( int lr=0; lr<=eps_ab; lr++)
              {
                for (int nr=0; nr<=(eps_ab-lr)/2; nr++ )
                {
                  int Lam_min = std::abs( Lf-lr ); // triangle condition
                  Lam_min += ( eps_ab+Lam_min+lr )%2; //  Moshinsky trans. conserves parity
                  int Lam_max = std::min( Lf+lr, std::min( eps_ab-lr-2*nr, eps_cd) ); // triangle condition and energy conservation
                  for (int Lam=Lam_min; Lam<=Lam_max; Lam+=2)
                  {
                    int Ncom = (eps_ab - 2*nr - lr - Lam)/2;
                    double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,Lf); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                    
                    for (int npr=0; npr<=(eps_cd-2*Ncom-Lam); npr++)
                    {
                        int lpr = eps_cd-2*Ncom-Lam-2*npr;
                        if (  (lpr+lr)%2 !=0 ) continue;
                        if ( (std::abs(Lam-lpr)>Li)  or ( (Lam+lpr)<Li) ) continue;

                        double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nc,lc,nd,ld,Li); // " " " "
                        double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nd,ld,nc,lc,Li); // (anti-symmetric part)
                        double integral = GetM0nuIntegral(e2max,nr,lr,npr,lpr,hw,transition,rho,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the NDBD class
                        double factMT = modelspace.phase(Li + lr + Lam)*GetM0nuT6j(Lf,lr,Lam,lpr,Li,T6jList); // factor from TMT
//                        
                        double tempfactMTSH = factMT*listSH[PairFN(lr,lpr)]; // just for the lines below, listSH from the spherical harmonics of order 2
                        sumMT += Df*Di*tempfactMTSH*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*tempfactMTSH*integral; // (anti-symmetric part)

                      } // end of for-loop over: lpr
                    } // end of for-loop over: Ncom
                  } // end of for-loop over: nr
                } // end of for-loop over: lr

              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
            } // end of for-loop over: Li
          } // end of for-loop over: Lf
          double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          // std::cout<<"na ="<<na<<"  nb = "<<nb<<"  nc = "<<nc<<"  nd = "<<nd<<std::endl;
          // std::cout<<"la ="<<la<<"  lb = "<<lb<<"  lc = "<<lc<<"  ld = "<<ld<<std::endl;
          M0nuT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuT_TBME.profiler.timer["M0nuT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    M0nuT_TBME.profiler.timer["M0nuT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuT_TBME.profiler.timer["M0nuT _Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuT_TBME;
  }
}
