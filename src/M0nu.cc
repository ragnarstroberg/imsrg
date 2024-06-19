#include "M0nu.hh"
#include "AngMom.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

using namespace PhysConst;

/// Namespace for functions related to neutrinoless double beta decay
namespace M0nu
{

  /// Converts (a,b,c,d) in base (maxa+1,maxb+1,maxc+1,maxd+1) to an ordered decimal integer
  /// eg: for (maxa,maxb,maxc,maxd) = (1,1,1,1) decimalgen would convert (0,1,1,0) to 6, ie: a binary to decimal converter
  /// NOTE: to make comparisons between such decimals, keep (maxa,maxb,maxc,maxd) consistent
  /// PS: I tots made this thing up, did some testing and it seemed to work...
  /// ... hopefully it's good, but could be a source of weird bugs if anything comes up with missing integrals wrt dq, see GetIntegral(...)
  int decimalgen(int a, int b, int c, int d, int maxb, int maxc, int maxd)
  {
    int coeff1 = maxd + 1;
    int coeff2 = (maxc + 1)*coeff1;
    int coeff3 = (maxb + 1)*coeff2;
    return coeff3*a + coeff2*b + coeff1*c + d; // eg: (0,1,1,0) -> (2*2*2)*0 + (2*2)*1 + (2)*1 + 0 = 6
  }

  //Return the phase due to a certain quantity, i.e. (-1)^x
  int phase(int x)
  {
    return x%2==0 ? 1 : -1;
  }

  //Radial wave function written in momentum space
  double HO_Radial_psi_mom(int n, int l, double hw, double p)
  {
     double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON*0.5));
     double x = p*b;
     double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / SQRTPI / gsl_sf_doublefact(2*n+2*l+1) * pow(b,3.0) );
     //double Norm = sqrt(2*gsl_sf_fact(n)/gsl_sf_gamma(n+l+0.5)/(2*n+2*l+1)*b*b*b);
     double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
     double psi = phase(n)*Norm * pow(x,l) * exp(-x*x*0.5) * L;
     return psi;
  }


  //The following four functions are used in the defintions of the neutrino potentials
  //They are implemented from Equation (4.14) of Charlie Payne's thesis. 
  //Note that they take the momentum square as arquments.
  double gv_func(double qsq)
  {
    double gV = NUCLEON_VECTOR_G/pow((1.0 + (qsq/(CUTV*CUTV))),2); // from Equation (4.14) of my thesis
    return gV;
  }

  double ga_func(double qsq)
  {
    double gA = NUCLEON_AXIAL_G/pow((1.0 + (qsq/(CUTA*CUTA))),2);
    return gA;
  }

  double gp_func(double qsq)
  {
    double gP = (2*M_PROTON*ga_func(qsq))/(qsq + M_PION_CHARGED*M_PION_CHARGED);
    return gP;
  }

  double gm_func(double qsq)
  {
    double gM = MAGMOM*gv_func(qsq);
    return gM;
  }


  /// Form factors of the neutrino potential of Gamow-Teller transition
  /// Implemented from Equation (4.12) of Charlie Payne's thesis
  double GTFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
    double Mprosq  = M_PROTON*M_PROTON; // the proton mass squared [MeV^2]
    double gA      = ga_func(qsq);
    double gP      = gp_func(qsq); 
    double gM      = gm_func(qsq); 
    double ff      = ((gA*gA) - ((gA*gP*qsq)/(3*M_PROTON)) + (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(6*Mprosq)))/(NUCLEON_AXIAL_G*NUCLEON_AXIAL_G); 
    return ff;
  }

  /// Form factors of the neutrino potential of Gamow-Teller transition
  /// Implemented from Equation (4.11) of Charlie Payne's thesis
  double FermiFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
    double gV      = gv_func(qsq);
    double ff      = (gV*gV)/(NUCLEON_VECTOR_G*NUCLEON_VECTOR_G); //Equation (4.11) of Charlie's thesis
    return ff;
  }

  /// Form factors of the neutrino potential of Gamow-Teller transition
  /// Implemented from Equation (4.13) of Charlie Payne's thesis
  double TensorFormFactor(double q)
  {
    double qsq     = q*q; // q squared [MeV^2]
    double Mprosq  = M_PROTON*M_PROTON; // the proton mass squared [MeV^2]
    double gA      = ga_func(qsq);
    double gP      = gp_func (qsq);
    double gM      = gm_func(qsq);
    double ff      = (((gA*gP*qsq)/(3*M_PROTON)) - (pow(gP*qsq,2)/(12*Mprosq)) + ((gM*gM*qsq)/(12*Mprosq)))/(NUCLEON_AXIAL_G*NUCLEON_AXIAL_G); //Equation (4.13) of Charlie's thesis
    return -ff;
  }


  //Abbreviation used in the definition of the partial waves expansion in
  //K. Erkelenz, R. Alzetta, and K. Holinde, Nucl. Phys. A 176, 413 (1971).
  //The factor of z^l from the paper is ignored since l=0 in all the needed cases.
  double A(double p, double pp, int J, double Eclosure, std::string transition, gsl_integration_glfixed_table * t, int norm, int size)
  {
    std::map<std::string, std::function<double(double)> > FFlist = {{"F",FermiFormFactor},{"GT",GTFormFactor},{"T",TensorFormFactor}}; // make a map from std::string to a function which takes in q and returns a form factor
    double A = 0;
    for (int i=0;i<size;i++)
    {
      double xi;
      double wi;
      gsl_integration_glfixed_point(-1,1,i,&xi,&wi,t);
      double q = sqrt(pp*pp+p*p-2*p*pp*xi);
      double ff = FFlist[transition](q);
      double temp = wi*ff/(q*q+q*Eclosure)*gsl_sf_legendre_Pl(J,xi);
      if (norm == 1)
      {
        temp *= HBARC*HBARC/(q*q);
      }
      A += temp;
    }
    return A*PI*HBARC*HBARC;//HBARC*HBARC factor to make A in fm^2
  }

  //Following functions are for the precalculation and the caching
  //of the A function for the partial wave decomposition
  //for the p values later used in the quadrature
  uint64_t AHash(int index_p,int index_pp, int J, int norm)
  { 
     return   (((uint64_t)(index_p)) << 30)
            + (((uint64_t)(index_pp)) << 20)
            + (((uint64_t)(J)) << 10)
            + ((uint64_t)(norm));
  }

  void AUnHash(uint64_t key, uint64_t& index_p, uint64_t& index_pp, uint64_t& J, uint64_t& norm)
  {
     index_p  = (key >> 30) & 0x3FFL;
     index_pp = (key >> 20) & 0x3FFL;
     J        = (key >> 10) & 0x3FFL;
     norm     = (key      ) & 0x3FFL;
  }

  //Precomputes the values of the A functions  and caches them for efficiency
  std::unordered_map<uint64_t,double> PrecalculateA(int e2max,double Eclosure, std::string transition, int npoints)
  {
    std::unordered_map<uint64_t,double> AList;
    int Jmax = e2max+1;
    std::vector<uint64_t> KEYS;
    int maxnorm = 0;
    if (transition=="T")
    {
      maxnorm = 1;
    }
    for (int norm = 0; norm <= maxnorm; norm++)
    {
        for (int J=0; J<=Jmax;J++)
        {
          for (int index_p = 0 ; index_p< npoints; index_p++)
          {
            for (int index_pp = index_p; index_pp<npoints; index_pp++)
            {
              uint64_t key = AHash(index_p,index_pp,J,norm);
              KEYS.push_back(key);
              AList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }//end of for loop over j 
          }//end of for loop over i
        }//end of for loop over J
    }//end of for loop over norm
    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(npoints);
    int size = 100;
    gsl_integration_glfixed_table * t2 = gsl_integration_glfixed_table_alloc(size);
    #pragma omp parallel for schedule(dynamic, 1)// this works as long as the gsl_function handle is within this for-loop
    for (size_t n=0; n<KEYS.size(); n++)
    {
      uint64_t key = KEYS[n];
      uint64_t index_p,index_pp,J,norm;
      AUnHash(key, index_p,index_pp,J,norm);
      double p,pp,wi,wj;
      gsl_integration_glfixed_point(0,25,index_p,&p,&wi,t);
      gsl_integration_glfixed_point(0,25,index_pp,&pp,&wj,t);
      AList[key] = A(p*HBARC,pp*HBARC,J,Eclosure,transition,t2,norm, size); 
    }
    gsl_integration_glfixed_table_free(t);
    gsl_integration_glfixed_table_free(t2);
    return AList;
  }

  //Fetch the cached value of a specific A functions 
  double GetA(int index_p, int index_pp, int J,int norm, std::unordered_map<uint64_t,double> &AList)
  {
    double A;
    if (index_p >index_pp) std::swap(index_p,index_pp);
    uint64_t key = AHash(index_p,index_pp,J,norm);
    auto it = AList.find(key);
    if (it != AList.end()) // return what we've found
    {
      A =it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl
                  << "   I shouldn't be here in GetA(" << index_p << ", " << index_pp << ", " << J << ", " << norm << "):    "
                  << "key = " << std::hex << key << std::dec << std::endl;
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetA(%d, %d, %d, %d):   key =%llx",index_p,index_pp,J,norm,key);
        exit(EXIT_FAILURE);
    }
    return A; 
  }


  //Function for the partial wave decomposition given by Eq 4.20 of 
  //K. Erkelenz, R. Alzetta, and K. Holinde, Nucl. Phys. A 176, 413 (1971).
  //We treat the fermi and gt as the same since we include the spin part for the gt when
  //computing the TBMEs. Since for all the possible cases used a value of J equal to l in
  //A(p,pp,J,l), we pass L rather than J. Furthermore, l = 0 by definition
  //for this case. Fetches the values of A from the precomputed list for efficiency.
   double W_fermi_gt(double p, double pp,int index_p, int index_pp, int l, int lp, int J, std::unordered_map<uint64_t,double>& AList)
  {
    double W = 0;
    W += 2*GetA(index_p,index_pp,l,0,AList);
    return W;
  }

  //Function for the partial wave decomposition given by Eq 4.24 of 
  //K. Erkelenz, R. Alzetta, and K. Holinde, Nucl. Phys. A 176, 413 (1971).
  //Only the cases for S=1 have been implemented since this decay can only happen with S=1
  //The operator is composed of 3x(tensor)-(spin) of the paper. Fetches the values of A from the precomputed list for efficiency.
  double W_tensor(double p, double pp,int index_p, int index_pp, int l, int lp, int J, std::unordered_map<uint64_t,double>& AList)
  {
    double W = 0;
    if (lp==J && J==l)
    {
       W += 6*((pp*pp+p*p)*GetA(index_p,index_pp,J,1,AList)-2*pp*p/(2*J+1)*(J*GetA(index_p,index_pp,J+1,1,AList)+(J+1)*GetA(index_p,index_pp,J-1,1,AList)))-2*GetA(index_p,index_pp,J,0,AList);
    }
    else if (lp ==J-1 && l==lp)
    {
      W += 6*((pp*pp+p*p)*GetA(index_p,index_pp,J-1,1,AList)-2*p*pp*GetA(index_p,index_pp,J,1,AList))/(2*J+1)-2*GetA(index_p,index_pp,J-1,0,AList);
    }
    else if (lp == J+1 && l==lp)
    {
      W += 6*(-(pp*pp+p*p)*GetA(index_p,index_pp,J+1,1,AList)+2*p*pp*GetA(index_p,index_pp,J,1,AList))/(2*J+1)-2*GetA(index_p,index_pp,J+1,0,AList);
    }
    else if (lp == J+1 && l==J-1)
    {
      W += -(12*sqrt(J*(J+1))/(2*J+1))*(p*p*GetA(index_p,index_pp,J+1,1,AList)+pp*pp*GetA(index_p,index_pp,J-1,1,AList)-2*p*pp*GetA(index_p,index_pp,J,1,AList));//Minus sign in front is to match with Takayuki's number but I'm not sure why it should be there...
    }
    else if (lp == J-1 && l==J+1)
    {
      W += -(12*sqrt(J*(J+1))/(2*J+1))*(p*p*GetA(index_p,index_pp,J-1,1,AList)+pp*pp*GetA(index_p,index_pp,J+1,1,AList)-2*p*pp*GetA(index_p,index_pp,J,1,AList));//"""
    }
    else
    {
      std::cout<<"Problem..."<<std::endl;
    }
    return W; 
  }


  //Integrand for the integral of the neutrino potential over mometum space. 
  double fq(double p, double pp,int index_p, int index_pp, int n, int l, int np, int lp, int J, double hw, std::string transition, double Eclosure, std::string src, std::unordered_map<uint64_t,double>& AList)
  { 
    
    std::map<std::string, std::function<double(double,double,int,int,int,int,int,std::unordered_map<uint64_t,double>& A)> > WList = {{"F",W_fermi_gt},{"GT",W_fermi_gt},{"T",W_tensor}};
    double W = WList[transition](p,pp,index_p,index_pp,l,lp,J,AList);
    return p*p*pp*pp*HO_Radial_psi_mom(n,l,hw,p)*W*HO_Radial_psi_mom(np,lp,hw,pp);
  }

  //Performs the integral over moemntum space using Gauss-Legendre quadrature quadrilateral rule
  //Uses the fact that fq(p,pp,...) = fq(pp,p,...) to reduce the number of points to evaluate
  //This is only valid however if l=lp and therefore we do the sum with all the points if 
  //that l !=lp
  double integrate_dq(int n, int l, int np, int lp, int J, double hw, std::string transition, double Eclosure, std::string src, int npoints, gsl_integration_glfixed_table * t, std::unordered_map<uint64_t,double>& AList)
  { 
    double I = 0;
    for (int i = 0 ; i< npoints; i++)
    {
      double xi;
      double wi;
      gsl_integration_glfixed_point(0,25,i,&xi,&wi,t);
      // if (l!=lp)
      // {
        for (int j = 0; j<npoints; j++)
        {
          double xj;
          double wj;
          gsl_integration_glfixed_point(0,25,j,&xj,&wj,t);
          I += wi*wj*fq(xi,xj,i,j,n,l,np,lp,J,hw,transition,Eclosure,src,AList); 
        }
      // }
      // else
      // {
      //   for (int j = i; j<=i; j++)
      //   {
      //     double xj;
      //     double wj;
      //     gsl_integration_glfixed_point(0,25,j,&xj,&wj,t);
      //     I += wi*wj*fq(xi,xj,i,j,n,l,np,lp,J,hw,transition,Eclosure,src,AList);
      //   }
      //   for (int j = i+1; j<npoints; j++)
      //   {
      //     double xj;
      //     double wj;
      //     gsl_integration_glfixed_point(0,25,j,&xj,&wj,t);
      //     I += 2*wi*wj*fq(xi,xj,i,j,n,l,np,lp,J,hw,transition,Eclosure,src,AList); 
      //   }
      // }
      
    }
    return I;
  }


  //Functions to precompute the integrals of the neutrino potentials for efficiency
  //when computing the TBMEs
  uint64_t IntHash(int n, int l, int np, int lp, int J)
  {
     return   (((uint64_t)(n)) << 40)
            + (((uint64_t)(l)) << 30)
            + (((uint64_t)(np)) << 20)
            + (((uint64_t)(lp)) << 10)
            + ((uint64_t)(J));
  }

  void IntUnHash(uint64_t key, uint64_t& n, uint64_t& l, uint64_t& np, uint64_t& lp, uint64_t& J)
  {
     n  = (key >> 40) & 0x3FFL;
     l  = (key >> 30) & 0x3FFL;
     np = (key >> 20) & 0x3FFL;
     lp = (key >> 10) & 0x3FFL;
     J  = (key      ) & 0x3FFL;
  }

  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals(int e2max, double hw, std::string transition, double Eclosure, std::string src)
  {
    IMSRGProfiler profiler;
    double t_start_pci = omp_get_wtime(); // profiling (s)
    std::unordered_map<uint64_t,double> IntList;
    int size=500;
    std::cout<<"calculating integrals wrt dp and dpp..."<<std::endl;
    std::unordered_map<uint64_t,double> AList = PrecalculateA(e2max,Eclosure,transition,size);
    std::cout<<"Done precomputing A's."<<std::endl;
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    std::vector<uint64_t> KEYS;
    if (transition == "F" or transition == "GT")
    {
      for (int S = 0; S<=1; S++)
      {
        for (int n=0; n<=maxn; n++)
        {
          for (int l=0; l<=maxl; l++)
          {
            int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
            for (int np=tempminnp; np<=maxnp; np++)
            {
              int minJ = abs(l-S);
              int tempmaxJ = l+S;
              for (int J = minJ; J<= tempmaxJ; J++)
              {
                uint64_t key = IntHash(n,l,np,l,J);
                KEYS.push_back(key);
                IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
              }
            }
          }
        }
      }
      
    }
    else if (transition == "T")
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=1; l<=maxl; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
          //int tempminnp = 0;
          for (int np=tempminnp; np<=maxnp; np++)
          {
            int tempminlp = (n==np ? l : 1); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
            int maxlp = std::min(l+2,maxl);
            for (int lp = tempminlp; lp<=maxlp; lp++)
            { 
              if ((abs(lp-l) != 2) and (abs(lp-l) != 0)) continue;
              int minJ = std::max(abs(l-1),abs(lp-1));
              int tempmaxJ = std::min(l+1,lp+1);
              for (int J = minJ; J<= tempmaxJ; J++)
              {
                uint64_t key = IntHash(n,l,np,lp,J);
                KEYS.push_back(key);
                IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
              }
            }
          }
        }
      }
    }
    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
    #pragma omp parallel for schedule(dynamic, 1)// this works as long as the gsl_function handle is within this for-loop
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp,J;
      IntUnHash(key, n,l,np,lp,J);
      IntList[key] = integrate_dq(n,l,np,lp,J,hw,transition,Eclosure,src,size,t,AList); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
      // Uncomment if you want to verify integrals values
      // std::stringstream intvalue;
      // intvalue<<n<<" "<<l<<" "<<np<<" "<<lp<<" "<<J<<" "<<IntList[key]<<std::endl;
      // std::cout<<intvalue.str();
      
    }
    gsl_integration_glfixed_table_free(t);
    
    std::cout<<"...done calculating the integrals"<<std::endl;
    std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
    return IntList;
  }

  //Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp,int J, double hw, std::string transition, double Eclosure, std::string src, std::unordered_map<uint64_t,double> &IntList)
  {
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
    //std::cout<<"n ="<<n<<", l ="<<l<<", np = "<<np<<", S = "<<S<<", J = "<<J<<std::endl;
    uint64_t key = IntHash(n,l,np,lp,J);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double integral;
      int size = 500;
      gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
      std::unordered_map<uint64_t,double> AList = PrecalculateA(e2max,Eclosure,transition,size);
      integral = integrate_dq(n, l, np, lp,J,hw, transition, Eclosure, src,size, t,AList);
      gsl_integration_glfixed_table_free(t);
      if (omp_get_num_threads() >= 2)
      {

        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl
                  << "   I shouldn't be here in GetIntegral(" << n << ", " << l << ", " << np << ", " << lp << ", " <<J << "):    "
                  << "key = " << std::hex << key << std::dec
                  << "   integral = " << integral << std::endl;
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d, %d):   key =%llx   integral=%f\n",n,l,np,lp,J,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;
      
      return integral;
    }
  }


  /// Gamow Teller operator for neutrinoless double beta decay. Opertor is written in momentum space and takes the form
  /// \f{equation} 
  ///       O_{GT}(\bold{q'}) = \frac{R}{2\pi^2}\frac{h_{GT}(q')}{q'(q'+E_c)}(\boldsymbol{\sigma_1} \cdot \boldsymbol{\sigma_2}) \tau_{1,+}\tau_{2,+}
  /// \f}
  /// Where \f$ h_{GT} \f$ is the neutrino potenital  impletmented above and \f$ E_c \f$ is the closure energy.
  /// Operator is then evaluated in the lab frame oscialltor basis.
  /// More detail on how to obatin the form of the operator can be found in https://drive.google.com/file/d/1QaMfuvQ7I3NM5h_ppjyIPsap3SWmCe2o/view?usp=sharing
  /// and on how to evaluate in lab frame in https://drive.google.com/file/d/1C6E2HnzSJ1bzMoIKWfH1GaZKjqaAEluG/view?usp=sharing
  Operator GamowTeller(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "GT";
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    M0nuGT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
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
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat; // set below based on "reduced" variable
      if (reduced == false)
      {
        Jhat = 1.0; // for non-reduced elements, to compare with JE
      }
      else //if (reduced == "true")
      { 
        Jhat = sqrt(2*J+1); // the hat factor of J
      }
      t_start_omp = omp_get_wtime(); // profiling (s)
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
          int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
          int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
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
            int L_min = std::max(abs(la-lb),abs(lc-ld));
            int L_max = std::min(la+lb,lc+ld);
            for (int L = L_min; L <= L_max; L++) // ...and sum over orbital angular momentum
            {
              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)              
              double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
              double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
              double nNJab = normab*AngMom::NineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
              double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
              double nNJcd = normcd*AngMom::NineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
              double nNJdc = normcd*AngMom::NineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
              double bulk = Seval*nNJab*nNJcd; // bulk product of the above
              double bulkas = Seval*nNJab*nNJdc; // (anti-symmetric part)
              int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
              for (int nr = 0; nr <= tempmaxnr; nr++)
              {
                double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                {
                  int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                      double integral = 0;
                      double normJrel;
                      double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                      double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                      double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                      int minJrel= abs(lr-S);
                      int maxJrel = lr+S;
                      for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                      {
                        normJrel = sqrt((2*Jrel+1)*(2*L+1))*phase(L+lr+S+J)*AngMom::SixJ(Lam,lr,L,S,J,Jrel);
                        integral += normJrel*normJrel*GetM0nuIntegral(e2max,nr,lr,npr,lr,Jrel,hw,transition,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                      }//end of for-loop over Jrel
                      sumMT += Df*Di*integral; // perform the Moshinsky transformation
                      sumMTas += Df*asDi*integral; // (anti-symmetric part)
                    } // end of for-loop over: lr
                  } // end of for-loop over: Ncom
                } // end of if: npr \in \Nat_0
              } // end of for-loop over: nr
              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
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

  /// Fermi operator for neutrinoless double beta decay. Opertor is written in momentum space and takes the form
  /// \f[ 
  ///       O_{F}(\bold{q'}) = \frac{R}{2\pi^2}\frac{h_{F}(q')}{q'(q'+E_c)} \tau_{1,+}\tau_{2,+}
  /// \f]
  /// Where \f$ h_{F} \f$ is the neutrino potenital  impletmented above and \f$ E_c \f$ is the closure energy.
   /// Operator is then evaluated in the lab frame oscialltor basis.
  /// More detail on how to obatin the form of the operator can be found in https://drive.google.com/file/d/1QaMfuvQ7I3NM5h_ppjyIPsap3SWmCe2o/view?usp=sharing
  /// and on how to evaluate in lab frame in https://drive.google.com/file/d/1C6E2HnzSJ1bzMoIKWfH1GaZKjqaAEluG/view?usp=sharing
  Operator Fermi(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "F";
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuF_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    std::cout<<"     reduced            =  "<<reduced<<std::endl;
    M0nuF_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuF_TBME.profiler.timer["M0nuF_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuF_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
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
          int eps_ab = 2*na + la + 2*nb + lb; // for conservation of energy in the Moshinsky brackets
          int eps_cd = 2*nc + lc + 2*nd + ld; // for conservation of energy in the Moshinsky brackets
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          for (int S=0; S<=1; S++) // sum over total spin...
          {
            int L_min = std::max(std::max(abs(la-lb),abs(lc-ld)), abs(J-S));
            int L_max = std::min(std::min(la+lb,lc+ld),J+S);
            for (int L = L_min; L<=L_max;L++)
            {
              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)
              double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
              double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
              double nNJab = normab*AngMom::NineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
              double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
              double nNJcd = normcd*AngMom::NineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
              double nNJdc = normcd*AngMom::NineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
              double bulk = nNJab*nNJcd; // bulk product of the above
              double bulkas = nNJab*nNJdc; // (anti-symmetric part)
              int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
              for (int nr = 0; nr <= tempmaxnr; nr++)
              {
                double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                if ((npr >= 0) and (npr == floor(npr)))
                {
                  int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                      double integral = 0;
                      double normJrel;
                      double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                      double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                      double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                      int minJrel= abs(lr-S);
                      int maxJrel = lr+S;
                      for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                      {
                        normJrel = sqrt((2*Jrel+1)*(2*L+1))*phase(L+lr+S+J)*AngMom::SixJ(Lam,lr,L,S,J,Jrel);
                        integral += normJrel*normJrel*GetM0nuIntegral(e2max,nr,lr,npr,lr,Jrel,hw,transition,Eclosure,src,IntList); // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                      }//end of for-loop over Jrel
                      sumMT += Df*Di*integral; // perform the Moshinsky transformation
                      sumMTas += Df*asDi*integral; // (anti-symmetric part)
                    } // end of for-loop over: lr
                  } // end of for-loop over: Ncom
                } // end of if: npr \in \Nat_0
              } // end of for-loop over: nr
              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
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




  /// Tensor operator for neutrinoless double beta decay. Opertor is written in momentum space and takes the form
  /// \f[ 
  ///        O_{T}(\bold{q'}) = -\frac{R}{2\pi^2}\frac{h_{T}(q')}{q'(q'+E_c)}[3(\boldsymbol{\sigma_1}\cdot\boldsymbol{\hat{q'}})(\boldsymbol{\sigma_2}\cdot\boldsymbol{\hat{q'}})-(\boldsymbol{\sigma_1} \cdot \boldsymbol{\sigma_2})] \tau_{1,+}\tau_{2,+}
  /// \f]
  /// Where \f$ h_{T} \f$ is the neutrino potenital  impletmented above and \f$ E_c \f$ is the closure energy. The minus factor up front as been included in the neutrino potential for simplicity.
   /// Operator is then evaluated in the lab frame oscialltor basis.
  /// More detail on how to obatin the form of the operator can be found in https://drive.google.com/file/d/1QaMfuvQ7I3NM5h_ppjyIPsap3SWmCe2o/view?usp=sharing
  /// and on how to evaluate in lab frame in https://drive.google.com/file/d/1C6E2HnzSJ1bzMoIKWfH1GaZKjqaAEluG/view?usp=sharing
  Operator Tensor(ModelSpace& modelspace, double Eclosure, std::string src)
  {
    bool reduced = true;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    std::string transition = "T";
    // run through the initial set-up routine
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    std::cout<<"     reduced            =  "<<reduced<<std::endl;
    M0nuT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
    const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals(e2max,hw,transition, Eclosure, src); // pre-calculate the needed integrals over dq and dr, for efficiency
    M0nuT_TBME.profiler.timer["M0nuT_1_sur"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout<<"calculating M0nu TBMEs..."<<std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto& itmat : M0nuT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
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
      #pragma omp parallel for schedule(dynamic,1) // need to do: PreCalculateMoshinsky(), PreCalcT6j, and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
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
          double sumLS = 0; // for the wave functions decomposition
          double sumLSas = 0; // (anti-symmetric part)
          int S = 1;
          int Lf_min = std::max(std::abs(la-lb), std::abs(J-S));
          int Lf_max = std::min(la+lb, J+S);
          for (int Lf = Lf_min; Lf<=Lf_max; Lf++) // sum over angular momentum coupled to l_a and l_b
          {
            double normab = sqrt((2*Lf+1)*(2*S+1)*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
            double nNJab = normab*AngMom::NineJ(la,lb,Lf,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
            int Li_min = std::max(std::abs(lc-ld), std::max(std::abs(J-S), std::abs(Lf-2)));
            int Li_max = std::min( lc+ld, std::min( J+S, Lf+2) );
            for (int Li = Li_min; Li <= Li_max; Li++) // sum over angular momentum coupled to l_c and l_d
            { 
              double sumMT = 0; // for the Moshinsky transformation
              double sumMTas = 0; // (anti-symmetric part)
              double normcd = sqrt((2*Li+1)*(2*S+1)*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
              double nNJcd = normcd*AngMom::NineJ(lc,ld,Li,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
              double nNJdc = normcd*AngMom::NineJ(ld,lc,Li,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)

              double bulk = nNJab*nNJcd; // bulk product of the above
              double bulkas = nNJab*nNJdc; // (anti-symmetric part)
              for ( int lr=1; lr<=eps_ab; lr++)
              {
                for (int nr=0; nr<=(eps_ab-lr)/2; nr++)
                {
                  int tempmaxNcom = std::min((eps_ab-2*nr-lr)/2, eps_cd);
                  for (int Ncom=0; Ncom<=tempmaxNcom; Ncom++ )
                  {
                    int Lam = eps_ab - 2*nr - lr - 2*Ncom;
                    if ( (Lam+2*Ncom) > eps_cd ) continue;
                    if ( (std::abs(Lam-lr)>Lf)  or ( (Lam+lr)<Lf) ) continue;
                    if ( (lr+Lam+eps_ab)%2>0 ) continue;
                    for (int npr=0; npr<=(eps_cd-2*Ncom-Lam)/2; npr++)
                    {
                        int lpr = eps_cd-2*Ncom-Lam-2*npr;
                        if (  (lpr+lr)%2 >0 ) continue;
                        if (lpr<1) continue;
                        if ( (std::abs(Lam-lpr)>Li)  or ( (Lam+lpr)<Li) ) continue;
                        double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,Lf); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                        double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nc,lc,nd,ld,Li); // " " " "
                        double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lpr,nd,ld,nc,lc,Li);// (anti-symmetric part)
                        double integral = 0;
                        double normJrel, normJrelp;
                        int minJrel = std::max(abs(lr-S),abs(lpr-S));
                        int maxJrel = std::min(lr+S,lpr+S);
                        for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                        {
                          if ( (std::abs(J-Jrel)>Lam)  or ( (Jrel+J)<Lam) ) continue;
                          normJrel  = sqrt((2*Jrel+1)*(2*Lf+1))*phase(Lf+lr+J+S)*AngMom::SixJ(Lam,lr,Lf,S,J,Jrel);
                          normJrelp = sqrt((2*Jrel+1)*(2*Li+1))*phase(Li+lpr+J+S)*AngMom::SixJ(Lam,lpr,Li,S,J,Jrel);
                          integral += normJrel*normJrelp*GetM0nuIntegral(e2max,nr,lr,npr,lpr,Jrel,hw,transition,Eclosure,src,IntList);
                        }
                        sumMT += Df*Di*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*integral; // (anti-symmetric part)
                      } // end of for-loop over: lpr
                    } // end of for-loop over: Ncom
                  } // end of for-loop over: nr
                } // end of for-loop over: lr
              sumLS += bulk*sumMT; // perform the LS-coupling sum
              sumLSas += bulkas*sumMTas; // (anti-symmetric part)
            } // end of for-loop over: Li
          } // end of for-loop over: Lf        
          // double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*sumLS; // compute the final matrix element, anti-symmetrize          
          double Mtbme = asNorm(ia,ib)*asNorm(ic,id)*prefact*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          M0nuT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
      M0nuT_TBME.profiler.timer["M0nuT_3_omp"] += omp_get_wtime() - t_start_omp; // profiling (r)
    } // end of for-loop over: auto
    M0nuT_TBME.profiler.timer["M0nuT_2_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    M0nuT_TBME.profiler.timer["M0nuT _Op"] += omp_get_wtime() - t_start; // profiling (r)
    return M0nuT_TBME;
  }

  Operator DGT_Op(ModelSpace& modelspace)
  {
    Operator DGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    DGT_TBME.SetHermitian(); // it should be Hermitian
    modelspace.PreCalculateMoshinsky();
    // create the TBMEs of DGT
    // auto loops over the TBME channels and such
    std::cout<<"calculating DGT TBMEs..."<<std::endl;
    for (auto& itmat : DGT_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0]; // grab the channel count from auto
      int chket = itmat.first[1]; // " " " " " "
      TwoBodyChannel& tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel& tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets(); // get the number of bras
      int nkets = tbc_ket.GetNumberKets(); // get the number of kets
      int J = tbc_bra.J; // NOTE: by construction, J := J_ab == J_cd := J'
      //double Jhat = 1;
      double Jhat = sqrt(2*J + 1); // the hat factor of J
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
          int eps_ab = 2*na+la+2*nb+lb;
          int eps_cd = 2*nc+lc+2*nd+ld;
          double ja = oa.j2/2.0;
          double jb = ob.j2/2.0;
          double jc = oc.j2/2.0;
          double jd = od.j2/2.0; // ...for convenience
          double sumLS = 0; // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0; // (anti-symmetric part)
          if (la+lb == lc+ld)
          {
            for (int S=0; S<=1; S++) // sum over total spin...
            {
              int Seval;
              Seval = 2*S*(S + 1) - 3;
              int L_min = std::max(abs(la-lb),abs(lc-ld));
              int L_max = std::min(la+lb,lc+ld);
              for (int L = L_min; L <= L_max; L++) // ...and sum over orbital angular momentum
              {
                double sumMT = 0; // for the Moshinsky transformation
                double sumMTas = 0; // (anti-symmetric part)              
                double tempLS = (2*L + 1)*(2*S + 1); // just for efficiency, only used in the three lines below
                double normab = sqrt(tempLS*(2*ja + 1)*(2*jb + 1)); // normalization factor for the 9j-symbol out front
                double nNJab = normab*AngMom::NineJ(la,lb,L,0.5,0.5,S,ja,jb,J); // the normalized 9j-symbol out front
                double normcd = sqrt(tempLS*(2*jc + 1)*(2*jd + 1)); // normalization factor for the second 9j-symbol
                double nNJcd = normcd*AngMom::NineJ(lc,ld,L,0.5,0.5,S,jc,jd,J); // the second normalized 9j-symbol
                double nNJdc = normcd*AngMom::NineJ(ld,lc,L,0.5,0.5,S,jd,jc,J); // (anti-symmetric part)
                double bulk = Seval*nNJab*nNJcd; // bulk product of the above
                double bulkas = Seval*nNJab*nNJdc; // (anti-symmetric part)
                int tempmaxnr = floor((eps_ab - L)/2.0); // just for the limits below
                for (int nr = 0; nr <= tempmaxnr; nr++)
                {
                  double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of my thesis
                  if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of my thesis
                        double integral = 0;
                        double normJrel;
                        double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                        double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                        double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                        int minJrel= abs(lr-S);
                        int maxJrel = lr+S;
                        for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                        {
                          normJrel = sqrt((2*Jrel+1)*(2*L+1))*phase(Lam+lr+S+J)*AngMom::SixJ(Lam,lr,L,S,J,Jrel);
                          integral += normJrel*normJrel*1; // grab the pre-calculated integral wrt dq and dr from the IntList of the modelspace class
                        }//end of for-loop over Jrel
                        sumMT += Df*Di*integral; // perform the Moshinsky transformation
                        sumMTas += Df*asDi*integral; // (anti-symmetric part)
                      } // end of for-loop over: lr
                    } // end of for-loop over: Ncom
                  } // end of if: npr \in \Nat_0
                } // end of for-loop over: nr
                sumLS += bulk*sumMT; // perform the LS-coupling sum
                sumLSas += bulkas*sumMTas; // (anti-symmetric part)
              } // end of for-loop over: L
            } // end of for-loop over: S
          }//end of if la+lb==lc+ld
          
          double Mtbme = 2*asNorm(ia,ib)*asNorm(ic,id)*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas); // compute the final matrix element, anti-symmetrize
          DGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      }// end of for-loop over: ibra
    } // end of for-loop over: auto
    std::cout<<"...done calculating DGT TBMEs"<<std::endl;
    return DGT_TBME;
  }

}// end namespace M0nu


