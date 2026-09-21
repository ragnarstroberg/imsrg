#include "M0nu.hh"
#include "AngMom.hh"
#include "Pwd.hh"
// #include "Helicity.hh"
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

  //Two-body radial wave function written in momentum space
  double HO_Radial_psi_mom(int n, int l, double hw, double p)
  {
     double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON*0.5));
     double x = p*b;
     double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / SQRTPI / gsl_sf_doublefact(2*n+2*l+1) * pow(b,3.0) );
     double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
     double psi = phase(n)*Norm * pow(x,l) * exp(-x*x*0.5) * L;
     return psi;
  }

  /// Different form component required for the form factors
  /// separated by their vector, vector-axial, induced pseudo-scalar...
  /// parts as defined in JHEP12(2018)097, eq A.7 , A.8 and A.9
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

  double hF_VV(double qsq)
  {
    return gv_func(qsq)*gv_func(qsq)/(NUCLEON_VECTOR_G*NUCLEON_VECTOR_G);
  }

  double hGT_AA(double qsq)
  {
    return ga_func(qsq)*ga_func(qsq) / (NUCLEON_AXIAL_G * NUCLEON_AXIAL_G);
  }

  double hGT_AP(double qsq)
  {
    double gA = ga_func(qsq);
    double gP = gp_func(qsq);
    return -(gA * gP * qsq) / (3 * M_NUCLEON * NUCLEON_AXIAL_G * NUCLEON_AXIAL_G);
  }

  double hGT_PP(double qsq)
  {
    double Mnucsq = M_NUCLEON * M_NUCLEON; // the proton mass squared [MeV^2]
    double gP = gp_func(qsq);
    return (pow(gP * qsq, 2) / (12 * Mnucsq * NUCLEON_AXIAL_G * NUCLEON_AXIAL_G));
  }

  double hGT_MM(double qsq)
  {
    double Mnucsq = M_NUCLEON * M_NUCLEON; // the proton mass squared [MeV^2]
    double gM = gm_func(qsq);
    return ((gM * gM * qsq) / (6 * Mnucsq* NUCLEON_AXIAL_G * NUCLEON_AXIAL_G));
  }

  double hT_AA(double qsq)
  {
    return hGT_AA(qsq);
  }

  double hT_AP(double qsq)
  {
    //extra factor of -1 from fourier trans r -> p space
    return hGT_AP(qsq);
  }

  double hT_PP(double qsq)
  {
    //extra factor of -1 from fourier trans r -> p space
    return hGT_PP(qsq);
  }

  double hT_MM(double qsq)
  {
    // extra factor of -1 from fourier trans r -> p space
    return -hGT_MM(qsq)/2;
  }


  /// Form factors of the neutrino potential of Gamow-Teller transition
  double GTFormFactor(double qsq)
  {
    //  qsq = q squared [MeV^2]
    double ff = hGT_AA(qsq)+hGT_AP(qsq)+hGT_PP(qsq)+hGT_MM(qsq); 
    return ff;
  }

  /// Form factors of the neutrino potential of Fermi transition
  double FermiFormFactor(double qsq)
  {
    //  qsq = q squared [MeV^2]
    double ff = hF_VV(qsq) ;
    return ff;
  }

  /// Form factors of the neutrino potential of tensor transition
  double TensorFormFactor(double qsq)
  {
    //  qsq = q squared [MeV^2]
    double ff = hT_AP(qsq)+hT_PP(qsq)+hT_MM(qsq);
    return ff;
  }

  double potential_closure(double p, double pp, double z, double Eclosure, std::function<double(double)> formfactor)
  {
    double  q = sqrt(p*p+pp*pp-2*p*pp*z)*HBARC;
    return HBARC*HBARC*formfactor(q*q)/(q*q+Eclosure*q);
  }

  /// Potential in position space is given as
  ///\f { equation }
  ///       H_{\alpha}(r) = \frac{2R}{2\pi m^2}\int_0^\infty q^2 h_{\alpha}(q) j_{\lambda}(q*r)
  ///\f
  /// which becomes
  ///\f { equation }
  ///       H_{\alpha}(q) = (-1))^\lambda\frac{R}{2\pi^2 m^2} h_{\alpha}(q)
  ///\f
  /// after doing the fourrier transform. For simplicity, we put the prefactor in the operator since lambda = 0 for
  /// GT and F and 2 for T. The mass to use is somewhat ambiguous. In Phys. Rev. C 91, 024613, they use m_e*m_p but
  /// in JHEP12(2018)097, they use m_pi^2. In the end it is only introduced to make unit consitent so
  /// it should simply be matched to the other factor in the formula.
  double potential_heavy_neutrino(double p, double pp, double z, std::function<double(double)> formfactor)
  {
    double qsq = (p * p + pp * pp - 2 * p * pp * z) * HBARC*HBARC;
    return  formfactor(qsq);
  }

  double integrate_dq(int n, int l, int np, int lp, int S, int J, double hw, PWD &pwd)
  { 
    int momentum_mesh_size = pwd.getMomentumMeshSize();
    gsl_integration_glfixed_table *t = pwd.getMomentumMesh();
    int max_momentum = pwd.getMaxMomentum();
    double I = 0;
    double wf, wfp;
    for (int i = 0 ; i<momentum_mesh_size ; i++)
    {
      double xi;
      double wi;
      gsl_integration_glfixed_point(0,max_momentum,i,&xi,&wi,t);
      wf = HO_Radial_psi_mom(n, l, hw, xi);
      for (int j = 0; j <momentum_mesh_size; j++)
      {
        double xj;
        double wj;
        gsl_integration_glfixed_point(0,max_momentum,j,&xj,&wj,t);
        wfp = HO_Radial_psi_mom(np, lp, hw, xj);
        double W = pwd.getW(xi, xj, i, j, S, l, lp, J) ;
        I += wi * wj * xi * xi * xj * xj * wf * wfp * W;
        // std::stringstream intvalue;
        // intvalue << S << ", " << l << ", " << lp << ", " << J << ", " << xi*HBARC << ", " << xj*HBARC <<", " << 2*pi*hW<< std::endl;
        // std::cout<<intvalue.str();
      }
    }
    return I;
  }



  //Functions to precompute the integrals of the neutrino potentials for efficiency
  //when computing the TBMEs
  uint64_t IntHash(int n, int l, int np, int lp, int S, int J)
  {
    return   (((uint64_t)(n)) << 50) 
           + (((uint64_t)(l)) << 40) 
           + (((uint64_t)(np)) << 30) 
           + (((uint64_t)(lp)) << 20) 
           + ((uint64_t)(S) << 10) 
           + ((uint64_t)(J));
  }

  void IntUnHash(uint64_t key, uint64_t &n, uint64_t &l, uint64_t &np, uint64_t &lp, uint64_t &S, uint64_t &J)
  {
    n = (key >> 50) & 0x3FFL;
    l = (key >> 40) & 0x3FFL;
    np = (key >> 30) & 0x3FFL;
    lp = (key >> 20) & 0x3FFL;
    S = (key >> 10) & 0x3FFL;
    J = (key)&0x3FFL;
  }

  std::unordered_map<uint64_t, double> PreCalculateMomentumSpaceIntegrals(int e2max, int Smin, int Smax, double hw, int Lrank, PWD &pwd)
  {
    std::unordered_map<uint64_t, double> IntList;
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    std::vector<uint64_t> KEYS;
    for (int S = Smin; S<=Smax; S++)
    {
      for (int n=0; n<=maxn; n++)
      {
        for (int l=0; l<=maxl-2*n; l++)
        {
          int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
          for (int np=tempminnp; np<=maxnp; np++)
          {
            if (Lrank == 0)
            {
              for (int np = tempminnp; np <= maxnp; np++)
              {
                int minJ = abs(l - S);
                int tempmaxJ = l + S;
                for (int J = minJ; J<=tempmaxJ; J++ )
                {
                  uint64_t key = IntHash(n, l, np, l, S, J);
                  KEYS.push_back(key);
                  IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
                }
              }
            }
            else
            {
              int tempminlp = (n == np ? l : 0); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
              int maxlp = std::min(l + Lrank, maxl);
              for (int lp = tempminlp; lp <= maxlp; lp++)
              {
                if ((abs(lp - l) != 2) and (abs(lp - l) != 0)) continue;
                int minJ = std::max(abs(l - S), abs(lp - S));
                int tempmaxJ = std::min(l + S, lp + S);
                for (int J = minJ; J <= tempmaxJ; J++)
                {
                  uint64_t key = IntHash(n, l, np, lp, S, J);
                  KEYS.push_back(key);
                  IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
                }
              }
            }
          }
        }
      }
    }
    sort(KEYS.begin(), KEYS.end());
    KEYS.erase(unique(KEYS.begin(), KEYS.end()), KEYS.end());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp,S,J;
      IntUnHash(key, n,l,np,lp,S,J);
   
      double integral = integrate_dq(n,l, np,lp, S , J, hw, pwd);
      IntList[key] = integral; // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
    }
    std::cout << "IntList has " << IntList.bucket_count() << " buckets and a load factor " << IntList.load_factor()
              << ", estimated storage ~= " << ((IntList.bucket_count() + IntList.size()) * (sizeof(size_t) + sizeof(void *))) / (1024.0 * 1024.0 * 1024.0) << " GB" << std::endl; // copied from (RS)
    // profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci;                                                                                                 // profiling (r)
    return IntList;
  }


  // Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetM0nuIntegral(int e2max, int n, int l, int np, int lp, int J, int S, std::unordered_map<uint64_t, double> &IntList)
  {
    int maxl = e2max;
    int maxnp = e2max / 2;
    int maxlp = e2max;
    int order1 = decimalgen(n, l, np, lp, maxl, maxnp, maxlp);
    int order2 = decimalgen(np, lp, n, l, maxl, maxnp, maxlp); // notice I was careful here with the order of maxl,maxnp,maxlp to make proper comparison
    if (order1 > order2)
    {
      std::swap(n, np); // using symmetry IntHash(n,l,np,lp) = IntHash(np,lp,n,l)
      std::swap(l, lp); // " " " " "
    }
    uint64_t key = IntHash(n, l, np, lp, S, J);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it->second;
    }
    else 
    {
      std::stringstream intvalue;
      intvalue<<"Integral value not precalculated..."<<std::endl;
      intvalue<<n<<", "<<l<<", "<<np<<", "<<lp<<", "<<S<<", "<<J<<std::endl;
      intvalue<<"Exiting program"<<std::endl;
      std::cout<<intvalue.str();
      exit(EXIT_FAILURE);
    }
  }

  /// Talmi-Moshinsky trasnform. Takes in the operator in relative frame and return them in the lab frame. Update value for both the 
  /// normal sum and the antisymmetric part.
  void TalmiMoshinkyTransform(ModelSpace &modelspace, double& sumTMT, double&sumTMTas, int na, int la, int nb, int lb, int nc, int lc, int nd, int ld, 
    int Li, int Lf, int S, int J, int rank, int e2max, std::unordered_map<uint64_t, double>& RelativeFrameOPList)
  {
    
    int eps_ab = 2*na+la+2*nb+lb;
    int eps_cd = 2*nc+lc+2*nd+ld;
    if ((abs(Lf-Li) > rank) or (Lf+Li < rank))
    {
      return;
    }
    int maxnr = floor((eps_ab-Lf)/2.0);
    for (int nr = 0; nr <= maxnr; nr++)
    {
      int maxNcom = maxnr - nr;
      for (int Ncom = 0; Ncom <= maxNcom; Ncom++)
      {
        int minlr = ceil((eps_ab-Lf)/2) - Ncom - nr;
        int maxlr = floor((eps_ab + Lf) / 2) - Ncom - nr;
        for (int lr = minlr; lr <= maxlr; lr++)
        {
          int Lam = eps_ab - 2*(Ncom+nr) - lr;
          if ((Lam + 2 * Ncom) > eps_cd) continue;
          if ((std::abs(Lam - lr) > Lf) or ((Lam + lr) < Lf)) continue;
          if ((lr + Lam + eps_ab) % 2 > 0) continue;
          int minlpr = std::max(std::abs(Lam-Li), std::abs(lr-rank));
          int maxlpr = std::min(Lam+Li, lr+rank);
          for (int lpr = minlpr; lpr <= maxlpr; lpr++)
          {
            if (abs(lr-lpr)%2 >0) continue;
            if ((eps_cd - eps_ab+lr-lpr)%2 > 0) continue;
            int npr = (eps_cd - eps_ab+lr-lpr)/2 + nr;
            if (npr < 0) continue;
            double Df = modelspace.GetMoshinsky(Ncom, Lam, nr, lr, na, la, nb, lb, Lf);     // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
            double Di = modelspace.GetMoshinsky(Ncom, Lam, npr, lpr, nc, lc, nd, ld, Li);   // " " " "
            double asDi = modelspace.GetMoshinsky(Ncom, Lam, npr, lpr, nd, ld, nc, lc, Li); // (anti-symmetric part)
            double integral = 0;
            double normJrel, normJrelp;
            int minJrel = std::max(abs(lr - S), abs(lpr - S));
            int maxJrel = std::min(lr + S, lpr + S);
            for (int Jrel = minJrel; Jrel <= maxJrel; Jrel++)
            {
              normJrel = sqrt((2 * Jrel + 1) * (2 * Lf + 1)) * phase(Lf + lr + J + S) * AngMom::SixJ(Lam, lr, Lf, S, J, Jrel);
              normJrelp = sqrt((2 * Jrel + 1) * (2 * Li + 1)) * phase(Li + lpr + J + S) * AngMom::SixJ(Lam, lpr, Li, S, J, Jrel);
              double intval = GetM0nuIntegral(e2max, nr, lr, npr, lpr, Jrel, S, RelativeFrameOPList);
              if (std::abs(intval) < 1e-16) continue;
              integral += normJrel * normJrelp * intval;
            }
            sumTMT += Df * Di * integral;     // perform the Moshinsky transformation
            sumTMTas += Df * asDi * integral; // (anti-symmetric part)
            
          }
        }
      }
    }
    return;
  }

  /// General format for a two-body scalar operator. This is based on the partial wave-decomposition from
  /// https://doi.org/10.1016/0375-9474(71)90279-X. Requires to pass a PWD class with the angular and momentum mesh 
  /// as well as the potential for each 6 possible cases for a scalar operator in momentum space predefined. The A integrals should also be 
  /// precalculated so that the maximal order of the Legendre polynomials to be evaluated can be set for efficiency. If it isn't, the A integrals
  /// will be calculated for all possible orders.
  Operator TwoBody_Scalar_operator(ModelSpace &modelspace, PWD &pwd, int Jrank, int Lrank, int Smin = 0, int Smax = 1)
  {
    double t_start, t_start_tbme;     // profiling (v)
    t_start = omp_get_wtime();                     // profiling (s)
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    Operator Op_TBME(modelspace, Jrank, 2, 0, 2);      // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    Op_TBME.SetHermitian();                        // it should be Hermitian
    Op_TBME.is_reduced = true;
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    if (pwd.getAsize() == 0)
    {
      std::cout << "Angular integrals for partial wave decomposition are not calcualted..." << std::endl;
      std::cout << "Computing angular integrals." << std::endl;
      pwd.calcA(e2max, 2);
      std::cout << "Done precomputing A integrals." << std::endl;
    }
    std::unordered_map<uint64_t, double> IntList = PreCalculateMomentumSpaceIntegrals(e2max, Smin, Smax, hw, Lrank, pwd); // pre-calculate the needed integrals over dp and dpp
    std::cout << "Done precomputing integrals." << std::endl;
    pwd.clearA();
    pwd.freeMomentumMesh();
    Op_TBME.profiler.timer["Op_PrecalculateMomentumIntegrals"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout << "calculating TBMEs..." << std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto &itmat : Op_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0];                                    // grab the channel count from auto
      int chket = itmat.first[1];                                    // " " " " " "
      TwoBodyChannel &tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel &tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets();                           // get the number of bras
      int nkets = tbc_ket.GetNumberKets();                           // get the number of kets
      int J = tbc_bra.J;                                             // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat = sqrt(2 * J + 1);                                 //We are computing reduced matrix elements
      #pragma omp parallel for schedule(dynamic, 1) // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);     // get the final state = <ab|
        int ia = bra.p;                      // get the integer label a
        int ib = bra.q;                      // get the integer label b
        Orbit &oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit &ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        int na = oa.n;                       // this is just...
        int nb = ob.n;
        int la = oa.l;
        int lb = ob.l;
        int eps_ab = 2 * na + la + 2 * nb + lb; // for conservation of energy in the Moshinsky brackets
        double ja = oa.j2 / 2.0;
        double jb = ob.j2 / 2.0;
        for (int iket = 0; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);     // get the initial state = |cd>
          int ic = ket.p;                      // get the integer label c
          int id = ket.q;                      // get the integer label d
          Orbit &oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit &od = modelspace.GetOrbit(id); // get the |d> state orbit
          int nc = oc.n;
          int nd = od.n;
          int lc = oc.l;
          int ld = od.l;
          int eps_cd = 2 * nc + lc + 2 * nd + ld; // for conservation of energy in the Moshinsky brackets
          double jc = oc.j2 / 2.0;
          double jd = od.j2 / 2.0;     // ...for convenience
          double sumLS = 0;            // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0;          // (anti-symmetric part)
          for (int S = Smin; S <= Smax; S++) // sum over total spin...
          {
            int Lf_min = std::max(std::abs(la - lb), std::abs(J - S));
            int Lf_max = std::min(la + lb, J + S);
            for (int Lf = Lf_min; Lf <= Lf_max; Lf++) // sum over angular momentum coupled to l_a and l_b
            {
              double normab = sqrt((2 * Lf + 1) * (2 * S + 1) * (2 * ja + 1) * (2 * jb + 1)); // normalization factor for the 9j-symbol out front
              double nNJab = normab * AngMom::NineJ(la, lb, Lf, 0.5, 0.5, S, ja, jb, J);      // the normalized 9j-symbol out front
              int Li_min = std::max(std::abs(lc - ld), std::max(std::abs(J - S), std::abs(Lf - Lrank)));
              int Li_max = std::min(lc + ld, std::min(J + S, Lf + Lrank));
              for (int Li = Li_min; Li <= Li_max; Li++) // sum over angular momentum coupled to l_c and l_d
              {

                double sumMT = 0;                                                               // for the Moshinsky transformation
                double sumMTas = 0;                                                             // (anti-symmetric part)
                double normcd = sqrt((2 * Li + 1) * (2 * S + 1) * (2 * jc + 1) * (2 * jd + 1)); // normalization factor for the second 9j-symbol
                double nNJcd = normcd * AngMom::NineJ(lc, ld, Li, 0.5, 0.5, S, jc, jd, J);      // the second normalized 9j-symbol
                double nNJdc = normcd * AngMom::NineJ(ld, lc, Li, 0.5, 0.5, S, jd, jc, J);      // (anti-symmetric part)
                double bulk = nNJab * nNJcd;   // bulk product of the above
                double bulkas = nNJab * nNJdc; // (anti-symmetric part)
                TalmiMoshinkyTransform(modelspace, sumMT, sumMTas, na, la, nb, lb, nc, lc, nd, ld,
                                       Li, Lf, S, J, Lrank, e2max, IntList);
                sumLS += bulk * sumMT;                                                                                         // perform the LS-coupling sum
                sumLSas += bulkas * sumMTas;                                                                                   // (anti-symmetric part)
              }                                                                                                                // end of for-loop over: Li
            }                                                                                                                  // end of for-loop over: Lf
          }                                                                                                                    // end of for-loop over: S                                                                                                                                                                                                                // end of for-loop over: S
          double Mtbme = asNorm(ia, ib) * asNorm(ic, id) * Jhat * (sumLS - phase(jc + jd - J) * sumLSas);           // compute the final matrix element, anti-symmetrize
          Op_TBME.TwoBody.SetTBME(chbra, chket, ibra, iket, Mtbme);                                                            // set the two-body matrix elements (TBME) to Mtbme
        }                                                                                                                      // end of for-loop over: iket
      }                                                                                                                        // end of for-loop over: ibra
    }                                                                                                                          // end of for-loop over: auto
    std::cout << "...done calculating M0nu TBMEs" << std::endl;
    Op_TBME.profiler.timer["Op_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    return Op_TBME;  
  }

  /// Gamow Teller operator for neutrinoless double beta decay. Opertor is written in momentum space and takes the form
  /// \f{equation}
  ///       O_{GT}(\bold{q'}) = \frac{R}{2\pi^2}\frac{h_{GT}(q')}{q'(q'+E_c)}(\boldsymbol{\sigma_1} \cdot \boldsymbol{\sigma_2}) \tau_{1,+}\tau_{2,+}
  /// \f}
  /// Where \f$ h_{GT} \f$ is the neutrino potenital  impletmented above and \f$ E_c \f$ is the closure energy.
  /// Operator is then evaluated in the lab frame oscialltor basis.
  /// More detail on how to obatin the form of the operator can be found in https://drive.google.com/file/d/1QaMfuvQ7I3NM5h_ppjyIPsap3SWmCe2o/view?usp=sharing
  /// and on how to evaluate in lab frame in https://drive.google.com/file/d/1C6E2HnzSJ1bzMoIKWfH1GaZKjqaAEluG/view?usp=sharing
  Operator GamowTeller(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [fm]
    const double prefact = Rnuc / (PI * PI);       // factor in-front of M0nu TBME, extra global factor of 2 since we use <p|\tau|n> = sqrt(2) [fm]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;                                       // Class for the partial wave decomposition
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp,  double z){return prefact*potential_closure(p,pp,z,Eclosure,formfactor);}, "spin-spin");
    pwd.calcA(e2max,0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuGT_TBME =  TwoBody_Scalar_operator(modelspace, pwd, 0, 0);
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
  Operator Fermi(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor)
  {
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp,  double z){return prefact*potential_closure(p,pp,z,Eclosure,formfactor);}, "central");
    pwd.calcA(e2max,0);
    pwd.freeAngularMesh();
    std::cout<<"Done precomputing A's."<<std::endl;
    Operator M0nuF_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 0);
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
    Operator Tensor(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor)
    {
      double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
      int e2max = modelspace.GetE2max(); // 2*emax
      int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
      const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [MeV^-1]
      const double prefact = Rnuc/(PI*PI); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
      modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
      PWD pwd;
      pwd.initializeAngularMesh(100);
      pwd.initializeMomentumMesh(500);
      pwd.setMaxMomentum(25);
      pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp, double z)
                       {
                         double qsq = p*p+pp*pp-2*p*pp*z;
                         return prefact*(3/qsq)*potential_closure(p, pp, z, Eclosure, formfactor); 
                        },
                        "tensor");
      pwd.setPotential([prefact,Eclosure, formfactor](double p, double pp, double z)
                       {
                         return -prefact*potential_closure(p, pp, z, Eclosure, formfactor); 
                        },
                       "spin-spin");
      pwd.calcA(e2max, 0);
      pwd.freeAngularMesh();
      std::cout << "Done precomputing A's." << std::endl;
      Operator M0nuT_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 2, 1);
      return M0nuT_TBME;
    }

  /// Contact operator for neutrinoless double beta decay. Opertor is written in momentum space. Only factor missing is the gvv/gA^2 which is multiplied to the NME
  /// at the end since gvv depends on the interactions. Values of gvv can be found in arXiv:2105.05415.
  /// The non-local  contact oprator in momentum space takes the form
  /// \f[
  ///        O_{C}(p,p') = - 4  \frac{R}{pi} exp((\frac{p}{\Lambda})^{2n})exp((\frac{p'}{\Lambda})^{2n})
  /// \f]
  /// Operator is then evaluated in the lab frame oscialltor basis.
  Operator Contact(ModelSpace& modelspace, double regulator_cutoff, int regulator_power)
  {
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double scale_ct = (M_NUCLEON*NUCLEON_AXIAL_G*NUCLEON_AXIAL_G*HBARC)/(4*F_PI*F_PI);
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = -Rnuc / (16 * SQRT2 * PI * PI * PI) * scale_ct * scale_ct; // factor in-front of M0nu TBME, includes the -2 from the -2gvv as well as the 1/2 from the convention used by Cirgilano and all when computing gvv. Only thing missing is gvv/gA^2
    std::string reg_type = "nonlocal";
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setRegulator(regulator_cutoff,regulator_power,reg_type);
    pwd.setPotential([prefact](double p, double pp, double z)
                      {
                        return prefact;//
                      },
                      "central");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuC_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 0, 0, 0);
    return M0nuC_TBME;
  }

  Operator GamowTellerHeavy(ModelSpace &modelspace,  std::string src, std::function<double(double)> formfactor)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [fm]
    const double prefact = Rnuc * HBARC*HBARC / (PI * PI* M_PROTON* M_ELECTRON);   // factor in-front of M0nu TBME, extra global factor of 2 since we use <p|\tau|n> = sqrt(2) [fm]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;                                       // Class for the partial wave decomposition
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setPotential([prefact, formfactor](double p, double pp, double z)
                     { return prefact * potential_heavy_neutrino(p, pp, z, formfactor); },
                     "spin-spin");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuGT_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 0);
    return M0nuGT_TBME;
  }

  Operator FermiHeavy(ModelSpace &modelspace, std::string src, std::function<double(double)> formfactor)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [fm]
    const double prefact = Rnuc*HBARC*HBARC/ (PI * PI * M_PROTON * M_ELECTRON); // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setPotential([prefact, formfactor](double p, double pp, double z)
                     { return prefact * potential_heavy_neutrino(p, pp, z, formfactor); },
                     "central");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuF_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 0);
    return M0nuF_TBME;
  }

  Operator TensorHeavy(ModelSpace &modelspace,  std::string src, std::function<double(double)> formfactor)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [MeV^-1]
    const double prefact = Rnuc *HBARC*HBARC / (PI * PI * M_PROTON * M_ELECTRON);       // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(500);
    pwd.setMaxMomentum(25);
    pwd.setPotential([prefact, formfactor](double p, double pp, double z)
                     {
                         double qsq = p*p+pp*pp-2*p*pp*z;
                         return prefact*(3/qsq)*potential_heavy_neutrino(p, pp, z, formfactor); },
                     "tensor");
    pwd.setPotential([prefact,  formfactor](double p, double pp, double z)
                     { return -prefact * potential_heavy_neutrino(p, pp, z, formfactor); },
                     "spin-spin");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuT_TBME = TwoBody_Scalar_operator(modelspace, pwd, 0, 2, 1);
    return M0nuT_TBME;
  }

  /// Double Gamow-Teller operator. It is simply the 0vbb GT operator with the neutrino potential set to 1.
  /// More infromation on the process can be found in 10.1103/PhysRevLett.120.142502.
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
                    // if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                    if (npr==nr)
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
                            integral += normJrel*normJrel*1;
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



  /// The following functions are to obtain radial distributions of the neutrinoless double beta decay matrix elements.


  double HO_Radial_psi(int n, int l, double hw, double r)
  {
    double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON*0.5) );
    double x = r/b;
    double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / SQRTPI / gsl_sf_doublefact(2*n+2*l+1) / pow(b,3.0) );
    double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
    double psi = Norm * pow(x,l) * exp(-x*x*0.5) * L;
    return psi;
  }

  double fq_radial_GT(double q, double Eclosure, double r12)
  {
    return gsl_sf_bessel_j0(q*r12)*q*GTFormFactor(q*HBARC)/(q+Eclosure/HBARC);
  }

  double integrate_dq_radial_GT(double Eclosure, double r12,  int npoints, gsl_integration_glfixed_table * t)
  {
    double I = 0;
    for (int i = 0 ; i< npoints; i++)
    {
      double xi;
      double wi;
      gsl_integration_glfixed_point(0,30,i,&xi,&wi,t);
      I += wi*fq_radial_GT(xi,Eclosure,r12);
    }
    return I;
  }

  std::unordered_map<uint64_t,double> PreCalculateM0nuIntegrals_R(int e2max, double hw, double Eclosure, double r12)
  {
    IMSRGProfiler profiler;
    double t_start_pci = omp_get_wtime(); // profiling (s)
    std::unordered_map<uint64_t,double> IntList;
    int size=1000;
    std::cout<<"calculating integrals wrt dq..."<<std::endl;
    int maxn = e2max/2;
    int maxl = e2max;
    int maxnp = e2max/2;
    std::vector<uint64_t> KEYS;
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
              uint64_t key = IntHash(n,l,np,l,S,J);
              KEYS.push_back(key);
              IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
            }
          }
        }
      }
    }

    gsl_integration_glfixed_table * t = gsl_integration_glfixed_table_alloc(size);
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp,S,J;
      IntUnHash(key, n,l,np,lp,S,J);
      IntList[key] = r12*r12*HO_Radial_psi(n, l, hw, r12)*HO_Radial_psi(np, lp, hw, r12)*integrate_dq_radial_GT(Eclosure,r12,size,t); // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
    }
    gsl_integration_glfixed_table_free(t);

    std::cout<<"...done calculating the integrals"<<std::endl;
    std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
    return IntList;
  }

  //Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetM0nuIntegral_R(int e2max, int n, int l, int np, int lp, int S, int J, double hw, double Eclosure, double r12, std::unordered_map<uint64_t,double> &IntList)
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
    // long int key = IntHash(n,l,np,lp); // if I ever get that version working...
    // std::cout<<"n ="<<n<<", l ="<<l<<", np = "<<np<<", S = "<<S<<", J = "<<J<<std::endl;
    uint64_t key = IntHash(n,l,np,lp,S,J);
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
      integral = r12*r12*HO_Radial_psi(n, l, hw, r12)*HO_Radial_psi(np, lp, hw, r12)*integrate_dq_radial_GT(Eclosure,r12,size,t);
      gsl_integration_glfixed_table_free(t);
      if (omp_get_num_threads() >= 2)
      {
        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl;
        std::cout << "   I shouldn't be here in GetIntegral(" << n << " , " << l << " , " << np << " , " << lp << " , " << S << " , " << J << "):   key = " << key << "   integral = " << integral << std::endl; 
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetIntegral(%d, %d, %d, %d, %d, %d):   key =%lx   integral=%f\n",n,l,np,lp,S,J,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;

      return integral;
    }
  }

  

  Operator GamowTeller_R(ModelSpace& modelspace, double Eclosure, double r12)
  {
    bool reduced = true;
    r12 =  r12*SQRT2;
    double t_start, t_start_tbme, t_start_omp; // profiling (v)
    t_start = omp_get_wtime(); // profiling (s)
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max(); // 2*emax
    Operator M0nuGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    M0nuGT_TBME.SetHermitian(); // it should be Hermitian
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]
    const double prefact = 4*Rnuc/(PI); // factor in-front of M0nu TBME, extra global 2 from <p|\tau|n> = sqrt(2) [fm]
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    std::unordered_map<uint64_t,double> IntList = PreCalculateM0nuIntegrals_R(e2max, hw, Eclosure, r12); // pre-calculate the needed integrals over dq and dr, for efficiency
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
                double npr = ((eps_cd - eps_ab)/2.0) + nr; // via Equation (4.73) of Charlie's thesis
                if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                {
                  int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                  for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                  {
                    int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                    int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                    for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                    {
                      int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of Charlie's thesis
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
                        integral += normJrel*normJrel*GetM0nuIntegral_R(e2max,nr,lr,npr,lr,S,Jrel,hw,Eclosure,r12,IntList); // grab the pre-calculated integral wrt dq from the IntList of the modelspace class
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

  Operator DGT_R(ModelSpace& modelspace, double r12)
  {
    r12 = r12*SQRT2;
    Operator DGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    DGT_TBME.SetHermitian(); // it should be Hermitian
    modelspace.PreCalculateMoshinsky();
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
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
      // double Jhat = 1;
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
                  // if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                  if (npr==nr)
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of Charlie's thesis
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
                          integral += normJrel*normJrel*r12*r12*HO_Radial_psi(nr, lr, hw, r12)*HO_Radial_psi(npr, lr, hw, r12);
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

          double Mtbme = 2*asNorm(ia,ib)*asNorm(ic,id)*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas)/sqrt(3); // compute the final matrix element, anti-symmetrize
          DGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      }// end of for-loop over: ibra
    } // end of for-loop over: auto
    std::cout<<"...done calculating DGT TBMEs"<<std::endl;
    return DGT_TBME;
  }

  /// Calculate B coefficient for Talmi integral. Formula given in Brody and Moshinsky
  /// "Tables of Transformation Brackets for Nuclear Shell-Model Calculations"
  long double TalmiB(int na, int la, int nb, int lb, int p)
  {
    if ( (la+lb)%2>0 ) return 0;

    int q = (la+lb)/2;

    if ( std::max(na+la+p, nb+lb+p) < 10 )
    {
      return AngMom::TalmiB( na, la, nb, lb, p);
    }
    long double logB1 = (lgamma(2*p+2)-lgamma(p+1) +0.5*(  lgamma(na+1)+lgamma(nb+1)-lgamma(na+la+1)-lgamma(nb+lb+1) + lgamma(2*na+2*la+2) + lgamma(2*nb+2*lb+2)) - (na+nb)*LOG2   ) ;
    long double B2 = 0;
    int kmin = std::max(0, p-q-nb);
    int kmax = std::min(na, p-q);
    for (int k=kmin;k<=kmax;++k)
    {
      B2  += exp(logB1+lgamma(la+k+1)-lgamma(k+1) +lgamma(p-(la-lb)/2-k+1) -lgamma(2*p-la+lb-2*k+2)
                - lgamma(2*la+2*k+2) -lgamma(na-k+1)
                - lgamma(nb - p + q + k+1) - lgamma(p-q-k+1) );
    }

    return AngMom::phase(p-q) *  B2;
  }

  double RadialIntegral_Gauss( int na, int la, int nb, int lb, double sigma )
  {
    long double I = 0;
    int pmin = (la+lb)/2;
    int pmax = pmin + na + nb;
    double kappa = 1.0 / (2*sigma*sigma); // Gaussian exp( -x^2 / (2 sigma^2)) =>  exp(-kappa x^2)

      for (int p=pmin;p<=pmax;++p)
      {
        I += TalmiB(na,la,nb,lb,p) * pow(1+kappa, -(p+1.5) )  ; // Talmi integral is analytic
      }
    return I;
  }

  double integrateRcom(int Ncom, int  Lambda, double  hw, double Rnucl)
  {
    double I = 1-RadialIntegral_Gauss(Ncom,Lambda,Ncom,Lambda,Rnucl);
    return I;
  }

  std::unordered_map<uint64_t,double> PreCalculateDGTRComIntegrals(int e2max, double hw, double Rnucl)
  {
    IMSRGProfiler profiler;
    double t_start_pci = omp_get_wtime(); // profiling (s)
    std::unordered_map<uint64_t,double> IntList;
    int size=500;
    std::cout<<"calculating integrals wrt dRcom..."<<std::endl;
    int maxn = e2max/2;
    int maxl = e2max;
    std::vector<uint64_t> KEYS;
    for (int n=0; n<=maxn; n++)
    {
      for (int l=0; l<=maxl-2*n; l++)
      {
        uint64_t key = IntHash(n,l,n,l,0,0);
        KEYS.push_back(key);
        IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
      }
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i=0; i<KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n,l,np,lp,S,J;
      IntUnHash(key, n,l,np,lp,S,J);
      IntList[key] = integrateRcom(n,l,hw,Rnucl);// these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)

    }

    std::cout<<"...done calculating the integrals"<<std::endl;
    std::cout<<"IntList has "<<IntList.bucket_count()<<" buckets and a load factor "<<IntList.load_factor()
      <<", estimated storage ~= "<<((IntList.bucket_count() + IntList.size())*(sizeof(size_t) + sizeof(void*)))/(1024.0*1024.0*1024.0)<<" GB"<<std::endl; // copied from (RS)
    profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci; // profiling (r)
    return IntList;
  }

  //Get an integral from the IntList cache or calculate it (parallelization dependent)
  double GetDGTRcomIntegral(int Ncom, int Lam, double hw, double Rnucl, std::unordered_map<uint64_t,double> &IntList)
  {
    uint64_t key = IntHash(Ncom,Lam,Ncom,Lam,0,0);
    auto it = IntList.find(key);
    if (it != IntList.end()) // return what we've found
    {
      return it -> second;
    }
    else // if we didn't find it, calculate it and add it to the list!
    {
      double integral;
      int size = 500;
      integral = integrateRcom(Ncom,Lam,hw,Rnucl);
      if (omp_get_num_threads() >= 2)
      {
        std::cout << "DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!" << std::endl;;
        std::cout << "   I shouldn't be here in " << __func__ << " ( " << Ncom << " , " << Lam << ")    key =" << key << "  integral = " << integral << std::endl;
//        printf("DANGER!!!!!!!  Updating IntList inside a parellel loop breaks thread safety!\n");
//        printf("   I shouldn't be here in GetIntegral(%d, %d):   key =%lx   integral=%f\n",Ncom,Lam,key,integral);
        exit(EXIT_FAILURE);
      }
      IntList[key] = integral;

      return integral;
    }
  }

  Operator DGT_R_SurfaceLocalization(ModelSpace& modelspace, double r12)
  {
    r12 = r12*SQRT2;
    Operator DGT_TBME(modelspace,0,2,0,2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    DGT_TBME.SetHermitian(); // it should be Hermitian
    modelspace.PreCalculateMoshinsky();
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int Anuc = modelspace.GetTargetMass();
    const double Rnuc = R0*pow(Anuc,1.0/3.0); // the nuclear radius [fm]const
    int e2max = modelspace.GetE2max(); // 2*emax
    std::unordered_map<uint64_t,double> IntList = PreCalculateDGTRComIntegrals(e2max,hw,Rnuc);
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
      // double Jhat = 1;
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
                  // if ((npr >= 0) and ((eps_cd-eps_ab)%2 == 0))
                  if (npr==nr)
                  {
                    int tempmaxNcom = tempmaxnr - nr; // just for the limits below
                    for (int Ncom = 0; Ncom <= tempmaxNcom; Ncom++)
                    {
                      int tempminlr = ceil((eps_ab - L)/2.0) - (nr + Ncom); // just for the limits below
                      int tempmaxlr = floor((eps_ab + L)/2.0) - (nr + Ncom); // " " " " "
                      for (int lr = tempminlr; lr <= tempmaxlr; lr++)
                      {
                        int Lam = eps_ab - 2*(nr + Ncom) - lr; // via Equation (4.73) of Charlie's thesis
                        double integral = 0;
                        double normJrel;
                        double Rcomintegral;
                        double Df = modelspace.GetMoshinsky(Ncom,Lam,nr,lr,na,la,nb,lb,L); // Ragnar has -- double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);
                        double Di = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nc,lc,nd,ld,L); // " " " "
                        double asDi = modelspace.GetMoshinsky(Ncom,Lam,npr,lr,nd,ld,nc,lc,L); // (anti-symmetric part)
                        int minJrel= abs(lr-S);
                        int maxJrel = lr+S;
                        for (int Jrel = minJrel; Jrel<=maxJrel; Jrel++)
                        {
                          normJrel = sqrt((2*Jrel+1)*(2*L+1))*phase(Lam+lr+S+J)*AngMom::SixJ(Lam,lr,L,S,J,Jrel);
                          Rcomintegral = GetDGTRcomIntegral(Ncom, Lam, hw, Rnuc, IntList);
                          integral += normJrel*normJrel*r12*r12*HO_Radial_psi(nr, lr, hw, r12)*HO_Radial_psi(npr, lr, hw, r12)*exp(-r12*r12*0.5)*Rcomintegral;
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

          double Mtbme = 2*asNorm(ia,ib)*asNorm(ic,id)*Jhat*(sumLS - modelspace.phase(jc + jd - J)*sumLSas)/sqrt(3); // compute the final matrix element, anti-symmetrize
          DGT_TBME.TwoBody.SetTBME(chbra,chket,ibra,iket,Mtbme); // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      }// end of for-loop over: ibra
    } // end of for-loop over: auto
    std::cout<<"...done calculating DGT TBMEs"<<std::endl;
    return DGT_TBME;
  }


  //Functions to get the operator integrands as a function of p and p'.

  double dq_pp(int n, int l, int np, int lp, int S, int J, double hw, PWD &pwd, int index_p, int index_pp)
  {
    int momentum_mesh_size = pwd.getMomentumMeshSize();
    gsl_integration_glfixed_table *t = pwd.getMomentumMesh();
    int max_momentum = pwd.getMaxMomentum();
    double I;
    double wf, wfp;
    double xi;
    double wi;
    gsl_integration_glfixed_point(0, max_momentum, index_p, &xi, &wi, t);
    wf = HO_Radial_psi_mom(n, l, hw, xi);
    double xj;
    double wj;
    gsl_integration_glfixed_point(0, max_momentum, index_pp, &xj, &wj, t);
    wfp = HO_Radial_psi_mom(np, lp, hw, xj);
    double W = pwd.getW(xi, xj, index_p, index_pp, S, l, lp, J);
    I = wi * wj * xi * xi * xj * xj * wf * wfp * W;
    if (n==0 && l==0 && np ==0 and lp ==0)
    {
      std::cout << "p = " << xi * HBARC << ", pp" << xj*HBARC<<std::endl;
    }
    return I;
  }

    std::unordered_map<uint64_t, double> PreCalculateMomentumSpaceIntegrands(int e2max, int Smin, int Smax, double hw, int Lrank, PWD &pwd, int index_p, int index_pp)
    {
      std::unordered_map<uint64_t, double> IntList;
      int maxn = e2max / 2;
      int maxl = e2max;
      int maxnp = e2max / 2;
      std::vector<uint64_t> KEYS;
      for (int S = Smin; S <= Smax; S++)
      {
        for (int n = 0; n <= maxn; n++)
        {
          for (int l = 0; l <= maxl - 2 * n; l++)
          {
            int tempminnp = n; // NOTE: need not start from 'int np=0' since IntHash(n,l,np,l) = IntHash(np,l,n,l), by construction
            for (int np = tempminnp; np <= maxnp; np++)
            {
              if (Lrank == 0)
              {
                for (int np = tempminnp; np <= maxnp; np++)
                {
                  int minJ = abs(l - S);
                  int tempmaxJ = l + S;
                  for (int J = minJ; J <= tempmaxJ; J++)
                  {
                    uint64_t key = IntHash(n, l, np, l, S, J);
                    KEYS.push_back(key);
                    IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
                  }
                }
              }
              else
              {
                int tempminlp = (n == np ? l : 0); // NOTE: need not start from 'int lp=0' since IntHash(n,l,np,lp) = IntHash(np,lp,n,l), by construction
                int maxlp = std::min(l + Lrank, maxl);
                for (int lp = tempminlp; lp <= maxlp; lp++)
                {
                  if ((abs(lp - l) != 2) and (abs(lp - l) != 0))
                    continue;
                  int minJ = std::max(abs(l - S), abs(lp - S));
                  int tempmaxJ = std::min(l + S, lp + S);
                  for (int J = minJ; J <= tempmaxJ; J++)
                  {
                    uint64_t key = IntHash(n, l, np, lp, S, J);
                    KEYS.push_back(key);
                    IntList[key] = 0.0; // "Make sure eveything's in there to avoid a rehash in the parallel loop" (RS)
                  }
                }
              }
            }
          }
        }
      }
      sort(KEYS.begin(), KEYS.end());
      KEYS.erase(unique(KEYS.begin(), KEYS.end()), KEYS.end());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < KEYS.size(); i++)
    {
      uint64_t key = KEYS[i];
      uint64_t n, l, np, lp, S, J;
      IntUnHash(key, n, l, np, lp, S, J);

      double integral = dq_pp(n, l, np, lp, S, J, hw, pwd, index_p, index_pp);
      IntList[key] = integral; // these have been ordered by the above loops such that we take the "lowest" value of decimalgen(n,l,np,lp,maxl,maxnp,maxlp), see GetIntegral(...)
    }
    std::cout << "IntList has " << IntList.bucket_count() << " buckets and a load factor " << IntList.load_factor()
              << ", estimated storage ~= " << ((IntList.bucket_count() + IntList.size()) * (sizeof(size_t) + sizeof(void *))) / (1024.0 * 1024.0 * 1024.0) << " GB" << std::endl; // copied from (RS)
    // profiler.timer["PreCalculateM0nuIntegrals"] += omp_get_wtime() - t_start_pci;                                                                                                 // profiling (r)
    return IntList;
  }

  Operator TwoBody_Scalar_integrand(ModelSpace &modelspace, PWD &pwd, int Jrank, int Lrank, int index_p, int index_pp, int Smin = 0, int Smax = 1)
  {
    double t_start, t_start_tbme;                 // profiling (v)
    t_start = omp_get_wtime();                    // profiling (s)
    double hw = modelspace.GetHbarOmega();        // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();            // 2*emax
    Operator Op_TBME(modelspace, Jrank, 2, 0, 2); // NOTE: from the constructor -- Operator::Operator(ModelSpace& ms, int Jrank, int Trank, int p, int part_rank)
    Op_TBME.SetHermitian();                       // it should be Hermitian
    modelspace.PreCalculateMoshinsky();           // pre-calculate the needed Moshinsky brackets, for efficiency
    if (pwd.getAsize() == 0)
    {
      std::cout << "Angular integrals for partial wave decomposition are not calcualted..." << std::endl;
      std::cout << "Computing angular integrals." << std::endl;
      pwd.calcA(e2max, 2);
      std::cout << "Done precomputing A integrals." << std::endl;
    }
    std::unordered_map<uint64_t, double> IntList = PreCalculateMomentumSpaceIntegrands(e2max, Smin, Smax, hw, Lrank, pwd, index_p, index_pp); // pre-calculate the needed integrals over dp and dpp
    std::cout << "Done precomputing integrals." << std::endl;
    pwd.clearA();
    pwd.freeMomentumMesh();
    Op_TBME.profiler.timer["Op_PrecalculateMomentumIntegrals"] += omp_get_wtime() - t_start; // profiling (r)
    // create the TBMEs of M0nu
    // auto loops over the TBME channels and such
    std::cout << "calculating TBMEs..." << std::endl;
    t_start_tbme = omp_get_wtime(); // profiling (s)
    for (auto &itmat : Op_TBME.TwoBody.MatEl)
    {
      int chbra = itmat.first[0];                                    // grab the channel count from auto
      int chket = itmat.first[1];                                    // " " " " " "
      TwoBodyChannel &tbc_bra = modelspace.GetTwoBodyChannel(chbra); // grab the two-body channel
      TwoBodyChannel &tbc_ket = modelspace.GetTwoBodyChannel(chket); // " " " " "
      int nbras = tbc_bra.GetNumberKets();                           // get the number of bras
      int nkets = tbc_ket.GetNumberKets();                           // get the number of kets
      int J = tbc_bra.J;                                             // NOTE: by construction, J := J_ab == J_cd := J'
      double Jhat = sqrt(2 * J + 1);                                 // We are computing reduced matrix elements
#pragma omp parallel for schedule(dynamic, 1)                        // need to do: PreCalculateMoshinsky() and PreCalcIntegrals() [above] and then "#pragma omp critical" [below]
      for (int ibra = 0; ibra < nbras; ibra++)
      {
        Ket &bra = tbc_bra.GetKet(ibra);     // get the final state = <ab|
        int ia = bra.p;                      // get the integer label a
        int ib = bra.q;                      // get the integer label b
        Orbit &oa = modelspace.GetOrbit(ia); // get the <a| state orbit
        Orbit &ob = modelspace.GetOrbit(ib); // get the <b| state prbit
        int na = oa.n;                       // this is just...
        int nb = ob.n;
        int la = oa.l;
        int lb = ob.l;
        int eps_ab = 2 * na + la + 2 * nb + lb; // for conservation of energy in the Moshinsky brackets
        double ja = oa.j2 / 2.0;
        double jb = ob.j2 / 2.0;
        for (int iket = 0; iket < nkets; iket++)
        {
          Ket &ket = tbc_ket.GetKet(iket);     // get the initial state = |cd>
          int ic = ket.p;                      // get the integer label c
          int id = ket.q;                      // get the integer label d
          Orbit &oc = modelspace.GetOrbit(ic); // get the |c> state orbit
          Orbit &od = modelspace.GetOrbit(id); // get the |d> state orbit
          int nc = oc.n;
          int nd = od.n;
          int lc = oc.l;
          int ld = od.l;
          int eps_cd = 2 * nc + lc + 2 * nd + ld; // for conservation of energy in the Moshinsky brackets
          double jc = oc.j2 / 2.0;
          double jd = od.j2 / 2.0;           // ...for convenience
          double sumLS = 0;                  // for the Bessel's Matrix Elemets (BMEs)
          double sumLSas = 0;                // (anti-symmetric part)
          for (int S = Smin; S <= Smax; S++) // sum over total spin...
          {
            int Lf_min = std::max(std::abs(la - lb), std::abs(J - S));
            int Lf_max = std::min(la + lb, J + S);
            for (int Lf = Lf_min; Lf <= Lf_max; Lf++) // sum over angular momentum coupled to l_a and l_b
            {
              double normab = sqrt((2 * Lf + 1) * (2 * S + 1) * (2 * ja + 1) * (2 * jb + 1)); // normalization factor for the 9j-symbol out front
              double nNJab = normab * AngMom::NineJ(la, lb, Lf, 0.5, 0.5, S, ja, jb, J);      // the normalized 9j-symbol out front
              int Li_min = std::max(std::abs(lc - ld), std::max(std::abs(J - S), std::abs(Lf - Lrank)));
              int Li_max = std::min(lc + ld, std::min(J + S, Lf + Lrank));
              for (int Li = Li_min; Li <= Li_max; Li++) // sum over angular momentum coupled to l_c and l_d
              {

                double sumMT = 0;                                                               // for the Moshinsky transformation
                double sumMTas = 0;                                                             // (anti-symmetric part)
                double normcd = sqrt((2 * Li + 1) * (2 * S + 1) * (2 * jc + 1) * (2 * jd + 1)); // normalization factor for the second 9j-symbol
                double nNJcd = normcd * AngMom::NineJ(lc, ld, Li, 0.5, 0.5, S, jc, jd, J);      // the second normalized 9j-symbol
                double nNJdc = normcd * AngMom::NineJ(ld, lc, Li, 0.5, 0.5, S, jd, jc, J);      // (anti-symmetric part)
                double bulk = nNJab * nNJcd;                                                    // bulk product of the above
                double bulkas = nNJab * nNJdc;                                                  // (anti-symmetric part)
                TalmiMoshinkyTransform(modelspace, sumMT, sumMTas, na, la, nb, lb, nc, lc, nd, ld,
                                       Li, Lf, S, J, Lrank, e2max, IntList);
                sumLS += bulk * sumMT;       // perform the LS-coupling sum
                sumLSas += bulkas * sumMTas; // (anti-symmetric part)
              } // end of for-loop over: Li
            } // end of for-loop over: Lf
          } // end of for-loop over: S                                                                                                                                                                                                                // end of for-loop over: S
          double Mtbme = asNorm(ia, ib) * asNorm(ic, id) * Jhat * (sumLS - phase(jc + jd - J) * sumLSas); // compute the final matrix element, anti-symmetrize
          Op_TBME.TwoBody.SetTBME(chbra, chket, ibra, iket, Mtbme);                                       // set the two-body matrix elements (TBME) to Mtbme
        } // end of for-loop over: iket
      } // end of for-loop over: ibra
    } // end of for-loop over: auto
    std::cout << "...done calculating M0nu TBMEs" << std::endl;
    Op_TBME.profiler.timer["Op_tbme"] += omp_get_wtime() - t_start_tbme; // profiling (r)
    return Op_TBME;
  }

  Operator GamowTeller_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [fm]
    const double prefact = Rnuc / (PI * PI);       // factor in-front of M0nu TBME, extra global factor of 2 since we use <p|\tau|n> = sqrt(2) [fm]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;                                       // Class for the partial wave decomposition
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(10);
    pwd.setMaxMomentum(2);
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp,  double z){return prefact*potential_closure(p,pp,z,Eclosure,formfactor);}, "spin-spin");
    pwd.calcA(e2max,0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuGT_TBME =  TwoBody_Scalar_integrand(modelspace, pwd, 0, 0, index_p, index_pp);
    return M0nuGT_TBME;
  }

  Operator Fermi_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [fm]
    const double prefact = Rnuc / (PI * PI);       // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [fm]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(10);
    pwd.setMaxMomentum(2);
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp, double z)
                     { return prefact * potential_closure(p, pp, z, Eclosure, formfactor); }, "central");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuF_TBME = TwoBody_Scalar_integrand(modelspace, pwd, 0, 0, index_p, index_pp);
    return M0nuF_TBME;
  }

  Operator Tensor_integrand(ModelSpace &modelspace, double Eclosure, std::string src, std::function<double(double)> formfactor, int index_p, int index_pp)
  {
    double hw = modelspace.GetHbarOmega();         // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();             // 2*emax
    int Anuc = modelspace.GetTargetMass();         // the mass number for the desired nucleus
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0); // the nuclear radius [MeV^-1]
    const double prefact = Rnuc / (PI * PI);       // factor in-front of M0nu TBME, extra global 2 for nutbar (as confirmed by benchmarking with Ca48 NMEs) [MeV^-1]
    modelspace.PreCalculateMoshinsky();            // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(10);
    pwd.setMaxMomentum(2);
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp, double z)
                     {
                         double qsq = p*p+pp*pp-2*p*pp*z;
                         return prefact*(3/qsq)*potential_closure(p, pp, z, Eclosure, formfactor); },
                     "tensor");
    pwd.setPotential([prefact, Eclosure, formfactor](double p, double pp, double z)
                     { return -prefact * potential_closure(p, pp, z, Eclosure, formfactor); },
                     "spin-spin");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuT_TBME = TwoBody_Scalar_integrand(modelspace, pwd, 0, 2, index_p, index_pp, 1);
    return M0nuT_TBME;
  }

  Operator Contact_integrand(ModelSpace &modelspace, double regulator_cutoff, int regulator_power, int index_p, int index_pp)
  {
    double hw = modelspace.GetHbarOmega(); // oscillator basis frequency [MeV]
    int e2max = modelspace.GetE2max();     // 2*emax
    int Anuc = modelspace.GetTargetMass(); // the mass number for the desired nucleus
    const double scale_ct = (M_NUCLEON * NUCLEON_AXIAL_G * NUCLEON_AXIAL_G * HBARC) / (4 * F_PI * F_PI);
    const double Rnuc = R0 * pow(Anuc, 1.0 / 3.0);                                    // the nuclear radius [fm]
    const double prefact = -Rnuc / (16 * SQRT2 * PI * PI * PI) * scale_ct * scale_ct; // factor in-front of M0nu TBME, includes the -2 from the -2gvv as well as the 1/2 from the convention used by Cirgilano and all when computing gvv. Only thing missing is gvv/gA^2
    std::string reg_type = "nonlocal";
    modelspace.PreCalculateMoshinsky(); // pre-calculate the needed Moshinsky brackets, for efficiency
    PWD pwd;
    pwd.initializeAngularMesh(100);
    pwd.initializeMomentumMesh(10);
    pwd.setMaxMomentum(2);
    pwd.setRegulator(regulator_cutoff, regulator_power, reg_type);
    pwd.setPotential([prefact](double p, double pp, double z)
                     {
                       return prefact; //
                     },
                     "central");
    pwd.calcA(e2max, 0);
    pwd.freeAngularMesh();
    std::cout << "Done precomputing A's." << std::endl;
    Operator M0nuC_TBME = TwoBody_Scalar_integrand(modelspace, pwd, 0, 0, index_p, index_pp, 0, 0);
    return M0nuC_TBME;
  }

  } // end namespace M0nu
