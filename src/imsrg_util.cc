
#include "imsrg_util.hh"
#include "AngMom.hh"

using namespace AngMom;

/// imsrg_util namespace. Used to define some helpful functions.
namespace imsrg_util
{

 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2)
 {
   Operator NumOp = Operator(modelspace);
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

double HO_Radial_psi(int n, int l, double hw, double r)
{
   double b = sqrt( (HBARC*HBARC) / (hw * M_NUCLEON) );
   double x = r/b;
   double Norm = 2*sqrt( gsl_sf_fact(n) * pow(2,n+l) / M_SQRTPI / gsl_sf_doublefact(2*n+2*l+1) / pow(b,3.0) );
   double L = gsl_sf_laguerre_n(n,l+0.5,x*x);
   double psi = Norm * pow(x,l) * exp(-x*x*0.5) * L;
   return psi;
}

 // Just do the HF transformation
 vector<double> GetOccupationsHF(HartreeFock& hf)
 {
    ModelSpace* modelspace = hf.Hbare.modelspace;
    int norb = modelspace->GetNumberOrbits();
    vector<double> occupation(norb);

    for (int i=0; i<norb; ++i)
    {
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
      // Transform it to the normal-ordered HF basis
      Operator N_NO = hf.TransformToHFBasis(N_bare).DoNormalOrdering();
      occupation[i] = N_NO.ZeroBody;
      cout << oi.n << " " << oi.l << " " << oi.j2 << "/2 " << occupation[i] << endl;
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
      Orbit & oi = modelspace->GetOrbit(i);
      // Get the number operator for orbit i
      Operator N_bare = NumberOp(*modelspace,oi.n,oi.l,oi.j2,oi.tz2);
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
       Orbit & oi = modelspace.GetOrbit(i);
       for (int ir=0;ir<nr_steps;++ir)
       {
          dens[ir] += HO_density(oi.n,oi.l,hw,R[ir]) * occupation[i];
       }
     }
     return dens;
 }

 Operator Single_Ref_1B_Density_Matrix(ModelSpace& modelspace)
 {
    Operator DM(modelspace,0,0,0,2);
    for (index_t a : modelspace.holes)
    {
       DM.OneBody(a,a) = 1.0;
    }
    return DM;
 }

 double Get_Charge_Density( Operator& DM, double r)
 {
   ModelSpace* modelspace = DM.GetModelSpace();
   double hw = modelspace->GetHbarOmega();
   double rho=0;
   for (index_t i : modelspace->proton_orbits )
   {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( index_t j : DM.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        if (abs(DM.OneBody(i,j))<1e-7) continue;
        Orbit& oj = modelspace->GetOrbit(j);
        rho += DM.OneBody(i,j) * HO_Radial_psi(oi.n, oi.l, hw, r) * HO_Radial_psi(oj.n, oj.l, hw, r);
      }
   }
   return rho;
 }

/// Relative kinetic energy, includingthe hw/A factor
/// \f[
/// T_{rel} = T - T_{CM}
/// \f]
 Operator Trel_Op(ModelSpace& modelspace)
 {
   Operator Trel(modelspace);
   int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      Trel.OneBody(a,a) = 0.5 * hw * (2*oa.n + oa.l +3./2); 
      for ( int b : Trel.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
         if (b<=a) continue;
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.n == ob.n+1)
            Trel.OneBody(a,b) = 0.5 * hw * sqrt( (oa.n)*(oa.n + oa.l +1./2));
         else if (oa.n == ob.n-1)
            Trel.OneBody(a,b) = 0.5 * hw * sqrt( (ob.n)*(ob.n + ob.l +1./2));
         Trel.OneBody(b,a) = Trel.OneBody(a,b);
      }
   }
   Trel -= TCM_Op(modelspace);
   return Trel;
 }


/// Center of mass kinetic energy, including the hw/A factor
/// \f[
/// T = \frac{\hbar\omega}{A}\sum_{ij} t_{ij} a^{\dagger}_{i} a_{j} + \frac{\hbar\omega}{A}\frac{1}{4}\sum_{ijkl} t_{ijkl} a^{\dagger}_{i}a^{\dagger}_{j}a_{l}a_{k}
/// \f]
/// with a one-body piece
/// \f[
/// t_{ij} = \frac{1}{\hbar\omega} \left\langle i | T_{12} | j \right\rangle = \frac{1}{2}(2n_i+\ell_i+3/2) \delta_{ij} + \frac{1}{2}\sqrt{n_j(n_j+\ell_j+\frac{1}{2})} \delta_{n_i,n_j-1}\delta_{k_i k_j}
/// \f]
/// where \f$k\f$ labels all quantum numbers other than \f$n\f$ and a two-body piece
/// \f[  
/// t_{ijkl} = \frac{1}{\hbar\omega} \left\langle ij | (T^{CM}_{12} - T^{rel}_{12}) | kl \right\rangle
/// \f]
 Operator TCM_Op(ModelSpace& modelspace)
 {
   double t_start = omp_get_wtime();
   int N2max = modelspace.GetN2max();
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   Operator TcmOp = Operator(modelspace);
   TcmOp.SetHermitian();
   // One body piece = p**2/(2mA)
   int norb = modelspace.GetNumberOrbits();
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (int j : modelspace.OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         Orbit & oj = modelspace.GetOrbit(j);
         if (j<i) continue;
         double tij = 0;
         if (oi.n == oj.n) tij = 0.5*(2*oi.n+oi.l + 1.5) * hw/A;
         else if (oi.n == oj.n-1) tij = 0.5*sqrt(oj.n*(oj.n+oj.l + 0.5)) * hw/A;
         TcmOp.OneBody(i,j) = tij;
         TcmOp.OneBody(j,i) = tij;
      }
   }

   // Two body piece = 2*p1*p2/(2mA) = (Tcm-Trel)/A
   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1)  // In order to make this parallel, need to precompute the Moshinsky brackets.
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         Orbit & oi = modelspace.GetOrbit(bra.p);
         Orbit & oj = modelspace.GetOrbit(bra.q);
         if ( 2*(oi.n+oj.n)+oi.l+oj.l > N2max) continue;
         for (int iket=ibra;iket<nkets;++iket)
         {
            
            Ket & ket = tbc.GetKet(iket);
            Orbit & ok = modelspace.GetOrbit(ket.p);
            Orbit & ol = modelspace.GetOrbit(ket.q);
            if ( 2*(ok.n+ol.n)+ok.l+ol.l > N2max) continue;
            double p1p2 = Calculate_p1p2(modelspace,bra,ket,tbc.J) * hw/A;
            if (abs(p1p2)>1e-7)
            {
              TcmOp.TwoBody.SetTBME(ch,ibra,iket,p1p2);
            }
         }
      }
   }
   TcmOp.profiler.timer["TCM_Op"] += omp_get_wtime() - t_start;
   return TcmOp;
 }


 // evaluate <bra| p1*p2 | ket> , omitting the prefactor  m * hbar_omega
/// This returns the antisymmetrized J-coupled two body matrix element of \f$ \vec{p}_1 \cdot \vec{p}_2 / (m\hbar\omega) \f$.
/// The formula is
/// \f{eqnarray*}{
/// \frac{1}{m\hbar\omega}
/// \left \langle a b \right| \vec{p}_1 \cdot \vec{p}_2 \left| c d \right \rangle_J
/// = \frac{1}{\sqrt{(1+\delta_{ab})(1+\delta_{cd})}} &\sum\limits_{LS}
/// \left[ \begin{array}{ccc}
///  \ell_a & s_a & j_a \\
///  \ell_b & s_b & j_b \\
///  L      & S   & J
/// \end{array} \right]
/// \left[ \begin{array}{ccc}
///  \ell_c & s_c & j_c \\
///  \ell_d & s_d & j_d \\
///  L      & S   & J
/// \end{array} \right] & \\
/// & \times \sum\limits_{\substack{N_{ab}N_{cd} \Lambda \\ n_{ab} n_{cd} \lambda}} \mathcal{A}_{abcd}^{\lambda S} \times
/// \left\langle N_{ab}\Lambda n_{ab} \lambda | n_{a} \ell_{a} n_{b} \ell_{b} \right\rangle_{L} &
/// \left\langle N_{cd}\Lambda n_{cd} \lambda | n_{c} \ell_{c} n_{d} \ell_{d} \right\rangle_{L} \\ & \times
/// \left( \left\langle N_{ab}\Lambda | t_{cm} | N_{cd} \Lambda \right \rangle
/// -\left\langle n_{ab}\lambda | t_{rel} | n_{cd} \lambda \right \rangle \right)
/// \f}
/// The antisymmetrization factor \f$ \mathcal{A}_{abcd}^{\lambda S} \f$ ensures that
/// the relative wave function is antisymmetrized.
/// It is given by \f$ \mathcal{A}_{abcd}^{\lambda S} = \left|t_{za}+t_{zc}\right| + \left| t_{za} + t_{zd} \right| (-1)^{\lambda + S + \left|T_z\right|}\f$ .
///
/// The center-of-mass and relative kinetic energies can be found by the same equation as used in the one-body piece of TCM_Op()
 double Calculate_p1p2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
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
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   // p1*p2 only connects kets with delta N = 0,1 ==> delta E = 0,2
   if (abs(fab-fcd)>2 or abs(fab-fcd)%2 >0 ) return 0; 

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double p1p2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

                double tcm = 0;
                double trel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   tcm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) tcm = sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) tcm = sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   trel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) trel = sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) trel = sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                p1p2 += (tcm-trel) * prefactor ;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize. The 0.5 comes from t ~ 0.5 * (N+3/2) hw
   p1p2 *= 0.5 / sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return p1p2 ;

 }



////////////////////////////////////////////////////////////////////////////
/////////////  IN PROGRESS...  doesn't work yet, and for now it's slower.  /
////////////////////////////////////////////////////////////////////////////

 void Calculate_p1p2_all(Operator& OpIn)
 {
   ModelSpace* modelspace = OpIn.GetModelSpace();
//   modelspace->PreCalculateMoshinsky();
   for ( int ch : modelspace->SortedTwoBodyChannels )
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int parity = tbc.parity;
//      int Tz = tbc.Tz;
      arma::mat& MatJJ = OpIn.TwoBody.GetMatrix(ch);
      int nkets_JJ = tbc.GetNumberKets();

      // Find the maximum oscillator energy for this channel
      Ket& ketlast = tbc.GetKet( tbc.GetNumberKets()-1 );
      int emax_ket = 2*ketlast.op->n + 2*ketlast.oq->n + ketlast.op->l + ketlast.oq->l;

      vector<array<int,6>> JacobiBasis;  // L,S,N,Lambda,n,lambda
      for (int L=max(J-1,0); L<=J+1; ++L)
      {
       for ( int S=abs(J-L); S<=1; ++S)
       {
        for ( int N=0; N<=emax_ket/2; ++N )
        {
         for ( int Lambda=0; Lambda<=(emax_ket-2*N); ++Lambda)
         {
          for ( int lambda=abs(L-Lambda)+(L+parity)%2; lambda<=min(Lambda+L,emax_ket-2*N-Lambda); lambda+=2)
          {
           for ( int n =0; n<=(emax_ket-2*N-Lambda-lambda)/2; ++n)
           {
             JacobiBasis.push_back({L,S,N,Lambda,n,lambda});
           }
          }
         }
        }
       }
      }


      int nkets_Jacobi = JacobiBasis.size();
      arma::mat MatJacobi(nkets_Jacobi,nkets_Jacobi);
      arma::mat Trans(nkets_Jacobi,nkets_JJ);
      int n_nonzero = 0;
      for (int iJJ=0; iJJ<nkets_JJ; ++iJJ)
      {
        Ket & ket = tbc.GetKet(iJJ);
        int la = ket.op->l;
        int lb = ket.oq->l;
        int na = ket.op->n;
        int nb = ket.oq->n;
        double ja = ket.op->j2*0.5;
        double jb = ket.oq->j2*0.5;
//        int ta = ket.op->n;
//        int tb = ket.oq->n;
        for (int iJac=0; iJac<nkets_Jacobi; ++iJac)
        {
          int L      = JacobiBasis[iJac][0];
          int S      = JacobiBasis[iJac][1];
          int N      = JacobiBasis[iJac][2];
          int Lambda = JacobiBasis[iJac][3];
          int n      = JacobiBasis[iJac][4];
          int lambda = JacobiBasis[iJac][5];

//          int Asym = 1; // Fix this...
          double ninej = NormNineJ(la,0.5,ja,lb,0.5,jb,L,S,J);
          if (abs(ninej)<1e-6) continue;
          double mosh = modelspace->GetMoshinsky(N,Lambda,n,lambda,na,la,nb,lb,L);
          if (abs(mosh)<1e-6) continue;
          Trans(iJac,iJJ) = ninej * mosh;
          n_nonzero += 1;
           
        }
      }
      for (int i=0; i<nkets_Jacobi; ++i)
      {
        int Li      = JacobiBasis[i][0];
        int Si      = JacobiBasis[i][1];
        int Ni      = JacobiBasis[i][2];
        int Lambdai = JacobiBasis[i][3];
        int ni      = JacobiBasis[i][4];
        int lambdai = JacobiBasis[i][5];
        for (int j=0; j<nkets_Jacobi; ++j)
        {
          int Lj      = JacobiBasis[j][0];
          int Sj      = JacobiBasis[j][1];
          int Nj      = JacobiBasis[j][2];
          int Lambdaj = JacobiBasis[j][3];
          int nj      = JacobiBasis[j][4];
          int lambdaj = JacobiBasis[j][5];
          if(Li!=Lj or Si!=Sj or Lambdai!=Lambdaj or lambdai!=lambdaj) continue;
          double tcm = 0;
          double trel = 0;
          if (ni == nj)
          {
            if      (Ni == Nj)   tcm = (2*Ni+Lambdai+1.5);
            else if (Ni == Nj+1) tcm = sqrt(Ni*( Ni+Lambdai+0.5));
            else if (Ni == Nj-1) tcm = sqrt(Nj*( Nj+Lambdai+0.5));
          }
          if (Ni == Nj)
          {
            if      (ni == nj)   tcm = (2*ni+lambdai+1.5);
            else if (ni == nj+1) tcm = sqrt(ni*( ni+lambdai+0.5));
            else if (ni == nj-1) tcm = sqrt(nj*( nj+lambdai+0.5));
          }
          MatJacobi(i,j) = tcm - trel;
        }
      }


      MatJJ = Trans.t() * MatJacobi * Trans;
      cout << "ch = " << ch << "   size of JJ basis = " << nkets_JJ << "  size of Jacobi Basis = " << nkets_Jacobi << "   nonzero matrix elements = " << n_nonzero << endl;

   }
 }

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////





// Center of mass R^2, with the hw/A factor
/// Returns
/// \f[ 
/// R^{2}_{CM} = \left( \frac{1}{A}\sum_{i}\vec{r}_{i}\right)^2 =
/// \frac{1}{A^2} \left( \sum_{i}r_{i}^{2} + \sum_{i\neq j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.
 Operator R2CM_Op(ModelSpace& modelspace)
 {
   Operator R2cmOp = Operator(modelspace);

   unsigned int norb = modelspace.GetNumberOrbits();
   for (unsigned int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (auto j : modelspace.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         if (j<i) continue;
         Orbit & oj = modelspace.GetOrbit(j);
         double rij = 0;
         if (oi.n == oj.n)        rij = (2*oi.n+oi.l + 1.5);
         else if (oi.n == oj.n-1) rij = -sqrt(oj.n*(oj.n+oj.l + 0.5));
         R2cmOp.OneBody(i,j) = rij;
         R2cmOp.OneBody(j,i) = rij;
      }
   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J); 
             
            R2cmOp.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            R2cmOp.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   return R2cmOp * (HBARC*HBARC/M_NUCLEON/hw)/(A*A);
 }



// Center of mass R^2, with the hw/A factor
/// Returns
/// \f[ 
/// R^{2}_{CM} = \left( \frac{1}{A}\sum_{i}\vec{r}_{i}\right)^2 =
/// \frac{1}{A^2} \left( \sum_{i}r_{i}^{2} + 2\sum_{i<j}\vec{r}_i\cdot\vec{r}_j  \right)
/// \f]
/// evaluated in the oscillator basis.
 Operator Rp2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   return R2CM_Op(modelspace) + (A-2.0)/(A*Z)*R2_1body_Op(modelspace,"proton") - 2./(A*Z)*R2_2body_Op(modelspace,"proton");
 }

 Operator Rn2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   return R2CM_Op(modelspace) + (A-2.0)/(A*N)*R2_1body_Op(modelspace,"neutron") - 2./(A*N)*R2_2body_Op(modelspace,"neutron");
 }

 Operator Rm2_corrected_Op(ModelSpace& modelspace, int A, int Z)
 {
   return R2CM_Op(modelspace) + (A-2.0)/(A*Z)*R2_p1_Op(modelspace) - 2./(A*Z)*R2_p2_Op(modelspace);
 }



 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
/// Returns the normalized, anti-symmetrized, J-coupled, two-body matrix element of \f$ \frac{m\omega^2}{\hbar \omega} \vec{r}_1\cdot\vec{r}_2 \f$.
/// Calculational details are similar to Calculate_p1p2().
 double Calculate_r1r2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
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
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double r1r2=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              if (Lab<abs(Lam_ab-lam_ab) or Lab>(Lam_ab+lam_ab) ) continue;
              // factor to account for antisymmetrization

              int asymm_factor = (abs(bra.op->tz2+ket.op->tz2) + abs(bra.op->tz2+ket.oq->tz2)*modelspace.phase( lam_ab + Sab ))/ 2;
              if ( asymm_factor ==0 ) continue;

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
              {
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd and N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

                double r2cm = 0;
                double r2rel = 0;
                if (n_ab == n_cd)
                {
                  if      (N_ab == N_cd)   r2cm = (2*N_ab+Lam_ab+1.5);
                  else if (N_ab == N_cd+1) r2cm = -sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) r2cm = -sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   r2rel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) r2rel = -sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) r2rel = -sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                r1r2 += (r2cm-r2rel) * prefactor;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize.
   r1r2 *=  sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1r2 ;

 }




/// Center of mass Hamiltonian
/// \f[
/// \begin{align}
/// H_{CM} &= T_{CM} + \frac{1}{2} Am\omega^2 R^2 \\
///        &= T_{CM} + \frac{1}{2b^2} AR^2 \hbar\omega
/// \end{align}
/// \f]
 Operator HCM_Op(ModelSpace& modelspace)
 {
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   double oscillator_b2 = HBARC*HBARC/M_NUCLEON/hw;
   Operator HcmOp = TCM_Op(modelspace) + R2CM_Op(modelspace) * (0.5*A * hw / oscillator_b2);
   return HcmOp;
 }




/// Returns
/// \f[ r^2 = \sum_{i} r_{i}^2 \f]
///
Operator RSquaredOp(ModelSpace& modelspace)
{
   Operator r2 = Operator(modelspace);
   r2.OneBody.zeros();
   unsigned int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (unsigned int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
   r2.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   return r2;
}


 Operator R2_p1_Op(ModelSpace& modelspace)
 {
   return R2_1body_Op(modelspace,"proton");
 }

/// One-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p1} = \sum_{i} e_{i}{r}_i^2
/// \f]
 Operator R2_1body_Op(ModelSpace& modelspace,string option)
 {
   Operator r2(modelspace);
   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   auto orbitlist = modelspace.proton_orbits;
   if (option == "neutron") orbitlist = modelspace.neutron_orbits;
   else if (option == "matter")  orbitlist.insert(orbitlist.end(),modelspace.neutron_orbits.begin(),modelspace.neutron_orbits.end());
   else if (option != "proton") cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p1_Op !!!" << endl;
 
   for (unsigned int a : orbitlist )
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for ( unsigned int b : r2.OneBodyChannels.at({oa.l, oa.j2, oa.tz2}) )
      {
        if ( b < a ) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              r2.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              r2.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           r2.OneBody(b,a) = r2.OneBody(a,b);
        }
      }
   }
   r2.OneBody *= oscillator_b;
   return r2;
 }

 Operator R2_p2_Op(ModelSpace& modelspace)
 {
   return R2_2body_Op(modelspace,"proton");
 }

/// Two-body part of the proton charge radius operator.
/// Returns
/// \f[ 
/// \hat{R}^{2}_{p2} = \sum_{i\neq j} e_{i}\vec{r}_i\cdot\vec{r}_j 
/// \f]
/// evaluated in the oscillator basis.
 Operator R2_2body_Op(ModelSpace& modelspace,string option)
 {
   Operator Rp2Op(modelspace,0,0,0,2);
   double oscillator_b = (HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega());

   int nchan = modelspace.GetNumberTwoBodyChannels();
//   modelspace.PreCalculateMoshinsky();
   #pragma omp parallel for schedule(dynamic,1) 
   for (int ch=0; ch<nchan; ++ch)
   {
      TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
      int nkets = tbc.GetNumberKets();
      for (int ibra=0;ibra<nkets;++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         if (option == "proton" and bra.op->tz2>0) continue;
         else if (option == "neutron" and bra.op->tz2<0) continue;
         else if (option != "matter") cout << "!!! WARNING. BAD OPTION "  << option << " FOR imsrg_util::R2_p2_Op !!!" << endl;
         for (int iket=ibra;iket<nkets;++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            double mat_el = Calculate_r1r2(modelspace,bra,ket,tbc.J) * oscillator_b ; 
            Rp2Op.TwoBody.SetTBME(ch,ibra,iket,mat_el);
            Rp2Op.TwoBody.SetTBME(ch,iket,ibra,mat_el);
         }
      }
   }
   return Rp2Op;
 }




// Electric monopole operator
/// Returns
/// \f[ r_{e}^2 = \sum_{i} e_{i} r_{i}^2 \f]
///
Operator E0Op(ModelSpace& modelspace)
{
   Operator e0(modelspace);
   e0.EraseZeroBody();
   e0.OneBody.zeros();
//   unsigned int norbits = modelspace.GetNumberOrbits();
   double hw = modelspace.GetHbarOmega();
   for (unsigned int a : modelspace.proton_orbits)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      e0.OneBody(a,a) = (2*oa.n + oa.l +1.5); 
      for (unsigned int b : e0.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
      {
        if (b<=a) continue;
        Orbit & ob = modelspace.GetOrbit(b);
        {
           if (oa.n == ob.n+1)
              e0.OneBody(a,b) = -sqrt( (oa.n)*(oa.n + oa.l +0.5));
           else if (oa.n == ob.n-1)
              e0.OneBody(a,b) = -sqrt( (ob.n)*(ob.n + ob.l +0.5));
           e0.OneBody(b,a) = e0.OneBody(a,b);
        }
      }
   }
   e0.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   return e0;
}





 // Evaluate <bra | hcom | ket>, omitting the factor (hbar * omega) /(m * omega^2)
 double Calculate_hcom(ModelSpace& modelspace, Ket & bra, Ket & ket, int J)
 {
   Orbit & oa = modelspace.GetOrbit(bra.p);
   Orbit & ob = modelspace.GetOrbit(bra.q);
   Orbit & oc = modelspace.GetOrbit(ket.p);
   Orbit & od = modelspace.GetOrbit(ket.q);

   int na = oa.n;
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
   double jd = od.j2/2.0;

   int fab = 2*na + 2*nb + la + lb;
   int fcd = 2*nc + 2*nd + lc + ld;
   if (abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

   double sa,sb,sc,sd;
   sa=sb=sc=sd=0.5;

   double hcom=0;

   // First, transform to LS coupling using 9j coefficients
   for (int Lab=abs(la-lb); Lab<= la+lb; ++Lab)
   {
     for (int Sab=0; Sab<=1; ++Sab)
     {
       if ( abs(Lab-Sab)>J or Lab+Sab<J) continue;

       double njab = NormNineJ(la,sa,ja, lb,sb,jb, Lab,Sab,J);
       if (njab == 0) continue;
       int Scd = Sab;
       int Lcd = Lab;
       double njcd = NormNineJ(lc,sc,jc, ld,sd,jd, Lcd,Scd,J);
       if (njcd == 0) continue;

       // Next, transform to rel / com coordinates with Moshinsky tranformation
       for (int N_ab=0; N_ab<=fab/2; ++N_ab)  // N_ab = CoM n for a,b
       {
         for (int Lam_ab=0; Lam_ab<= fab-2*N_ab; ++Lam_ab) // Lam_ab = CoM l for a,b
         {
           int Lam_cd = Lam_ab; // tcm and trel conserve lam and Lam, ie relative and com orbital angular momentum
           for (int lam_ab=(fab-2*N_ab-Lam_ab)%2; lam_ab<= (fab-2*N_ab-Lam_ab); lam_ab+=2) // lam_ab = relative l for a,b
           {
              // factor to account for antisymmetrization
              int asymm_factor = 1;

//              if (ket.Tz!=0)
              if (ket.op->tz2 == ket.oq->tz2)
              {
                if ((lam_ab+Sab)%2>0) continue; // Pauli rule for identical particles
                asymm_factor = 2 ;
              }
              else if ( (bra.p + ket.p)%2 >0) // if we have pnnp or nppn, then pick up a phase
              {
                asymm_factor = modelspace.phase( lam_ab + Sab );
              }

              int lam_cd = lam_ab; // tcm and trel conserve lam and Lam
              int n_ab = (fab - 2*N_ab-Lam_ab-lam_ab)/2; // n_ab is determined by energy conservation

              double mosh_ab = modelspace.GetMoshinsky(N_ab,Lam_ab,n_ab,lam_ab,na,la,nb,lb,Lab);

              if (abs(mosh_ab)<1e-8) continue;

//              for (int N_cd=max(0,N_ab-1); N_cd<=N_ab+1; ++N_cd) // N_cd = CoM n for c,d
//              {
                int N_cd = N_ab;
                int n_cd = (fcd - 2*N_cd-Lam_cd-lam_cd)/2; // n_cd is determined by energy conservation
                if (n_cd < 0) continue;
                if  (n_ab != n_cd or N_ab != N_cd) continue;

                double mosh_cd = modelspace.GetMoshinsky(N_cd,Lam_cd,n_cd,lam_cd,nc,lc,nd,ld,Lcd);
                if (abs(mosh_cd)<1e-8) continue;

                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                hcom += (2*N_ab + Lam_ab - 2*n_ab - lam_ab) * prefactor;

//              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize
   hcom /= sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return hcom ;

 }


/// Returns the \f$ T^{2} \f$ operator
 Operator Isospin2_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.75);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     arma::mat& TB = T2.TwoBody.GetMatrix(ch);
     // pp,nn:  2<t2.t1> = 1/(2(1+delta_ab)) along diagonal
     if (abs(tbc.Tz) == 1)
     {
        TB.diag().fill(0.5); // pp,nn TBME's
        for (int ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           if (bra.p == bra.q)
           {
             TB(ibra,ibra) /= 2.;
           }
        }
     }
     else if (tbc.Tz == 0)
     {
        for (int ibra=0;ibra<tbc.GetNumberKets(); ++ibra)
        {
           Ket& bra = tbc.GetKet(ibra);
           Orbit& oa = modelspace.GetOrbit(bra.p);
           Orbit& ob = modelspace.GetOrbit(bra.q);
           for (int iket=ibra;iket<tbc.GetNumberKets(); ++iket)
           {
             Ket& ket = tbc.GetKet(iket);
             Orbit& oc = modelspace.GetOrbit(ket.p);
             Orbit& od = modelspace.GetOrbit(ket.q);
             if (oa.j2==oc.j2 and oa.n==oc.n and oa.l==oc.l
               and ob.j2==od.j2 and ob.n==od.n and ob.l==od.l )
             {
               // tz1 tz2 case
               if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
               {
                 TB(ibra,iket) -= 0.5;
               }
               // t+ t- case
               if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
               {
                 TB(ibra,iket) += 1.0;
               }
               // if a==b==c==d, we need to consider the exchange term
               if (oa.j2==ob.j2 and oa.n==ob.n and oa.l==ob.l)
               {
                  int phase = bra.Phase(tbc.J);
                  // tz1 tz2 case
                  if( oa.tz2 == oc.tz2 and ob.tz2==od.tz2)
                  {
                    TB(ibra,iket) += phase * 1.0;
                  }
                  // t+ t- case
                  if( oa.tz2 == od.tz2 and ob.tz2==oc.tz2)
                  {
                    TB(ibra,iket) -= phase * 0.5;
                  }
               }
               TB(iket,ibra) = TB(ibra,iket); // hermitian
             }

           }
        }
     }
   }
   return T2;
 }



  /// Returns a reduced electric multipole operator with units \f$ e\f$ fm\f$^{\lambda} \f$
  Operator ElectricMultipoleOp(ModelSpace& modelspace, int L)
  {
    Operator EL(modelspace, L,0,L%2,2);
    for (int i : modelspace.proton_orbits)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for ( int j : EL.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        // multiply radial integra by b^L = (hbar/mw)^L/2
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L) * pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*L) ;
        EL.OneBody(i,j) = modelspace.phase((oi.j2+1)/2) * sqrt( (oi.j2+1)*(oj.j2+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(oi.j2/2.0, L, oj.j2/2.0, 0.5,0, -0.5) * r2int;
        EL.OneBody(j,i) = modelspace.phase(oi.j2-oj.j2) * EL.OneBody(i,j);
      }
    }
    return EL;
  }

  /// Returns a reduced magnetic multipole operator with units \f$ \mu_{N}\f$ fm\f$ ^{\lambda-1} \f$
  Operator MagneticMultipoleOp(ModelSpace& modelspace, int L)
  {
    double gl[2] = {1.,0.};
    double gs[2] = {5.586, -3.826};
    Operator ML(modelspace, L,0,(L+1)%2,2);
    if (L<1)
    {
      cout << "A magnetic monopole operator??? Setting it to zero." << endl;
      return ML;
    }
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for ( int j : ML.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
        if (j<i) continue;
        Orbit& oj = modelspace.GetOrbit(j);
        // multiply radial integral by b^L = (hbar/mw)^L/2
        double r2int = RadialIntegral(oi.n,oi.l,oj.n,oj.l,L-1) * pow( HBARC*HBARC/M_NUCLEON/modelspace.GetHbarOmega(),0.5*(L-1));
        int kappa = ( modelspace.phase(oi.l+(oi.j2+1)/2) * (oi.j2+1) + modelspace.phase(oj.l+(oj.j2+1)/2) * (oj.j2+1) )/2;
        ML.OneBody(i,j) = modelspace.phase((oi.j2+1)/2) * sqrt( (oi.j2+1)*(oj.j2+1)*(2*L+1)/4./3.1415926) * AngMom::ThreeJ(oi.j2/2.0, L, oj.j2/2.0, 0.5,0, -0.5)
                        * (L - kappa) *(gl[(oi.tz2+1)/2]*(1+kappa/(L+1))-0.5*gs[(oi.tz2+1)/2] )  * r2int;
        ML.OneBody(j,i) = modelspace.phase(oi.j2-oj.j2) * ML.OneBody(i,j);
      }
    }
    return ML;
  }


/// Evaluate the radial integral \f[
/// \tilde{\mathcal{R}}^{\lambda}_{ab} = \int_{0}^{\infty} dx \tilde{g}_{n_a\ell_a}(x)x^{\lambda+2}\tilde{g}_{n_b\ell_b}(x)
/// \f]
/// where \f$ \tilde{g}(x) \f$ is the radial part of the harmonic oscillator wave function with unit oscillator length \f$ b=1 \f$
/// and \f$ x = r/b \f$.
/// To obtain the radial integral for some other oscillator length, multiply by \f$ b^{\lambda} \f$.
/// This implementation uses eq (6.41) from Suhonen.
/// Note this is only valid for \f$ \ell_a+\ell_b+\lambda\f$ = even.
  double RadialIntegral(int na, int la, int nb, int lb, int L)
  {
    int tau_a = max((lb-la+L)/2,0);
    int tau_b = max((la-lb+L)/2,0);
    int sigma_min = max(max(na-tau_a,nb-tau_b),0);
    int sigma_max = min(na,nb);
  
    double term1 = AngMom::phase(na+nb) * gsl_sf_fact(tau_a)*gsl_sf_fact(tau_b) * sqrt(gsl_sf_fact(na)*gsl_sf_fact(nb) / (gsl_sf_gamma(na+la+1.5)*gsl_sf_gamma(nb+lb+1.5) ) );
    double term2 = 0;
    for (int sigma=sigma_min; sigma<=sigma_max; ++sigma)
    {
      term2 += gsl_sf_gamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
    }
    return term1*term2;
  
  }





  Operator AllowedFermi_Op(ModelSpace& modelspace)
  {
    Operator Fermi(modelspace,0,1,0,2);
    Fermi.SetHermitian();
    const double M_fermi = 1.0;
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : Fermi.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.tz2 == oj.tz2) continue;
        Fermi.OneBody(i,j) = M_fermi;
      }
    }
    return Fermi;
  }

  Operator AllowedGamowTeller_Op(ModelSpace& modelspace)
  {
    Operator GT(modelspace,1,1,0,2);
    GT.SetHermitian();
    int norbits = modelspace.GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace.GetOrbit(i);
      for (int j : GT.OneBodyChannels[{oi.l,oi.j2,oi.tz2}] )
      {
        Orbit& oj = modelspace.GetOrbit(j);
        if (oi.n!=oj.n or oi.l != oj.l or oi.tz2==oj.tz2) continue;
        double sixj = modelspace.GetSixJ(0.5,0.5,1.0,oj.j2/2.,oi.j2/2.,oi.l);
        double M_gt = 2 * modelspace.phase(oi.l+oi.j2/2.0+1.5) * sqrt((oi.j2+1)*(oj.j2+1)) * sqrt(1.5) * sixj;
        GT.OneBody(i,j) = M_gt;
      }
    }
    return GT;
  }




  void Reduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    int norbits = modelspace->GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( int j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) *= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      int ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second *= sqrt(2*J+1.);
    }
  }

  void UnReduce(Operator& X)
  {
    ModelSpace* modelspace = X.GetModelSpace();
    int norbits = modelspace->GetNumberOrbits();
    for (int i=0; i<norbits; ++i)
    {
      Orbit& oi = modelspace->GetOrbit(i);
      for ( int j : X.OneBodyChannels.at({oi.l, oi.j2, oi.tz2}) )
      {
         X.OneBody(i,j) /= sqrt(oi.j2+1.);
      }
    }

    for ( auto& itmat : X.TwoBody.MatEl )
    {
      int ch_bra = itmat.first[0];
      int J = modelspace->GetTwoBodyChannel(ch_bra).J;
      itmat.second /= sqrt(2*J+1.);
    }

  }





  void CommutatorTest(Operator& X, Operator& Y)
  {
    Operator Zscalar(X);
    if ( (X.IsHermitian() and Y.IsHermitian()) or (X.IsAntiHermitian() and Y.IsAntiHermitian()) ) Zscalar.SetAntiHermitian();
    if ( (X.IsHermitian() and Y.IsAntiHermitian()) or (X.IsAntiHermitian() and Y.IsHermitian()) ) Zscalar.SetHermitian();
    Zscalar.Erase();
    Operator Ztensor(Zscalar);
    Operator Yred = Y;
    Reduce(Yred);

//    X.comm111ss(Y,Zscalar);
//    X.comm111st(Yred,Ztensor);
    Zscalar.comm111ss(X,Y);
    Ztensor.comm111st(X,Yred);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm111 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm111 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm121ss(Y,Zscalar);
    X.comm121st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm121 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm121 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm122ss(Y,Zscalar);
    X.comm122st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm122 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm122 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_pp_hh_221ss(Y,Zscalar);
    X.comm222_pp_hh_221st(Yred,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_pp_hh_221 norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_pp_hh_221 diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;

    Zscalar.Erase();
    Ztensor.Erase();
    X.comm222_phss(Y,Zscalar);
    Reduce(Y); // Not sure why I can't use Yred...
    X.comm222_phst(Y,Ztensor);
    Zscalar.Symmetrize();
    Ztensor.Symmetrize();
    UnReduce(Ztensor);
    cout << "comm222_ph norm = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << ",   " << Ztensor.OneBodyNorm() << " " << Ztensor.TwoBodyNorm() << endl;
    Zscalar -= Ztensor;
    cout << "comm222_ph diff = " << Zscalar.OneBodyNorm() << " " << Zscalar.TwoBodyNorm() << endl;



  }

}// namespace imsrg_util

