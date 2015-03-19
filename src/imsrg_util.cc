/** @file */


#include "imsrg_util.hh"
#include "AngMom.hh"

using namespace AngMom;

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


// Operator TCM_Op(ModelSpace& modelspace)
// {
//   return TCM_Op(modelspace, modelspace.GetN2max());
// }
// Center of mass kinetic energy, including the hw/A factor
// Operator TCM_Op(ModelSpace& modelspace, int N2max)
 Operator TCM_Op(ModelSpace& modelspace)
 {
   int N2max = modelspace.GetN2max();
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   Operator TcmOp = Operator(modelspace);

   // One body piece = p**2/(2mA) = (N+3/2)hw/A
   int norb = modelspace.GetNumberOrbits();
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (int j=i; j<norb; ++j)
      {
         Orbit & oj = modelspace.GetOrbit(j);
//         if ( 2*(oi.n+oj.n)+oi.l+oj.l > N2max) continue;
         if (oi.l != oj.l or oi.j2 != oj.j2 or oi.tz2 != oj.tz2) continue;
         double tij = 0;
         if (oi.n == oj.n) tij = 0.5*(2*oi.n+oi.l + 1.5) * hw/A;
         else if (oi.n == oj.n-1) tij = 0.5*sqrt(oj.n*(oj.n+oj.l + 0.5)) * hw/A;
         TcmOp.OneBody(i,j) = tij;
         TcmOp.OneBody(j,i) = tij;
      }
   }

   // Two body piece = 2*p1*p2/(2mA) = (Tcm-Trel)/A
   int nchan = modelspace.GetNumberTwoBodyChannels();
   #pragma omp parallel for schedule(dynamic,5) 
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
//            #pragma omp critical
            if (abs(p1p2)>1e-7)
            {
              TcmOp.TwoBody.SetTBME(ch,ibra,iket,p1p2);
              TcmOp.TwoBody.SetTBME(ch,iket,ibra,p1p2);
            }
         }
      }
   }
   return TcmOp;
 }


 // evaluate <bra| p1*p2 | ket> , omitting the prefactor  m * hbar_omega
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
   if (abs(fab-fcd)%2 >0) return 0; // p1*p2 only connects kets with delta N = 0,1
   if (abs(fab-fcd)>2) return 0; // p1*p2 only connects kets with delta N = 0,1

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
              // factor to account for antisymmetrization
              int asymm_factor = 1;

              if (ket.Tz!=0)
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
                p1p2 += (tcm-trel) * prefactor *0.5;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize
   p1p2 /= sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return p1p2 ;

 }

// Center of mass kinetic energy, with the hw/A factor
 Operator VCM_Op(ModelSpace& modelspace)
 {
   Operator VcmOp = Operator(modelspace);

   int norb = modelspace.GetNumberOrbits();
//   double one_body_prefactor = 1./(
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      for (int j=i; j<norb; ++j)
      {
         Orbit & oj = modelspace.GetOrbit(j);
         if (oi.l != oj.l or oi.j2 != oj.j2 or oi.tz2 != oj.tz2) continue;
         double tij = 0;
         if (oi.n == oj.n) tij = 0.5*(2*oi.n+oi.l + 1.5);
         else if (oi.n == oj.n-1) tij = -0.5*sqrt(oj.n*(oj.n+oj.l + 0.5));
         VcmOp.OneBody(i,j) = tij;
         VcmOp.OneBody(j,i) = tij;
      }
   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   #pragma omp parallel for schedule(dynamic,5) 
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
//            #pragma omp critical
            {
              VcmOp.TwoBody.SetTBME(ch,ibra,iket,mat_el);
              VcmOp.TwoBody.SetTBME(ch,iket,ibra,mat_el);
            }
         }
      }
   }
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   return VcmOp *hw / ( 2*A );
 }

 // Evaluate <bra | r1*r2 | ket>, omitting the factor (hbar * omega) /(m * omega^2)
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
              // factor to account for antisymmetrization
              int asymm_factor = 1;

              if (ket.Tz!=0)
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
                  else if (N_ab == N_cd+1) tcm = -sqrt(N_ab*( N_ab+Lam_ab+0.5));
                  else if (N_ab == N_cd-1) tcm = -sqrt(N_cd*( N_cd+Lam_ab+0.5));
                }
                if (N_ab == N_cd)
                {
                  if      (n_ab == n_cd)   trel = (2*n_ab+lam_ab+1.5);
                  else if (n_ab == n_cd+1) trel = -sqrt(n_ab*( n_ab+lam_ab+0.5));
                  else if (n_ab == n_cd-1) trel = -sqrt(n_cd*( n_cd+lam_cd+0.5));
                }
                double prefactor = njab * njcd * mosh_ab * mosh_cd * asymm_factor;
                r1r2 += (tcm-trel) * prefactor;

              } // N_cd
           } // lam_ab
         } // Lam_ab
       } // N_ab

     } // Sab
   } // Lab

   // normalize
   r1r2 /= sqrt((1.0+bra.delta_pq())*(1.0+ket.delta_pq()));
   return r1r2 ;

 }




// Center of mass kinetic energy, with the hw/A factor
 Operator HCM_Op(ModelSpace& modelspace)
 {
   Operator HcmOp = Operator(modelspace);

   int norb = modelspace.GetNumberOrbits();
   for (int i=0; i<norb; ++i)
   {
      Orbit & oi = modelspace.GetOrbit(i);
      HcmOp.OneBody(i,i) = (2*oi.n+oi.l + 1.5);
   }

   int nchan = modelspace.GetNumberTwoBodyChannels();
   #pragma omp parallel for schedule(dynamic,5) 
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
            mat_el += Calculate_p1p2(modelspace,bra,ket,tbc.J); // added this. not sure if I know what I'm doing...
//            #pragma omp critical
            {
              HcmOp.TwoBody.SetTBME(ch,ibra,iket,mat_el);
              HcmOp.TwoBody.SetTBME(ch,iket,ibra,mat_el);
            }
         }
      }
   }
   double hw = modelspace.GetHbarOmega();
   int A = modelspace.GetTargetMass();
   return HcmOp *hw / ( 2*A );
 }

Operator RSquaredOp(ModelSpace& modelspace)
{
   Operator r2 = Operator(modelspace);
   r2.OneBody.zeros();
   int norbits = modelspace.GetNumberOrbits();
   int A = modelspace.GetTargetMass();
   double hw = modelspace.GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      r2.OneBody(a,a) = (2*oa.n + oa.l +3./2); 
      for (int b=a+1;b<norbits;++b) 
      {
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.l == ob.l and oa.j2 == ob.j2 and oa.tz2 == ob.tz2)
         {
            if (oa.n == ob.n+1)
               r2.OneBody(a,b) = sqrt( (oa.n)*(oa.n + oa.l +1./2));
            else if (oa.n == ob.n-1)
               r2.OneBody(a,b) = sqrt( (ob.n)*(ob.n + ob.l +1./2));
            r2.OneBody(b,a) = r2.OneBody(a,b);
         }
      }
   }
   r2.OneBody *= (HBARC*HBARC/M_NUCLEON/hw);
   return r2;
}

// Electric monopole operator
Operator E0Op(ModelSpace& modelspace)
{
   Operator e0 = Operator(modelspace);
   e0.EraseZeroBody();
   int norbits = modelspace.GetNumberOrbits();
   int A = modelspace.GetTargetMass();
   double hw = modelspace.GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit & oa = modelspace.GetOrbit(a);
      if (oa.tz2 > 0 ) continue; // Only want to count the protons
      e0.OneBody(a,a) = (2*oa.n + oa.l +3./2); 
      for (int b=a+1;b<norbits;++b) 
      {
         Orbit & ob = modelspace.GetOrbit(b);
         if (oa.l == ob.l and oa.j2 == ob.j2 and oa.tz2 == ob.tz2)
         {
            if (oa.n == ob.n+1)
               e0.OneBody(a,b) = sqrt( (oa.n)*(oa.n + oa.l +1./2));
            else if (oa.n == ob.n-1)
               e0.OneBody(a,b) = sqrt( (ob.n)*(ob.n + ob.l +1./2));
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

              if (ket.Tz!=0)
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



 Operator Isospin2_Op(ModelSpace& modelspace)
 {
   Operator T2 = Operator(modelspace,0,0,0,2);
   T2.OneBody.diag().fill(0.75);

   for (int ch=0; ch<T2.nChannels; ++ch)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
//     arma::mat& TB = T2.TwoBody.at(ch).at(ch);
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

}


