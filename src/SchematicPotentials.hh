#ifndef SchematicPotentials_hh
#define SchematicPotentials_hh

#include "AngMom.hh"
#include "ModelSpace.hh"
#include "Operator.hh"
#include "imsrg_util.hh"


namespace SchematicPotentials
{

 // Surface delta interaction
 // see section 8.2 of Suhonen
 //
 Operator SurfaceDelta( ModelSpace& modelspace, float R0, float V0 )
 {
   Operator VSDI = Operator(modelspace, 0,0,0,2);
 
   int nch = modelspace.GetNumberTwoBodyChannels();
   double hw = modelspace.GetHbarOmega();
 
   for ( int ch=0; ch<nch; ch++)
   {
     TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
     int J = tbc.J;
     int nkets = tbc.GetNumberKets();
     for (int ibra=0; ibra<nkets; ibra++)
     {
       Ket& bra = tbc.GetKet(ibra);
       Orbit& oa = modelspace.GetOrbit(bra.p);
       Orbit& ob = modelspace.GetOrbit(bra.q);
       double ja = 0.5*oa.j2;
       double jb = 0.5*ob.j2;
 
       for (int iket=ibra; iket<nkets; iket++)
       {
         Ket& ket = tbc.GetKet(iket);
         Orbit& oc = modelspace.GetOrbit(ket.p);
         Orbit& od = modelspace.GetOrbit(ket.q);
         double jc = 0.5*oc.j2;
         double jd = 0.5*od.j2;
 
         double kappa_ac = R0 * imsrg_util::HO_Radial_psi( oa.n, oa.l, hw, R0)
                              * imsrg_util::HO_Radial_psi( oc.n, oc.l, hw, R0)  ;
         double kappa_bd = R0 * imsrg_util::HO_Radial_psi( ob.n, ob.l, hw, R0)
                              * imsrg_util::HO_Radial_psi( od.n, od.l, hw, R0)  ;
         double Kabcd = -V0 * kappa_ac * kappa_bd / (16.0 * PI);
 
         double Vabcd = 0;
 
         if (tbc.Tz==0)  // Suhonen 8.65
         {  // ignore the 1 + (-1)^(la+lb+lc+ld). This just ensures parity conservation
           Vabcd = Kabcd * 2 * sqrt( (oa.j2+1.)*(ob.j2+1.)*(oc.j2+1.)*(od.j2+1.) )
                    * (    AngMom::ThreeJ(ja,jb,J,0.5,0.5,-1) * AngMom::ThreeJ(jc,jd,J,0.5,0.5,-1)
                       - AngMom::phase(oa.l+oc.l+jb+jd)
                         * AngMom::ThreeJ(ja,jb,J,0.5,-0.5,0) * AngMom::ThreeJ(jc,jd,J,0.5,-0.5,0));
         }
         else  if (J%2 == tbc.parity) // Suhonen 8.66
         {
           Vabcd = -Kabcd * AngMom::phase(oa.l + oc.l + jb+jd) * 2 * 2
                     * sqrt( (oa.j2+1.)*(ob.j2+1.)*(oc.j2+1.)*(od.j2+1.) )
                     * AngMom::ThreeJ(ja,jb,J,0.5,-0.5,0) * AngMom::ThreeJ(jc,jd,J,0.5,-0.5,0);
         }
 
         VSDI.TwoBody.SetTBME( ch, ch, bra.p, bra.q, ket.p, ket.q, Vabcd );
 
       }
     }
 
   }
   return VSDI;
 }



}



#endif
