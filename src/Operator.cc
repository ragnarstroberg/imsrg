
#include "Operator.hh"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

Operator::Operator()
{}

Operator::Operator(ModelSpace* ms) // Create a zero-valued operator in a given model space
{
  modelspace = ms;
  hermitian = true;
  antihermitian = false;
  cross_coupled = false;
  ZeroBody = 0;
  int nOneBody = modelspace->GetNumberOrbits();
  int nKets = modelspace->GetNumberKets();
  OneBody = arma::mat(nOneBody,nOneBody,arma::fill::zeros);
  nChannels = modelspace->GetNumberTwoBodyChannels();

  for (int ch=0;ch<nChannels;++ch)
  {
      int npq = modelspace->GetTwoBodyChannel(ch).GetNumberKets();
      TwoBody[ch] = arma::mat(npq,npq,arma::fill::zeros);
  }

  // Not sure if this stuff is used...
  P_hole_onebody = arma::mat(nOneBody,nOneBody,arma::fill::zeros);
  P_particle_onebody = arma::mat(nOneBody,nOneBody,arma::fill::eye);
  for (int i=0;i<nOneBody;i++)
  {
     if (modelspace->GetOrbit(i)->hvq == 0)
     {
        P_hole_onebody(i,i) = 1;
        P_particle_onebody(i,i) = 0;
     }
  }
  
}

void Operator::Copy(const Operator& op)
{
   modelspace    = op.modelspace;
   nChannels     = modelspace->GetNumberTwoBodyChannels();
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
   cross_coupled = op.cross_coupled;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBody[ch] = op.TwoBody[ch];
   }
   P_hole_onebody = op.P_hole_onebody;
   P_particle_onebody = op.P_particle_onebody;
}


void Operator::PrintTwoBody()
{
   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      cout << "J,p,t = " << tbc.J << " " << tbc.parity << " " << tbc.Tz << endl;
      for (int ik=0;ik<tbc.GetNumberKets();++ik)
      {
         Ket * ket = tbc.GetKet(ik);
         cout << "| " << ket->p << " " << ket->q << " >  ";
      }
      cout << endl;
      TwoBody[ch].print();
      cout << endl;
   }
}

void Operator::PrintOut() 
{
/*
   //for (int j=0; j<JMAX; j++)
   for (int j=0; j<TwoBodyJmax; j++)
   {
     for (int p=0; p<=1; p++)
     {
        for (int t=-1; t<=1; t++)
        {
           TwoBodyChannel *tbc = GetTwoBodyChannel(j,p,t);
//           cout << tbc.J << " " << tbc.parity << " "
//                     << tbc.Tz << "  ===> " << tbc.GetNumberKets() << endl;
           cout << tbc->J << " " << tbc->parity << " "
                     << tbc->Tz << "  ===> " << tbc->GetNumberKets() << endl;
           if (tbc->GetNumberKets()>20) continue;
          // if (tbc.GetNumberKets()>20) continue;
           for (int i=0;i<tbc->GetNumberKets();i++)
           //for (int i=0;i<tbc.GetNumberKets();i++)
           {
             //for (int ii=0;ii<tbc.GetNumberKets();ii++)
             for (int ii=0;ii<tbc->GetNumberKets();ii++)
             {
                //cout << tbc.TBME(i,ii) << " ";
                cout << tbc->TBME(i,ii) << " ";
             }
                cout << endl;
           }
        }
     }
   }
*/
}


void Operator::Eye()
{
   ZeroBody = 1;
   OneBody.eye();
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBody[ch].eye();
   }
}



double Operator::GetTBME(int ch, int a, int b, int c, int d) const
{
   auto& tbc = cross_coupled ? modelspace->GetTwoBodyChannel_CC(ch) : modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc.GetNumberKets() or ket_ind > tbc.GetNumberKets() )
     return 0;
   Ket * bra = tbc.GetKet(bra_ind);
   Ket * ket = tbc.GetKet(ket_ind);

   double phase = 1;
   if (a>b) phase *= bra->Phase(tbc.J) ;
   if (c>d) phase *= ket->Phase(tbc.J);
   if (a==b) phase *= sqrt(1.0+bra->delta_pq());
   if (c==d) phase *= sqrt(1.0+ket->delta_pq());
   return phase * TwoBody[ch](bra_ind, ket_ind);
}

double Operator::GetTBME_NoPhase(int ch, int a, int b, int c, int d) const
{
   auto& tbc = cross_coupled ? modelspace->GetTwoBodyChannel_CC(ch) : modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc.GetNumberKets() or ket_ind > tbc.GetNumberKets() )
     return 0;
   Ket * bra = tbc.GetKet(bra_ind);
   Ket * ket = tbc.GetKet(ket_ind);

   double phase = 1;
   //if (a>b) phase *= bra->Phase(tbc.J) ;
   if (a>b) bra_ind += tbc.GetNumberKets() ;
   if (c>d) ket_ind += tbc.GetNumberKets() ;
   if (a==b) phase *= sqrt(1.0+bra->delta_pq());
   if (c==d) phase *= sqrt(1.0+ket->delta_pq());
   return phase * TwoBody[ch](bra_ind, ket_ind);
}


void Operator::SetTBME(int ch, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
//   cout << "in SetTBME  ch = " << ch << " abcd = " << a << b << c << d << " bra_ind = " << bra_ind << " ket_ind = " << ket_ind << endl;
   double phase = 1;
   if (a>b) phase *= tbc.GetKet(bra_ind)->Phase(tbc.J);
   if (c>d) phase *= tbc.GetKet(ket_ind)->Phase(tbc.J);
   TwoBody[ch](bra_ind,ket_ind) = phase * tbme;
   if (hermitian) TwoBody[ch](ket_ind,bra_ind) = phase * tbme;
   if (antihermitian) TwoBody[ch](ket_ind,bra_ind) = - phase * tbme;
}


double Operator::GetTBME(int ch, Ket *bra, Ket *ket) const
{
   return GetTBME(ch,bra->p,bra->q,ket->p,ket->q);
}

void Operator::SetTBME(int ch, Ket* ket, Ket* bra, double tbme)
{
   SetTBME(ch, bra->p,bra->q,ket->p,ket->q,tbme);
}

double Operator::GetTBME(int j, int p, int t, int a, int b, int c, int d) const
{
   //int ch = (t+1)*2*JMAX + p*JMAX + j;
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,a,b,c,d);
}

void Operator::SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   //int ch = (t+1)*2*JMAX + p*JMAX + j; // This should be altered to refer to the modelspace method.
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,a,b,c,d,tbme);
}

double Operator::GetTBME(int j, int p, int t, Ket* bra, Ket* ket) const
{
   //int ch = (t+1)*2*JMAX + p*JMAX + j;
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,bra,ket);
}

void Operator::SetTBME(int j, int p, int t, Ket* bra, Ket* ket, double tbme)
{
   //int ch = (t+1)*2*JMAX + p*JMAX + j;
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,bra,ket,tbme);
}


double Operator::GetTBMEmonopole(int a, int b, int c, int d)
{
   double mon = 0;
   Orbit *oa = modelspace->GetOrbit(a);
   Orbit *ob = modelspace->GetOrbit(b);
   Orbit *oc = modelspace->GetOrbit(c);
   Orbit *od = modelspace->GetOrbit(d);
   int Tzab = (oa->tz2 + ob->tz2)/2;
   int parityab = (oa->l + ob->l)%2;
   int Tzcd = (oc->tz2 + od->tz2)/2;
   int paritycd = (oc->l + od->l)%2;

   if (Tzab != Tzcd or parityab != paritycd) return 0;

   int jmin = abs(oa->j2 - ob->j2)/2;
   int jmax = (oa->j2 + ob->j2)/2;
   
   for (int J=jmin;J<=jmax;++J)
   {
      mon += (2*J+1) * GetTBME(J,parityab,Tzab,a,b,c,d);
   }
   mon /= (oa->j2 +1)*(ob->j2+1);
//   if (a==b) mon /= sqrt(2);
//   if (c==d) mon /= sqrt(2);
   return mon;
}

double Operator::GetTBMEmonopole(Ket * bra, Ket * ket)
{
   return GetTBMEmonopole(bra->p,bra->q,ket->p,ket->q);
}


// multiply operator by a scalar
Operator& Operator::operator*=(const double rhs)
{
   ZeroBody *= rhs;
   OneBody *= rhs;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBody[ch] *= rhs;
   }
   return *this;
}

Operator Operator::operator*(const double rhs)
{
   Operator opout = Operator(*this);
   opout *= rhs;
   return opout;
}

// Add non-member operator so we can multiply an operator
// by a scalar from the lhs, i.e. s*O = O*s
inline Operator operator*(const double lhs, Operator& rhs)
{
   return rhs * lhs;
}

Operator& Operator::operator+=(const Operator& rhs)
{
   ZeroBody += rhs.ZeroBody;
   OneBody += rhs.OneBody;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBody[ch] += rhs.TwoBody[ch];
   }
   return *this;
}

Operator Operator::operator+(const Operator& rhs)
{
   return ( Operator(*this) += rhs );
}


Operator Operator::DoNormalOrdering()
{
   Operator opNO = Operator(*this);
   opNO.ZeroBody = ZeroBody;
   opNO.OneBody = OneBody;
   int norbits = modelspace->GetNumberOrbits();

   //for (int k=0;k<norbits;++k)
   for (int &k : modelspace->holes) // C++11 syntax
   {
//      Orbit * orb = modelspace->GetOrbit(k);
//      if (orb->hvq > 0) continue;
      opNO.ZeroBody += (modelspace->GetOrbit(k)->j2+1) * OneBody(k,k);
   }
//   cout << "after first loop, zero body part is" << opNO.ZeroBody << endl;


   //////////////////////// v--This works, but it's the slow way (optimize if needed) /////////////////////
   for (int ch=0;ch<nChannels;++ch)
   {
      //TwoBodyChannel *tbc = GetTwoBodyChannel(ch);
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
      //opNO.SetTwoBody(ch, TwoBody[ch] );
      opNO.TwoBody[ch] = TwoBody[ch];
      int npq = tbc.GetNumberKets();
      int J = tbc.J;

      int ibra,iket;
      for (int a=0;a<norbits;++a)
      {
         Orbit *orba = modelspace->GetOrbit(a);
         for (int b=a;b<norbits;++b)
         {
            Orbit *orbb = modelspace->GetOrbit(b);
            if (orbb->j2 != orba->j2 or orbb->tz2 != orba->tz2 or orbb->l != orba->l) continue;
            //for (int c=0;c<norbits;++c)
            for (int &h : modelspace->holes)  // C++11 syntax
            {
              // if ( modelspace->GetOrbit(c)->hvq >0 ) continue;
               if ( (ibra = tbc.GetLocalIndex(min(a,h),max(a,h)))<0) continue;
               if ( (iket = tbc.GetLocalIndex(min(b,h),max(b,h)))<0) continue;
               Ket * bra = tbc.GetKet(ibra);
               Ket * ket = tbc.GetKet(iket);
//               double tbme = (2*J+1.0)/(orba->j2+1) * GetTBME(ch,bra,ket);
               //double tbme = (2*J+1.0)/(orba->j2+1) * sqrt( (1+bra->delta_pq())*(1+ket->delta_pq()) ) * GetTBME(ch,bra,ket);
               double tbme = (2*J+1.0)/(orba->j2+1)  * GetTBME(ch,bra,ket);

               if (a<=h and b<=h) opNO.OneBody(a,b) += tbme;
               if (a>h and b<=h) opNO.OneBody(a,b) += bra->Phase(J)*tbme;
               if (a<=h and b>h) opNO.OneBody(a,b) += ket->Phase(J)*tbme;
               if (a>h and b>h) opNO.OneBody(a,b) += bra->Phase(J)*ket->Phase(J)*tbme;
            }
            if (hermitian) opNO.OneBody(b,a) = opNO.OneBody(a,b);
            if (antihermitian) opNO.OneBody(b,a) = -opNO.OneBody(a,b);
         }
      }
   //////////////////////// ^-- This works, but it's the slow way (I think) /////////////////////


      for (int k=0;k<npq;++k) // loop over kets in this channel
      {
         Ket * ket = tbc.GetKet(k);
         Orbit *oa = modelspace->GetOrbit(ket->p);
         Orbit *ob = modelspace->GetOrbit(ket->q);

//         if (oa->hvq > 0 and ob->hvq > 0) continue; // need at least one orbit in the core

         // Sum over diagonal two-body terms with both orbits in the core
         if ( oa->hvq == 0 and ob->hvq == 0)
         {
//            opNO.ZeroBody += TwoBody[ch](k,k) * (2*tbc.J+1) / (1.0+ket->delta_pq());  // <ab|V|ab>  (a,b in core)
            opNO.ZeroBody += TwoBody[ch](k,k) * (2*tbc.J+1);  // <ab|V|ab>  (a,b in core)
         }

      } // for k
   } // for ch
//   cout << "after second loop, zero body part is" << opNO.ZeroBody << endl;
   return opNO;
}



void Operator::EraseTwoBody()
{
   for (int ch=0;ch<nChannels;++ch)
   {
       TwoBody[ch].zeros();
   }
}




//*****************************************************************************************
//   Transform to a cross-coupled basis, for use in the 2body commutator relation
//    using a Pandya tranformation 
//                                               ____________________________________________________________________________
//                                              |                                                                            |
//   |a    |b  <ab|_J'    \a   /c  <ac|_J       |  <ac|V|bd>_J = sum_J' (2J'+1) (-1)^(ja+jd+J+J'+1) { ja jc J } <ab|V|cd>_J' |
//   |     |               \  /                 |                                                   { jd jb J'}              |
//   |_____|                \/_____             |____________________________________________________________________________|
//   |  V  |      ===>         V  /\ 
//   |     |                     /  \
//   |c    |d  |cd>_J'          /b   \d   |bd>_J
//
//   STANDARD                 CROSS  
//   COUPLING                 COUPLED  
//                                      
void Operator::UpdateCrossCoupled()
{
   // loop over cross-coupled channels
   int norbits = modelspace->GetNumberOrbits();
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
//      TwoBody_CC_left[ch_cc] = arma::mat(nKets_cc,nKets_cc,arma::fill::zeros);
//      TwoBody_CC_right[ch_cc] = arma::mat(nKets_cc,nKets_cc,arma::fill::zeros);

//// These matrices don't actually need to be square. Save some space by fixing that...
//      TwoBody_CC_left[ch_cc] = arma::mat(2*nKets_cc,2*nKets_cc,arma::fill::zeros);
//      TwoBody_CC_right[ch_cc] = arma::mat(2*nKets_cc,2*nKets_cc,arma::fill::zeros);
      TwoBody_CC_left[ch_cc] = arma::mat(2*nKets_cc,2*nKets_cc,arma::fill::zeros);
      TwoBody_CC_right[ch_cc] = arma::mat(2*nKets_cc,2*nKets_cc,arma::fill::zeros);

      // loop over cross-coupled bras <ac| in this channel
      //for (int ibra_cc=0; ibra_cc<nKets_cc; ++ibra_cc)
      //for (int a=0; a<norbits; ++a)
      for (int& a : modelspace->holes)
      {
        //for (int c=0; c<norbits; ++c)
        for (int& c : modelspace->particles)
        {
         int ibra_cc = tbc_cc.GetLocalIndex(min(a,c),max(a,c));
         if (ibra_cc < 0) continue;
         Ket * bra_cc = tbc_cc.GetKet(ibra_cc);
//         if (a>c) ibra_cc += nKets_cc;
         //Ket * bra_cc = tbc_cc.GetKet(ibra_cc);
//         int a = bra_cc->p;
//         int c = bra_cc->q;
         Orbit * oa = modelspace->GetOrbit(a);
         Orbit * oc = modelspace->GetOrbit(c);
         double ja = oa->j2/2.0;
         double jc = oc->j2/2.0;

         // loop over cross-coupled kets |bd> in this channel
         //for (int iket_cc=0; iket_cc<nKets_cc; ++iket_cc)
         for (int b=0; b<norbits; ++b)
         {
         for (int d=0; d<norbits; ++d)
         {
           int iket_cc = tbc_cc.GetLocalIndex(min(b,d),max(b,d));
           if (iket_cc < 0) continue;
           Ket * ket_cc = tbc_cc.GetKet(iket_cc);
           if (b>d) iket_cc += nKets_cc;
       //     Ket * ket_cc = tbc_cc.GetKet(iket_cc);
       //     int b = ket_cc->p;
       //     int d = ket_cc->q;
            Orbit * ob = modelspace->GetOrbit(b);
            Orbit * od = modelspace->GetOrbit(d);
            double jb = ob->j2/2.0;
            double jd = od->j2/2.0;

            // Get Tz,parity and range of J for <ab || cd > coupling
            int Tz_std = (oa->tz2 + ob->tz2)/2;
            int parity_std = (oa->l + ob->l)%2;
//            if ( (oc->tz2 + od->tz2)/2 != Tz_std) continue;
//            if ( (oc->l + od->l)%2 != parity_std) continue;

            //int jmin = max(abs(oa->j2-ob->j2),abs(oc->j2-od->j2))/2;
            //int jmax = min(oa->j2+ob->j2,oc->j2+od->j2)/2;
            int jmin = max(abs(ja-jb),abs(jc-jd));
            int jmax = min(ja+jb,jc+jd);
//               if (ch_cc==0 )
//                  cout << endl << "a,b,c,d = " << a << " " << b << " " << c << " " << d << endl;


            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               //double sixj = modelspace->GetSixJ(oa->j2/2.,oc->j2/2.,tbc_cc.J,od->j2/2.,ob->j2/2.,J_std);
               double sixj = modelspace->GetSixJ(ja,jc,tbc_cc.J,jd,jb,J_std);
               if (sixj == 0) continue;
               //int phase = 1-2*(int(oa->j2+od->j2)/2+tbc_cc.J+J_std+1)%2;
               //int phase = 1-2*((int(oa->j2+od->j2)/2+tbc_cc.J+J_std)%2);
               int phase = modelspace->phase(J_std);
               double tbme = GetTBME(J_std,parity_std,Tz_std,a,b,c,d);
               sm += (2*J_std+1) * phase * sixj * tbme * sqrt(2*tbc_cc.J+1); // added in for convenience...
               //sm += (2*J_std+1) * phase * sixj * tbme;
               //if ( (ch_cc==120 and ibra_cc==10 and iket_cc==10) or (ch_cc==150 and ibra_cc==0 and iket_cc==0) )
               //if ( (ch_cc==151 and ibra_cc==0 and iket_cc==0) or (ch_cc==150 and ibra_cc==0 and iket_cc==0) )
//             if ( (ch_cc==151) )
//             {
//                cout << " Jcc = " << tbc_cc.J << " pcc = " << tbc_cc.parity << " Tz_cc = " << tbc_cc.Tz
//                     << " ibra_cc = " << ibra_cc << " iket_cc = " << iket_cc << endl;
//                cout << "  J= " << J_std << " p= " << parity_std << " Tz = " << Tz_std << " phase = " << phase << " sixj = " << sixj 
//                     << " < " << a << " " << b << " | V | " << c << " " << d  << " > = " << tbme << endl;
//                cout << "abcd = " << a << "," << b << "," << c << "," << d << endl;
//                cout << "    { " << ja << "  " << jc << "  " << tbc_cc.J << " }" << endl;
//                cout << "    { " << jd << "  " << jb << "  " << J_std    << " }" << endl;
//                cout << "   TwoBody_CC_left: " << sm << endl;
//                cout << endl;
//             }
            }
            //TwoBody_CC_left[ch_cc](ibra_cc,iket_cc) = sm;
            //TwoBody_CC_left[ch_cc](iket_cc,ibra_cc) = sm * modelspace->phase( (oa->j2+od->j2)/2+tbc_cc.J );
            TwoBody_CC_left[ch_cc](iket_cc,ibra_cc) = sm * modelspace->phase( ja+jd+tbc_cc.J );

            // Get Tz,parity and range of J for <cb || ad > coupling
            Tz_std = (oa->tz2 + od->tz2)/2;
            parity_std = (oa->l + od->l)%2;
            //jmin = max(abs(oc->j2-ob->j2),abs(oa->j2-od->j2))/2;
            //jmax = min(oc->j2+ob->j2,oa->j2+od->j2)/2;
            jmin = max(abs(jc-jb),abs(ja-jd));
            jmax = min(jc+jb,ja+jd);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               //double sixj = modelspace->GetSixJ(oa->j2/2.,oc->j2/2.,tbc_cc.J,ob->j2/2.,od->j2/2.,J_std);
               double sixj = modelspace->GetSixJ(ja,jc,tbc_cc.J,jb,jd,J_std);
               if (sixj == 0) continue;
               //int phase = 1-2*(int(oa->j2+od->j2)/2+tbc_cc.J+J_std+1)%2;
               //int phase = 1-2*((int(oa->j2+ob->j2)/2+tbc_cc.J+J_std)%2);
               int phase = modelspace->phase(J_std);
               //double tbme = GetTBME(J_std,parity_std,Tz_std,a,b,c,d);
               double tbme = GetTBME(J_std,parity_std,Tz_std,c,b,a,d);
               //sm += (2*J_std+1) * phase * sixj * tbme;
               sm += (2*J_std+1) * phase * sixj * tbme * sqrt(2*tbc_cc.J+1); // sqrt(2J+1) added for convenience?
               //if (ch_cc==120 and ibra_cc==10 and iket_cc==10)
               //if ( (ch_cc==120 and ibra_cc==10 and iket_cc==10) or (ch_cc==150 and ibra_cc==0 and iket_cc==0) )
               //if ( (ch_cc==151 and ibra_cc==0 and iket_cc==0) or (ch_cc==150 and ibra_cc==0 and iket_cc==0) )
//             if ( (ch_cc==151) )
//             {
//                cout << " Jcc = " << tbc_cc.J << " pcc = " << tbc_cc.parity << " Tz_cc = " << tbc_cc.Tz
//                     << " ibra_cc = " << ibra_cc << " iket_cc = " << iket_cc << endl;
//                cout << " J= " << J_std << " p= " << parity_std << " Tz = " << Tz_std << " phase = " << phase << " sixj = " << sixj
//                     << " < " << c << " " << b << " | V | " << a << " " << d  << " > = " << tbme << endl;
//                cout << "cbad = " << c << "," << b << "," << a << "," << d << endl;
//                cout << "    { " << ja << "  " << jc << "  " << tbc_cc.J << " }" << endl;
//                cout << "    { " << jb << "  " << jd << "  " << J_std    << " }" << endl;
//                cout << "   TwoBody_CC_right: " << sm << endl;
//                cout << endl;
//             }
            }
            //TwoBody_CC_right[ch_cc](ibra_cc,iket_cc) = sm * modelspace->phase(ja+jb+tbc_cc.J);
            TwoBody_CC_right[ch_cc](ibra_cc,iket_cc) = sm * modelspace->phase(ja+jd+tbc_cc.J+1);
//            TwoBody_CC_right[ch_cc](iket_cc,ibra_cc) = sm;


         }
      }
      
   }
   }
   }
}


//*****************************************************************************************
//  OpOut = exp(Omega) Op exp(-Omega)
//  Eventually, add some checking for convergence.
Operator Operator::BCH_Transform( Operator &Omega)
{
   Operator OpOut = *this;
   Operator OpNested = *this;
   double prefactor = 1.0;
   for (int i=1; i<6; ++i)
   {
      prefactor /= i;
      OpNested = Omega.Commutator(OpNested);
      OpOut += OpNested * prefactor;
   }
   return OpOut;
}


//*****************************************************************************************
//  OpOut = Omega_N+1, where
//  exp(Omega_N+1) = exp(dOmega) * exp(Omega_N)
//  Omega_N+1 = Omega_N + dOmega + 1/2[dOmega, Omega_N] + 1/12 [dOmega,[dOmega,Omega_N]] + 1/12 [Omega_N,[Omega_N,dOmega]] + ...
//  Eventually, add some checking for convergence.
Operator Operator::BCH_Product( Operator &dOmega)
{
   Operator OpOut = *this;
   //Operator OpNested = dOmega.Commutator(*this);
   Operator OpNested = Commutator(dOmega);
   OpOut += dOmega;
   OpOut += OpNested*(-1./2);
   //OpOut += ( dOmega.Commutator(OpNested) + OpNested.Commutator(*this) )*(1./12);
   Operator OpNested2 = Commutator(OpNested);
   OpOut += OpNested2*(1/12.);
   OpOut += dOmega.Commutator(OpNested2)*(-1/720);
   return OpOut;
}

// Frobenius norm of the operator
double Operator::Norm()
{
   //double nrm = ZeroBody*ZeroBody;
   double nrm = 0;
   double n1 = arma::norm(OneBody,"fro");
//   double n1 = arma::accu(OneBody.t()*OneBody);
   nrm += n1*n1;
//   nrm += n1;
//   cout << "Onebody =  " << endl; OneBody.print();
//   cout << "One body norm = " << n1 << endl;
//   cout << "One body norm = " << sqrt(n1) << endl;
   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      double n2 = arma::norm(TwoBody[ch],"fro");
//      double n2 = arma::accu(TwoBody[ch].t()*TwoBody[ch]);
//      cout << "in Norm. ch = " << ch << " n2 = " << n2*n2 << endl;
      nrm += n2*n2;
//      nrm += n2;
//      nrm += pow(arma::norm(TwoBody[ch],"fro"),2);
   }
   return sqrt(nrm);
}



Operator Operator::Commutator(Operator& opright)
{
   Operator out = opright;
//   out.EraseZeroBody();
//   out.EraseOneBody();
   out.EraseTwoBody();

   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   out.ZeroBody  = comm110(opright) + comm220(opright) ;

   out.OneBody  = comm111(opright) + comm121(opright) + comm221(opright);

   UpdateCrossCoupled();
   opright.UpdateCrossCoupled();

   comm122(opright,out);
   comm222_pp_hh(opright,out);
   comm222_ph(opright,out);

   return out;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
///////////////////////////////////////////////////////////////////////////////////////////

//*****************************************************************************************
//                ____Y    __         ____X
//          X ___(_)             Y___(_) 
//
//  [X1,Y1](0) = Sum_ab (2j_a+1) x_ab y_ba  (n_a-n_b) 
//             = Sum_a  (2j_a+1)  (xy-yx)_aa n_a
//
// -- AGREES WITH NATHAN'S RESULTS
double Operator::comm110(Operator& opright)
{
  if (IsHermitian() and opright.IsHermitian()) return 0; // I think this is the case
  if (IsAntiHermitian() and opright.IsAntiHermitian()) return 0; // I think this is the case

   double comm = 0;
   int norbits = modelspace->GetNumberOrbits();
   arma::mat xyyx = OneBody*opright.OneBody - opright.OneBody*OneBody;
   for ( int& a : modelspace->holes) // C++11 range-for syntax: loop over all elements of the vector <particle>
   {
      comm += (modelspace->GetOrbit(a)->j2+1)*xyyx(a,a);
   }
   return comm;
}


//*****************************************************************************************
//         __Y__       __X__
//        ()_ _()  -  ()_ _()
//           X           Y
//
//  [ X^(2), Y^(2) ]^(0) = 1/2 Sum_abcd  Sum_J (2J+1) x_abcd y_cdab (n_a n_b nbar_c nbar_d)
//                       = 1/2 Sum_J (2J+1) Sum_abcd x_abcd y_cdab (n_a n_b nbar_c nbar_d)  
//                       = 1/2 Sum_J (2J+1) Sum_ab  (X*P_pp*Y)_abab  P_hh
//
//  -- AGREES WITH NATHAN'S RESULTS (within < 1%)
double Operator::comm220( Operator& opright)
{
   if (IsHermitian() and opright.IsHermitian()) return 0; // I think this is the case
   if (IsAntiHermitian() and opright.IsAntiHermitian()) return 0; // I think this is the case
   double comm = 0;
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      comm += 2 * (2*tbc.J+1) * arma::trace( tbc.Proj_hh * TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] );

   }
   return comm;
}

//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |
//
// -- AGREES WITH NATHAN'S RESULTS
arma::mat Operator::comm111(Operator & opright)
{
   return OneBody*opright.OneBody - opright.OneBody*OneBody;
}

//*****************************************************************************************
//                                       |
//      i |              i |             |
//        |    ___.Y       |__X__        |
//        |___(_)    _     |   (_)__.    |  [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y    |        * sum_J (2J+1) x_biaj^(J)  
//                                       |      
//---------------------------------------*        = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
// -- AGREES WITH NATHAN'S RESULTS 
arma::mat Operator::comm121(Operator& opright)
{
   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(norbits,norbits,arma::fill::zeros);
   for (int i=0;i<norbits;++i)
   {
      Orbit *oi = modelspace->GetOrbit(i);
      for (int j=0;j<norbits; ++j) // Later make this j=i;j<norbits... and worry about Hermitian vs anti-Hermitian
      {
          Orbit *oj = modelspace->GetOrbit(j);
          if (oi->j2 != oj->j2 or oi->l != oj->l or oi->tz2 != oj->tz2) continue; // At some point, make a OneBodyChannel class...
          for (int &a : modelspace->holes)  // C++11 syntax
          {
             Orbit *oa = modelspace->GetOrbit(a);
             for (int b=0;b<norbits;++b)
             {
                Orbit *ob = modelspace->GetOrbit(b);
                comm(i,j) += (ob->j2+1) *  OneBody(a,b) * opright.GetTBMEmonopole(b,i,a,j) ;
                comm(i,j) -= (oa->j2+1) *  OneBody(b,a) * opright.GetTBMEmonopole(a,i,b,j) ;

                // comm211 part
                comm(i,j) -= (ob->j2+1) *  opright.OneBody(a,b) * GetTBMEmonopole(b,i,a,j) ;
                comm(i,j) += (oa->j2+1) *  opright.OneBody(b,a) * GetTBMEmonopole(a,i,b,j) ;
             }
          }
      }
   }
   return comm;
}


// Same as above, just reversed
//arma::mat Operator::comm211(Operator& opright)
//{
//   return -opright.comm121(*this);
//}


//*****************************************************************************************
//
//      i |              i |            [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//        |__Y__           |__X__           * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//        |    /\          |    /\          * (x_ciab y_abcj - y_ciab xabcj)
//        |   (  )   _     |   (  )                                                                                      
//        |____\/          |____\/       = 1/(2(2j+1)) sum_J (2J+1)
//      j | X            j |  Y            *  sum_c ( Pp*X*Phh*Y*Pp - Pp*Y*Phh*X*Pp)  - (Ph*X*Ppp*Y*Ph - Ph*Y*Ppp*X*Ph)_cicj
//                                     
//
// -- AGREES WITH NATHAN'S RESULTS 
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
arma::mat Operator::comm221(Operator& opright)
{

   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(norbits,norbits,arma::fill::zeros);
   Operator Mpp = opright;
   Operator Mhh = opright;

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      
      Mpp.TwoBody[ch] = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
      Mhh.TwoBody[ch] = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);

      // If commutator is hermitian or antihermitian, we only
      // need to do half the sum. Add this.
      for (int i=0;i<norbits;++i)
      {
         Orbit *oi = modelspace->GetOrbit(i);
         for (int j=0;j<norbits;++j)
         {
            Orbit *oj = modelspace->GetOrbit(j);
            if (oi->j2 != oj->j2 or oi->l != oj->l or oi->tz2 != oj->tz2) continue;
            double cijJ = 0;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->holes)
            {
               cijJ +=   Mpp.GetTBME(ch,i,c,j,c);
            // Sum c over particles and include the n_a * n_b terms
            }
            for (int &c : modelspace->particles)
            {
               cijJ +=  Mhh.GetTBME(ch,i,c,j,c);
            }
            comm(i,j) += (2*tbc.J+1.0)/(oi->j2 +1.0) * cijJ;
         } // for j
      } // for i
   } //for ch
   return comm;
}


//*****************************************************************************************
//
//    |     |               |      |           [X2,Y1](2) = sum_a ( Y_ia X_ajkl + Y_ja X_iakl - Y_ak X_ijal - Y_al X_ijka )
//    |     |___.Y          |__X___|         
//    |     |         _     |      |          
//    |_____|               |      |_____.Y        
//    |  X  |               |      |            
//
// -- AGREES WITH NATHAN'S RESULTS
void Operator::comm122(Operator& opright, Operator& opout )
{
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      int norbits = modelspace->GetNumberOrbits();
      for (int ibra = 0;ibra<npq; ++ibra)
      {
         Ket * bra = tbc.GetKet(ibra);
         int i = bra->p;
         int j = bra->q;
         for (int iket=0;iket<npq; ++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            int k = ket->p;
            int l = ket->q;
            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
               cijkl += OneBody(i,a) * opright.GetTBME(ch,a,j,k,l);
               cijkl += OneBody(j,a) * opright.GetTBME(ch,i,a,k,l);
               cijkl -= OneBody(a,k) * opright.GetTBME(ch,i,j,a,l);
               cijkl -= OneBody(a,l) * opright.GetTBME(ch,i,j,k,a);

               // comm212 terms
               cijkl -= opright.OneBody(i,a) * GetTBME(ch,a,j,k,l);
               cijkl -= opright.OneBody(j,a) * GetTBME(ch,i,a,k,l);
               cijkl += opright.OneBody(a,k) * GetTBME(ch,i,j,a,l);
               cijkl += opright.OneBody(a,l) * GetTBME(ch,i,j,k,a);
            }
            opout.TwoBody[ch](ibra,iket) += cijkl / sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()) );
         }
      }
   }
}



//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_pp(hh) = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
// -- AGREES WITH NATHAN'S RESULTS
//   No factor of 1/2 because the matrix multiplication corresponds to a restricted sum (a<=b) 
void Operator::comm222_pp_hh(Operator& opright, Operator& opout )
{
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat Mpp = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
      arma::mat Mhh = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);
      opout.TwoBody[ch] += Mpp - Mhh;
   }
}


void Operator::comm222_ph_slow(Operator& opright, Operator& opout )
{
   for (int ch=0;ch<nChannels;++ch)
   {
//      if (ch>2) break;
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int parity = tbc.parity;
      int Tz = tbc.Tz;
      int nkets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nkets; ++ibra)
      {
         Ket * bra = tbc.GetKet(ibra);
         int i = bra->p;
         int j = bra->q;
         Orbit * oi = modelspace->GetOrbit(i);
         Orbit * oj = modelspace->GetOrbit(j);
         double ji = oi->j2/2.;
         double jj = oj->j2/2.;
//         cout << "ji = " << ji << "  jj = " << jj << endl;
         for (int iket=0; iket<nkets; ++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            int k = ket->p;
            int l = ket->q;
            Orbit * ok = modelspace->GetOrbit(k);
            Orbit * ol = modelspace->GetOrbit(l);
            double jk = (ok->j2/2.);
            double jl = (ol->j2/2.);
            double commijkl = 0;
//            int phij = 1-2*((oi->j2+oj->j2)%4);
//            int phjk = 1-2*((oj->j2+oj->j2)%4);
//            int phil = 1-2*((oi->j2+ok->j2)%4);
//            int phjl = 1-2*((oj->j2+ok->j2)%4);
            for (int& a : modelspace->holes)
            {
               Orbit * oa = modelspace->GetOrbit(k);
               double ja = (oa->j2/2.);
               for (int & b : modelspace->particles)
               {
                  Orbit * ob = modelspace->GetOrbit(k);
                  double jb = (ob->j2/2.);
                  int j1min = 0;
                  int j1max = 6;  // you can do better than this ...
                  int j2min = 0;
                  int j2max = 6;
                  for (int J1=j1min; J1<=j1max; ++J1)
                  {
                     int ch1 = modelspace->GetTwoBodyChannelIndex(J1,parity,Tz); 
                     for (int J2=j2min; J2<=j2max; ++J2)
                     {
                        int ch2 = modelspace->GetTwoBodyChannelIndex(J2,parity,Tz); 
                        //double pref = 1-2*((J1+J2+J)%2) * (2*J1+1)*(2*J2+1);
                        double pref = (1-2*((int(ji+jk)+J1+J2)%2)) * (2*J1+1)*(2*J2+1);
                        int pref2 = 1-2*(J%2);
                        double nj1 = modelspace->GetNineJ( ja, jl, J1,   ji, J, jj,   J2, jk, jb );
                        double nj2 = modelspace->GetNineJ( ja, jl, J1,   jj, J, ji,   J2, jk, jb );
                        double nj3 = modelspace->GetNineJ( jb, jl, J1,   ji, J, jj,   J2, jk, ja );
                        double nj4 = modelspace->GetNineJ( jb, jl, J1,   jj, J, ji,   J2, jk, ja );

                        double me1 =         GetTBME(ch1,b,j,a,l) * opright.GetTBME(ch2,a,i,b,k) - opright.GetTBME(ch1,b,j,a,l) *         GetTBME(ch2,a,i,b,k);
                        double me2 = opright.GetTBME(ch1,b,i,a,l) *         GetTBME(ch2,a,j,b,k) -         GetTBME(ch1,b,i,a,l) * opright.GetTBME(ch2,a,j,b,k);
                        double me3 =         GetTBME(ch1,a,j,b,l) * opright.GetTBME(ch2,b,i,a,k) - opright.GetTBME(ch1,a,j,b,l) *         GetTBME(ch2,b,i,a,k);
                        double me4 = opright.GetTBME(ch1,a,i,b,l) *         GetTBME(ch2,b,j,a,k) -         GetTBME(ch1,a,i,b,l) * opright.GetTBME(ch2,b,j,a,k);

                        commijkl += pref * (nj1 * me1 * pref2 - nj2 * me2 );
                        commijkl -= pref * (nj3 * me3 * pref2 - nj4 * me4 );
                     }
                  }
               }
            }
            opout.TwoBody[ch](ibra,iket) = commijkl;
         }
      }
   }
}
//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.             [X2,Y2](2)_ph = -sum_ab (na-nb) sum_J' (-1)^(J'+ja+jb)
//
//                                                                             v----|    v----|                                    v----|    v----|
//   |          |      |          |              * [ (-1)^(ji+jl) { jk jl J } <bi|X|aj> <ak|Y|bl>   - (-1)^(jj+jl-J') { jk jl J } <bj|X|ai> <ak|Y|bl>
//   |     __Y__|      |     __X__|                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//   |    /\    |      |    /\    |
//   |   (  )   |  _   |   (  )   |                                            v----|    v----|                                    v----|    v----|
//   |____\/    |      |____\/    |            -  (-1)^(ji+jk-J') { jk jl J } <bi|X|aj> <al|Y|bk>   + (-1)^(jj+jk)    { jk jl J } <bj|X|ai> <al|Y|bk>  ]
//   |  X       |      |  Y       |                               { jj ji J'}   ^----|    ^----|                      { jj ji J'}   ^----|    ^----|
//           
//            
// -- This appears to agree with Nathan's results
//
void Operator::comm222_ph(Operator& opright, Operator& opout )
{

   int herm = -1;
   if ( (IsHermitian() and opout.IsAntiHermitian()) or ( IsAntiHermitian() and opright.IsHermitian() ) ) herm = 1;
   // Construct the intermediate matrices N1 and N2
   // These are implemented as operators in order to use
   // the TBME structure and accessor functions
   arma::mat N1[nChannels];
   arma::mat N2[nChannels];
   for (int ch=0;ch<nChannels;++ch)
   {
      // maybe implicitly include the Proj_ph_cc part into the right hand operators?
      TwoBodyChannel_CC& tbc_CC = modelspace->GetTwoBodyChannel_CC(ch);
      N1[ch] =         TwoBody_CC_left[ch] *  tbc_CC.Proj_ph_cc  *  opright.TwoBody_CC_right[ch];
      N2[ch] = opright.TwoBody_CC_left[ch] *  tbc_CC.Proj_ph_cc  *          TwoBody_CC_right[ch];

   }

   // Now evaluate the commutator for each channel (standard coupling)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int nKets = tbc.GetNumberKets();
      for (int ibra=0; ibra<nKets; ++ibra)
      {
         Ket * bra = tbc.GetKet(ibra);
         int i = bra->p;
         int j = bra->q;
         Orbit* oi = modelspace->GetOrbit(i);
         Orbit* oj = modelspace->GetOrbit(j);
         double ji = oi->j2/2.;
         double jj = oj->j2/2.;

         //for (int iket=0; iket<nKets; ++iket)
         for (int iket=ibra; iket<nKets; ++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            int k = ket->p;
            int l = ket->q;
            Orbit* ok = modelspace->GetOrbit(k);
            Orbit* ol = modelspace->GetOrbit(l);
            double jk = ok->j2/2.;
            double jl = ol->j2/2.;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi->l+ok->l)%2;
            int Tz_cc = abs(oi->tz2+ok->tz2)/2;
            //int jmin = max(abs(int(ji-jj)),abs(int(jk-jl)));
            //int jmax = min(int(ji+jj),int(jk-jl));
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            double comm = 0;
//            cout << "Begin J loop at ch = " << ch << ", nkets(120) = " << modelspace->GetTwoBodyChannel(120).GetNumberKets() << endl;
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+tbc.J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k) indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l) indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
                 double me1 = N1[ch_cc](indx_jl,indx_ik);
                 double me2 = N2[ch_cc](indx_jl,indx_ik);
                 double me3 = N2[ch_cc](indx_ik,indx_jl);
                 double me4 = N1[ch_cc](indx_ik,indx_jl);
               comm -= (1) * phase * sixj * (me1-me2-me3+me4);
//               if (ch==modelspace->GetTwoBodyChannelIndex(1,0,1) and ((ibra == 0 and iket == 3) or (ibra==3 and iket ==0)) )
//               {
//                  cout << "Jprime = " << Jprime << " parity_cc = " << parity_cc << " Tz_cc = " << Tz_cc 
//                       << "  ijkl = " << i << "," << j << "," << k << "," << l
//                       << " me1 = " << me1 << " me2 = " << me2 
//                       << " me3 = " << me3 << " me4 = " << me4
//                       << " sixj("<<jk<<","<<jl<<","<<tbc.J<<","<<jj<<","<<ji<<","<<Jprime<<") = " << sixj << " phase = " << phase
//                       <<  "  comm(" << ibra << "," << iket << ") = " << comm
//                       << endl;
//               }
            }

            parity_cc = (oi->l+ol->l)%2;
            Tz_cc = abs(oi->tz2+ol->tz2)/2;
            jmin = max(abs(int(ji-jl)),abs(int(jk-jj)));
            jmax = min(int(ji+jl),int(jk+jj));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {

               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_il = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,l),max(i,l));
               int indx_jk = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,k),max(j,k));
               if (i>l) indx_il += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>k) indx_jk += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               double me2 = N2[ch_cc](indx_il,indx_jk);
               double me3 = N2[ch_cc](indx_jk,indx_il);
               double me4 = N1[ch_cc](indx_jk,indx_il);
               comm -= (1) * phase * sixj * (me1-me2-me3+me4);
//               if (ch==modelspace->GetTwoBodyChannelIndex(1,0,1) and ((ibra == 0 and iket == 3) or (ibra==3 and iket ==0)) )
//               {
//                  cout << "~Jprime = " << Jprime << " parity_cc = " << parity_cc << " Tz_cc = " << Tz_cc
//                       << "  ijkl = " << i << "," << j << "," << k << "," << l
//                       << " me1 = " << me1 << " me2 = " << me2 
//                       << " me3 = " << me3 << " me4 = " << me4
//                       << " sixj("<<jk<<","<<jl<<","<<tbc.J<<","<<ji<<","<<jj<<","<<Jprime<<") = " << sixj << " phase = " << phase
//                       <<  "  comm(" << ibra << "," << iket << ") = " << comm
//                       << endl;
//               }
            }
            comm /= sqrt((1.0+bra->delta_pq())*(1.0+ket->delta_pq()));
            opout.TwoBody[ch](ibra,iket) += comm;
            if (iket != ibra)
            {
               opout.TwoBody[ch](iket,ibra) += herm*comm;
            }
         }
      }
   }
}

