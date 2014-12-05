
#include "Operator.hh"
#include "AngMom.hh"
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

double  Operator::bch_transform_threshold = 1e-6;
double  Operator::bch_product_threshold = 1e-4;

Operator::Operator(ModelSpace* ms) // Create a zero-valued operator in a given model space
{
  modelspace = ms;
  hermitian = true;
  antihermitian = false;
//  cross_coupled = false;
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

}

void Operator::Copy(const Operator& op)
{
   modelspace    = op.modelspace;
   nChannels     = modelspace->GetNumberTwoBodyChannels();
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
//   cross_coupled = op.cross_coupled;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBody[ch] = op.TwoBody[ch];
   }
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
   TwoBodyChannel& tbc =  modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc.GetNumberKets() or ket_ind > tbc.GetNumberKets() )
     return 0;
   Ket * bra = tbc.GetKet(bra_ind);
   Ket * ket = tbc.GetKet(ket_ind);

   double phase = 1;
   if (a>b) phase *= bra->Phase(tbc.J) ;
   if (c>d) phase *= ket->Phase(tbc.J);
   if (a==b) phase *= sqrt(2.);
   if (c==d) phase *= sqrt(2.);
   return phase * TwoBody[ch](bra_ind, ket_ind);
}


void Operator::SetTBME(int ch, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
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
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,a,b,c,d);
}

void Operator::SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   SetTBME(ch,a,b,c,d,tbme);
}

double Operator::GetTBME(int j, int p, int t, Ket* bra, Ket* ket) const
{
   int ch = modelspace->GetTwoBodyChannelIndex(j,p,t);
   return GetTBME(ch,bra,ket);
}

void Operator::SetTBME(int j, int p, int t, Ket* bra, Ket* ket, double tbme)
{
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


// divide operator by a scalar
Operator& Operator::operator/=(const double rhs)
{
   return *this *=(1.0/rhs);
}

Operator Operator::operator/(const double rhs)
{
   Operator opout = Operator(*this);
   opout /= rhs;
   return opout;
}

// Add operators
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

// Subtract operators
Operator& Operator::operator-=(const Operator& rhs)
{
   ZeroBody -= rhs.ZeroBody;
   OneBody -= rhs.OneBody;
   for (int ch=0; ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      TwoBody[ch] -= rhs.TwoBody[ch];
   }
   return *this;
}

Operator Operator::operator-(const Operator& rhs)
{
   return ( Operator(*this) -= rhs );
}

//********************************************************

Operator Operator::DoNormalOrdering()
{
   Operator opNO = Operator(*this);

   // Trivial parts
   opNO.ZeroBody = ZeroBody;
   opNO.OneBody = OneBody;
   for (int ch=0;ch<nChannels;++ch)
   {
      opNO.TwoBody[ch] = TwoBody[ch];
   }



   for (int &k : modelspace->holes) // loop over hole orbits
   {
      opNO.ZeroBody += (modelspace->GetOrbit(k)->j2+1) * OneBody(k,k);
   }

   int norbits = modelspace->GetNumberOrbits();

   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;

      // Zero body part
      for (int& iket : tbc.KetIndex_hh) // loop over hole-hole kets in this channel
      {
        opNO.ZeroBody += TwoBody[ch](iket,iket) * (2*J+1);  // <ab|V|ab>  (a,b in core)
      }

      // One body part
      int ibra,iket;
      for (int a=0;a<norbits;++a)
      {
         Orbit *orba = modelspace->GetOrbit(a);
         int bstart = IsNonHermitian() ? 0 : a; // If it's neither hermitian or anti, we need to do the full sum
         for (int b=bstart;b<norbits;++b)
         {
            Orbit *orbb = modelspace->GetOrbit(b);
            if (orbb->j2 != orba->j2 or orbb->tz2 != orba->tz2 or orbb->l != orba->l) continue;
            for (int &h : modelspace->holes)  // C++11 syntax
            {
               if ( (ibra = tbc.GetLocalIndex(min(a,h),max(a,h)))<0) continue;
               if ( (iket = tbc.GetLocalIndex(min(b,h),max(b,h)))<0) continue;
               Ket * bra = tbc.GetKet(ibra);
               Ket * ket = tbc.GetKet(iket);
               double tbme = (2*J+1.0)/(orba->j2+1)  * GetTBME(ch,bra,ket);
               int phase = 1;
               if (a>h) phase *= bra->Phase(J);
               if (b>h) phase *= ket->Phase(J);

               opNO.OneBody(a,b) += phase * tbme;  // <ah|V|bh>
            }
            if (hermitian) opNO.OneBody(b,a) = opNO.OneBody(a,b);
            if (antihermitian) opNO.OneBody(b,a) = -opNO.OneBody(a,b);
         }
      }


   } // loop over channels
   return opNO;
}



void Operator::EraseTwoBody()
{
   for (int ch=0;ch<nChannels;++ch)
   {
       TwoBody[ch].zeros();
   }
}

void Operator::CalculateKineticEnergy()
{
   OneBody.zeros();
   int norbits = modelspace->GetNumberOrbits();
   int A = modelspace->GetTargetMass();
   float hw = modelspace->GetHbarOmega();
   for (int a=0;a<norbits;++a)
   {
      Orbit * orba = modelspace->GetOrbit(a);
      OneBody(a,a) = 0.5 * hw * (2*orba->n + orba->l +3./2); 
      for (int b=a+1;b<norbits;++b)  // make this better once OneBodyChannel is implemented
      {
         Orbit * orbb = modelspace->GetOrbit(b);
         if (orba->l == orbb->l and orba->j2 == orbb->j2 and orba->tz2 == orbb->tz2)
         {
            if (orba->n == orbb->n+1)
               OneBody(a,b) = 0.5 * hw * sqrt( (orba->n)*(orba->n + orba->l +1./2));
            else if (orba->n == orbb->n-1)
               OneBody(a,b) = 0.5 * hw * sqrt( (orbb->n)*(orbb->n + orbb->l +1./2));
            OneBody(b,a) = OneBody(a,b);
         }
      }
   }
//   OneBody *= (1-1./A);

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
   #pragma omp parallel for
   for (int ch_cc=0; ch_cc<nChannels; ++ch_cc)
   {
      TwoBodyChannel& tbc_cc = modelspace->GetTwoBodyChannel_CC(ch_cc);
      int nKets_cc = tbc_cc.GetNumberKets();
      int nph_kets = tbc_cc.KetIndex_ph.size();
      int J_cc = tbc_cc.J;

//   These matrices don't actually need to be square, since we only care about
//   the particle-hole kets, which we sum over via matrix multiplication:
//   [ ]                  [     ]
//   |L|  x  [  R   ]  =  |  M  |
//   [ ]                  [     ]
//
      TwoBody_CC_left[ch_cc]  = arma::mat(2*nKets_cc, nph_kets,   arma::fill::zeros);
      TwoBody_CC_right[ch_cc] = arma::mat(nph_kets,   2*nKets_cc, arma::fill::zeros);

      // loop over cross-coupled ph bras <ac| in this channel
      for (int i_ph=0; i_ph<nph_kets; ++i_ph)
      {
         Ket * bra_cc = tbc_cc.GetKet( tbc_cc.KetIndex_ph[i_ph] );
         int a = bra_cc->p;
         int c = bra_cc->q;
         Orbit * oa = modelspace->GetOrbit(a);
         Orbit * oc = modelspace->GetOrbit(c);
         double ja = oa->j2/2.0;
         double jc = oc->j2/2.0;

         // loop over cross-coupled kets |bd> in this channel
         // we go to 2*nKets to include |bd> and |db>
         for (int iket_cc=0; iket_cc<2*nKets_cc; ++iket_cc)
         {
            Ket * ket_cc = tbc_cc.GetKet(iket_cc%nKets_cc);
            int b = iket_cc < nKets_cc ? ket_cc->p : ket_cc->q;
            int d = iket_cc < nKets_cc ? ket_cc->q : ket_cc->p;
            Orbit * ob = modelspace->GetOrbit(b);
            Orbit * od = modelspace->GetOrbit(d);
            double jb = ob->j2/2.0;
            double jd = od->j2/2.0;

            int phase_ad = modelspace->phase(ja+jd);

            // Get Tz,parity and range of J for <ab || cd > coupling
            int Tz_std = (oa->tz2 + ob->tz2)/2;
            int parity_std = (oa->l + ob->l)%2;
            int jmin = max(abs(ja-jb),abs(jc-jd));
            int jmax = min(ja+jb,jc+jd);
            double sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jc,J_cc,jd,jb,J_std);
               if (sixj == 0) continue;
               int phase = modelspace->phase(J_std);
               double tbme = GetTBME(J_std,parity_std,Tz_std,a,b,c,d);
               sm += (2*J_std+1) * phase * sixj * tbme ; 
            }
            TwoBody_CC_left[ch_cc](iket_cc,i_ph) = sm * phase_ad;


            // Get Tz,parity and range of J for <cb || ad > coupling
            Tz_std = (oa->tz2 + od->tz2)/2;
            parity_std = (oa->l + od->l)%2;
            jmin = max(abs(jc-jb),abs(ja-jd));
            jmax = min(jc+jb,ja+jd);
            sm = 0;
            for (int J_std=jmin; J_std<=jmax; ++J_std)
            {
               double sixj = modelspace->GetSixJ(ja,jc,J_cc,jb,jd,J_std);
               if (sixj == 0) continue;
               int phase = modelspace->phase(J_std);
               double tbme = GetTBME(J_std,parity_std,Tz_std,c,b,a,d);
               sm += (2*J_std+1) * phase * sixj * tbme ;
            }
            TwoBody_CC_right[ch_cc](i_ph,iket_cc) = - sm * phase_ad;

         }
      }
   }
}


//*****************************************************************************************
//  OpOut = exp(Omega) Op exp(-Omega)
Operator Operator::BCH_Transform( Operator &Omega)
{
//   double bch_transform_threshold = 1e-6;
   int max_iter = 12;
   double nx = Norm();
   double ny = Omega.Norm();
   Operator OpOut = *this;
   Operator OpNested = *this;
   for (int i=1; i<max_iter; ++i)
   {
      OpNested = Omega.Commutator(OpNested) / i;
      
      OpOut += OpNested;

      if (OpNested.Norm() < (nx+ny)*bch_transform_threshold) return OpOut;

   }
   cout << "Warning: BCH_Transform didn't coverge after "<< max_iter << " nested commutators" << endl;
   return OpOut;
}


//*****************************************************************************************
// Baker-Campbell-Hausdorff formula
//  returns Z, where
//  exp(Z) = exp(X) * exp(Y).
//  Z = X + Y + 1/2[X, Y]
//     + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]]
//     - 1/24 [Y,[X,[X,Y]]]
//     - 1/720 [Y,[Y,[Y,[Y,X]]]] - 1/720 [X,[X,[X,[X,Y]]]]
//     + ...
Operator Operator::BCH_Product( Operator &Y)
{
   Operator& X = *this;
//   double bch_product_threshold = 1e-4;

   Operator Z = X + Y; 

   Operator XY = X.Commutator(Y);
   Z += XY*(1./2);    // [X,Y]
   double nx = X.Norm();
   double ny = Y.Norm();
   double nc1 = XY.Norm();

   if ( nc1/2 < (nx+ny)*bch_product_threshold ) return Z;

   Operator YYX = XY.Commutator(Y); // [[X,Y],Y] = [Y,[Y,X]]
   double nc2 = YYX.Norm();
   Z += YYX * (1/12.);      // [Y,[Y,X]]

   if ( nc2/12 < (nx+ny)*bch_product_threshold ) return Z;

   Operator XXY = X.Commutator(XY); // [X,[X,Y]]
   double nc3 = XXY.Norm();
   Z += XXY * (1/12.);      // [X,[X,Y]]

   if ( nc3/12 < (nx+ny)*bch_product_threshold ) return Z;

   cout << "Warning: BCH product expansion not converged after 3 nested commutators!" << endl;

   Operator YXXY = Y.Commutator(XXY); // [Y,[X,[X,Y]]]
   double nc4 = YXXY.Norm();
   Z += YXXY*(-1./24);

   Operator YYYX = Y.Commutator(YYX) ;
   Operator YYYYX = Y.Commutator(YYYX) ;
   double nc5 = YYYYX.Norm();
   Operator XXXY =  X.Commutator(XXY) ; // [X,[X,[X,[X,Y]]]]
   Operator XXXXY = X.Commutator(XXXY) ; // [X,[X,[X,[X,Y]]]]
   double nc6 = XXXXY.Norm();
   Z += (YYYYX + XXXXY)*(-1./720);

   if ( nc6/720. < (nx+ny)*bch_product_threshold ) return Z;
   cout << "Warning: BCH product expansion not converged after 5 nested commutators!" << endl;

   return Z;
}

// Frobenius norm of the operator
double Operator::Norm()
{
   double nrm = 0;
   double n1 = OneBodyNorm();
   double n2 = TwoBodyNorm();
   return sqrt(n1*n1+n2*n2);
}

double Operator::OneBodyNorm()
{
   return arma::norm(OneBody,"fro");
}

double Operator::TwoBodyNorm()
{
   double nrm = 0;
   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      double n2 = arma::norm(TwoBody[ch],"fro");
      nrm += n2*n2;
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

   out.OneBody  = comm111(opright) + comm121(opright);// + comm221(opright);

   UpdateCrossCoupled();
   opright.UpdateCrossCoupled();

   comm122(opright,out);
//   comm222_pp_hh(opright,out);
   comm222_pp_hh_221(opright, out);
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
   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      double c = 2 * (2*tbc.J+1) * arma::trace( tbc.Proj_hh * TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] );
      #pragma omp critical
      comm += c;

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
   #pragma omp parallel for
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

   #pragma omp parallel for schedule(dynamic,5)
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
            #pragma omp critical
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
// Right now, this is the slowest one...
void Operator::comm122(Operator& opright, Operator& opout )
{
   int herm = opout.IsHermitian() ? 1 : -1;

   #pragma omp parallel for schedule(dynamic,5)
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
         int indx_ij = ibra;
         double pre_ij = i==j ? sqrt(2) : 1;
         for (int iket=ibra;iket<npq; ++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            int k = ket->p;
            int l = ket->q;
            int indx_kl = iket;
            double pre_kl = k==l ? sqrt(2) : 1;

            double cij = 0;
            double ckl = 0;
            double cijkl = 0;
            for (int a=0;a<norbits;++a)
            {
               int indx_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
               int indx_ia = tbc.GetLocalIndex(min(i,a),max(i,a));
               int indx_al = tbc.GetLocalIndex(min(a,l),max(a,l));
               int indx_ka = tbc.GetLocalIndex(min(k,a),max(k,a));

               if (indx_aj >= 0)
               {
                  double pre_aj = a>j ? tbc.GetKet(indx_aj)->Phase(tbc.J) : 1;
                  if (a==j) pre_aj *= sqrt(2);
                  ckl += pre_aj  * OneBody(i,a) * opright.TwoBody[ch](indx_aj,indx_kl);
                  ckl -= pre_aj  * opright.OneBody(i,a) * TwoBody[ch](indx_aj,indx_kl);
               }

               if (indx_ia >= 0)
               {
                  double pre_ia = i>a ? tbc.GetKet(indx_ia)->Phase(tbc.J) : 1;
                  if (i==a) pre_ia *= sqrt(2);
                  ckl += pre_ia * OneBody(j,a) * opright.TwoBody[ch](indx_ia,indx_kl);
                  ckl -= pre_ia * opright.OneBody(j,a) * TwoBody[ch](indx_ia,indx_kl);
               }

               if (indx_al >= 0)
               {
                  double pre_al = a>l ? tbc.GetKet(indx_al)->Phase(tbc.J) : 1;
                  if (a==l) pre_al *= sqrt(2);
                  cij -= pre_al * OneBody(a,k) * opright.TwoBody[ch](indx_ij,indx_al);
                  cij += pre_al * opright.OneBody(a,k) * TwoBody[ch](indx_ij,indx_al);
               }

               if (indx_ka >= 0)
               {
                  double pre_ka = k>a ? tbc.GetKet(indx_ka)->Phase(tbc.J) : 1;
                  if (k==a) pre_ka *= sqrt(2);
                  cij -= pre_ka * OneBody(a,l) * opright.TwoBody[ch](indx_ij,indx_ka);
                  cij += pre_ka * opright.OneBody(a,l) * TwoBody[ch](indx_ij,indx_ka);
               }

            }
            cijkl = (ckl*pre_kl + cij*pre_ij) / sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()) );
            #pragma omp critical
            {
            opout.TwoBody[ch](ibra,iket) += cijkl;
            if (ibra != iket)
               opout.TwoBody[ch](iket,ibra) += herm * cijkl;
            }
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
   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0; ch<nChannels; ++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      arma::mat Mpp = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
      arma::mat Mhh = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);
      opout.TwoBody[ch] += Mpp - Mhh;
   }
}







// Since comm222_pp_hh and comm211 both require the construction of 
// the intermediate matrices Mpp and Mhh, we can combine them and
// only calculate the intermediates once.
void Operator::comm222_pp_hh_221(Operator& opright, Operator& opout )
{

   int herm = opout.IsHermitian() ? 1 : -1;
   int norbits = modelspace->GetNumberOrbits();
   Operator Mpp = opright;
   Operator Mhh = opright;

   #pragma omp parallel for schedule(dynamic,5)
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      
//      Mpp.TwoBody[ch] = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
//      Mhh.TwoBody[ch] = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);
      Mpp.TwoBody[ch] = TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch];
      Mhh.TwoBody[ch] = TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch];
      if (opout.IsHermitian())
      {
         Mpp.TwoBody[ch] += Mpp.TwoBody[ch].t();
         Mhh.TwoBody[ch] += Mhh.TwoBody[ch].t();
      }
      else if (opout.IsAntiHermitian())
      {
         Mpp.TwoBody[ch] -= Mpp.TwoBody[ch].t();
         Mhh.TwoBody[ch] -= Mhh.TwoBody[ch].t();
      }
      else
      {
        Mpp.TwoBody[ch] -= opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch];
        Mhh.TwoBody[ch] -= opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch];
      }


      // The two body part
      opout.TwoBody[ch] += Mpp.TwoBody[ch] - Mhh.TwoBody[ch];

//      #pragma omp parallel for schedule(dynamic,3)
      for (int i=0;i<norbits;++i)
      {
         Orbit *oi = modelspace->GetOrbit(i);
         int jmin = opout.IsNonHermitian() ? 0 : i;
         //for (int j=0;j<norbits;++j)
         for (int j=jmin;j<norbits;++j)
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
            //comm(i,j) += (2*tbc.J+1.0)/(oi->j2 +1.0) * cijJ;
            cijJ *= (2*tbc.J+1.0)/(oi->j2+1.0);
            opout.OneBody(i,j) += cijJ;
            if (! opout.IsNonHermitian() and i!=j)
            {
               opout.OneBody(j,i) += herm * cijJ;
            }
         } // for j
      } // for i
   } //for ch
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
//                      double nj1 = modelspace->GetNineJ( ja, jl, J1,   ji, J, jj,   J2, jk, jb );
//                      double nj2 = modelspace->GetNineJ( ja, jl, J1,   jj, J, ji,   J2, jk, jb );
//                      double nj3 = modelspace->GetNineJ( jb, jl, J1,   ji, J, jj,   J2, jk, ja );
//                      double nj4 = modelspace->GetNineJ( jb, jl, J1,   jj, J, ji,   J2, jk, ja );
                        double nj1 = AngMom::NineJ( ja, jl, J1,   ji, J, jj,   J2, jk, jb );
                        double nj2 = AngMom::NineJ( ja, jl, J1,   jj, J, ji,   J2, jk, jb );
                        double nj3 = AngMom::NineJ( jb, jl, J1,   ji, J, jj,   J2, jk, ja );
                        double nj4 = AngMom::NineJ( jb, jl, J1,   jj, J, ji,   J2, jk, ja );

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

   int herm = opout.IsHermitian() ? 1 : -1;

   // Construct the intermediate matrices N1 and N2
   arma::mat N1[nChannels];
   arma::mat N2[nChannels];
 // probably better not to parallelize, since the matrix muliplication can use OMP
   for (int ch=0;ch<nChannels;++ch)
   {
      N1[ch] =         TwoBody_CC_left[ch] *    opright.TwoBody_CC_right[ch];
      N2[ch] = opright.TwoBody_CC_left[ch] *            TwoBody_CC_right[ch];
   }

   // Now evaluate the commutator for each channel (standard coupling)
   #pragma omp parallel for schedule(dynamic,5)
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

         for (int iket=ibra; iket<nKets; ++iket)
         {
            Ket * ket = tbc.GetKet(iket);
            int k = ket->p;
            int l = ket->q;
            Orbit* ok = modelspace->GetOrbit(k);
            Orbit* ol = modelspace->GetOrbit(l);
            double jk = ok->j2/2.;
            double jl = ol->j2/2.;


            double comm = 0;

            // now loop over the cross coupled TBME's
            int parity_cc = (oi->l+ok->l)%2;
            int Tz_cc = abs(oi->tz2+ok->tz2)/2;
            int jmin = max(abs(int(jj-jl)),abs(int(jk-ji)));
            int jmax = min(int(jj+jl),int(jk+ji));
            for (int Jprime=jmin; Jprime<=jmax; ++Jprime)
            {
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,jj,ji,Jprime);
               int phase = modelspace->phase(ji+jl+tbc.J);
               int ch_cc = modelspace->GetTwoBodyChannelIndex(Jprime,parity_cc,Tz_cc);
               int indx_ik = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(i,k),max(i,k));
               int indx_jl = modelspace->GetTwoBodyChannel_CC(ch_cc).GetLocalIndex(min(j,l),max(j,l));
               if (i>k)
                 indx_ik += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>l)
                 indx_jl += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double me1 = N1[ch_cc](indx_jl,indx_ik);
               double me2 = N2[ch_cc](indx_jl,indx_ik);
               double me3 = N2[ch_cc](indx_ik,indx_jl);
               double me4 = N1[ch_cc](indx_ik,indx_jl);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
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
               if (i>l)
                 indx_il += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               if (j>k)
                 indx_jk += modelspace->GetTwoBodyChannel_CC(ch_cc).GetNumberKets();
               double sixj = modelspace->GetSixJ(jk,jl,tbc.J,ji,jj,Jprime);
               int phase = modelspace->phase(ji+jl);
               double me1 = N1[ch_cc](indx_il,indx_jk);
               double me2 = N2[ch_cc](indx_il,indx_jk);
               double me3 = N2[ch_cc](indx_jk,indx_il);
               double me4 = N1[ch_cc](indx_jk,indx_il);
               comm -= (2*Jprime+1) * phase * sixj * (me1-me2-me3+me4);
            }

            comm /= sqrt((1.0+bra->delta_pq())*(1.0+ket->delta_pq()));
            #pragma omp critical
            {
            opout.TwoBody[ch](ibra,iket) += comm;
            if (iket != ibra)
            {
               opout.TwoBody[ch](iket,ibra) += herm*comm;
            }
            }
         }
      }
   }
}

