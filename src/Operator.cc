
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
  cout << "In Operator the modelspace has norbits = " << modelspace->GetNumberOrbits() << endl;
  hermitian = true;
  antihermitian = false;
  ZeroBody = 0;
  int nOneBody = modelspace->GetNumberOrbits();
  int nKets = modelspace->GetNumberKets();
  cout << "OneBody size = " << nOneBody*nOneBody << endl;
  cout << "TwoBody size = " << nKets*(nKets+1)/2 << endl;
  OneBody = arma::mat(nOneBody,nOneBody,arma::fill::zeros);
//  TwoBodyJmax = 0;
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



double Operator::GetTBME(int ch, int a, int b, int c, int d) const
{
   TwoBodyChannel &tbc = modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc.GetNumberKets() or ket_ind > tbc.GetNumberKets() )
   {
     return 0;
   }
   double phase = 1;
   Ket * bra = tbc.GetKet(bra_ind);
   Ket * ket = tbc.GetKet(ket_ind);
   if (a>b) phase *= bra->Phase(tbc.J) ;
   if (c>d) phase *= ket->Phase(tbc.J);
   if (a==b) phase *= sqrt(1.0+bra->delta_pq());
   if (c==d) phase *= sqrt(1.0+ket->delta_pq());
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
   int ch = (t+1)*2*JMAX + p*JMAX + j;
   return GetTBME(ch,a,b,c,d);
}

void Operator::SetTBME(int j, int p, int t, int a, int b, int c, int d, double tbme)
{
   int ch = (t+1)*2*JMAX + p*JMAX + j;
   SetTBME(ch,a,b,c,d,tbme);
}

double Operator::GetTBME(int j, int p, int t, Ket* bra, Ket* ket) const
{
   int ch = (t+1)*2*JMAX + p*JMAX + j;
   return GetTBME(ch,bra,ket);
}

void Operator::SetTBME(int j, int p, int t, Ket* bra, Ket* ket, double tbme)
{
   int ch = (t+1)*2*JMAX + p*JMAX + j;
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
   for (int &k : modelspace->hole) // C++11 syntax
   {
//      Orbit * orb = modelspace->GetOrbit(k);
//      if (orb->hvq > 0) continue;
      opNO.ZeroBody += (modelspace->GetOrbit(k)->j2+1) * OneBody(k,k);
   }
   cout << "after first loop, zero body part is" << opNO.ZeroBody << endl;


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
            for (int &h : modelspace->hole)  // C++11 syntax
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
   cout << "after second loop, zero body part is" << opNO.ZeroBody << endl;
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
   Operator OpNested = dOmega.Commutator(*this);
   Operator OpOut = *this + dOmega + OpNested*(1./2) + ( dOmega.Commutator(OpNested) + OpNested.Commutator(*this) )*(1./12);
   return OpOut;
}

// Frobenius norm of the operator
double Operator::Norm()
{
   double nrm = ZeroBody*ZeroBody;
   double n1 = arma::norm(OneBody,"fro");
//   double n1 = arma::accu(OneBody.t()*OneBody);
   nrm += n1*n1;
//   nrm += n1;
   cout << "Onebody =  " << endl; OneBody.print();
   cout << "One body norm = " << n1 << endl;
//   cout << "One body norm = " << sqrt(n1) << endl;
   for (int ch=0;ch<modelspace->GetNumberTwoBodyChannels();++ch)
   {
      double n2 = arma::norm(TwoBody[ch],"fro");
//      double n2 = arma::accu(TwoBody[ch].t()*TwoBody[ch]);
      nrm += n2*n2;
//      nrm += n2;
//      nrm += pow(arma::norm(TwoBody[ch],"fro"),2);
   }
   return sqrt(nrm);
}



Operator Operator::Commutator(Operator& opright)
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();
   if ( (IsHermitian() and opright.IsHermitian()) or (IsAntiHermitian() and opright.IsAntiHermitian()) ) out.SetAntiHermitian();
   else if ( (IsHermitian() and opright.IsAntiHermitian()) or (IsAntiHermitian() and opright.IsHermitian()) ) out.SetHermitian();
   else out.SetNonHermitian();

   out.ZeroBody  = comm110(opright) + comm220(opright) ;

   out.OneBody  = comm111(opright) + comm121(opright) + comm211(opright) + comm221(opright);

   for (int ch=0;ch<nChannels; ++ch)
   {
      out.TwoBody[ch] = comm122(opright,ch) + comm212(opright,ch) + comm222_ph(opright,ch) + comm222_pphh(opright,ch);
   }

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
   for ( int& a : modelspace->hole) // C++11 range-for syntax: loop over all elements of the vector <particle>
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
          for (int &a : modelspace->hole)  // C++11 syntax
          {
             Orbit *oa = modelspace->GetOrbit(a);
             for (int b=0;b<norbits;++b)
             {
                Orbit *ob = modelspace->GetOrbit(b);
                comm(i,j) += (ob->j2+1) *  OneBody(a,b) * opright.GetTBMEmonopole(b,i,a,j) ;
                comm(i,j) -= (oa->j2+1) *  OneBody(b,a) * opright.GetTBMEmonopole(a,i,b,j) ;
             }
          }
      }
   }
   return comm;
}


// Same as above, just reversed
arma::mat Operator::comm211(Operator& opright)
{
   return -opright.comm121(*this);
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
//                                     
//
arma::mat Operator::comm221(Operator& opright)
{

   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(norbits,norbits,arma::fill::zeros);
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();

      arma::mat Mpph = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
      arma::mat Mhhp = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);
//      cout << "Mpph:" << endl; Mpph.print();
//      cout << "Mhhp:" << endl; Mhhp.print();

      for (int i=0;i<norbits;++i)
      {
         Orbit *oi = modelspace->GetOrbit(i);
         double angmomfactor = (2*tbc.J+1.0)/(oi->j2 +1.0);
         for (int j=0;j<norbits;++j)
         {
            Orbit *oj = modelspace->GetOrbit(j);
            if (oi->j2 != oj->j2 or oi->l != oj->l or oi->tz2 != oj->tz2) continue;
            // Sum c over holes and include the nbar_a * nbar_b terms
            for (int &c : modelspace->hole)
            {
               int ibra = tbc.GetLocalIndex(c,i);
               int iket = tbc.GetLocalIndex(c,j);
               if (ibra<0 or iket < 0) continue;
               //if (ch==0)
               if (i==0 and j==0)
               {
                  //cout << "ibra = " << ibra << "  iket = " << iket << "  Mpph(ibra,iket) = " << Mpph(ibra,iket) << endl;
                  cout << "i,j,c = " << i << " " << j << " " << c << " J = " << tbc.J <<  "  Mpph(ibra,iket) = " << Mpph(ibra,iket) << endl;
               }
               comm(i,j) +=  angmomfactor * Mpph(ibra,iket);
            }
            // Sum c over particles and include the n_a * n_b terms
            for (int &c : modelspace->particles)
            {
               int ibra = tbc.GetLocalIndex(c,i);
               int iket = tbc.GetLocalIndex(c,j);
               if (ibra<0 or iket < 0) continue;
               //if (ch==0)
               if (i==0 and j==0)
               {
                  //cout << "ibra = " << ibra << "  iket = " << iket << "  Mhhp(ibra,iket) = " << Mhhp(ibra,iket) << endl;
                  cout << "i,j,c = " << i << " " << j << " " << c << " J = " << tbc.J << "  Mhhp(ibra,iket) = " << Mhhp(ibra,iket) << endl;
               }
               comm(i,j) +=  angmomfactor * Mhhp(ibra,iket);
            }

         } // for j
      } // for i
   } //for ch
   
   return 0.5*comm;
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
arma::mat Operator::comm122(Operator& opright, int ch )
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   int npq = tbc.GetNumberKets();
   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(npq,npq,arma::fill::zeros);
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
         for (int a=0;a<norbits;++a)
         {
            comm(ibra,iket) += OneBody(i,a) * opright.GetTBME(ch,a,j,k,l);
            comm(ibra,iket) += OneBody(j,a) * opright.GetTBME(ch,i,a,k,l);
            comm(ibra,iket) -= OneBody(a,k) * opright.GetTBME(ch,i,j,a,l);
            comm(ibra,iket) -= OneBody(a,l) * opright.GetTBME(ch,i,j,k,a);
         }
         comm(ibra,iket) /= sqrt( (1.0+bra->delta_pq())*(1.0+ket->delta_pq()) );
      }
   }
   return comm;
}

// Same as above, just reversed
arma::mat Operator::comm212(Operator& opright, int ch )
{
   return -opright.comm122(*this, ch);
}



//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_ph = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
// -- APPEARS TO AGREE WITH NATHAN'S RESULTS
//   although I'm not sure why there's a missing 1/2... 
arma::mat Operator::comm222_ph(Operator& opright, int ch )
{
   TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
   arma::mat P = tbc.Proj_pp - tbc.Proj_hh;
   arma::mat comm  =  TwoBody[ch] * P * opright.TwoBody[ch]  - opright.TwoBody[ch] * P * TwoBody[ch];

   //return 0.5*comm;
   return 1*comm;
}

//*****************************************************************************************
//
//  THIS IS THE BIG UGLY ONE.
//
//  |          |      |          |
//  |     __Y__|      |     __X__|
//  |    /\    |      |    /\    |
//  |   (  )   |  _   |   (  )   |
//  |____\/    |      |____\/    |
//  |  X       |      |  Y       |
//
arma::mat Operator::comm222_pphh(Operator& opright, int ch )
{
   arma::mat comm = TwoBody[ch];
   return comm;
}

