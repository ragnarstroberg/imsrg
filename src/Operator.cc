
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
   nChannels     = op.nChannels;
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
   TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
   int bra_ind = tbc.GetLocalIndex(min(a,b),max(a,b));
   int ket_ind = tbc.GetLocalIndex(min(c,d),max(c,d));
   if (bra_ind < 0 or ket_ind < 0 or bra_ind > tbc.GetNumberKets() or ket_ind > tbc.GetNumberKets() )
   {
     return 0;
   }
   double phase = 1;
   if (a>b) phase *= tbc.GetKet(bra_ind)->Phase(tbc.J);
   if (c>d) phase *= tbc.GetKet(ket_ind)->Phase(tbc.J);
   return phase * TwoBody[ch](bra_ind, ket_ind);
}

void Operator::SetTBME(int ch, int a, int b, int c, int d, double tbme)
{
   TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
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
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
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
               //double tbme = (2*J+1.0)/(orba->j2+1)*tbc.GetTBME(bra,ket);
               double tbme = (2*J+1.0)/(orba->j2+1) * GetTBME(ch,bra,ket);

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
            //opNO.ZeroBody += tbc->TBME(k,k) * (2*tbc->J+1) / (1.0+ket->delta_pq());  // <ab|V|ab>  (a,b in core)
            opNO.ZeroBody += TwoBody[ch](k,k) * (2*tbc.J+1) / (1.0+ket->delta_pq());  // <ab|V|ab>  (a,b in core)
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
       //GetTwoBodyChannel(ch)->TBME.zeros();
   }
}



//Operator Operator::Commutator(const Operator& opright)
Operator Operator::Commutator(Operator& opright)
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   out.ZeroBody  = comm110(opright) + comm220(opright) ;

   out.OneBody  = comm111(opright) + comm211(opright) - opright.comm211(*this)  + comm221(opright);

   for (int ch=0;ch<nChannels; ++ch)
   {
      out.TwoBody[ch] = comm212(opright,ch) - opright.comm212(*this,ch) + comm222_ph(opright,ch) + comm222_pphh(opright,ch);
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
double Operator::comm110(Operator& opright)
{
  if (IsHermitian() and opright.IsHermitian()) return 0; // I think this is the case
  if (IsAntiHermitian() and opright.IsAntiHermitian()) return 0; // I think this is the case

   double comm = 0;
   int norbits = modelspace->GetNumberOrbits();
   arma::mat xyyx = OneBody*opright.OneBody - opright.OneBody*OneBody;
   for ( int& a : modelspace->particles) // C++11 range-for syntax: loop over all elements of the vector <particle>
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
//                       = 1/2 Sum_J (2J+1) Sum_ab  P_hh*(X*P_pp*Y)_abab  
//
//
double Operator::comm220( Operator& opright)
{
   double comm = 0;
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      //comm += 0.5*(2*tbleft->J+1) * arma::accu( tbleft->Proj_hh * tbleft->TBME * tbleft->Proj_pp * tbright->TBME);
      //comm += 0.5*(2*tbc.J+1) * arma::accu( tbc.Proj_hh * TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch]);
      comm += (2*tbc.J+1) * arma::accu( tbc.Proj_hh * TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch]);
   }
   return comm;
}

//*****************************************************************************************
//
//        |____. Y          |___.X
//        |        _        |
//  X .___|            Y.___|              [X1,Y1](1)  =  XY - YX
//        |                 |

arma::mat Operator::comm111(Operator & opright)
{
   return OneBody*opright.OneBody - opright.OneBody*OneBody;
}

//*****************************************************************************************
// 
//      i |              i |
//        |    ___.Y       |__X__
//        |___(_)    _     |   (_)__.     [X2,Y1](1)  =  1/(2j_i+1) sum_ab(n_a-n_b)y_ab 
//      j | X            j |        Y            * sum_J (2J+1) x_biaj^(J)  
//                                             
//                                               = 1/(2j+1) sum_a n_a sum_J (2J+1)
//                                                  * sum_b y_ab x_biaj - yba x_aibj
//
arma::mat Operator::comm211(Operator& opright)
{
//  cout << "In comm211" << endl;
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
             for (int ch=0;ch<nChannels;++ch)
             {
                TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
                if (tbc.GetNumberKets() < 1) continue;
                int ind_ai = tbc.GetLocalIndex(min(a,i),max(a,i));
                int ind_aj = tbc.GetLocalIndex(min(a,j),max(a,j));
                if ( ind_ai < 0 and ind_aj < 0) continue; 
                double commijJ = 0;
                for (int b=0;b<norbits; ++b)
                {
                   commijJ += opright.OneBody(a,b) * GetTBME(ch,b,i,a,j) - opright.OneBody(b,a) * GetTBME(ch,a,i,b,j);
                }
                comm(i,j) += commijJ*(2*tbc.J+1);
             }
          }
          comm(i,j) /= (modelspace->GetOrbit(i)->j2 + 1);
      }
   }
   return comm;
}

//*****************************************************************************************
//
//      i |             i |
//        |__Y___         |__X___
//        |____(_)  _     |____(_)     [X2,Y2](1)  =  1/(2(2j_i+1)) sum_J (2J+1) 
//      j | X           j |  Y            * sum_abc (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
//                                         * (x_ciab y_abcj - y_ciab xabcj)
//
//                                   = 1/(2(2j+1)) sum_J (2J+1)
//                                     *  sum_c ( Pp*X*Phh*Y - Pp*Y*Phh*X)  - (Ph*X*Ppp*Y - Ph*Y*Ppp*X)_cicj
//
arma::mat Operator::comm221(Operator& opright)
{

   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(norbits,norbits,arma::fill::zeros);
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();

      arma::mat Mpph = (TwoBody[ch] * tbc.Proj_pp * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_pp * TwoBody[ch]);
      arma::mat Mhhp = (TwoBody[ch] * tbc.Proj_hh * opright.TwoBody[ch] - opright.TwoBody[ch] * tbc.Proj_hh * TwoBody[ch]);

      for (int i=0;i<norbits;++i)
      {
         Orbit *oi = modelspace->GetOrbit(i);
         double angmomfactor = (2*tbc.J+1)/(2*(oi->j2 +1));
         for (int j=0;j<norbits;++j)
         {
            Orbit *oj = modelspace->GetOrbit(j);
            if (oi->j2 != oj->j2 or oi->l != oj->l or oi->tz2 != oj->tz2) continue;
            for (int &c : modelspace->hole)
            {
               int ibra = tbc.GetLocalIndex(c,i);
               int iket = tbc.GetLocalIndex(c,j);
               if (ibra<0 or iket < 0) continue;
               comm(i,j) +=  angmomfactor * Mpph(ibra,iket);
            }
            for (int &c : modelspace->particles)
            {
               int ibra = tbc.GetLocalIndex(c,i);
               int iket = tbc.GetLocalIndex(c,j);
               if (ibra<0 or iket < 0) continue;
               comm(i,j) +=  angmomfactor * Mhhp(ibra,iket);
            }
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
arma::mat Operator::comm212(Operator& opright, int ch )
{
   TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
   int npq = tbc.GetNumberKets();
   int norbits = modelspace->GetNumberOrbits();
   arma::mat comm = arma::mat(npq,npq,arma::fill::zeros);
   for (int a=0;a<norbits;++a)
   {
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
            comm(ibra,iket) += opright.OneBody(i,a) * GetTBME(ch,a,j,k,l);
            comm(ibra,iket) += opright.OneBody(j,a) * GetTBME(ch,i,a,k,l);
            comm(ibra,iket) -= opright.OneBody(a,k) * GetTBME(ch,i,j,a,l);
            comm(ibra,iket) -= opright.OneBody(a,l) * GetTBME(ch,i,j,k,a);
         }
      }
   }
   return comm;
}

//*****************************************************************************************
//
//  |     |      |     |   
//  |__Y__|      |__x__|   [X2,Y2](2)_ph = 1/2 sum_ab (X_ijab Y_abkl - Y_ijab X_abkl)(1 - n_a - n_b)
//  |     |  _   |     |                = 1/2 [ X*(P_pp-P_hh)*Y - Y*(P_pp-P_hh)*X ]
//  |__X__|      |__Y__|   
//  |     |      |     |   
//
arma::mat Operator::comm222_ph(Operator& opright, int ch )
{
   TwoBodyChannel tbc = modelspace->GetTwoBodyChannel(ch);
   arma::mat comm  = (        TwoBody[ch] * (tbc.Proj_pp - tbc.Proj_hh) * opright.TwoBody[ch])
                   - (opright.TwoBody[ch] * (tbc.Proj_pp - tbc.Proj_hh) *         TwoBody[ch]);

   return 0.5*comm;
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

