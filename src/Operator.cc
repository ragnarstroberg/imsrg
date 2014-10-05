
#include "Operator.hh"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

//TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace * ms)

TwoBodyChannel::TwoBodyChannel()
{
}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, Operator * op)
//TwoBodyChannel::TwoBodyChannel(int j, int p, int t, bool Herm=false, bool AHerm=false)
{
   J = j;
   parity = p;
   Tz = t;
   hermitian = op->IsHermitian();
   antihermitian = op->IsAntiHermitian();
   modelspace = op->GetModelSpace();
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1);
   for (int i=0;i<nk;i++)
   {
      Ket *ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      //if (CheckChannel_pq(ket->p, ket->q))
      {
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         NumberKets++;
      }
   }
   TBME = arma::mat(NumberKets, NumberKets, arma::fill::zeros);

}


TwoBodyChannel::TwoBodyChannel(const TwoBodyChannel& rhs)
{
   J                 = rhs.J;
   parity            = rhs.parity;
   Tz                = rhs.Tz;
   hermitian         = rhs.hermitian;
   antihermitian     = rhs.antihermitian;
   modelspace        = rhs.modelspace;
   NumberKets = rhs.NumberKets;
   TBME              = rhs.TBME;
   KetMap            = rhs.KetMap;
   KetList           = rhs.KetList;
}


TwoBodyChannel& TwoBodyChannel::operator=(const TwoBodyChannel& rhs)
{
   J                 = rhs.J;
   parity            = rhs.parity;
   Tz                = rhs.Tz;
   hermitian         = rhs.hermitian;
   antihermitian     = rhs.antihermitian;
   modelspace        = rhs.modelspace;
   NumberKets = rhs.NumberKets;
   TBME              = rhs.TBME;
   KetMap            = rhs.KetMap;
   KetList           = rhs.KetList;
   return *this;
}


TwoBodyChannel& TwoBodyChannel::operator+=(const TwoBodyChannel& rhs)
{
   if (J == rhs.J and parity == rhs.parity and Tz == rhs.Tz)
   {
      hermitian         = (hermitian and rhs.hermitian);
      antihermitian     = (antihermitian and rhs.antihermitian);
      TBME              += rhs.TBME;
   }
   return *this;
}

TwoBodyChannel& TwoBodyChannel::operator-=(const TwoBodyChannel& rhs)
{
   if (J == rhs.J and parity == rhs.parity and Tz == rhs.Tz)
   {
      hermitian         = (hermitian and rhs.hermitian);
      antihermitian     = (antihermitian and rhs.antihermitian);
      TBME              -= rhs.TBME;
   }
   return *this;
}

TwoBodyChannel& TwoBodyChannel::operator+(const TwoBodyChannel& rhs)
{
   return TwoBodyChannel( *this) += rhs;
}

TwoBodyChannel& TwoBodyChannel::operator-(const TwoBodyChannel& rhs)
{
   return TwoBodyChannel( *this) -= rhs;
}



float TwoBodyChannel::GetTBME(int bra, int ket) const 
{
   //int bra_ind = Getpqconfig(bra);
   //int ket_ind = Getpqconfig(ket);
   int bra_ind = GetLocalIndex(bra);
   int ket_ind = GetLocalIndex(ket);
   if (bra_ind < 0 or ket_ind < 0) return 0;
   if (bra_ind > NumberKets or ket_ind > NumberKets) return 0;
   return TBME(bra_ind, ket_ind);
}

float TwoBodyChannel::GetTBME(Ket *bra, Ket *ket) const
{
//   int bra_ind = Getpqconfig(bra->p,bra->q);
//   int ket_ind = Getpqconfig(ket->p,ket->q);
   int bra_ind = GetLocalIndex(bra->p,bra->q);
   int ket_ind = GetLocalIndex(ket->p,ket->q);
   if (bra_ind < 0 or ket_ind < 0) return 0;
   if (bra_ind > NumberKets or ket_ind > NumberKets) return 0;
   return TBME(bra_ind, ket_ind);
}

void TwoBodyChannel::SetTBME(int bra, int ket, float tbme)
{
   int bra_ind = GetLocalIndex(bra);
   int ket_ind = GetLocalIndex(ket);
   if ( (bra_ind < 0) or (ket_ind < 0) or
       (bra_ind > NumberKets) or (ket_ind > NumberKets) )
   {
     cout << "Bad assignment of tbme!"
         << " < " << bra << " | V | "
         << ket << " > "
         << "J,parity,Tz = " << J << "," << parity << "," << Tz
         << endl;
     return;
   }
   TBME(bra_ind, ket_ind) = tbme;
}

void TwoBodyChannel::SetTBME(Ket *bra, Ket *ket, float tbme)
{
   int bra_ind = GetLocalIndex(bra->p,bra->q);
   int ket_ind = GetLocalIndex(ket->p,ket->q);
   if ( (bra_ind < 0) or (ket_ind < 0) or
       (bra_ind > NumberKets) or (ket_ind > NumberKets) )
   {
            int ja,la,ta,jb,lb,tb;
   int a = bra->p;
   int b = bra->q;
   ja = modelspace->GetOrbit(a)->j2;
   la = modelspace->GetOrbit(a)->l;
   ta = modelspace->GetOrbit(a)->tz2;
   jb = modelspace->GetOrbit(b)->j2;
   lb = modelspace->GetOrbit(b)->l;
   tb = modelspace->GetOrbit(b)->tz2;
     cout << "Bad assignment of tbme!"
         << " < " << bra->p << " " << bra->q << " | V | "
         << ket->p << " " << ket->q << " > "
         << "J,parity,Tz = " << J << "," << parity << "," << Tz
         << "  bra_ind = " << bra_ind << "  ket_ind = " << ket_ind
         << endl
         << "    modelspace ket index = " << modelspace->GetKetIndex(ket)
         << "    modelspace bra index = " << modelspace->GetKetIndex(bra)
         << " out of " << modelspace->GetNumberKets()
         << endl
         << "  check_ket_channel = " << CheckChannel_ket(bra) << endl
         << " check parity: " << (la+lb)%2 << " = " << parity << "?" << endl
         << " check Tz: " << ta+tb << " = " << 2*Tz << "?" << endl
         << " check J: " << ja+jb << " >= " << 2*J << "?" << endl
         << " check J: " << abs(ja-jb) << " <= " << 2*J << "?" << endl

         << endl;
     return;
   }
   TBME(bra_ind, ket_ind) = tbme;
   //if (op->IsHermitian() )
   if ( hermitian )
      TBME(ket_ind, bra_ind) = tbme; 
   //if (op->IsAntiHermitian() )
   if ( antihermitian )
      TBME(ket_ind, bra_ind) = -tbme; 
}


bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
{
   if ((p==q) and (J%2 != 0)) return false; // Pauli principle
   Orbit * op = modelspace->GetOrbit(p);
   Orbit * oq = modelspace->GetOrbit(q);

   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

//===================================================================================
//===================================================================================
//  START IMPLEMENTATION OF OPERATOR METHODS
//===================================================================================
//===================================================================================

Operator::Operator()
{}

Operator::Operator(ModelSpace* ms)
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
//  OneBody.resize( nOneBody*nOneBody, 0.0);
  TwoBodyJmax = 0;
  nChannels = 2*3*JMAX;

  for (int J=0;J<JMAX;J++)
  {
    for (int par=0;par<=1;par++)
    {
      for (int tz=-1;tz<=1;tz++)
      {
         //TwoBody[J][par][(tz+1)] = new TwoBodyChannel(J,par,tz,this);
         //if (J>TwoBodyJmax && TwoBody[J][par][tz+1]->GetNumberKets() > 0 ) TwoBodyJmax = J;
         TwoBody[J][par][(tz+1)] = TwoBodyChannel(J,par,tz,this);
         if (J>TwoBodyJmax && TwoBody[J][par][tz+1].GetNumberKets() > 0 ) TwoBodyJmax = J;
      }
    }
  }

  for (int i=0;i<nOneBody;i++)
  {
     OneBody(i,i) = modelspace->GetOrbit(i)->spe;
  }

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

Operator::Operator(const Operator& op)
{
   modelspace = op.modelspace;
   hermitian = op.hermitian;
   antihermitian = op.antihermitian;
   ZeroBody = op.ZeroBody;
   OneBody = op.OneBody;
   for (int J=0;J<JMAX;J++)
   {
     for (int par=0;par<=1;par++)
     {
       for (int tz=0;tz<=2;tz++)
       {
          TwoBody[J][par][(tz+1)] = TwoBodyChannel(J,par,tz,this);
          TwoBody[J][par][(tz+1)].TBME = op.TwoBody[J][par][tz].TBME;
       }
     }
   }
   TwoBodyJmax = op.TwoBodyJmax;
   nChannels = op.nChannels;
}


Operator Operator::operator=(const Operator& rhs)
{
   return Operator(rhs);
}



void Operator::PrintOut() 
{
   //for (int j=0; j<JMAX; j++)
   for (int j=0; j<TwoBodyJmax; j++)
   {
     for (int p=0; p<=1; p++)
     {
        for (int t=-1; t<=1; t++)
        {
           TwoBodyChannel tbc = GetTwoBodyChannel(j,p,t);
           cout << tbc.J << " " << tbc.parity << " "
                     << tbc.Tz << "  ===> " << tbc.GetNumberKets() << endl;
           if (tbc.GetNumberKets()>20) continue;
           for (int i=0;i<tbc.GetNumberKets();i++)
           {
             for (int ii=0;ii<tbc.GetNumberKets();ii++)
             {
                cout << tbc.TBME(i,ii) << " ";
             }
                cout << endl;
           }
        }
     }
   }
}




Operator Operator::DoNormalOrdering()
{
   Operator opNO = Operator(modelspace);
   opNO.ZeroBody = ZeroBody;
   opNO.OneBody = OneBody;
   int norbits = modelspace->GetNumberOrbits();
   int nkets = modelspace->GetNumberKets();

   for (int k=0;k<norbits;++k)
   {
      Orbit * orb = modelspace->GetOrbit(k);
      if (orb->hvq > 0) continue;
      opNO.ZeroBody += (orb->j2+1) * OneBody(k,k);
   }

   int nchan = this->GetNumberTwoBodyChannels();
   for (int i=0;i<nchan;++i)
   {
      TwoBodyChannel tbc = GetTwoBodyChannel(i);
      opNO.SetTwoBodyChannel(i, tbc );
      int npq = tbc.GetNumberKets();
      for (int k=0;k<npq;++k)
      {
         Ket * ket = tbc.GetKet(k);
         Orbit *oa = modelspace->GetOrbit(ket->p);
         Orbit *ob = modelspace->GetOrbit(ket->q);
         if ( oa->hvq == 0 and ob->hvq == 0)
         {
            opNO.ZeroBody += tbc.TBME(k,k) * (2*tbc.J+1) / (1.0+ket->delta_pq());
         }

         for (int m=k;m<npq;++m)
         {
            Ket * bra = tbc.GetKet(m);
            Orbit *oc = modelspace->GetOrbit(bra->p);
            Orbit *od = modelspace->GetOrbit(bra->q);
            int J = tbc.J;
            if (bra->q == ket->q and ob->hvq==0)
            {
               opNO.OneBody(bra->p,ket->p) += (2*J+1)/(oa->j2+1) * tbc.TBME(m,k);
               opNO.OneBody(ket->p,bra->p) += (2*J+1)/(oa->j2+1) * tbc.TBME(m,k);
            }
            if (bra->p == ket->p and oa->hvq==0)
            {
               opNO.OneBody(bra->q,ket->q) += (2*J+1)/(ob->j2+1) * tbc.TBME(m,k);
            }
            if (bra->p == ket->q and oa->hvq==0)
            {
               opNO.OneBody(bra->q,ket->p) += (2*J+1)/(ob->j2+1) * bra->Phase(J) * tbc.TBME(m,k);
            }
         }

      }
   }
   return opNO;
}



void Operator::EraseTwoBody()
{
   for (int J=0;J<JMAX;J++)
   {
      for (int p=0;p<2;p++)
      {
         for (int t=0;t<3;t++)
         {
             TwoBody[J][p][t].TBME.zeros();
         }
      }
   }
}

/*
void Operator::WriteOneBody(const char* filename)
{
   ofstream obfile;
   obfile.open(filename, ofstream::out);
   int norbits = modelspace->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      for (int j=0;j<norbits;j++)
      {
         if (abs(OneBody(i,j)) > 1e-6)
         {
            obfile << i << "    " << j << "       " << OneBody(i,j) << endl;
         }
      }
   }
   obfile.close();
}
*/

/*
void Operator::WriteTwoBody(const char* filename)
{
   ofstream tbfile;
   tbfile.open(filename, ofstream::out);
   for (int Tz=-1;Tz<=1;++Tz)
   {
      for (int p=0;p<=1;++p)
      {
         for (int J=0;J<JMAX;++J)
         {
            int npq = TwoBody[J][p][Tz+1].GetNumberKets();
            if (npq<1) continue;
            for (int i=0;i<npq;++i)
            {
               int iket = TwoBody[J][p][Tz+1].GetKetIndex(i);
               Ket *bra = modelspace->GetKet(iket);
               Orbit *oa = modelspace->GetOrbit(bra->p);
               Orbit *ob = modelspace->GetOrbit(bra->q);
               for (int j=i;j<npq;++j)
               {
                  int jket = TwoBody[J][p][Tz+1].GetKetIndex(j);
                  Ket *ket = modelspace->GetKet(jket);
                  //double tbme = TwoBody[J][p][Tz+1].TBME(i,j);
                  double tbme = TwoBody[J][p][Tz+1].GetTBME(bra,ket);
                  if ( abs(tbme)<1e-4 ) continue;
                  Orbit *oc = modelspace->GetOrbit(ket->p);
                  Orbit *od = modelspace->GetOrbit(ket->q);
                  int wint = 4;
                  int wfloat = 12;
//                  tbfile 
//                         << setw(wint)   << oa->n  << setw(wint) << oa->l  << setw(wint)<< oa->j2 << setw(wint) << oa->tz2 
//                         << setw(wint+2) << ob->n  << setw(wint) << ob->l  << setw(wint)<< ob->j2 << setw(wint) << ob->tz2 
//                         << setw(wint+2) << oc->n  << setw(wint) << oc->l  << setw(wint)<< oc->j2 << setw(wint) << oc->tz2 
//                         << setw(wint+2) << od->n  << setw(wint) << od->l  << setw(wint)<< od->j2 << setw(wint) << od->tz2 
//                         << setw(wint+3) << J << setw(wfloat) << std::fixed << tbme// << endl;
//                         << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
//                         //<< setw(wint+2) << p   << setw(wint) << Tz << setw(wint) << J << setw(wfloat) << std::fixed << tbme << endl;
//
                  tbfile 
                         << setw(wint) << bra->p
                         << setw(wint) << bra->q
                         << setw(wint) << ket->p
                         << setw(wint) << ket->q
                         << setw(wint+3) << J << setw(wfloat) << std::fixed << tbme// << endl;
                         << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
                         //<< setw(wint+2) << p   << setw(wint) << Tz << setw(wint) << J << setw(wfloat) << std::fixed << tbme << endl;
               }
            }
         }
      }
   }
   tbfile.close();
}

*/

//Operator Operator::Commutator(const Operator& opleft, const Operator& opright)
Operator Operator::Commutator(const Operator& opright)
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   out.ZeroBody  = comm110(OneBody, opright.OneBody);
   out.OneBody  = comm111(OneBody, opright.OneBody);

   //for (int J=0;J<JMAX;J++)
   for (int J=0;J<TwoBodyJmax;J++)
   {
      for (int p=0;p<=1;p++)
      {
         for (int t=0;t<=2;t++)
         {
            out.ZeroBody += comm220(TwoBody[J][p][t], opright.TwoBody[J][p][t]);

            out.OneBody += comm211(TwoBody[J][p][t], opright.OneBody) - comm211(opright.TwoBody[J][p][t],OneBody);
            out.OneBody += comm221(TwoBody[J][p][t], opright.TwoBody[J][p][t]);

            out.TwoBody[J][p][t]  = comm112(OneBody, opright.OneBody);
            out.TwoBody[J][p][t] += comm212(TwoBody[J][p][t], opright.OneBody) - comm212(opright.TwoBody[J][p][t],OneBody);
            out.TwoBody[J][p][t] += comm222(TwoBody[J][p][t], opright.TwoBody[J][p][t]);
         }
      }
   }

   return out;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
///////////////////////////////////////////////////////////////////////////////////////////

float Operator::comm110(const arma::mat& obleft, const arma::mat& obright)
{
   float comm = 0;
   comm = trace(P_hole_onebody*(obleft*obright - obright*obleft));
   
   return comm;
}


//*****************************************************************************************

float Operator::comm220(const TwoBodyChannel& tbleft, const TwoBodyChannel& tbright)
{
   float comm = 0;
   int npq = tbleft.GetNumberKets();
   // Projector onto hole-hole states
   arma::mat P_hh = arma::mat(npq,npq,arma::fill::zeros);
   // Projector onto partile-particle states
   arma::mat P_pp = arma::mat(npq,npq,arma::fill::zeros);
   for(int i=0;i<npq;i++)
   {
      int ketindex = tbleft.GetKetIndex(i);
      int a = modelspace->GetKet(ketindex)->p;
      int b = modelspace->GetKet(ketindex)->q;
      if (modelspace->GetOrbit(a)->hvq == 0 and modelspace->GetOrbit(b)->hvq==0) P_hh(i,i)=1;
      if (modelspace->GetOrbit(a)->hvq > 0 and modelspace->GetOrbit(b)->hvq>0) P_pp(i,i)=1;
   }
   comm = arma::trace(P_hh*(tbleft.TBME * P_pp * tbright.TBME - tbright.TBME * P_pp * tbleft.TBME));

   return comm;
}

//*****************************************************************************************

arma::mat Operator::comm111(const arma::mat& obleft, const arma::mat& obright)
{
   arma::mat comm = obleft*obright - obright*obleft;
   return comm;
}

//*****************************************************************************************

// NOT YET FINISHED...
arma::mat Operator::comm211(const TwoBodyChannel& tbleft, const arma::mat& obright)
{
   arma::mat comm = obright;
   return comm;
}

//*****************************************************************************************

// NOT YET FINISHED...
arma::mat Operator::comm221(const TwoBodyChannel& tbleft, const TwoBodyChannel& tbright)
{
   int npq = tbleft.GetNumberKets();
   // Projector onto hole-hole states, etc
   arma::mat P_hh = arma::mat(npq,npq,arma::fill::zeros);
   arma::mat P_pp = arma::mat(npq,npq,arma::fill::zeros);
   arma::mat P_h = arma::mat(npq,npq,arma::fill::zeros);
   arma::mat P_p = arma::mat(npq,npq,arma::fill::zeros);
   for(int i=0;i<npq;i++)
   {
      int ketindex = tbleft.GetKetIndex(i);
      int a = modelspace->GetKet(ketindex)->p;
      int b = modelspace->GetKet(ketindex)->q;
      if (modelspace->GetOrbit(a)->hvq == 0 and modelspace->GetOrbit(b)->hvq==0) P_hh(i,i)=1;
      if (modelspace->GetOrbit(a)->hvq > 0 and modelspace->GetOrbit(b)->hvq>0) P_pp(i,i)=1;
      if (modelspace->GetOrbit(a)->hvq == 0 ) P_h(i,i)=1;
      if (modelspace->GetOrbit(a)->hvq > 0 ) P_p(i,i)=1;
   }

   arma::mat AB = P_h*(tbleft.TBME * P_pp * tbright.TBME - tbright.TBME * P_pp * tbleft.TBME);
             AB += P_p * (tbleft.TBME * P_hh * tbright.TBME - tbright.TBME * P_hh * tbleft.TBME);
   for (int a=0;a<npq;a++)
   {
      for (int b=0;b<npq;b++)
      {
         
      }
   }
   
   arma::mat comm = arma::mat(modelspace->GetNumberOrbits(),modelspace->GetNumberOrbits(),arma::fill::zeros);
   return comm;
}

//*****************************************************************************************

TwoBodyChannel Operator::comm112(const arma::mat& obleft, const arma::mat& obright)
{
   TwoBodyChannel comm = TwoBodyChannel();
   return comm;
}

//*****************************************************************************************

TwoBodyChannel Operator::comm212(const TwoBodyChannel& tbleft, const arma::mat& obright)
{
   TwoBodyChannel comm = tbleft;
   return comm;
}

//*****************************************************************************************

TwoBodyChannel Operator::comm222(const TwoBodyChannel& tbleft, const TwoBodyChannel& tbright)
{
   TwoBodyChannel comm = tbleft;
   return comm;
}

//*****************************************************************************************


