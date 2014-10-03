
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
   Number_pq_configs = 0;
   int nk = modelspace->GetNumberKets();
   pqindex.resize(nk,-1);
   for (int i=0;i<nk;i++)
   {
      Ket *ket = modelspace->GetKet(i);
      if (CheckChannel_pq(ket->p, ket->q))
      {
         pqindex[i] = Number_pq_configs;
         KetList.push_back(i);
         Number_pq_configs++;
      }
   }
   TBME = arma::mat(Number_pq_configs, Number_pq_configs, arma::fill::zeros);

}


TwoBodyChannel::TwoBodyChannel(const TwoBodyChannel& rhs)
{
   J                 = rhs.J;
   parity            = rhs.parity;
   Tz                = rhs.Tz;
   hermitian         = rhs.hermitian;
   antihermitian     = rhs.antihermitian;
   modelspace        = rhs.modelspace;
   Number_pq_configs = rhs.Number_pq_configs;
   TBME              = rhs.TBME;
   pqindex           = rhs.pqindex;
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
   Number_pq_configs = rhs.Number_pq_configs;
   TBME              = rhs.TBME;
   pqindex           = rhs.pqindex;
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



float TwoBodyChannel::GetTBME(int bra, int ket)
{
   int bra_ind = Getpqconfig(bra);
   int ket_ind = Getpqconfig(ket);
   if (bra_ind < 0 || ket_ind < 0) return 0;
   if (bra_ind > Number_pq_configs || ket_ind > Number_pq_configs) return 0;
   return TBME(bra_ind, ket_ind);
}

float TwoBodyChannel::GetTBME(Ket *bra, Ket *ket)
{
   int bra_ind = Getpqconfig(bra->p,bra->q);
   int ket_ind = Getpqconfig(ket->p,ket->q);
   if (bra_ind < 0 || ket_ind < 0) return 0;
   if (bra_ind > Number_pq_configs || ket_ind > Number_pq_configs) return 0;
   return TBME(bra_ind, ket_ind);
}

void TwoBodyChannel::SetTBME(int bra, int ket, float tbme)
{
   int bra_ind = Getpqconfig(bra);
   int ket_ind = Getpqconfig(ket);
   if ( (bra_ind < 0) or (ket_ind < 0) or
       (bra_ind > Number_pq_configs) or (ket_ind > Number_pq_configs) )
   {
     cerr << "Bad assignment of tbme!"
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
   int bra_ind = Getpqconfig(bra->p,bra->q);
   int ket_ind = Getpqconfig(ket->p,ket->q);
   if ( (bra_ind < 0) or (ket_ind < 0) or
       (bra_ind > Number_pq_configs) or (ket_ind > Number_pq_configs) )
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
     //cerr << "Bad assignment of tbme!"
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
         << "  check_pq_channel = " << CheckChannel_pq(bra->p,bra->q) << endl
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


bool TwoBodyChannel::CheckChannel_pq(int p, int q)
{
   int jp,lp,tp,jq,lq,tq;
   jp = modelspace->GetOrbit(p)->j2;
   lp = modelspace->GetOrbit(p)->l;
   tp = modelspace->GetOrbit(p)->tz2;
   jq = modelspace->GetOrbit(q)->j2;
   lq = modelspace->GetOrbit(q)->l;
   tq = modelspace->GetOrbit(q)->tz2;

   if ((lp+lq)%2 != parity) return false;
   if ((tp+tq) != 2*Tz) return false;
   if (jp+jq < 2*J) return false;
   if (abs(jp-jq) > 2*J) return false;
   if ((p==q) and (J%2!=0)) return false; // Pauli principle
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

  for (int J=0;J<JMAX;J++)
  {
    for (int par=0;par<=1;par++)
    {
      for (int tz=-1;tz<=1;tz++)
      {
         //TwoBody[J][par][(tz+1)] = new TwoBodyChannel(J,par,tz,this);
         //if (J>TwoBodyJmax && TwoBody[J][par][tz+1]->Number_pq_configs > 0 ) TwoBodyJmax = J;
         TwoBody[J][par][(tz+1)] = TwoBodyChannel(J,par,tz,this);
         if (J>TwoBodyJmax && TwoBody[J][par][tz+1].Number_pq_configs > 0 ) TwoBodyJmax = J;
      }
    }
  }
  GetSPEFromModelSpace();
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
}


Operator Operator::operator=(const Operator& rhs)
{
   return Operator(rhs);
}


void Operator::GetSPEFromModelSpace()
{
   int norbits = modelspace->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      OneBody(i,i) = modelspace->GetOrbit(i)->spe;
   }
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
           TwoBodyChannel *tbc = GetTwoBodyChannel(j,p,t);
           cout << tbc->J << " " << tbc->parity << " "
                     << tbc->Tz << "  ===> " << tbc->Number_pq_configs << endl;
           if (tbc->Number_pq_configs>20) continue;
           for (int i=0;i<tbc->Number_pq_configs;i++)
           {
             for (int ii=0;ii<tbc->Number_pq_configs;ii++)
             {
                cout << tbc->TBME(i,ii) << " ";
             }
                cout << endl;
           }
        }
     }
   }
}


Operator Operator::DoNormalOrdering()
{
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
            int npq = TwoBody[J][p][Tz+1].Number_pq_configs;
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
/*                  tbfile 
                         << setw(wint)   << oa->n  << setw(wint) << oa->l  << setw(wint)<< oa->j2 << setw(wint) << oa->tz2 
                         << setw(wint+2) << ob->n  << setw(wint) << ob->l  << setw(wint)<< ob->j2 << setw(wint) << ob->tz2 
                         << setw(wint+2) << oc->n  << setw(wint) << oc->l  << setw(wint)<< oc->j2 << setw(wint) << oc->tz2 
                         << setw(wint+2) << od->n  << setw(wint) << od->l  << setw(wint)<< od->j2 << setw(wint) << od->tz2 
                         << setw(wint+3) << J << setw(wfloat) << std::fixed << tbme// << endl;
                         << "    < " << bra->p << " " << bra->q << " | V | " << ket->p << " " << ket->q << " >" << endl;
                         //<< setw(wint+2) << p   << setw(wint) << Tz << setw(wint) << J << setw(wfloat) << std::fixed << tbme << endl;
*/
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

Operator Operator::Commutator(const Operator& opleft, const Operator& opright)
{
   Operator out = opright;
   out.EraseZeroBody();
   out.EraseOneBody();
   out.EraseTwoBody();

   out.ZeroBody  = comm110(opleft.OneBody, opright.OneBody);
   out.OneBody  = comm111(opleft.OneBody, opright.OneBody);

   for (int J=0;J<JMAX;J++)
   {
      for (int p=0;p<2;p++)
      {
         for (int t=0;t<3;t++)
         {
//            out.ZeroBody += comm210(opleft.TwoBody[J][p][t], opright.OneBody) - comm210(opright.TwoBody[J][p][t],opleft.OneBody);
            out.ZeroBody += comm220(opleft.TwoBody[J][p][t], opright.TwoBody[J][p][t]);

            out.OneBody += comm211(opleft.TwoBody[J][p][t], opright.OneBody) - comm211(opright.TwoBody[J][p][t],opleft.OneBody);
            out.OneBody += comm221(opleft.TwoBody[J][p][t], opright.TwoBody[J][p][t]);

            out.TwoBody[J][p][t]  = comm112(opleft.OneBody, opright.OneBody);
            out.TwoBody[J][p][t] += comm212(opleft.TwoBody[J][p][t], opright.OneBody) - comm212(opright.TwoBody[J][p][t],opleft.OneBody);
            out.TwoBody[J][p][t] += comm222(opleft.TwoBody[J][p][t], opright.TwoBody[J][p][t]);
         }
      }
   }

   return out;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
// Below is the implementation of the commutators in the various channels
//

float Operator::comm110(arma::mat obleft, arma::mat obright)
{
   float comm = 0;
   comm = trace(P_hole_onebody*(obleft*obright - obright*obleft));
   
   return comm;
}


//*****************************************************************************************

float Operator::comm220(TwoBodyChannel tbleft, TwoBodyChannel tbright)
{
   float comm = 0;
   int npq = tbleft.Number_pq_configs;
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

arma::mat Operator::comm111(arma::mat obleft, arma::mat obright)
{
   arma::mat comm = obleft*obright - obright*obleft;
   return comm;
}

//*****************************************************************************************

arma::mat Operator::comm211(TwoBodyChannel tbleft, arma::mat obright)
{
   arma::mat comm = obright;
   return comm;
}

//*****************************************************************************************

// NOT YET FINISHED...
arma::mat Operator::comm221(TwoBodyChannel tbleft, TwoBodyChannel tbright)
{
   int npq = tbleft.Number_pq_configs;
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

TwoBodyChannel Operator::comm112(arma::mat obleft, arma::mat obright)
{
   TwoBodyChannel comm = TwoBodyChannel();
   return comm;
}

//*****************************************************************************************

TwoBodyChannel Operator::comm212(TwoBodyChannel tbleft, arma::mat obright)
{
   TwoBodyChannel comm = tbleft;
   return comm;
}

//*****************************************************************************************

TwoBodyChannel Operator::comm222(TwoBodyChannel tbleft, TwoBodyChannel tbright)
{
   TwoBodyChannel comm = tbleft;
   return comm;
}

//*****************************************************************************************


