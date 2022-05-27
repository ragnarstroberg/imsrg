
#include "ModelSpace.hh"


TwoBodyChannel::~TwoBodyChannel()
{
//  std::cout << "In TwoBodyChannel destructor" << std::endl;
}

TwoBodyChannel::TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace *ms)
  :  J(j), parity(p), Tz(t) , modelspace(ms)
{
//  Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
//  Initialize(j,p,t, ms);
  Initialize();
}

TwoBodyChannel::TwoBodyChannel(int ch, ModelSpace *ms)
 : modelspace(ms)
{
//   Initialize(ch,ms);
//   int j,p,t;
   ms->UnpackTwoBodyChannelIndex(ch, J,parity,Tz);
//   Initialize(J,parity,Tz, ms);
   Initialize();
//   Initialize(ch,ms);
}

//void TwoBodyChannel::Initialize(int ch, ModelSpace *ms)
//void TwoBodyChannel::Initialize(int J, int parity, int Tz, ModelSpace *ms)
void TwoBodyChannel::Initialize()
{
//   modelspace = ms;
//   modelspace->UnpackTwoBodyChannelIndex(ch,  J,parity,Tz);
//   int tbjmax = modelspace->TwoBodyJmax;

   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      Ket &ket = modelspace->GetKet(i);
//      std::cout << "J p Tz = " << J << " " << parity << " " << Tz << "   checking ket " << i << " -> " << ket.p << " , " << ket.q << std::endl;
      if ( CheckChannel_ket(ket) )
      {
//         std::cout << "       yes " << std::endl;
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         NumberKets++;
      }
   }
   KetIndex_pp = GetKetIndexFromList(modelspace->KetIndex_pp);
   KetIndex_hh = GetKetIndexFromList(modelspace->KetIndex_hh);
   KetIndex_ph = GetKetIndexFromList(modelspace->KetIndex_ph);
   KetIndex_cc = GetKetIndexFromList(modelspace->KetIndex_cc);
   KetIndex_vc = GetKetIndexFromList(modelspace->KetIndex_vc);
   KetIndex_qc = GetKetIndexFromList(modelspace->KetIndex_qc);
   KetIndex_vv = GetKetIndexFromList(modelspace->KetIndex_vv);
   KetIndex_qv = GetKetIndexFromList(modelspace->KetIndex_qv);
   KetIndex_qq = GetKetIndexFromList(modelspace->KetIndex_qq);
   std::vector<double> occvec;
   std::vector<double> unoccvec;
   for (index_t i=0;i<modelspace->KetIndex_hh.size();++i)
   {
      if (CheckChannel_ket(modelspace->GetKet(modelspace->KetIndex_hh[i])))
      {
        occvec.push_back( modelspace->Ket_occ_hh[i]);
        unoccvec.push_back( modelspace->Ket_unocc_hh[i]);
      }
   }
   Ket_occ_hh = arma::vec(occvec);
   Ket_unocc_hh = arma::vec(unoccvec);
   occvec.clear();
   unoccvec.clear();
   for (index_t i=0;i<modelspace->KetIndex_ph.size();++i)
   {
      if (CheckChannel_ket(modelspace->GetKet(modelspace->KetIndex_ph[i])))
      {
        occvec.push_back( modelspace->Ket_occ_ph[i]);
        unoccvec.push_back( modelspace->Ket_unocc_ph[i]);
      }
   }
   Ket_occ_ph = arma::vec(occvec);
   Ket_unocc_ph = arma::vec(unoccvec);
}


//int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 
size_t TwoBodyChannel::GetLocalIndex(int p, int q) const
{
 if (p<=q)
   return KetMap[modelspace->GetKetIndex(p,q)];
 else
   return KetMap[modelspace->GetKetIndex(q,p)] + NumberKets;
} 

// get pointer to ket using local index
const Ket & TwoBodyChannel::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; 
Ket & TwoBodyChannel::GetKet(int i) { return modelspace->GetKet(KetList[i]);}; 


//bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
bool TwoBodyChannel::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->index==oq->index) and (J%2 != 0)) return false; // Pauli principle
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (std::abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

const arma::uvec& TwoBodyChannel::GetKetIndex_pp() const { return KetIndex_pp;};
const arma::uvec& TwoBodyChannel::GetKetIndex_hh() const { return KetIndex_hh;};
const arma::uvec& TwoBodyChannel::GetKetIndex_ph() const { return KetIndex_ph;};
const arma::uvec& TwoBodyChannel::GetKetIndex_cc() const { return KetIndex_cc;};
const arma::uvec& TwoBodyChannel::GetKetIndex_vc() const { return KetIndex_vc;};
const arma::uvec& TwoBodyChannel::GetKetIndex_qc() const { return KetIndex_qc;};
const arma::uvec& TwoBodyChannel::GetKetIndex_vv() const { return KetIndex_vv;};
const arma::uvec& TwoBodyChannel::GetKetIndex_qv() const { return KetIndex_qv;};
const arma::uvec& TwoBodyChannel::GetKetIndex_qq() const { return KetIndex_qq;};



arma::uvec TwoBodyChannel::GetKetIndexFromList(std::vector<index_t>& vec_in)
{
   std::vector<index_t> index_list (std::min(vec_in.size(),KetList.size()));
   auto it = set_intersection(KetList.begin(),KetList.end(),vec_in.begin(),vec_in.end(),index_list.begin());
   index_list.resize(it-index_list.begin());
   for (auto& x : index_list)
   {
     x = KetMap[x];
   }
   return arma::uvec(index_list);
}

//************************************************************************
//************************************************************************
//************************************************************************

TwoBodyChannel_CC::~TwoBodyChannel_CC()
{
//   std::cout << "In TwoBodyChannel_CC destructor" << std::endl;
}

TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
//  : J(j), parity(p) , Tz(t) , modelspace(ms)
  : TwoBodyChannel(j,p,t,ms)
{
//  Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
  Initialize();
//  Initialize(j,p,t, ms);
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
// : modelspace(ms)
  : TwoBodyChannel(N,ms)
{
//   int j,p,t;
//   ms->UnpackTwoBodyChannelIndex_CC(N, j,p,t);
   ms->UnpackTwoBodyChannelIndex_CC(N, J,parity,Tz);
//   Initialize(j,p,t, ms);
   Initialize();
//   Initialize(N,ms);
}


// Check if orbits pq participate in this cross-coupled two-body channel
// Difference from regular channels:
// no Pauli rule, <pp||nn> is allowed. But |Tz| is still conserved, <-- only for rankT=0 operators.
// i.e. <pp||pn> is not allowed. So we use |Tz| rather than Tz,
// and don't use Tz=-1.
// Another way of formulating it is that "isospin-parity", i.e. |Tz|%2, is conserved the same way parity is.
bool TwoBodyChannel_CC::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->l + oq->l)%2 != parity)    return false;
   if (op->j2 + oq->j2 < 2*J)          return false;
   if (std::abs(op->j2 - oq->j2) > 2*J)     return false;
//   if (modelspace->single_species)
//   {
//     if (std::abs(op->tz2 + oq->tz2) != 2*std::abs(Tz)) return false;
//   }
//   else
   if (not modelspace->single_species)
   {
//     if (std::abs(op->tz2 + oq->tz2) != 2*Tz) return false;
     if (std::abs(op->tz2 - oq->tz2) != 2*Tz) return false;
   }

   return true;
}


