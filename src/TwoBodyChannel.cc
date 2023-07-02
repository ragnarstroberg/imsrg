
#include "ModelSpace.hh"


TwoBodyChannel_base::~TwoBodyChannel_base()
{}

TwoBodyChannel_base::TwoBodyChannel_base()
{}

TwoBodyChannel_base::TwoBodyChannel_base(int j, int p, int t, ModelSpace *ms)
  :  J(j), parity(p), Tz(t) , modelspace(ms)
{}

TwoBodyChannel_base::TwoBodyChannel_base(int ch, ModelSpace *ms)
 : modelspace(ms)
{}

void TwoBodyChannel_base::Initialize()
{
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      Ket &ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
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


size_t TwoBodyChannel_base::GetLocalIndex(int p, int q) const
{
 if (p<=q)
   return KetMap[modelspace->GetKetIndex(p,q)];
 else
   return KetMap[modelspace->GetKetIndex(q,p)] + NumberKets;
} 

// get pointer to ket using local index
const Ket & TwoBodyChannel_base::GetKet(int i) const { return modelspace->GetKet(KetList[i]);}; 
Ket & TwoBodyChannel_base::GetKet(int i) { return modelspace->GetKet(KetList[i]);}; 


bool TwoBodyChannel_base::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   std::cout << "====////////========/////// SHOULDNT BE HERE LINE " <<__LINE__ << "=====///////=======" << std::endl;
   return false;
}



const arma::uvec& TwoBodyChannel_base::GetKetIndex_pp() const { return KetIndex_pp;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_hh() const { return KetIndex_hh;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_ph() const { return KetIndex_ph;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_cc() const { return KetIndex_cc;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_vc() const { return KetIndex_vc;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_qc() const { return KetIndex_qc;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_vv() const { return KetIndex_vv;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_qv() const { return KetIndex_qv;};
const arma::uvec& TwoBodyChannel_base::GetKetIndex_qq() const { return KetIndex_qq;};



// Take in a vector of indices corresponding to the global list of kets contained in ModelSpace
// select the indices that belong to this TwoBodyChannel and replace them with their
// corresponding local indices (relevant for making dense matrices), and cast this
// to an armadillo uvec, suitable for indexing armadillo matrices.
arma::uvec TwoBodyChannel_base::GetKetIndexFromList(std::vector<index_t>& vec_in)
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
//****** Implementation of pp channels (i.e. not Pandya tansformed) ******
//************************************************************************


TwoBodyChannel::TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace *ms)
 : TwoBodyChannel_base( j,p,t,ms )
{
   Initialize();
}

TwoBodyChannel::TwoBodyChannel(int N, ModelSpace *ms)
{
   ms->UnpackTwoBodyChannelIndex_CC(N, J,parity,Tz);
   Initialize();
}


// Check if orbits pq participate in this two-body channel
bool TwoBodyChannel::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->index==oq->index) and (J%2 != 0)) return false; // Pauli principle
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (std::abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}




//************************************************************************
//************ Implementation of Pandya tansformed  channels *************
//************************************************************************


TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
 : TwoBodyChannel_base( j,p,t,ms )
{
   Initialize();
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
{
   ms->UnpackTwoBodyChannelIndex_CC(N, J,parity,Tz);
   Initialize();
}

// Check if orbits pq participate in this cross-coupled two-body channel
// Difference from regular channels:
// no Pauli rule, and isospin difference isconserved, rather than isospin sum
// Since we don't distinguish between orderings ab and ba, we just use |tza - tzb| and don't use Tz=-1.
bool TwoBodyChannel_CC::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->l + oq->l)%2 != parity)    return false;
   if (op->j2 + oq->j2 < 2*J)          return false;
   if (std::abs(op->j2 - oq->j2) > 2*J)     return false;
   if (std::abs(op->tz2 - oq->tz2) != 2*Tz) return false;

   return true;
}



