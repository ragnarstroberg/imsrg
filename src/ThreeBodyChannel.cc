
#include "ModelSpace.hh"

ThreeBodyChannel::ThreeBodyChannel(int twoj, int p, int twotz, ModelSpace* ms)
 : twoJ(twoj), parity(p), twoTz(twotz), modelspace(ms)
{
  Initialize( );
}


void ThreeBodyChannel::Initialize()
{
  // Figure out which 3 body kets belong to this channel
  KetList.resize(0);
  KetMap.resize(0);
  for ( size_t iket3=0; iket3<modelspace->Kets3.size(); iket3++ )
  {
    auto& ket3 = modelspace->Kets3[iket3];
    if (CheckChannel_ket(ket3) )
    {
//      if (twoJ==1 and parity==1 and twoTz==-3)
//      {
//        std::cout << "@@@@@@  passed the test for iket3 = " << iket3 << "    " << ket3.p << " " << ket3.q << " " << ket3.r << " " << ket3.Jpq << std::endl;
//      }
      KetList.push_back(iket3); // list of indices[ 1,  3, 4,  8,...] of locations in ModelSpace::Kets3 that participate in this channel
      if ( iket3 >= KetMap.size()) KetMap.resize(2*iket3+1,-1);
      KetMap[iket3] = KetList.size()-1; // a mapping of global ket index to where it sits in this channels KetList (-1 means it doesn't participate)
    }
  }
//  std::cout << "Initializing ThreeBodyChannel with JPT " << twoJ << " " << parity << " " << twoTz << " and KetMap looks like" << std::endl;
//  for ( auto k : KetMap ) std::cout << " " << k;
//  std::cout << std::endl;
}


int ThreeBodyChannel::GetNumber3bKets()
{
  return KetList.size();
}

// We will need to deal with requests of indices with the wrong ordering...
size_t ThreeBodyChannel::GetLocalIndex( int p, int q, int r, int Jpq )
{
//   std::cout << "|||||In " << __func__ << " pqr Jpq = " << p << " " << q << " " << r << "  " << Jpq << std::endl;
//   std::cout <<  "       Ket3Index = " << modelspace->GetKet3Index(p,q,r,Jpq) << std::endl;
//   std::cout << "|||||In " << __func__ << " Ket3Index = " << modelspace->GetKet3Index(p,q,r,Jpq) << std::endl;
//   std::cout << "pqr,Jpq = " << p << " " << q << " " << r << " " << Jpq << "    and this channel has JPT = " << twoJ << " " << parity << " " << twoTz << std::endl;
//   return KetMap[ modelspace->GetKet3Index(p,q,r,Jpq)];
   return KetMap.at( modelspace->GetKet3Index(p,q,r,Jpq));
}


bool ThreeBodyChannel::CheckChannel_ket( Ket3& ket3)
{
  if ( (ket3.op->l + ket3.oq->l + ket3.oR->l)%2 != parity ) return false;  // need to match parity
  if ( (ket3.op->tz2 + ket3.oq->tz2 + ket3.oR->tz2) != twoTz ) return false;  // need to match Tz
  if ( std::abs( 2*ket3.Jpq - ket3.oR->j2)>twoJ or (2*ket3.Jpq+ket3.oR->j2)<twoJ ) return false; // triangle condition
// TODO Put these back in. This is only for testing...
  if ( twoJ == (ket3.op->j2 + ket3.oq->j2 + ket3.oR->j2) )
  {
    if ( ket3.p==ket3.q   or  ket3.p==ket3.r  or ket3.q==ket3.r ) return false;
  }
  if ( twoJ==(ket3.op->j2+ket3.oq->j2+ket3.oR->j2)  and (ket3.p==ket3.q or ket3.q==ket3.r) ) return false; // fully stretched state can't have p==q
  if ( (ket3.p==ket3.q) and std::abs(twoJ-ket3.oR->j2)==2*ket3.op->j2 ) return false; // fully stretched state can't have p==q
  if ( (ket3.q==ket3.r) and std::abs(twoJ-ket3.op->j2)==2*ket3.oq->j2 ) return false; // fully stretched state can't have p==q
  if ( (ket3.p==ket3.q) and (ket3.p==ket3.r) and twoJ>(3*ket3.op->j2 - 3) ) return false;
  return true;
}


Ket3& ThreeBodyChannel::GetKet(size_t iket)
{
 return modelspace->GetKet3(KetList[iket]);
};


size_t ThreeBodyChannel::GetNumberKets(){ return KetList.size();};

