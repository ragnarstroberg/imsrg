
#include "ModelSpace.hh"



//ThreeBodyChannel(int twoj, int p, int twotz, ModelSpace* ms)
ThreeBodyChannel::ThreeBodyChannel(int twoj, int p, int twotz, ModelSpace* ms)
 : twoJ(twoj), parity(p), twoTz(twotz), modelspace(ms)
{
//  std::cout << "IN " << __func__ << std::endl;
  Initialize( );
}


void ThreeBodyChannel::Initialize()
{
//  std::cout << "IN " << __func__ << " JPT = " << twoJ << " " << parity << " " << twoTz << std::endl;
  // Figure out which 3 body kets belong to this channel
  KetList.resize(0);
  KetMap.resize(0);
  for ( size_t iket3=0; iket3<modelspace->Kets3.size(); iket3++ )
  {
//    std::cout << "iket3 = " << iket3 << std::endl;
    auto& ket3 = modelspace->Kets3[iket3];
    if (CheckChannel_ket(ket3) )
    {
//      std::cout << "pushing back  " << iket3 << "    " << ket3.p << " " << ket3.q << " " << ket3.r << " " << ket3.Jpq << std::endl;
      KetList.push_back(iket3); // list of indices[ 1,  3, 4,  8,...] of locations in ModelSpace::Kets3 that participate in this channel
//      std::cout << "maybe resizing? iket3 =" << iket3 << " KetMap size = "<<  KetMap.size() << std::endl;
      if ( iket3 >= KetMap.size()) KetMap.resize(2*iket3+1,-1);
//      std::cout << "mapping Ketlist size = " << KetList.size() << " KetMap size = " << KetMap.size() << std::endl;
      KetMap[iket3] = KetList.size()-1; // a mapping of global ket index to where it sits in this channels KetList (-1 means it doesn't participate)
//      std::cout << "done mapping.  " << iket3 << " -> " << KetMap[iket3] << std::endl;
    }
  }
//  std::cout << "done with " << __func__ << std::endl;
}


int ThreeBodyChannel::GetNumber3bKets()
{
  return KetList.size();
}

// We will need to deal with requests of indices with the wrong ordering...
size_t ThreeBodyChannel::GetLocalIndex( int p, int q, int r, int Jpq )
{
//   std::cout << "IN " << __func__ << std::endl;
//   std::cout << "jpt = " << twoJ << " " << parity << " " << twoTz << std::endl;
//   std::cout << "size of KetMap = " << KetMap.size() << std::endl;
//   std::cout << " looking for " << modelspace->GetKet3Index(p,q,r,Jpq) << std::endl;
//   std::cout << "pqr Jpq = " << p << " " << q << " " << r << " " << Jpq << std::endl;
//   std::cout << "by the way, minus one looks like " << size_t(-1) << std::endl;
   return KetMap[ modelspace->GetKet3Index(p,q,r,Jpq)];
}


bool ThreeBodyChannel::CheckChannel_ket( Ket3& ket3)
{
  if ( (ket3.op->l + ket3.oq->l + ket3.oR->l)%2 != parity ) return false; 
  if ( (ket3.op->tz2 + ket3.oq->tz2 + ket3.oR->tz2) != twoTz ) return false; 
  if ( std::abs( 2*ket3.Jpq - ket3.oR->j2)>twoJ or (2*ket3.Jpq+ket3.oR->j2)<twoJ ) return false;
  return true;
}


