
#include "ThreeBodyMEpn.hh"
#include "AngMom.hh"
#include <set>
#include <tuple>
#include <sstream>



ThreeBodyMEpn::ThreeBodyMEpn()
: herm(1), PN_mode(false) {}

ThreeBodyMEpn::ThreeBodyMEpn(ModelSpace* ms)
 : modelspace(ms), isospin3BME(ms), herm(1), PN_mode(false)
{
//  Allocate();
}

ThreeBodyMEpn::ThreeBodyMEpn(ModelSpace* ms, int e3max)
: modelspace(ms), isospin3BME(ms,e3max), E3max(e3max), PN_mode(false)
{
//  Allocate();
}


ThreeBodyMEpn::ThreeBodyMEpn(const ThreeBodyMEpn& tbme)
 : modelspace(tbme.modelspace), MatEl(tbme.MatEl), isospin3BME(tbme.isospin3BME), emax(tbme.emax), E3max(tbme.E3max),
   herm(tbme.herm), PN_mode(tbme.PN_mode)
{
//  Allocate();
}


void ThreeBodyMEpn::Allocate()
{
  if (PN_mode) Allocate_PN();
  else Allocate_Isospin();
}


void ThreeBodyMEpn::Allocate_Isospin()
{
  isospin3BME.Allocate();
}


// This will need to be more elaborate if we want to use tensor 3-body.
void ThreeBodyMEpn::Allocate_PN()
{
  total_dimension = 0;
  size_t nch = modelspace->GetNumberThreeBodyChannels();
  for (size_t ch=0; ch<nch; ch++)
//  for (auto Tbc : modelspace->ThreeBodyChannels )
  {
    ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel( ch );
    size_t channel_dim = Tbc.GetNumber3bKets(); // Number of kets in this 3body J,p,Tz channel
    MatEl[{ch,ch}] = SymmMatrix<ME_type>( channel_dim );
    total_dimension += channel_dim * (channel_dim+1)/2;
  }
}



/// interface methods. When calling these, the user shouldn't need to care whether
/// we're storing the matrix elements in isospin or PN formalism.

//access in PN formalism (regardless of how it's stored)

ThreeBodyMEpn::ME_type ThreeBodyMEpn::GetME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n) const
{
  if (PN_mode) return GetME_pn_PN( Jab_in, Jde_in, J2, i,j,k,l,m,n);
  else return isospin3BME.GetME_pn(Jab_in, Jde_in, J2, i,j,k,l,m,n);
}

void ThreeBodyMEpn::SetME_pn(  int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V)
{
  if (PN_mode)  SetME_pn_PN( Jab_in, Jde_in, J2, i,j,k,l,m,n,V);
//  else isospin3BME.SetME_pn( Jab_in, Jde_in, J2, i,j,k,l,m,n,V);
}

void ThreeBodyMEpn::AddToME_pn(int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V)
{
  if (PN_mode)  AddToME_pn_PN( Jab_in, Jde_in, J2, i,j,k,l,m,n,V);
  else isospin3BME.AddToME_pn( Jab_in, Jde_in, J2, i,j,k,l,m,n,V);
}

//access in isospin formalism (regardless of how it's stored)

ThreeBodyMEpn::ME_type ThreeBodyMEpn::GetME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n)
{
  if (PN_mode) return GetME_PN( Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n);
  else return isospin3BME.GetME(Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n);
}

void ThreeBodyMEpn::SetME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V)
{
  if (PN_mode) SetME_PN( Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n,V);
  else isospin3BME.SetME(Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n,V);
}

void ThreeBodyMEpn::AddToME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V)
{
  if (PN_mode) AddToME_PN( Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n,V);
  else isospin3BME.AddToME(Jab_in, Jde_in, J2, tab_in,tde_in,twoT,i,j,k,l,m,n,V);
}







/// These are the ones that eventually get called, but typically the other methods will
/// be more convenient to call.
ThreeBodyMEpn::ME_type ThreeBodyMEpn::GetME_pn_PN_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket) const
{
//  std::cout << "IN " << __func__ << std::endl;
//  std::cout << "ch_bra ,ch_ket " << ch_bra << " " << ch_ket << std::endl;
//  std::cout << std::endl << MatEl.at({ch_bra,ch_ket}).FullMatrix() <<std::endl << std::endl;
//  return  MatEl[{ch_bra,ch_ket}].Get(ibra,iket) ;
  auto matrix = MatEl.at({ch_bra,ch_ket});
  size_t dim = matrix.size();
  if (ibra>=dim or iket>=dim)
  {
//    return 0;
    std::ostringstream oss;
    oss << __func__ << " ch_bra,ch_ket " << ch_bra << " " << ch_ket << "  ibra,iket " << ibra << " " << iket;
//    throw std::domain_error( "GetME_pn_ch, ibra or iket > dim");
    throw std::domain_error( oss.str() );
  }
  return  matrix.Get(ibra,iket) ;
}

// We have this here in case we want to set a matrix element, but we store it in a different
// coupling order. In that case, we need to add to multiple matrix elements with the appropriate
// recoupling coefficients included (see below).
void ThreeBodyMEpn::AddToME_pn_PN_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyMEpn::ME_type matel)
{
//  std::cout << "IN " << __func__ << std::endl;
//  if (ch_bra == 12 and ch_ket==12 and ibra==6 and iket==6)  std::cout << "IN " << __func__ << "  adding " << matel << std::endl;
  auto& symm = MatEl.at({ch_bra,ch_ket});
  size_t dim = symm.size();
  if (ibra>=dim or iket>=dim)
  {
//    return;
    std::cout << "ibra,iket matel " << ibra << " " << iket << " " << matel << std::endl;
    throw std::domain_error( "AddToME_pn_ch, ibra or iket > dim");
  }
//  std::cout << "about to Get ibra,iket = " << ibra << " " << iket << "  dim = " << symm.size() << std::endl;
  auto val = symm.Get(ibra,iket);
  symm.Put(ibra,iket, val + matel);
}

void ThreeBodyMEpn::SetME_pn_PN_ch(size_t ch_bra, size_t ch_ket, size_t ibra, size_t iket, ThreeBodyMEpn::ME_type matel)
{
//  std::cout << "IN " << __func__ << std::endl;
  auto& matrix = MatEl.at({ch_bra,ch_ket});
  size_t dim = matrix.size();
  if ( ibra>=dim or iket>=dim)
  {
//    return;
    throw std::domain_error( "SetME_pn_ch, ibra or iket > dim");
  }
  MatEl.at({ch_bra,ch_ket}).Put(ibra,iket, matel);
}







void ThreeBodyMEpn::AddToME_pn_PN(  int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyMEpn::ME_type me_add )
{

  std::vector<double>  recouple_bra, recouple_ket;
  std::vector<size_t>  ibra, iket;

  if (a==b and Jab%2>0) return;
  if (d==e and Jde%2>0) return;

//  std::cout << "Do the recoupling" << std::endl;
//  size_t ch_bra = GetKetIndex_withRecoupling( twoJ, Jab, a,b,c, ibra, recouple_bra );
//  size_t ch_ket = GetKetIndex_withRecoupling( twoJ, Jde, d,e,f, iket, recouple_ket );
  size_t ch_bra = GetKetIndex_withRecoupling( Jab, twoJ, a,b,c, ibra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( Jde, twoJ, d,e,f, iket, recouple_ket );

  if ( ibra.size()<1 or iket.size()<1) return;
//  std::cout << "IN " << __func__ << "  ch: " << ch_bra << " " << ch_ket << "   abcdef " << a << " " << b << " " << c << " " << d << " "<< e << " " << f << "    Jab Jde twoJ = " << Jab << " " << Jde << " " << twoJ << std::endl;
  if ( ch_bra != ch_ket) return ;
  ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch_bra);
//  std::cout << "channel " << ch_bra << " corresponds to " << Tbc.twoJ << " " << Tbc.parity << " " << Tbc.twoTz << std::endl;


  double overlap_bra_ket_in = 0;
  for ( int i=0; i<ibra.size(); i++)
  {
    for (int j=0; j<iket.size(); j++)
    {
      if (ibra[i] == iket[j])  overlap_bra_ket_in += recouple_bra[i] * recouple_ket[j] ;
    }
  }

//  std::cout << "looping" << std::endl;
//  std::cout << " just looking over ibra and iket real quck..." << std::endl;
//  std::cout << "sizes: " << ibra.size() << " " << iket.size() << std::endl;
//  for ( auto i : ibra ) std::cout << i << " ";
//  std::cout << std::endl;
//  for ( auto i : iket ) std::cout << i << " ";
//  std::cout << std::endl;
//  if ( Jab==1 and Jde==1 and twoJ==3 and a==5 and b==1 and c==0 and d==5 and e==1 and f==0) std::cout << "IN " << __func__ << "  adding " << me_add << std::endl;
//  if ( Jab==1 and Jde==1 and twoJ==3 ) std::cout << "IN " << __func__ << " " << a << " " << b<< " " << c << " " << d << " " << e << " " << f << "  adding " << me_add << std::endl;
//  std::cout << "IN " << __func__ << " " << a << " " << b<< " " << c << " " << d << " " << e << " " << f << "  adding " << me_add << std::endl;
  double me_out = 0;
//  double symmetry_factor = 2;
//  double symmetry_factor = 1;
//  if ( ibra[0]==iket[0] ) symmetry_factor = (1 + herm);
  double symmetry_factor = (ibra[0] == iket[0] and Jab==Jde) ? 0.5 : 1;
  double normalization_denom = 1 + herm * overlap_bra_ket_in * overlap_bra_ket_in;
//  std::cout << " normalization_denom = " << normalization_denom << std::endl;
  if ( std::abs(normalization_denom) < 1e-8 ) return;
  symmetry_factor = 1.0 / normalization_denom;

//  std::cout << "sizes : " << ibra.size() << "  " << iket.size() << std::endl;
//  herm = +1;
////  if (Jab==Jde and a==d and b==e and c==f) symmetry_factor = 1;
  for ( int i=0; i<ibra.size(); i++)
  {
    for (int j=0; j<iket.size(); j++)
    {
//     symmetry_factor = (ibra[i] == iket[j]) ? 1 + herm : 1;
//     double symmetry_factor2 = (ibra[i] == iket[j]) ? 0.5 : 1;
     double symmetry_factor2 = ( ibra[i] == iket[j]) ? 1+herm : 1;
//      if ( iket[j] > ibra[i] and std::find( ibra.begin(), ibra.end(), iket[j]) != ibra.end() ) continue;
//     if ( ibra[i] == iket[j] ) symmetry_factor2 +=herm;
//if ( Jab==1 and Jde==1 and twoJ==3 )      std::cout << " call AddToME_pn_ch ( " << ch_bra << ", " << ch_ket << ", " << ibra[i] << " " << iket[j] << std::endl;
//     std::cout << " call AddToME_pn_ch (  ij= " << i << " " << j << "   " << ch_bra << ", " << ch_ket << ", " << ibra[i] << " " << iket[j] << "   recouple : " << recouple_bra[i] << " " << recouple_ket[j] << "  symmetry = " << symmetry_factor << " " << symmetry_factor2 << std::endl;
       AddToME_pn_PN_ch( ch_bra, ch_ket, ibra[i], iket[j], recouple_bra[i] * recouple_ket[j] * me_add * symmetry_factor * symmetry_factor2   );
//       if ( iket[j] == ibra[i])
//       AddToME_pn_ch( ch_ket, ch_bra, iket[j], ibra[i], recouple_ket[j] * recouple_bra[i] * me_add * herm );
//       AddToME_pn_ch( ch_ket, ch_bra, iket[j], ibra[i], recouple_ket[j] * recouple_bra[i] * me_add * herm * symmetry );
    }
  }
//  std::cout << " done adding,  the matrix element is now " << GetME_pn(Jab,Jde,twoJ,a,b,c,d,e,f) << std::endl;
//  std::cout << "Done." << std::endl;
}



// This does the recoupling twice, so it could certainly be more efficient.
// But it's easier to code and understand this way, so I'll change it if need be.
void ThreeBodyMEpn::SetME_pn_PN(  int Jab_in, int Jde_in, int twoJ, int a, int b, int c, int d, int e, int f, ThreeBodyMEpn::ME_type me_set )
{
  double me_previous = GetME_pn_PN( Jab_in, Jde_in, twoJ, a,b,c,d,e,f);
//  std::cout << "IN " << __func__ << "  me_set = " << me_set << "  me_previous = " << me_previous << std::endl;
  AddToME_pn_PN( Jab_in, Jde_in, twoJ, a,b,c,d,e,f,  me_set-me_previous );
}



//// Let's not implement this until we actually need it...
////
//ME_type ThreeBodyMEpn::GetME_pn(size_t ch_bra, size_t ch_ket, Ket3& bra, Ket3& ket) const
//{
//  ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
//  ThreeBodyChannel& Tbc_ket = modelspace->GetThreeBodyChannel(ch_ket);
//  size_t ibra = Tbc_bra.GetLocalIndex( bra );
//  size_t iket = Tbc_bra.GetLocalIndex( ket );
//  return GetME_pn( chbra, chket, ibra, iket );
//}






ThreeBodyMEpn::ME_type ThreeBodyMEpn::GetME_pn_PN(int Jab, int Jde, int twoJ, int a, int b, int c, int d, int e, int f) const
{

//  std::cout << __func__ << " begin" << std::endl;
  std::vector<double> recouple_bra;
  std::vector<double> recouple_ket;
  std::vector<size_t> ibra;
  std::vector<size_t> iket;
  size_t ch_bra = GetKetIndex_withRecoupling( Jab, twoJ, a,b,c, ibra, recouple_bra );
  size_t ch_ket = GetKetIndex_withRecoupling( Jde, twoJ, d,e,f, iket, recouple_ket );
  if ( ch_bra != ch_ket) return 0;

//  std::cout << "Start loop" << std::endl;

  double me_out = 0;
  for ( int i=0; i<ibra.size(); i++)
  {
    for (int j=0; j<iket.size(); j++)
    {
//       std::cout << "i,j" << i << " " << j << "  ibra iket " << ibra[i] << " " << iket[j] << "  recouple  " << recouple_bra[i] << " " << recouple_ket[j] << std::endl;
//      if ( std::abs(recouple_bra[i])<1e-8 or std::abs(recouple_ket[j])<1e-8) continue;
//      if ( ibra[i]==size_t(-1) or iket[j]==size_t(-1) )
//      {
//        std::cout << __func__ << "  got a -1 index lookin up " << Jab << " " << Jde << " " << twoJ << "  " << a << " " << b << " " << c << " " << d << " " << e << " " << f << std::endl;
//      }
      me_out += recouple_bra[i] * recouple_ket[j] * GetME_pn_PN_ch( ch_bra, ch_ket, ibra[i], iket[j] );
//      std::cout << "ch_bra,ch_ket " << ch_bra << " " << ch_ket << "  ibra,iket " << ibra[i] << " " << iket[j]
//                << "  recouple bra,ket " << recouple_bra[i] << " " << recouple_ket[j] << "  me_out = " << me_out << std::endl;
    }
  }

  return me_out;
}



// Matrix element is given in isospin formalism. Convert to proton/neutron.
//  
//
void ThreeBodyMEpn::SetME_PN(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V) 
{
//  std::cout << " ENTER " << __func__ <<  "  J " << Jab_in << " " << Jde_in << " " << J2 << "  T " << tab_in << " " << tde_in << " " << twoT << "   ijklm " << i << " "<< j << " " << k << " " << l << " " << m << " " << n <<"    V " << V << std::endl;
  if (i==j and (Jab_in+tab_in)%2==0) return;
  if (l==m and (Jde_in+tde_in)%2==0) return;
  double me_current = GetME( Jab_in, Jde_in, J2, tab_in, tde_in, twoT, i,j,k,l,m,n);
  double me_shift = V - me_current;
//  if ( std::abs(me_shift)<1e-8) return;
  AddToME(Jab_in,Jde_in,J2, tab_in,tde_in,twoT, i,j,k,l,m,n, me_shift);
//  if (i==4 and j==0 and k==0 and l==4 and m==0 and n==0 and Jab_in==0 and Jde_in==0 and J2==1)
//  if (i==2 and j==0 and k==0 and l==2 and m==0 and n==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  {
//    std::cout << "IN " << __func__ << " and tab,tde,twoT=  " << tab_in << " " << tde_in << " " << twoT << "     V = " << V << "  me_current was " << me_current
//              << "  now Get is " <<  GetME( Jab_in, Jde_in, J2, tab_in, tde_in, twoT, i,j,k,l,m,n) << std::endl;
//  }
}

//void ThreeBodyMEpn::SetME(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V) 
void ThreeBodyMEpn::AddToME_PN(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ThreeBodyMEpn::ME_type V) 
{

//  std::cout << "IN " <<__func__ << std::endl;

//  if (i==4 and j==0 and k==0 and l==4 and m==0 and n==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  {
//    std::cout << std::endl << "IN " << __func__ << " and tab,tde,twoT=  " << tab_in << " " << tde_in << " " << twoT << "     V = " << V << std::endl;
//    std::cout << " ENTER " << __func__ <<  "  J " << Jab_in << " " << Jde_in << " " << J2 << "  T " << tab_in << " " << tde_in << " " << twoT << "   ijklm " << i << " "<< j << " " << k << " " << l << " " << m << " " << n <<"    V " << V << std::endl;
//  }

  typedef std::tuple<int,int,int,int,int,int,int,int,int,ThreeBodyMEpn::ME_type>  element_info;
  std::set< element_info > elements_to_set;


  if (i==j and (Jab_in + tab_in)%2==0) return;
  if (l==m and (Jde_in + tde_in)%2==0) return;

  for (int twoTz=-twoT; twoTz<=twoT; twoTz+=2)
  {
  for (int tz2i : {-1,1} )
  {
    for (int tz2j : {-1,1} )
    {
      double clebsch_ij = AngMom::CG(0.5,0.5*tz2i, 0.5,0.5*tz2j, tab_in, 0.5*(tz2i+tz2j) );
      if ( std::abs(clebsch_ij)<1e-7) continue;
      for (int tz2k : {-1,1} )
      {
//        double clebsch_ijk = AngMom::CG(tab_in, 0.5*(tz2i+tz2j), 0.5, 0.5*tz2k, 0.5*twoT, 0.5*(tz2i+tz2j+tz2k) );
        double clebsch_ijk = AngMom::CG(tab_in, 0.5*(tz2i+tz2j), 0.5, 0.5*tz2k, 0.5*twoT, 0.5*twoTz );
        if ( std::abs(clebsch_ijk)<1e-7) continue;
        for (int tz2l : {-1,1} )
        {
          for (int tz2m : {-1,1} )
          {
            int tz2n = tz2i + tz2j + tz2k - tz2l - tz2m;
            if (std::abs(tz2n) != 1) continue;
            double clebsch_lm = AngMom::CG(0.5,0.5*tz2l, 0.5,0.5*tz2m, tde_in, 0.5*(tz2l+tz2m) );
            if ( std::abs(clebsch_lm)<1e-7) continue;
//            double clebsch_lmn = AngMom::CG(tde_in, 0.5*(tz2l+tz2m), 0.5, 0.5*tz2n, 0.5*twoT, 0.5*(tz2l+tz2m+tz2n) );
            double clebsch_lmn = AngMom::CG(tde_in, 0.5*(tz2l+tz2m), 0.5, 0.5*tz2n, 0.5*twoT, 0.5*twoTz );
            if ( std::abs(clebsch_lmn)<1e-7) continue;
            size_t ipn = 2*(i/2) + (tz2i+1)/2;
            size_t jpn = 2*(j/2) + (tz2j+1)/2;
            size_t kpn = 2*(k/2) + (tz2k+1)/2;
            size_t lpn = 2*(l/2) + (tz2l+1)/2;
            size_t mpn = 2*(m/2) + (tz2m+1)/2;
            size_t npn = 2*(n/2) + (tz2n+1)/2;
//            if (i==j and jpn>ipn) continue;
//            if (i==k and kpn>ipn) continue;
//            if (j==k and kpn>jpn) continue;
//            if (l==m and mpn>lpn) continue;
//            if (l==n and npn>lpn) continue;
//            if (m==n and npn>mpn) continue;
//            if (jpn>ipn or kpn>jpn or kpn>ipn or mpn>lpn or npn>mpn or npn>lpn or lpn>ipn) continue;
//            if (ipn==lpn and jpn==mpn and kpn==npn and Jde_in > Jab_in) continue;

//  if (ipn==5 and jpn==1 and kpn==0 and lpn==5 and mpn==1 and npn==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  if (ipn==5 and jpn==1 and kpn==0 and lpn==5 and mpn==1 and npn==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  if (i==4 and j==0 and k==0 and l==4 and m==0 and n==0 and Jab_in==0 and Jde_in==0 and J2==1)
  {
//    std::cout << "Before calling AddToME_pn, the matrix element is " << GetME_pn(1,1,3,5,1,0,5,1,0) << std::endl;
//    std::cout << "ijklmn_iso : " << i << " " << j << " " << k << " " << l << " " << m << " " << n << std::endl;
//    std::cout << "ijklmn_pn : " << ipn << " " << jpn << " " << kpn << " " << lpn << " " << mpn << " " << npn << std::endl;
//    std::cout << " Before calling AddToME_pn, the matrix element is " << GetME_pn(Jab_in,Jde_in,J2,ipn,jpn,kpn,lpn,mpn,npn) << std::endl;
//    std::cout << " clebsch:  " << clebsch_ij << " " << clebsch_ijk << " " << clebsch_lm << " " << clebsch_lmn << "    V = " << V << std::endl;
  }
            double me_old = GetME_pn(Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn);
            ME_type me_new = me_old + clebsch_ij * clebsch_ijk * clebsch_lm * clebsch_lmn * V;
            elements_to_set.insert(std::make_tuple( Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn,  me_new ) );
//            AddToME_pn( Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn,  clebsch_ij * clebsch_ijk * clebsch_lm * clebsch_lmn * V);
//            AddToME_pn( Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn,  clebsch_ij * clebsch_ijk * clebsch_lm * clebsch_lmn * me_shift);
//  if (i==4 and j==0 and k==0 and l==4 and m==0 and n==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  if (ipn==5 and jpn==1 and kpn==0 and lpn==5 and mpn==1 and npn==0 and Jab_in==1 and Jde_in==1 and J2==3)
//  if (i==4 and j==0 and k==0 and l==4 and m==0 and n==0 and Jab_in==0 and Jde_in==0 and J2==1)
  {
//            std::cout << " Called AddToME_pn   " << Jab_in << " " << Jde_in << " " << J2 << " " << ipn << " " << jpn << " " << kpn << " " << lpn << " " << mpn << " " << npn << "   with clebsch  "<< clebsch_ij << " " << clebsch_ijk << " " << clebsch_lm << " " << clebsch_lmn << "     V = " << V << std::endl;
//            std::cout << "Afterwards, <510|V|510> (1,1,3) = " << GetME_pn(1,1,3,5,1,0,5,1,0) << std::endl;
//    std::cout << " Afterwards, the matrix element is " << GetME_pn(Jab_in,Jde_in,J2,ipn,jpn,kpn,lpn,mpn,npn) << std::endl;
  }
          }
        }
      }
    }
  }
  }//for twoTz
  for ( auto& elem : elements_to_set )
  {
    SetME_pn(std::get<0>(elem), std::get<1>(elem), std::get<2>(elem), std::get<3>(elem), std::get<4>(elem),
             std::get<5>(elem), std::get<6>(elem), std::get<7>(elem), std::get<8>(elem), std::get<9>(elem) );

//    std::cout << std::get<0>(elem) << " " << std::get<1>(elem) << " " <<  std::get<2>(elem) << " " << std::get<3>(elem) << " " << std::get<4>(elem)
//              << " "  << std::get<5>(elem) << " " << std::get<6>(elem) << " " <<  std::get<7>(elem) << " " << std::get<8>(elem) << " " << std::get<9>(elem)
//              << std::endl;
//    std::cout << " offending term: " << GetME_pn(1,1,3,3,1,0,3,1,0) << std::endl;
  }

}

// isospin version
ThreeBodyMEpn::ME_type ThreeBodyMEpn::GetME_PN(  int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n ) 
{

//  std::cout << std::endl << " **** IN " << __func__ << "  J,t, ijklmn " << Jab_in << " " << Jde_in << " " << J2 << "  " << tab_in << " " << tde_in << " " << twoT << "   " << i << " " << j << " " << k << " " << l << " " << m << " "<< n << std::endl;
  ThreeBodyMEpn::ME_type me_iso = 0;
  int twoTz = -twoT; // it should be independent of Tz, so we just pick one
  if (i==j and (Jab_in+tab_in)%2==0) return 0;
  if (l==m and (Jde_in+tde_in)%2==0) return 0;
  for (int tz2i : {-1,1} )
  {
    for (int tz2j : {-1,1} )
    {
      double clebsch_ij = AngMom::CG(0.5,0.5*tz2i, 0.5,0.5*tz2j, tab_in, 0.5*(tz2i+tz2j) );
      if ( std::abs(clebsch_ij)<1e-7) continue;
      for (int tz2k : {-1,1} )
      {
        if ( (tz2i + tz2j + tz2k) != twoTz) continue;
        double clebsch_ijk = AngMom::CG(tab_in, 0.5*(tz2i+tz2j), 0.5, 0.5*tz2k, 0.5*twoT, 0.5*twoTz );
        if ( std::abs(clebsch_ijk)<1e-7) continue;
        for (int tz2l : {-1,1} )
        {
          for (int tz2m : {-1,1} )
          {
            int tz2n = twoTz - tz2l - tz2m;
            if (std::abs(tz2n)!=1) continue;
            double clebsch_lm = AngMom::CG(0.5,0.5*tz2l, 0.5,0.5*tz2m, tde_in, 0.5*(tz2l+tz2m) );
            if ( std::abs(clebsch_lm)<1e-7) continue;
            double clebsch_lmn = AngMom::CG(tde_in, 0.5*(tz2l+tz2m), 0.5, 0.5*tz2n, 0.5*twoT, 0.5*twoTz );
            if ( std::abs(clebsch_lmn)<1e-7) continue;
            size_t ipn = 2*(i/2) + (tz2i+1)/2;
            size_t jpn = 2*(j/2) + (tz2j+1)/2;
            size_t kpn = 2*(k/2) + (tz2k+1)/2;
            size_t lpn = 2*(l/2) + (tz2l+1)/2;
            size_t mpn = 2*(m/2) + (tz2m+1)/2;
            size_t npn = 2*(n/2) + (tz2n+1)/2;
            double me_pn = GetME_pn( Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn);
            me_iso += me_pn  * clebsch_ij * clebsch_ijk * clebsch_lm * clebsch_lmn ;
//            std::cout << " ** IN " << __func__ << "   ijklmn_pn: " << ipn << " " << jpn << " " << kpn << " " << lpn << " " << mpn << " " << npn << "   Js " << Jab_in << " " << Jde_in << " " << J2
//                      << "   clebsch: " << clebsch_ij << " " << clebsch_ijk << " " << clebsch_lm << " " << clebsch_lmn << "  mp_pn = " << me_pn << "  me_iso = " << me_iso << std::endl;
          }
        }
      }
    }
  }
  return me_iso;
}




/*
void ThreeBodyMEpn::SetME_isospin5(  int Jab_in, int Jde_in, int J2, int i, int j, int k, int l, int m, int n, std::array<double,5>& isospin_5plet)
{
   // The 5plet stores (t,t,T) = (0,0,1), (0,1,1), (1,0,1), (1,1,1), (1,1,3)
  std::array<double,5> tab = {0,0,1,1,1};
  std::array<double,5> tde = {0,1,0,1,1};
  std::array<double,5> twoT  = {1,1,1,1,3};

//  std::cout << "ENTER" << __func__ <<  " with " << Jab_in << " " << Jde_in << " " << J2 << " " << i << " " << j << " " << k << " " << l << " " << m << " " << n << " ";
//  for (int i=0;i<5;i++) std::cout << " " << isospin_5plet[i];
//  std::cout  << std::endl;
 
  for (int tz2i : {-1,1} )
  {
   for (int tz2j : {-1,1} )
   {
    int tzab = (tz2i+tz2j)/2;
    for (int tz2k : {-1,1} )
    {
     int twoTz = tz2i + tz2j + tz2k;
     for (int tz2l : {-1,1} )
     {
      for (int tz2m : {-1,1} )
      {
       int tzde = (tz2l+tz2m)/2;
       for (int tz2n : {-1,1} )
       {
          if ( (tz2l+tz2m+tz2n) != twoTz ) continue;
            size_t ipn = 2*(i/2) + (tz2i+1)/2;
            size_t jpn = 2*(j/2) + (tz2j+1)/2;
            size_t kpn = 2*(k/2) + (tz2k+1)/2;
            size_t lpn = 2*(l/2) + (tz2l+1)/2;
            size_t mpn = 2*(m/2) + (tz2m+1)/2;
            size_t npn = 2*(n/2) + (tz2n+1)/2;
            double me_pn = 0;
            for (size_t index=0; index<5; index++)
            {
              if (std::abs(isospin_5plet[index])<1e-8) continue;
              if (std::abs(twoTz)>twoT[index]) continue;
              if (std::abs(tzab)>tab[index]) continue;
              if (std::abs(tzde)>tde[index]) continue;
//              std::cout << "( " << 0.5 << " " << 0.5*tz2i << " " << 0.5 << " " << 0.5*tz2j << " "<< tab[index] << " " << tzab << std::endl;
//              std::cout << "( " << 0.5 << " " << 0.5*tz2l << " " << 0.5 << " " << 0.5*tz2m << " "<< tde[index] << " " << tzde << std::endl;
//              std::cout << "( " << tzab << " " << tab[index] << " " << 0.5 << " " << 0.5*tz2k << " "<< twoT[index] << " " << 0.5*twoTz << std::endl;
//              std::cout << "( " << tzde << " " << tde[index] << " " << 0.5 << " " << 0.5*tz2n << " "<< twoT[index] << " " << 0.5*twoTz << std::endl;
              me_pn += AngMom::CG( 0.5,0.5*tz2i, 0.5,0.5*tz2j,  tab[index], tzab) * 
                       AngMom::CG( 0.5,0.5*tz2l, 0.5,0.5*tz2m,  tde[index], 0.5*(tz2l+tz2m)) *
                       AngMom::CG( tab[index],tzab,  0.5,0.5*tz2k, 0.5*twoT[index], 0.5*twoTz) *
                       AngMom::CG( tde[index],tzde,  0.5,0.5*tz2n, 0.5*twoT[index], 0.5*twoTz) * isospin_5plet[index] ;

            if (i==2 and j==2 and k==2 and l==2 and m==0 and n==0 and J2==3)
            {
            std::cout << "   adding " <<  Jab_in << " " << Jde_in << " " << J2 << "   " << ipn << " " << jpn << " " << kpn << " " << lpn << " " << mpn << " " << npn
                      << "  =  " <<  AngMom::CG( 0.5,0.5*tz2i, 0.5,0.5*tz2j,  tab[index], tzab)
                      << " * " << AngMom::CG( 0.5,0.5*tz2l, 0.5,0.5*tz2m,  tde[index], 0.5*(tz2l+tz2m))
                      << " * " << AngMom::CG( tab[index],tzab,  0.5,0.5*tz2k, 0.5*twoT[index], 0.5*twoTz)
                      << " * " << AngMom::CG( tde[index],tzde,  0.5,0.5*tz2n, 0.5*twoT[index], 0.5*twoTz)
                      << "   *   " << isospin_5plet[index]
                      <<  "   =>  " << me_pn << std::endl;
            }


            }

            SetME_pn( Jab_in, Jde_in, J2, ipn,jpn,kpn,lpn,mpn,npn, me_pn);

            if (i==2 and j==2 and k==2 and l==2 and m==0 and n==0 and J2==3)
            {
              std::cout << "   setting " <<  Jab_in << " " << Jde_in << " " << J2 << "   " << ipn << " " << jpn << " " << kpn << " " << lpn << " " << mpn << " " << npn
                        <<  "   =  " << me_pn << std::endl;
              std::cout << " Target ME: " << GetME_pn( 0,1,3,  3,2,2, 2,1,0) << std::endl;
            }

       }
      }
     }
    }
   }
  }

}

*/



void ThreeBodyMEpn::TransformToPN()
{
  std::cout << " " << __func__ << "   changing storage from isospin to proton/neutron" << std::endl;
  Allocate_PN();
  for ( auto& iter : MatEl )
  {
    size_t ch_bra = iter.first[0];
    size_t ch_ket = iter.first[1]; // ch_bra and ch_ket are presumably the same...
//    std::cout << "ch_bra " << ch_bra << std::endl;

    ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch_bra);
//    std::cout << "This channel has " << Tbc.twoJ << " " << Tbc.parity << " " << Tbc.twoTz << std::endl;
    int twoJ = Tbc.twoJ;
    size_t nkets = Tbc.GetNumber3bKets();
    for (size_t ibra=0; ibra<nkets; ibra++)
    {
//      std::cout << " ibra " << ibra << std::endl;
      Ket3& bra = Tbc.GetKet(ibra);
      for (size_t iket=0; iket<=ibra; iket++)
      {
//        std::cout << "  iket " << iket << std::endl;
        Ket3& ket = Tbc.GetKet(iket);
//        std::cout << "   calling GetME_pn  " << twoJ << " " << bra.Jpq << " " << ket.Jpq << " " << bra.p << " " << bra.q << " " << bra.r
//                  << " " << ket.p << " " << ket.q << " " << ket.r << std::endl;
        double me_pn = isospin3BME.GetME_pn( bra.Jpq, ket.Jpq,twoJ,  bra.p, bra.q, bra.r, ket.p, ket.q, ket.r );
//        std::cout << "   read out me_pn = " << me_pn << std::endl;
        SetME_pn_PN_ch( ch_bra, ch_ket, ibra, iket, me_pn);
//        std::cout << "    done setting. set value is " << GetME_pn_PN( bra.Jpq, ket.Jpq, twoJ, bra.p, bra.q, bra.r, ket.p, ket.q, ket.r) << std::endl;
      }
    }
//    std::cout << std::endl;
//    Print( ch_bra, ch_ket);
  }
  // hopefully free up memory?
  std::vector<ThreeBME_type>().swap( isospin3BME.MatEl );
  std::unordered_map<size_t, size_t>().swap( isospin3BME.OrbitIndexHash );
  PN_mode = true;
//  std::cout << "Done" << std::endl;
}




//size_t ThreeBodyMEpn::GetKetIndex_withRecoupling( int twoJ, int Jab_in, size_t a_in, size_t b_in, size_t c_in, std::vector<size_t>& iket , std::vector<double>& recouple) const
size_t ThreeBodyMEpn::GetKetIndex_withRecoupling( int Jab_in, int twoJ, size_t a_in, size_t b_in, size_t c_in, std::vector<size_t>& iket , std::vector<double>& recouple) const
{
//  std::cout << "IN " << __func__ << std::endl;

  int a,b,c;
  int recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);

  int permutation_phase = PermutationPhase( recoupling_case );

  Orbit& oa = modelspace->GetOrbit(a);
  Orbit& ob = modelspace->GetOrbit(b);
  Orbit& oc = modelspace->GetOrbit(c);
  int parity = ( oa.l + ob.l + oc.l )%2;
  int twoTz = ( oa.tz2 + ob.tz2 + oc.tz2 );
//  std::cout << "Call modelspace->GetThreeBodyChannelIndex " << twoJ << " " << parity << " " << twoTz << std::endl;
  int ch = modelspace->GetThreeBodyChannelIndex( twoJ, parity, twoTz );

//  std::cout << "Before staring the loop. ch = " << ch << std::endl;
  if (ch < 0 ) return ch;
//  auto Tbc = modelspace->GetThreeBodyChannel(ch);
//  std::cout << " this points to the channel with " << Tbc.twoJ << " " << Tbc.parity << " " << Tbc.twoTz << std::endl;
//  std::cout << "  which should be channel " << modelspace->GetThreeBodyChannelIndex( Tbc.twoJ, Tbc.parity, Tbc.twoTz) << std::endl;
//  std::cout << "    IN " << __func__ << "  recoupling_case = " << recoupling_case << std::endl;

  int Jab_min = std::max( std::abs(oa.j2-ob.j2), std::abs(oc.j2-twoJ) )/2;
  int Jab_max = std::min( oa.j2+ob.j2, oc.j2+twoJ)/2;

  if (  ( a_in==b_in and (Jab_in%2)>0 ) 
//     or ( a_in==b_in and a_in==c_in and oa.j2<3)
     or ( a_in==b_in and (Jab_in>oa.j2-1) )
     or ( a_in==b_in and a_in==c_in and ( twoJ > (3*oa.j2-3)) )
     or ( twoJ==(oa.j2+ob.j2+oc.j2) and (a==b or b==c) ) )
  {
    Jab_max = Jab_min-1;
  }

  double ja = oa.j2*0.5;
  double jb = ob.j2*0.5;
  double jc = oc.j2*0.5;
//  if (recoupling_case==ABC or recoupling_case==BAC)
//  {
//    Jab_min = Jab_max = Jab_in;
//  }

  // Loop over possible values of Jab with the new ordering of a,b,c and
  // fill a vector of the index of where each of those |a,b,c,Jab> states lives
  // as well as the recouplng coefficient.
//  std::cout << "Start loop Jab min,max = " << Jab_min << " " << Jab_max << std::endl;
//  std::cout << "================================================================" << std::endl;
  for (int Jab=Jab_min; Jab<=Jab_max; Jab++)
  {
//    std::cout << "|||| a b c Jab = " << a << " " << b <<" " << c << " " << Jab << std::endl;
//    std::cout << " ch = " << ch << std::endl;
//    std::cout << "modelspace has  "<< modelspace->GetNumberThreeBodyChannels() << "  3b channels" << std::endl;
//    std::cout << " Get3bchannel = " << std::endl;
//    modelspace->GetThreeBodyChannel(ch);
//    std::cout << "That was ok. Now try getting the local index" << std::endl;
//    size_t localindex = modelspace->GetThreeBodyChannel(ch).GetLocalIndex(a,b,c,Jab);
//    std::cout << "that was ok. I got a localindex of " << localindex << std::endl;
//    std::cout << " now pusihg back " << modelspace->GetThreeBodyChannel(ch).GetLocalIndex(a,b,c,Jab) << std::endl;
    if (a==b and (Jab%2)>0) continue;
//    std::cout << "  looking for index,  ch = " << ch << "  a,b,c,Jab = " << a << " " << b << " " << c << " " << Jab << std::endl;
//    ThreeBodyChannel& Tbc = modelspace->GetThreeBodyChannel(ch);
//    std::cout << " This Tbc has " << Tbc.twoJ << " " << Tbc.parity << " " << Tbc.twoTz << std::endl;
    size_t index = modelspace->GetThreeBodyChannel(ch).GetLocalIndex( a,b,c,Jab );
    if ( index == size_t(-1) ) continue;

//    std::cout << "  computing coefficient" << std::endl;
    double coefficient =  permutation_phase * RecouplingCoefficient( recoupling_case, ja,jb,jc, Jab_in, Jab, twoJ );
    if (std::abs(coefficient)<1e-10) continue;

    iket.push_back( index );
    recouple.push_back( coefficient );
  }
//  std::cout << "done" << std::endl;
//  std::cout << "returning vectors with sizes " << iket.size() << " " << recouple.size() << std::endl;

  return ch;

}






// Define some constants for the various permutations of three indices
// for use in RecouplingCoefficient and SortOrbits
const int ThreeBodyMEpn::ABC = 0;
const int ThreeBodyMEpn::BCA = 1;
const int ThreeBodyMEpn::CAB = 2;
const int ThreeBodyMEpn::ACB = 3;
const int ThreeBodyMEpn::BAC = 4;
const int ThreeBodyMEpn::CBA = 5;

/// If we make an odd permutation of indices, we get a fermionic minus sign.
int ThreeBodyMEpn::PermutationPhase( int recoupling_case ) const
{
  return  (recoupling_case < 3) ?  +1 : -1;
}


//*******************************************************************
/// Coefficients for recoupling three body matrix elements.
/// Note that this does not include the -1 factor for an odd
/// permutation of fermionic operators. That is handled in ThreeBodyMEpn::PermutationPhase.
/// Here, we just only with the angular momentum / isospin recoupling factors
//*******************************************************************
double ThreeBodyMEpn::RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const
{
//   std::cout << "IN " << __func__ << std::endl;
   if ( std::abs(int(ja-jb))>Jab  or int(ja+jb)<Jab) return 0;
   if ( std::abs(int(jc-twoJ/2.))>Jab  or int(jc+twoJ/2.)<Jab) return 0;
   switch (recoupling_case)
   {
    case ABC: return Jab==Jab_in ? 1 : 0;
    case BCA: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    case CAB: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case ACB: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, twoJ/2., Jab_in);
    case BAC: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case CBA: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, twoJ/2., Jab_in);
    default: return 0;
    }
}




//*******************************************************************
/// Rearrange orbits (abc) so that a>=b>=c
/// and return an int which reflects the required reshuffling
/// - 0: (abc)_in -> (abc)
/// - 1: (bca)_in -> (abc)
/// - 2: (cab)_in -> (abc)
/// - 3: (acb)_in -> (abc)  -- odd permutation
/// - 4: (bac)_in -> (abc)  -- odd permutation
/// - 5: (cba)_in -> (abc)  -- odd permutation
//*******************************************************************
int ThreeBodyMEpn::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
{
//   a_in -= a_in%2;
//   b_in -= b_in%2;
//   c_in -= c_in%2;
   a=a_in;
   b=b_in;
   c=c_in;
   if (a<b)  std::swap(a,b);
   if (b<c)  std::swap(b,c);
   if (a<b)  std::swap(a,b);

   int recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? ABC : ACB;
   else if (a_in==b)  recoupling_case = (b_in==a) ? BAC : BCA;
   else               recoupling_case = (b_in==a) ? CAB : CBA;

   return recoupling_case;
}
  
size_t ThreeBodyMEpn::size()
{
  size_t nelem = 0;
  for (auto iter : MatEl ) nelem += iter.second.size();
  return nelem;

}

void ThreeBodyMEpn::Erase()
{
  for (auto iter : MatEl ) iter.second.zeros();
}


double ThreeBodyMEpn::Norm() const
{
  double norm = 0;
  for ( auto iter : MatEl ) norm += iter.second.Norm();
  return norm;
}


void ThreeBodyMEpn::Print(size_t ch_bra, size_t ch_ket)
{
  ThreeBodyChannel& Tbc_bra = modelspace->GetThreeBodyChannel(ch_bra);
  std::cout << "Channel " << ch_bra << "  J p Tz = " << Tbc_bra.twoJ << " " << Tbc_bra.parity << " " << Tbc_bra.twoTz << std::endl;
  std::cout << "Kets: ";
  for ( size_t iket : Tbc_bra.KetList )
  {
    Ket3& ket = modelspace->GetKet3(iket);
    std:: cout << "(" << ket.p << "," << ket.q << "," << ket.r << ";" << ket.Jpq << ")  ";
  }
  std::cout << std::endl;
  std::cout << MatEl.at({ch_bra,ch_ket}).FullMatrix() << std::endl << std::endl;
}


void ThreeBodyMEpn::PrintAll()
{
  for (auto& iter : MatEl )
  {
    if ( iter.second.size() > 0 )
    Print( iter.first[0], iter.first[1] );
  }
}



ThreeBodyMEpn& ThreeBodyMEpn::operator*=(const double rhs)
{
  for (auto iter : MatEl ) iter.second *= rhs;
  return *this;
}

ThreeBodyMEpn& ThreeBodyMEpn::operator+=(const ThreeBodyMEpn& rhs)
{
  for (auto iter : MatEl ) iter.second += rhs.MatEl.at(iter.first);
  return *this;
}

ThreeBodyMEpn& ThreeBodyMEpn::operator-=(const ThreeBodyMEpn& rhs)
{
  for (auto iter : MatEl ) iter.second -= rhs.MatEl.at(iter.first);
  return *this;
}

