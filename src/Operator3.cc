
#include "Operator3.hh"


/////////////////// CONSTRUCTORS /////////////////////////////////////////
Operator3::Operator3()
{
   modelspace = NULL;
   nChannels = 0;
   hermitian = true;
   antihermitian = false;
}

Operator3::Operator3(ModelSpace& ms) // Create a zero-valued operator in a given model space
{
  modelspace = &ms;
  hermitian = true;
  antihermitian = false;
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

/////////// COPY METHOD //////////////////////////
void Operator3::Copy(const Operator3& op)
{
   modelspace    = op.modelspace;
   nChannels     = op.nChannels;
   hermitian     = op.hermitian;
   antihermitian = op.antihermitian;
   ZeroBody      = op.ZeroBody;
   OneBody       = op.OneBody;
   TwoBody       = op.TwoBody;
   ThreeBody     = op.ThreeBody;
}

void AllocateThreeBody()
{

}


//double GetThreeBodyME(int J, int Jprime, int K2, int Tz2, int parity, int i, int j, int k, int l, int m, int n)
double GetThreeBodyME(int Jab, int Jde, int J, int tab, int tde, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
   // reorder so a>=b>=c and d>=e>=f
   int a = max( max(a_in,b_in), c_in);
   int b = max( min(a_in,b_in), min(max(a_in,b_in),c_in) );
   int c = min( min(a_in,b_in), c_in);
   int d = max( max(d_in,e_in), f_in);
   int e = max( min(d_in,e_in), min(max(d_in,e_in),f_in) );
   int f = min( min(d_in,e_in), f_in);


   // need to figure out phase factor associated with the reordering
   // also need to map  a => a/2

   int isospin_index = tde + 2*tab + 2*(T2-1); // ranges from 0-4
   int J_index = J + J_max*Jab + J_max*Jab_max*Jde;
   int orbit_index_bra =  a*(a+1)*(a+2)/6 + b*(b+1)/2 + c;
   int orbit_index_ket =  d*(d+1)*(d+2)/6 + e*(e+1)/2 + f;
   long int orb_indx = (orbit_index_bra * orbit_index_bra+1)/2 + orbit_index_ket;
   // this index can get huge, and lots of the terms are zero. need to map this index
   // to a more densely packed index.


}



//
//   Gamma(2)^J_ijkl = V(2)^J_ijkl + Sum_a n_a  Sum_K (2K+1)/(2J+1) V(3)^JJK_ijakla
//
//
Operator Operator3::DoNormalOrdering3()
{
   opNO3 = Operator3(*modelspace);
//   #pragma omp parallel for
   for (int ch=0;ch<nChannels;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      for (int ibra=0; ibra<tbc.GetNumberKets(); ++ibra)
      {
         Ket & bra = tbc.GetKet(ibra);
         int i = bra.p;
         int j = bra.q;
         for (int iket=ibra; iket<tbc.GetNumberKets(); ++iket)
         {
            Ket & ket = tbc.GetKet(iket);
            int k = ket.p;
            int l = ket.q;
            for (int& a : modelspace->holes)
            {
               Orbit & oa = modelspace->GetOrbit(a);
               int Tz2 = 2*tbc.Tz + oa.tz2;
               int kmin2 = abs(2*J-oa.j2);
               int kmax2 = 2*J+oa.j2;
               int parity = tbc.parity * (oa.l%2);
               for (int K2=kmin; K2<=kmax; ++K2)
               {
                  #pragma omp critical
                  opNO3.TwoBody[ch](ibra,iket) += (K2+1) * GetThreeBodyME(J,Jprime,K2,Tz2,parity,i,j,a,k,l,a);
               }
            }
            opNO3.TwoBody[ch](ibra,iket) /= (2*J+1);
         }
      }
   }
   Operator opNO2 = opNO3.DoNormalOrdering();
   opNO2.ScaleZeroBody(1./3.);
   opNO2.ScaleOneBody(1./2.);
   opNO2 += DoNormalOrdering();
   return opNO2;

}


