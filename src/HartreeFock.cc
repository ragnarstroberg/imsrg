
#include "HartreeFock.hh"
#include "ModelSpace.hh"

#ifndef SQRT2
  #define SQRT2 1.4142135623730950488
#endif


using namespace std;

HartreeFock::HartreeFock(Operator& hbare) : Hbare(hbare)
{
   tolerance = 1e-10;
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C         = arma::mat(norbits,norbits,arma::fill::eye);
   Vij       = arma::mat(norbits,norbits,arma::fill::zeros);
   V3ij       = arma::mat(norbits,norbits,arma::fill::zeros);
   F         = arma::mat(norbits,norbits);
   Vmon      = arma::mat(nKets,nKets);
   Vmon_exch = arma::mat(nKets,nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);
   t = Hbare.OneBody;
   energies = t.diag();
   BuildMonopoleV();
   if (hbare.GetParticleRank()>2)
   {
      BuildMonopoleV3();
   }
   UpdateDensityMatrix();
   UpdateF();
}


//*********************************************************************
/// Diagonalize and update the Fock matrix until convergence.
/// Then, call ReorderCoefficients() to make sure the index
/// ordering and phases are preserved in the transformation
/// from the original basis to the Hatree-Fock basis.
//*********************************************************************
void HartreeFock::Solve()
{
   iterations = 0; // counter so we don't go on forever
   int maxiter = 1000;

   for (iterations=0; iterations<maxiter; ++iterations)
   {
      Diagonalize();  // Diagonalize the Fock matrix
      ReorderCoefficients(); // Reorder here so we can identify the hole orbits.
      UpdateDensityMatrix();  // 1 body density matrix, used in UpdateF()
      UpdateF();  // Update the Fock matrix
      if ( CheckConvergence() ) break;
   }

   if (iterations==maxiter)
   {
      cerr << "Warning: Hartree-Fock calculation didn't converge after " << maxiter << " iterations." << endl;
   }
   else
   {
      cout << "HF converged after " << iterations << " iterations. " << endl;
   }

   ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
   CalcEHF();
}


//*********************************************************************
/// Calculate the HF energy.
/// \f{eqnarray*} E_{HF} &=& \sum_{\alpha} t_{\alpha\alpha} 
///                    + \frac{1}{2}\sum_{\alpha\beta} V_{\alpha\beta\alpha\beta}
///                    + \frac{1}{6}\sum_{\alpha\beta\gamma} V_{\alpha\beta\gamma\alpha\beta\gamma} \\
///    &=& \sum_{ij} (2j_i+1) \rho_{ij} ( t_{ij} +\frac{1}{2}\tilde{V}^{(2)}_{ij} + \frac{1}{6}\tilde{V}^{(3)}_{ij} )
/// \f}
/// Where the matrices \f{eqnarray*}
///  \tilde{V}^{(2)}_{ij} &=& \sum_{ab} \rho_{ab}\bar{V}^{(2)}_{iajb} \\
///  \tilde{V}^{(2)}_{ij} &=& \sum_{abcd} \rho_{ab}\rho_{cd} \bar{V}^{(3)}_{iacjbd} \\
/// \f}
/// have already been calculated by UpdateF().
//*********************************************************************
void HartreeFock::CalcEHF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   EHF = 0;
   double e3hf = 0;
   int norbits = ms->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      Orbit& oi = ms->GetOrbit(i);
      int jfactor = oi.j2 +1;
      for (int j : ms->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}))
      {
         EHF +=  rho(i,j) * jfactor * (t(i,j)+0.5*Vij(i,j)+1./6*V3ij(i,j));
         e3hf +=  rho(i,j) * jfactor * (1./6*V3ij(i,j));
      }
   }
   cout << "e3hf = " << e3hf << endl;
}


//*********************************************************************
/// [See Suhonen eq. 4.85]
/// Diagonalize the fock matrix \f$ <a|F|b> \f$ and put the
/// eigenvectors in \f$C(i,\alpha) = <i|\alpha> \f$
/// and eigenvalues in the vector energies.
/// Save the last vector of energies to check
/// for convergence.
/// Different channels are diagonalized independently.
/// This guarantees that J,Tz, and \f$ \pi \f$ remain good. 
//*********************************************************************
void HartreeFock::Diagonalize()
{
   prev_energies = energies;
   for ( auto& it : Hbare.GetModelSpace()->OneBodyChannels )
   {
      arma::uvec orbvec(it.second);
      arma::mat F_ch = F.submat(orbvec,orbvec);
      arma::mat C_ch;
      arma::vec E_ch;
      bool success = false;
      int diag_tries = 0;
      while ( not success)
      {
         success = arma::eig_sym(E_ch, C_ch, F_ch);
         ++diag_tries;
         if (diag_tries > 5)
         {
           cout << "Hartree-Fock: Failed to diagonalize the submatrix " 
                << " on iteration # " << iterations << ". The submatrix looks like:" << endl;
           F_ch.print();
         }
      }
      // Update the full overlap matrix C and energy vector
      energies(orbvec) = E_ch;
      C.submat(orbvec,orbvec) = C_ch;

   }
}


//*********************************************************************
/// Construct an unnormalized two-body monopole interaction
/// \f[ \langle ab | \bar{V}^{(2)} | cd \rangle = \sqrt{(1+\delta_{ab})(1+\delta_{cd})} \sum_{J} (2J+1) \langle ab | V^{(2)} | cd \rangle_{J} \f]
/// This method utilizes the operator method  Operator::GetTBMEmonopole() 
///
//*********************************************************************
void HartreeFock::BuildMonopoleV()
{
   Vmon.zeros();
   Vmon_exch.zeros();
   int nKets = Hbare.GetModelSpace()->GetNumberKets();
   for (int ibra=0;ibra<nKets;++ibra)
   {
      Ket & bra = Hbare.GetModelSpace()->GetKet(ibra);
      int a = bra.p;
      int b = bra.q;
      Orbit & oa = Hbare.GetModelSpace()->GetOrbit(a);
      Orbit & ob = Hbare.GetModelSpace()->GetOrbit(b);
      double norm = (oa.j2+1)*(ob.j2+1);
      for (int iket=ibra;iket<nKets;++iket)
      {
         Ket & ket = Hbare.GetModelSpace()->GetKet(iket);
         int c = ket.p;
         int d = ket.q;
         Vmon(ibra,iket)       = Hbare.TwoBody.GetTBMEmonopole(a,b,c,d)*norm;
         Vmon_exch(ibra,iket)  = Hbare.TwoBody.GetTBMEmonopole(a,b,d,c)*norm;
      }
   }
   Vmon = arma::symmatu(Vmon);
   Vmon_exch = arma::symmatu(Vmon_exch);
}



//*********************************************************************
/// Construct an unnormalized three-body monopole interaction
/// \f[ \langle iab | \bar{V}^{(3)} | jcd \rangle =
///     \sum\limits_{J,J_{12}}\sum_{Tt_{12}}(2J+1)(2T+1) 
///       \langle (ia)J_{12}t_{12};b JT| V^{(3)} | (jc)J_{12}t_{12}; d JT\rangle \f]
///
//*********************************************************************
void HartreeFock::BuildMonopoleV3()
{
  // First, allocate
  ModelSpace * modelspace = Hbare.GetModelSpace();
  int norbits = modelspace->GetNumberOrbits();
  vector< pair<const array<int,6>,double>*> entries;
  for (int a=0; a<norbits; ++a)
  {
    Orbit& oa = modelspace->GetOrbit(a);
    int ea = 2*oa.n + oa.l;
    for (int b : modelspace->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
    {
      Orbit& ob = modelspace->GetOrbit(b);
      int eb = 2*ob.n + ob.l;
 
        for (int c=0; c<norbits; ++c)
        {
          Orbit& oc = modelspace->GetOrbit(c);
          int ec = 2*oc.n + oc.l;
          for (int d : modelspace->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
          {
            Orbit& od = modelspace->GetOrbit(d);
            int ed = 2*od.n + od.l;
 
            for (int i=0; i<norbits; ++i)
            {
              Orbit& oi = modelspace->GetOrbit(i);
              int ei = 2*oi.n + oi.l;
              if ( ea+ec+ei > Hbare.E3max ) continue;
              for (int j : modelspace->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
              {
                Orbit& oj = modelspace->GetOrbit(j);
                int ej = 2*oj.n + oj.l;
                if ( eb+ed+ej > Hbare.E3max ) continue;
                if ( (oi.l+oa.l+ob.l+oj.l+oc.l+od.l)%2 >0) continue;
                Vmon3[{a,c,i,b,d,j}] = 0;
                entries.push_back(& (*Vmon3.find({a,c,i,b,d,j})) );
              }
            }
          }
        }
      }
    }

   // the calculation takes longer, so parallelize that part
//   #pragma omp parallel for
   for (auto it=entries.begin(); it<entries.end(); ++it)
   {
      int a=(*it)->first[0];
      int c=(*it)->first[1];
      int i=(*it)->first[2];
      int b=(*it)->first[3];
      int d=(*it)->first[4];
      int j=(*it)->first[5];
      double& v = (*it)->second;
      int j2a = modelspace->GetOrbit(a).j2;
      int j2c = modelspace->GetOrbit(c).j2;
      int j2i = modelspace->GetOrbit(i).j2;
      int j2b = modelspace->GetOrbit(b).j2;
      int j2d = modelspace->GetOrbit(d).j2;
      int j2j = modelspace->GetOrbit(j).j2;
 
      int j2min = max( abs(j2a-j2c), abs(j2b-j2d) )/2;
      int j2max = min (j2a+j2c, j2b+j2d)/2;
      for (int j2=j2min; j2<=j2max; ++j2)
      {
        int Jmin = max( abs(2*j2-j2i), abs(2*j2-j2j) );
        int Jmax = 2*j2 + min(j2i, j2j);
        for (int J=Jmin; J<=Jmax; J+=2)
        {
           v += Hbare.GetThreeBodyME_pn(j2,j2,J,a,c,i,b,d,j) * (J+1);
        }
      }
   }
}



//*********************************************************************
/// one-body density matrix 
/// \f$ <i|\rho|j> = \sum\limits_{\beta} n_{\beta} <i|\beta> <\beta|j> \f$
/// where \f$n_{\beta} \f$ ensureas that beta runs over HF orbits in
/// the core (i.e. below the fermi surface)
//*********************************************************************
void HartreeFock::UpdateDensityMatrix()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   // Pcore is a projector onto core orbits.
   arma::mat Pcore = arma::mat(norbits,norbits,arma::fill::zeros);
   for (auto& beta : ms->holes)
   {
       Pcore(beta,beta) = 1;
   }

   rho = C * Pcore * C.t();
}


//*********************************************************************
///  [See Suhonen eq 4.85]
/// \f[ F_{ij} = t_{ij}  +  \frac{1}{2j_i+1}\sum_{ab} \rho_{ab} \bar{V}^{(2)}_{iajb}
///               + \frac{1}{2(2j_i+1)}\sum_{abcd}\rho_{ab} \rho_{cd} \bar{V}^{(3)}_{iacjbd}  \f]
/// * \f$ F \f$ is the Fock matrix, to be diagonalized
/// * \f$ t \f$ is the kinetic energy
/// * \f$\rho\f$ is the density matrix defined in UpdateDensityMatrix()
/// * \f$ \bar{V}^{(2)} \f$ is the monopole component of the 2-body interaction defined in BuildMonopoleV().
/// * \f$ \bar{V}^{(3)} \f$ is the monopole component of the 3-body interaction devined in BuildMonopoleV3().
//*********************************************************************
void HartreeFock::UpdateF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int bra, ket;
   Vij.zeros();
   V3ij.zeros();

//   #pragma omp parallel for
   for (int i=0;i<norbits;i++)
   {
      Orbit& oi = ms->GetOrbit(i);
      for (int j : ms->OneBodyChannels.at({oi.l,oi.j2,oi.tz2}) )
      {
         Orbit& oj = ms->GetOrbit(j);
         if (j<i) continue;
         for (int a=0;a<norbits;++a)
         {
            Orbit& oa = ms->GetOrbit(a);
            bra = ms->GetKetIndex(min(i,a),max(i,a));
            for (int b : ms->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}) )
            {
               Orbit& ob = ms->GetOrbit(b);
               ket = ms->GetKetIndex(min(j,b),max(j,b));
               // 2body term <ai|V|bj>
               if (a>i xor b>j)
                  Vij(i,j) += rho(a,b)*Vmon_exch(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
               else
                  Vij(i,j) += rho(a,b)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>



               if (Hbare.GetParticleRank()<3) continue;


               // 3body term  <aci|V|bdj> <a|rho|b> <c|rho|d>
               for (int c=0;c<norbits;++c)
               {
                 Orbit& oc = ms->GetOrbit(c);
                 if ( 2*(oi.n+oa.n+oc.n)+oi.l+oa.l+oc.l > Hbare.E3max ) continue;
                 for (int d : ms->OneBodyChannels.at({oc.l,oc.j2,oc.tz2}) )
                 {
                   Orbit& od = ms->GetOrbit(d);
                   if ( 2*(oj.n+ob.n+od.n)+oj.l+ob.l+od.l > Hbare.E3max ) continue;
                   if ( (oi.l+oa.l+oc.l+oj.l+ob.l+od.l)%2 >0) continue;

                   V3ij(i,j) += rho(a,b) * rho(c,d) * Vmon3[{a,c,i,b,d,j}];
                 }
               }
           }
         }
         Vij(i,j) /= (oi.j2+1);
//         Vij(j,i) = Vij(i,j);  // Hermitian & real => symmetric
         V3ij(i,j) /= (oi.j2+1); //
//         V3ij(j,i) = V3ij(i,j);  // Hermitian & real => symmetric
      }
   }
   Vij = arma::symmatu(Vij);
   V3ij = arma::symmatu(V3ij);
   F = t + Vij + 0.5*V3ij;
}



//********************************************************
/// Check for convergence using difference in s.p. energies
/// between iterations.
/// Converged when
/// \f[ \delta_{e} \equiv \sqrt{ \sum_{i}(e_{i}^{(n)}-e_{i}^{(n-1)})^2} < \textrm{tolerance} \f]
/// where \f$ e_{i}^{(n)} \f$ is the \f$ i \f$th eigenvalue of the Fock matrix after \f$ n \f$ iterations.
///
//********************************************************
bool HartreeFock::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   double ediff = sqrt(arma::dot(de,de)) / energies.size();
   return (ediff < tolerance);
}


//**********************************************************************
/// Eigenvectors/values come out of the diagonalization energy-ordered.
/// We want them ordered corresponding to the input ordering, i.e. we want
/// the matrix C to be maximal and positive along the diagonal.
/// For a 3x3 matrix this would be something like
/// \f[
/// \left(
/// \begin{array}{rrr}
///   -0.8 & 0.2 & -0.6 \\
///   -0.3 & 0.3 &  0.9 \\
///    0.2 & 0.9 & -0.4 \\
/// \end{array}\right)
/// \rightarrow
/// \left(\begin{array}{rrr} 
///    0.8 & -0.6 & 0.2  \\
///    0.3 &  0.9 & 0.3  \\
///   -0.2 & -0.4 & 0.9  \\
/// \end{array}\right)
/// \f]
//**********************************************************************
void HartreeFock::ReorderCoefficients()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   int nswaps = 10;

   // first, reorder them so we still know the quantum numbers
   while (nswaps>0) // loop until we don't have to make any more swaps
   {
     nswaps = 0;
     for (int i=0;i<C.n_rows;++i) // loop through rows -> original basis states
     {
        arma::rowvec row = C.row(i);
        arma::uword imax; 
        double maxval = abs(row).max(imax);
        if (imax == i) continue;
        if (maxval > abs(C(imax,imax) ))
        {
           C.swap_cols(i,imax);
           energies.swap_rows(i,imax);
           ++nswaps;
        }
     }
   }
   // Make sure the diagonal terms are positive. (For easier comparison).
   for (int i=0;i<C.n_rows;++i) // loop through original basis states
   {
      if ( C(i,i) < 0 )
      {
         C.col(i) *= -1;
      }
   }

}



//**************************************************************************
/// Takes in an operator expressed in the basis of the original Hamiltonian,
/// and returns that operator in the Hartree-Fock basis.
/// \f[ t_{HF} = C^{\dagger} t_{HO} C \f]
/// \f[ V_{HF}^{J} = D^{\dagger} V^{J}_{HO} D \f]
/// The matrix \f$ D \f$ is defined as
/// \f[ D_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }  C_{a\alpha} C_{b\beta} \f]
/// The factor in the square root is due to the fact that we're using normalized TBME's.
/// Since only kets with \f$ a\leq b\f$ are stored, we can use the antisymmetry of the TBME's and define
/// \f[ D(J)_{ab\alpha\beta} \equiv \sqrt{ \frac{1+\delta_{ab}} {1+\delta_{\alpha\beta}} }
///      \left( C_{a\alpha} C_{b\beta} -(1-\delta_{ab})(-1)^{j_a+j_b-J} C_{b\alpha}C_{a\beta}\right) \f]
///
//**************************************************************************
Operator HartreeFock::TransformToHFBasis( Operator& OpIn)
{
   Operator OpHF = OpIn;
   OpHF.Erase();

   if (OpIn.rank_J + OpIn.rank_T + OpIn.parity > 0)
   {
      cout << "Warning, HF transformation hasn't been checked for tensor operators" << endl;
   }

   cout << "Transform one body" << endl;
   // Easy part:
   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>
   // where |i> is the original basis and |a> is the HF basis
   OpHF.OneBody = C.t() * OpIn.OneBody * C;


   cout << "Transform two body" << endl;
   // Medium-difficulty part:
   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.
   int nchan = OpIn.GetModelSpace()->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = OpIn.GetModelSpace()->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;

      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>

      // loop over all possible original basis configurations <pq| in this J,p,Tz channel.
      // and all possible HF configurations |p'q'> in this J,p,Tz channel                                    
      // bra is in the original basis, ket is in the HF basis                                              
      // i and j are the indices of the matrix D for this channel                    
      for (int i=0; i<npq; ++i)    
      {
         Ket & bra = tbc.GetKet(i);   
         for (int j=0; j<npq; ++j)    
         {
            Ket & ket = tbc.GetKet(j); 
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q);
            if (bra.p!=bra.q)
            {
               D(i,j) += C(bra.q,ket.p) * C(bra.p,ket.q) * bra.Phase(tbc.J);
            }
            if (bra.p==bra.q)    D(i,j) *= SQRT2;
            if (ket.p==ket.q)    D(i,j) /= SQRT2;
         }
      }

//     auto& IN  =  OpIn.TwoBody[ch].at(ch);
//     auto& OUT =  OpHF.TwoBody[ch].at(ch);
     auto& IN  =  OpIn.TwoBody.GetMatrix(ch);
     auto& OUT =  OpHF.TwoBody.GetMatrix(ch);
     OUT  =    D.t() * IN * D;
   }

   return OpHF;
}


/// Returns the normal-ordered Hamiltonian in the Hartree-Fock basis
/// \f[ E_0 = E_{HF} \f]
/// \f[ f_{ij} = C^{\dagger} F_{ij} C \f]
/// \f[ \Gamma^{J}_{ijkl} = D^{\dagger} \left[V^{J}_{ijkl} + \sum_{ab}\sum_{J_3}(2J_{3}+1)\rho_{ab}V^{JJJ_{3}}_{ijaklb}\right] D \f]
/// Where the matrix \f$ D\f$ is the same as the one defined in TransformToHFBasis().
///
Operator HartreeFock::GetNormalOrderedH()
{
   ModelSpace *modelspace = Hbare.GetModelSpace();
   Operator HNO = Operator(*modelspace,0,0,0,2);
   HNO.ZeroBody = EHF;
   HNO.OneBody = C.t() * F * C;

   int nchan = modelspace->GetNumberTwoBodyChannels();
   int norb = modelspace->GetNumberOrbits();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      int J = tbc.J;
      int npq = tbc.GetNumberKets();

      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat V3NO  = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>

//      #pragma omp parallel for schedule(dynamic,1)
      for (int i=0; i<npq; ++i)    
      {
         Ket & bra = tbc.GetKet(i);
         for (int j=0; j<npq; ++j)
         {
            Ket & ket = tbc.GetKet(j); 
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q);
            if (bra.p!=bra.q)
            {
               D(i,j) += C(bra.q,ket.p) * C(bra.p,ket.q) * bra.Phase(J);
            }
            if (bra.p==bra.q)    D(i,j) *= SQRT2;
            if (ket.p==ket.q)    D(i,j) /= SQRT2;

            if (Hbare.GetParticleRank()<3) continue;
            if (i>j) continue;
            for (int a=0; a<norb; ++a)
            {
              Orbit & oa = modelspace->GetOrbit(a);
              if ( 2*oa.n+oa.l+bra.E2 > Hbare.GetE3max() ) continue;
              for (int b : modelspace->OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
              {
                Orbit & ob = modelspace->GetOrbit(b);
                if ( 2*ob.n+ob.l+ket.E2 > Hbare.GetE3max() ) continue;
                int J3min = abs(2*J-oa.j2);
                int J3max = 2*J + oa.j2;
                for (int J3=J3min; J3<=J3max; J3+=2)
                {
                  V3NO(i,j) += rho(a,b) * (J3+1) * Hbare.GetThreeBodyME_pn(J,J,J3,bra.p,bra.q,a,ket.p,ket.q,b);
                }
              }
            }
            V3NO(i,j) /= (2*J+1);
            V3NO(j,i) = V3NO(i,j);
         }
      }

     auto& V2  =  Hbare.TwoBody.GetMatrix(ch);
     auto& OUT =  HNO.TwoBody.GetMatrix(ch);
     OUT  =    D.t() * (V2 + V3NO) * D;
   }
   
   return HNO;

}




