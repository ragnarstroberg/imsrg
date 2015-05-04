#include "ThreeBodyME.hh"
#include "AngMom.hh"


ThreeBodyME::~ThreeBodyME()
{}

ThreeBodyME::ThreeBodyME()
: modelspace(NULL),E3max(0)
{
//   cout << "Default ThreeBodyME constructor" << endl;
// MatEl.resize(0);
}

ThreeBodyME::ThreeBodyME(ModelSpace* ms)
: modelspace(ms), E3max(ms->N3max)
{}

ThreeBodyME::ThreeBodyME(ModelSpace* ms, int e3max)
: modelspace(ms),E3max(e3max)
{}

// Confusing nomenclature: J2 means 2 times the total J of the three body system
// this should most likely be done with vectors?


void ThreeBodyME::Allocate()
{
  E3max = modelspace->GetN3max();
  cout << "Begin AllocateThreeBody() with E3max = " << E3max << endl;
  int norbits = modelspace->GetNumberOrbits();
  int nvectors = 0;
  int total_dimension = 0;
  int total_dimension_with_Jcut = 0;
  int lmax = 500*norbits; // maybe do something with this later...

  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   int ea = 2*oa.n+oa.l;
   if (ea>E3max) break;
   vector<vector<vector<vector<vector<vector<ThreeBME_type>>>>>> vecb;
   for (int b=0; b<=a; b+=2)
   {
     if (oa.l > lmax) break;
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;

     int Jab_min = abs(oa.j2-ob.j2)/2;
     int Jab_max = (oa.j2+ob.j2)/2;
     vector<vector<vector<vector<vector<ThreeBME_type>>>>> vecc;
     for (int c=0; c<=b; c+=2)
     {
       if (ob.l > lmax) break;
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       vector<vector<vector<vector<ThreeBME_type>>>> vecd;
       for (int d=0; d<=a; d+=2)
       {
         if (oc.l > lmax) break;
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
//         if (ed>E3max) break;
         vector<vector<vector<ThreeBME_type>>> vece;
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           if (od.l > lmax) break;
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
//           if ((ed+ee)>E3max) break;
           vector<vector<ThreeBME_type>> vecf;
           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             if (oe.l > lmax) break;
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 or of.l > lmax) 
             {
               vecf.push_back( vector<ThreeBME_type>(0,0.0) );
               continue;
             }
             int Jde_min = abs(od.j2-oe.j2)/2;
             int Jde_max = (od.j2+oe.j2)/2;

             int dim = 0;
             for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
             {
              for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
              {
                int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
                int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
                for (int J2=J2_min; J2<=J2_max; J2+=2)
                {
                  dim += 5; // 5 different isospin combinations
                  if (Jab==Jde) total_dimension_with_Jcut +=5;
                } //J2
              } //Jde
             } //Jab
             vecf.push_back( vector<ThreeBME_type>(dim,0.0) );
             ++nvectors;
             total_dimension += dim;
//              if (dim >0)
//              cout << "--- " << a << " " << b << " " << c << " " << d << " " << e << " " << f << "  :  " << dim << endl;

           } //f
           vece.push_back( vecf );
           ++nvectors;
         } //e
         vecd.push_back( vece );
         ++nvectors;
       } //d
       vecc.push_back( vecd );
       ++nvectors;
     } //c
     vecb.push_back( vecc );
     ++nvectors;
   } //b
   MatEl.push_back( vecb );
   ++nvectors;
  } //a
  cout << "Allocated " << total_dimension << " three body matrix elements (" <<  total_dimension * sizeof(ThreeBME_type)/1024./1024./1024. << " GB), "
       << nvectors << " vectors (" << nvectors * 3/1024./1024./1024. <<" GB)." << endl;
  cout << "Using only Jab = Jde elements, this would be " << total_dimension_with_Jcut << " (" << total_dimension_with_Jcut * sizeof(ThreeBME_type)/1024./1024./1024. << " GB)." << endl;

}


/*
void ThreeBodyME::Allocate()
{
  E3max = modelspace->GetN3max();
  cout << "Begin AllocateThreeBody() with E3max = " << E3max << endl;
  int norbits = modelspace->GetNumberOrbits();
  int nvectors = 0;
  int total_dimension = 0;

  // Allocate pointers to a
  int dim = 0;
  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   if ((2*oa.n + oa.l)>E3max) break;
   ++dim;
  }
  MatEl.resize(dim); // big ugly pointer
  nvectors += dim;

  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   int ea = 2*oa.n+oa.l;
   if (ea>E3max) break;
   dim = 0;
   for (int b=0; b<=a; b+=2)
   {
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;
     ++dim;
   }
             if (dim<1) continue;
   MatEl[a/2].resize(dim);
   nvectors += dim;

   for (int b=0; b<=a; b+=2)
   {
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;
     dim = 0;
     for (int c=0; c<=b; c+=2)
     {
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       ++dim;
     }
     if (dim<1) continue;
     MatEl[a/2][b/2].resize(dim);
     nvectors += dim;
     for (int c=0; c<=b; c+=2)
     {
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       dim = 0;
       for (int d=0; d<=a; d+=2)
       {
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
         if (ed>E3max) break;
         ++dim;
       }
             if (dim<1) continue;
       MatEl[a/2][b/2][c/2].resize(dim);
       nvectors += dim;

       for (int d=0; d<=a; d+=2)
       {
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
         if (ed>E3max) break;
         dim = 0;
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
           if ((ed+ee)>E3max) break;
           ++dim;
         }
             if (dim<1) continue;
         MatEl[a/2][b/2][c/2][d/2].resize(dim);
         nvectors += dim;
      
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
           if ((ed+ee)>E3max) break;
           dim = 0;
           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             ++dim;
           }
             if (dim<1) continue;
           MatEl[a/2][b/2][c/2][d/2][e/2].resize(dim);
           nvectors += dim;

           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             dim = 0;

             int Jab_min = abs(oa.j2-ob.j2)/2;
             int Jde_min = abs(od.j2-oe.j2)/2;
             int Jab_max = (oa.j2+ob.j2)/2;
             int Jde_max = (od.j2+oe.j2)/2;


             for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
             {
              for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
              {
                int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
                int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
                for (int J2=J2_min; J2<=J2_max; J2+=2)
                {
                  dim += 5; // 5 different isospin combinations
                } //J2
              } //Jde
             } //Jab
             if (dim<1) continue;
             MatEl[a/2][b/2][c/2][d/2][e/2][f/2].resize(dim,0.0);
             total_dimension += dim;
 
           } //f
         } //e
       } //d
     } //c
   } //b
  } //a
  cout << "Allocated " << total_dimension << " three body matrix elements (" <<  total_dimension * sizeof(ThreeBME_type)/1024./1024./1024. << " GB), "
       << nvectors << " vectors (" << nvectors * 24./1024./1024./1024. <<" GB)." << endl;

}
*/




/*
void ThreeBodyME::Allocate()
{
  E3max = modelspace->GetN3max();
  cout << "Begin AllocateThreeBody() with E3max = " << E3max << endl;
  vector<double> zerovector(5,0.0);
  int norbits = modelspace->GetNumberOrbits();
  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   if ((2*oa.n + oa.l)>E3max) continue;
   for (int b=0; b<=a; b+=2)
   {
    Orbit& ob = modelspace->GetOrbit(b);
    if ((2*oa.n+oa.l + 2*ob.n+ob.l) > E3max) continue;
    for (int c=0; c<=b; c+=2)
    {
     Orbit& oc = modelspace->GetOrbit(c);
     if ((2*oa.n+oa.l + 2*ob.n+ob.l + 2*oc.n+oc.l) > E3max) continue;

     // Begin loop over ket states
     for( int d=0; d<=a; d+=2)
     {
      Orbit& od = modelspace->GetOrbit(d);
      for (int e=0; e<= (d==a ? b : d); e+=2)
      {
       Orbit& oe = modelspace->GetOrbit(e);
       for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
       {
        Orbit& of = modelspace->GetOrbit(f);

        // conserve parity
        if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0) continue;


        int Jab_min = abs(oa.j2-ob.j2)/2;
        int Jde_min = abs(od.j2-oe.j2)/2;
        int Jab_max = (oa.j2+ob.j2)/2;
        int Jde_max = (od.j2+oe.j2)/2;


        for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
        {
         for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
         {
           int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
           int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
           for (int J2=J2_min; J2<=J2_max; J2+=2)
           {
             MatEl[{a,b,c,d,e,f,J2,Jab,Jde}] = {0.,0.,0.,0.,0.};
           } //J2
         } //Jde
        } //Jab
       } //f
      } //e
     } //d
    } //c
   } //b
  } //a
}
*/





//*******************************************************************
/// Get three body matrix element in proton-neutron formalism.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
ThreeBME_type ThreeBodyME::GetME_pn(int Jab_in, int Jde_in, int J2, int a, int b, int c, int d, int e, int f)
{

   double tza = modelspace->GetOrbit(a).tz2*0.5;
   double tzb = modelspace->GetOrbit(b).tz2*0.5;
   double tzc = modelspace->GetOrbit(c).tz2*0.5;
   double tzd = modelspace->GetOrbit(d).tz2*0.5;
   double tze = modelspace->GetOrbit(e).tz2*0.5;
   double tzf = modelspace->GetOrbit(f).tz2*0.5;

   double Vpn=0;
   int Tmin = min( abs(tza+tzb+tzc), abs(tzd+tze+tzf) );
   for (int tab=abs(tza+tzb); tab<=1; ++tab)
   {
      // CG calculates the Clebsch-Gordan coefficient
      double CG1 = AngMom::CG(0.5,tza, 0.5,tzb, tab, tza+tzb);
      for (int tde=abs(tzd+tze); tde<=1; ++tde)
      {
         double CG2 = AngMom::CG(0.5,tzd, 0.5,tze, tde, tzd+tze);
         if (CG1*CG2==0) continue;
         for (int T=Tmin; T<=3; ++T)
         {
           double CG3 = AngMom::CG(tab,tza+tzb, 0.5,tzc, T/2., tza+tzb+tzc);
           double CG4 = AngMom::CG(tde,tzd+tze, 0.5,tzf, T/2., tzd+tze+tzf);
           if (CG3*CG4==0) continue;
           Vpn += CG1*CG2*CG3*CG4*GetME(Jab_in,Jde_in,J2,tab,tde,T,a,b,c,d,e,f);

         }
      }
   }
   return Vpn;
}


//*******************************************************************
/// Get three body matrix element in isospin formalism
/// \f$ V_{abcdef}^{J_{ab}J_{de}Jt_{ab}t_{de}T} \f$
///  (which is how they're stored).
/// The elements are stored with the following restrictions: \f$ a\geq b \geq c\f$,
/// \f$ d\geq e \geq f\f$, \f$ a\geq d\f$. If \f$ a=d\f$ then \f$ b \geq e\f$,
/// and if \f$ b=e \f$ then \f$ c \geq f \f$.
/// Other orderings are obtained by recoupling on the fly.
//*******************************************************************
ThreeBME_type ThreeBodyME::GetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in)
{
//   cout << "===== GET =====" << endl;
   return AddToME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in,0.0);
}

//*******************************************************************
/// Set a three body matrix element. Since only a subset of orbit
/// orderings are stored, we need to recouple if the input ordering
/// is different.
//*******************************************************************
void ThreeBodyME::SetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, ThreeBME_type V)
{
   AddToME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in,V);
}

//*******************************************************************
/// Since setting and getting three body matrix elements requires
/// almost identical code, they are combined into one function
/// which adds \f$V_{in}\f$ to the matrix element and returns
/// its value \f$V_{out}\f$. To access a matrix element, we set
/// \f$V_{in}=0\f$. To set the matrix element, we simply
/// disregard $\V_{out}\f$.
//*******************************************************************

ThreeBME_type ThreeBodyME::AddToME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, ThreeBME_type V_in)
{

   // Re-order so that a>=b>=c, d>=e>=f
   int a,b,c,d,e,f;
   int abc_recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);
   int def_recoupling_case = SortOrbits(d_in,e_in,f_in,d,e,f);

   if (d>a or (d==a and e>b) or (d==a and e==b and f>c))
   {
      swap(a,d);
      swap(b,e);
      swap(c,f);
      swap(Jab_in,Jde_in);
      swap(tab_in,tde_in);
      swap(abc_recoupling_case, def_recoupling_case);
   }



   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);
   if (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) return 0;
   if (2*(od.n+oe.n+of.n)+od.l+oe.l+of.l > E3max) return 0;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;
   double je = oe.j2*0.5;
   double jf = of.j2*0.5;

   int Jab_min = abs(ja-jb);
   int Jde_min = abs(jd-je);
   int Jab_max = (ja+jb);
   int Jde_max = (jd+je);

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;



   auto& vj = MatEl.at(a/2).at(b/2).at(c/2).at(d/2).at(e/2).at(f/2);
   
//   cout << "    accessing " << a << " " << b << " " << c << " " << d << " " << e << " " << f << " size = " << vj.size() << endl;
//   cout << " size = " << vj.size() << endl;
   double V_out = 0;
   int J_index = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     double Cj_abc = RecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,J2);
     // Pick up a -1 for odd permutations
     if (abc_recoupling_case>2) Cj_abc *= -1;

     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_def = RecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,J2);
       // Pick up a -1 for odd permutations
       if (def_recoupling_case>2) Cj_def *= -1;

       int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
       int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
       if (J2_min>J2_max) continue;
       J_index += (J2-J2_min)/2*5;

       if (J2>=J2_min and J2<=J2_max)
       {
       for (int tab=tab_min; tab<=tab_max; ++tab)
       {
         double Ct_abc = RecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,T2);
         for (int tde=tde_min; tde<=tde_max; ++tde)
         {
           double Ct_def = RecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,T2);

           int Tindex = 2*tab + tde + (T2-1)/2;

           vj.at(J_index + Tindex) += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj.at(J_index+Tindex);


         }
       }
       }
       J_index += (J2_max-J2+2)/2*5;
     }
   }
   return V_out;
}




//*******************************************************************
/// Coefficients for recoupling three body matrix elements
//*******************************************************************
double ThreeBodyME::RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J)
{
   switch (recoupling_case)
   {
    case 0: return Jab==Jab_in ? 1 : 0;
    case 1: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    case 2: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 3: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 4: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case 5: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
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
int ThreeBodyME::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c)
{
   a_in -= a_in%2;
   b_in -= b_in%2;
   c_in -= c_in%2;
   a=a_in;
   b=b_in;
   c=c_in;
   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);

   int recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? 0 : 3;
   else if (a_in==b)  recoupling_case = (b_in==a) ? 4 : 1;
   else               recoupling_case = (b_in==a) ? 2 : 5;

   return recoupling_case;
}


void ThreeBodyME::Erase()
{
   MatEl.clear();
}




