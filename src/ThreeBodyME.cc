#include "ThreeBodyME.hh"
#include "AngMom.hh"


ThreeBodyME::ThreeBodyME()
{}

ThreeBodyME::ThreeBodyME(ModelSpace* ms)
: modelspace(ms)
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
//   cout << "===== SET =====" << endl;
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

//   int A=4;
//   int B=2;
//   int C=0;
//   int D=0;
//   int E=0;
//   int F=0;
//   cout << "  BEGIN Sizes: " << endl;
//   cout << "     "  << MatEl.size() << endl; 
//   cout << "     "  << MatEl[A/2].size() << endl; 
//   cout << "     "  << MatEl[A/2][B/2].size() << endl; 
//   cout << "     "  << MatEl[A/2][B/2][C/2].size() << endl; 
//   cout << "     "  << MatEl[A/2][B/2][C/2][D/2].size() << endl; 
//   cout << "     "  << MatEl[A/2][B/2][C/2][D/2][E/2].size() << endl; 
//   cout << "     "  << MatEl[A/2][B/2][C/2][D/2][E/2][F/2].size() << endl;


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

//   int Jab_min = max( abs(ja-jb), abs(J2*0.5-jc) );
//   int Jab_max = min( ja+jb, J2*0.5+jc );
//   int Jde_min = max( abs(jd-je), abs(J2*0.5-jf) );
//   int Jde_max = min( jd+je, J2*0.5+jf );

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;


//   cout << "After swapping, " << a << " " << b << " " << c << " " << d << " " << e << " " << f << endl;
//
//   cout << "    setting vj" << endl;
//   cout << "  Sizes: " << endl;
//   cout << "     "  << MatEl.size() << endl; 
//   cout << "     "  << MatEl[a/2].size() << endl; 
//   cout << "     "  << MatEl[a/2][b/2].size() << endl; 
//   cout << "     "  << MatEl[a/2][b/2][c/2].size() << endl; 
//   cout << "     "  << MatEl[a/2][b/2][c/2][d/2].size() << endl; 
//   cout << "     "  << MatEl[a/2][b/2][c/2][d/2][e/2].size() << endl; 
//   cout << "     "  << MatEl[a/2][b/2][c/2][d/2][e/2][f/2].size() << endl; 
   auto& vj = MatEl.at(a/2).at(b/2).at(c/2).at(d/2).at(e/2).at(f/2);
//   cout << "    done. size = " << vj.size() << endl;
   double V_out = 0;
   int J_index = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     double Cj_abc = RecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,J2);
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_def = RecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,J2);
       int J2_min = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2));
       int J2_max = min( 2*Jab+oc.j2, 2*Jde+of.j2);
       if (J2_min>J2_max) continue;
       J_index += (J2-J2_min)/2*5;

//       if (abs(Cj_abc*Cj_def) < 1e-6) continue;
//       if ( abs(V_in)<1e-6 )
//       {
//         cout << "Jab_in = " << Jab_in << "  Jab = " << Jab << "  "
//              << "Jde_in = " << Jde_in << "  Jde = " << Jde << "  "
//              << "Jab_min = " << Jab_min << "  Jde_min = "<< Jde_min << "  J2 = " << J2 << endl;
//       }
       if (J2>=J2_min and J2<=J2_max)
       {
       for (int tab=tab_min; tab<=tab_max; ++tab)
       {
         double Ct_abc = RecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,T2);
         for (int tde=tde_min; tde<=tde_max; ++tde)
         {
           double Ct_def = RecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,T2);

           int Tindex = 2*tab + tde + (T2-1)/2;

//           vj[Tindex] += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
//           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj[Tindex];
//           if (a==28 and b==4 and c==0 and d==26 and e==4 and f==0 )
//            if (abs(V_in)<1e-6)
//            {
//            cout << "ja = " << oa.j2 << "/2  jb = " << ob.j2 << "/2  Jab = " << Jab << " jc = " << oc.j2
//                 << "/2  jd = " << od.j2 << "/2  je = " << oe.j2 << "/2    Jde = " << Jde << " jf = " << of.j2
//                 << "/2  J2 = " << J2 << "/2  J_index = " << J_index << "  Tindex = " << Tindex << endl;
//           }
//           vj[J_index + Tindex] += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
//           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj[J_index+Tindex];
           vj.at(J_index + Tindex) += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj.at(J_index+Tindex);
//           if ( abs(V_in)>1e-6 )
//           if ( abs(V_in)<1e-6 )
//           {
//             cout << "<<<<< Setting J_index = " << J_index << ", Tindex = " << Tindex << "  vj = " << vj[J_index + Tindex] << endl;
//           }
//           if ( abs(V_in)<1e-6 and abs(vj[J_index+Tindex])>1e-6)
//           if ( abs(V_in)<1e-6 )
//           {
//              cout << "  tab = " << tab << "  tde = " << tde << "  case_abc " << abc_recoupling_case << "  case_def " << def_recoupling_case 
//                   << "   Cj: " << Cj_abc << " " << Cj_def << "   Ct: " << Ct_abc << " " << Ct_def << " vj = " << vj[J_index+Tindex] << "  V_out = " << V_out << endl;
//           }
//           if (abs(V_in)<1e-6)
//           {
//           cout << "J_index = " << J_index << "  Tindex = " << Tindex << "   vj = " << vj[J_index + Tindex] << "   V_out = " << V_out
//                << "   Cjs: " << Cj_abc << " " << Cj_def << "   Cts: " << Ct_abc << " " << Ct_def
//                << "  element = " << MatEl[a/2][b/2][c/2][d/2][e/2][f/2][J_index + Tindex]
//                << endl;
//           }

         }
       }
       }
       J_index += (J2_max-J2+2)/2*5;
     }
   }
//   if (abs(V_in)<1e-6)
//   cout << "Returning V_out = " << V_out << endl;
   return V_out;
}



/*
double ThreeBodyME::AddToME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, double V_in)
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


   int Jab_min = max( abs(ja-jb), abs(J2*0.5-jc) );
   int Jab_max = min( ja+jb, J2*0.5+jc );
   int Jde_min = max( abs(jd-je), abs(J2*0.5-jf) );
   int Jde_max = min( jd+je, J2*0.5+jf );

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;

   if (abc_recoupling_case==0 or abc_recoupling_case==2)
   {
     Jab_min = Jab_max = Jab_in;
     tab_min = tab_max = tab_in;
   }
   if (def_recoupling_case==0 or def_recoupling_case==2)
   {
      Jde_min = Jde_max = Jde_in;
      tde_min = tde_max = tde_in;
   }

   double V_out = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {
     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_abc = RecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,J2);
       double Cj_def = RecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,J2);

       array<double,5>& vj = MatEl.at({a,b,c,d,e,f,J2,Jab,Jde});
       for (int tab=tab_min; tab<=tab_max; ++tab)
       {
         for (int tde=tde_min; tde<=tde_max; ++tde)
         {
           double Ct_abc = RecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,T2);
           double Ct_def = RecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,T2);

           int Tindex = 2*tab + tde + (T2-1)/2;

           vj[Tindex] += Cj_abc * Cj_def * Ct_abc * Ct_def * V_in;
           V_out += Cj_abc * Cj_def * Ct_abc * Ct_def * vj[Tindex];

         }
       }
     }
   }
   return V_out;
}
*/







//*******************************************************************
/// Coefficients for recoupling three body matrix elements
//*******************************************************************
double ThreeBodyME::RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J)
{
   switch (recoupling_case)
   {
    case 0: return Jab==Jab_in ? 1 : 0;
    case 1: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 2: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case 3: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    case 4: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case 5: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    default: return 0;
    }
}


//*******************************************************************
/// Rearrange orbits (abc) so that a>=b>=c
/// and return an int which reflects the required reshuffling
/// - 0: (abc)_in -> (abc)
/// - 1: (acb)_in -> (abc)
/// - 2: (bac)_in -> (abc)
/// - 3: (bca)_in -> (abc)
/// - 4: (cab)_in -> (abc)
/// - 5: (cba)_in -> (abc)
//*******************************************************************
int ThreeBodyME::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c)
{
   a=a_in;
   b=b_in;
   c=c_in;
   if (a<b)  swap(a,b);
   if (b<c)  swap(b,c);
   if (a<b)  swap(a,b);

   int recoupling_case;
   if (a_in==a)       recoupling_case = (b_in==b) ? 0 : 1;
   else if (a_in==b)  recoupling_case = (b_in==a) ? 2 : 3;
   else               recoupling_case = (b_in==a) ? 4 : 5;

   a -= a%2;
   b -= b%2;
   c -= c%2;
   return recoupling_case;
}


void ThreeBodyME::Erase()
{
// for (auto& it_Orb: MatEl)
// {
//    it_Orb.second = {0.,0.,0.,0.,0.,};
// }
}




