#include "ThreeBodyME.hh"
#include "AngMom.hh"


ThreeBodyME::~ThreeBodyME()
{}

ThreeBodyME::ThreeBodyME()
: modelspace(NULL),E3max(0),total_dimension(0)
{
}

ThreeBodyME::ThreeBodyME(ModelSpace* ms)
: modelspace(ms), E3max(ms->E3max), total_dimension(0)
{}

ThreeBodyME::ThreeBodyME(ModelSpace* ms, int e3max)
: modelspace(ms),E3max(e3max), total_dimension(0)
{}

// Define some constants for the various permutations of three indices
// for use in RecouplingCoefficient and SortOrbits
const int ThreeBodyME::ABC = 0;
const int ThreeBodyME::BCA = 1;
const int ThreeBodyME::CAB = 2;
const int ThreeBodyME::ACB = 3;
const int ThreeBodyME::BAC = 4;
const int ThreeBodyME::CBA = 5;


/// Hash function to map six indices to a single long unsigned int.
/// Each index gets 10 bits, for a maximum of 1024 indices.
/// If we have good isospin, this means we're ok up to emax=43.
// In the future, if this isn't sufficient, we could change
// the key type to a bitset<128>, which would allow up to emax=1446.
size_t ThreeBodyME::KeyHash(size_t a,size_t b,size_t c,size_t d,size_t e,size_t f) const
{
 return (   (a/2)
         + ((d/2) << 10)
         + ((b/2) << 20)
         + ((e/2) << 30)
         + ((c/2) << 40)
         + ((f/2) << 50) );

}

/*
// Confusing nomenclature: J2 means 2 times the total J of the three body system
void ThreeBodyME::Allocate()
{
  MatEl.clear();
  OrbitIndex.clear();
  E3max = modelspace->GetE3max();
  int norbits = modelspace->GetNumberOrbits();
  std::cout << "Begin AllocateThreeBody() with E3max = " << E3max << " norbits = " << norbits << std::endl;
  int nvectors = 0;
  int lmax = 500*norbits; // maybe do something with this later...

  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   int ea = 2*oa.n+oa.l;
   if (ea>E3max) break;
   vector<vector<vector<vector<vector<size_t>>>>> vecb;
   for (int b=0; b<=a; b+=2)
   {
     if (oa.l > lmax) break;
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;

     int Jab_min = std::abs(oa.j2-ob.j2)/2;
     int Jab_max = (oa.j2+ob.j2)/2;
     vector<vector<vector<vector<size_t>>>> vecc;
     for (int c=0; c<=b; c+=2)
     {
       if (ob.l > lmax) break;
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       vector<vector<vector<size_t>>> vecd;
       for (int d=0; d<=a; d+=2)
       {
         if (oc.l > lmax) break;
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
         vector<vector<size_t>> vece;
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           if (od.l > lmax) break;
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
           vector<size_t> vecf;
           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             if (oe.l > lmax) break;
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 or of.l > lmax) 
             {
               vecf.push_back( -1 );
               continue;
             }
             vecf.push_back( total_dimension );
             int Jde_min = std::abs(od.j2-oe.j2)/2;
             int Jde_max = (od.j2+oe.j2)/2;

             for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
             {
              for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
              {
                int J2_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
                int J2_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
                for (int J2=J2_min; J2<=J2_max; J2+=2)
                {
//                  dim += 5; // 5 different isospin combinations
                  total_dimension += 5; // 5 different isospin combinations
                } //J2
              } //Jde
             } //Jab
             ++nvectors;

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
//   MatEl.push_back( vecb );
   OrbitIndex.push_back( vecb);
   ++nvectors;
  } //a
  MatEl.resize(total_dimension,0.0);
  MatEl.shrink_to_fit();
  std::cout << "Allocated " << total_dimension << " three body matrix elements (" <<  total_dimension * sizeof(ThreeBME_type)/1024./1024./1024. << " GB), "
       << nvectors << " vectors (" << nvectors * sizeof(vector<size_t>)/1024./1024./1024. <<" GB)." << std::endl;

}

*/



// Confusing nomenclature: J2 means 2 times the total J of the three body system
void ThreeBodyME::Allocate()
{
  MatEl.clear();
  OrbitIndexHash.clear();
  E3max = modelspace->GetE3max();
  int norbits = modelspace->GetNumberOrbits();
  std::cout << "Begin AllocateThreeBody() with E3max = " << E3max << " norbits = " << norbits << std::endl;
  int lmax = 500*norbits; // maybe do something with this later...

  for (int a=0; a<norbits; a+=2)
  {
   Orbit& oa = modelspace->GetOrbit(a);
   int ea = 2*oa.n+oa.l;
   if (ea>E3max) break;
   for (int b=0; b<=a; b+=2)
   {
     if (oa.l > lmax) break;
     Orbit& ob = modelspace->GetOrbit(b);
     int eb = 2*ob.n+ob.l;
     if ((ea+eb)>E3max) break;

     int Jab_min = std::abs(oa.j2-ob.j2)/2;
     int Jab_max = (oa.j2+ob.j2)/2;
     for (int c=0; c<=b; c+=2)
     {
       if (ob.l > lmax) break;
       Orbit& oc = modelspace->GetOrbit(c);
       int ec = 2*oc.n+oc.l;
       if ((ea+eb+ec)>E3max) break;
       for (int d=0; d<=a; d+=2)
       {
         if (oc.l > lmax) break;
         Orbit& od = modelspace->GetOrbit(d);
         int ed = 2*od.n+od.l;
         for (int e=0; e<= (d==a ? b : d); e+=2)
         {
           if (od.l > lmax) break;
           Orbit& oe = modelspace->GetOrbit(e);
           int ee = 2*oe.n+oe.l;
           for (int f=0; f<=((d==a and e==b) ? c : e); f+=2)
           {
             if (oe.l > lmax) break;
             Orbit& of = modelspace->GetOrbit(f);
             int ef = 2*of.n+of.l;
             if ((ed+ee+ef)>E3max) break;
             if ((oa.l+ob.l+oc.l+od.l+oe.l+of.l)%2>0 or of.l > lmax) 
             {
               continue;
             }
             OrbitIndexHash[ KeyHash(a,b,c,d,e,f) ] = total_dimension;
             int Jde_min = std::abs(od.j2-oe.j2)/2;
             int Jde_max = (od.j2+oe.j2)/2;

             for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
             {
              for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
              {
                int J2_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
                int J2_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
                for (int J2=J2_min; J2<=J2_max; J2+=2)
                {
                  total_dimension += 5; // 5 different isospin combinations
                } //J2
              } //Jde
             } //Jab
           } //f
         } //e
       } //d
     } //c
   } //b
  } //a
  MatEl.resize(total_dimension,0.0);
  std::cout << "Allocated " << total_dimension << " three body matrix elements (" <<  total_dimension * sizeof(ThreeBME_type)/1024./1024./1024. << " GB), "
       << std::endl << "  number of buckets in hash table: " << OrbitIndexHash.bucket_count() << "  and load factor = " << OrbitIndexHash.load_factor()
       << "  estimated storage ~ " << ((OrbitIndexHash.bucket_count()+OrbitIndexHash.size()) * (sizeof(size_t)+sizeof(void*))) / (1024.*1024.*1024.) << " GB"
       << std::endl;

}





//*******************************************************************
/// Get three body matrix element in proton-neutron formalism.
/// \f[
///  V_{abcdef}^{(pn)} = \sum_{t_{ab} t_{de} T} <t_a t_b | t_{ab}> <t_d t_e | t_{de}>
///  <t_{ab} t_c | T> <t_{de} t_f| T> V_{abcdef}^{t_{ab} t_{de} T}
/// \f]
//*******************************************************************
ThreeBME_type ThreeBodyME::GetME_pn(int Jab_in, int Jde_in, int J2, int a, int b, int c, int d, int e, int f) const
{

   double tza = modelspace->GetOrbit(a).tz2*0.5;
   double tzb = modelspace->GetOrbit(b).tz2*0.5;
   double tzc = modelspace->GetOrbit(c).tz2*0.5;
   double tzd = modelspace->GetOrbit(d).tz2*0.5;
   double tze = modelspace->GetOrbit(e).tz2*0.5;
   double tzf = modelspace->GetOrbit(f).tz2*0.5;

   double Vpn=0;
   int Tmin = std::min( std::abs(tza+tzb+tzc), std::abs(tzd+tze+tzf) );
   for (int tab=std::abs(tza+tzb); tab<=1; ++tab)
   {
      // CG calculates the Clebsch-Gordan coefficient
      double CG1 = AngMom::CG(0.5,tza, 0.5,tzb, tab, tza+tzb);
      for (int tde=std::abs(tzd+tze); tde<=1; ++tde)
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
ThreeBME_type ThreeBodyME::GetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const
{
   auto elements =  AccessME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in);
   double me = 0;
   for (auto elem : elements) me += MatEl.at(elem.first) * elem.second;
   return me;
}

//*******************************************************************
/// Set a three body matrix element. Since only a subset of orbit
/// orderings are stored, we need to recouple if the input ordering
/// is different.
//*******************************************************************
void ThreeBodyME::SetME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in, ThreeBME_type V)
{
   auto elements = AccessME(Jab_in,Jde_in,J2,tab_in,tde_in,T2,a_in,b_in,c_in,d_in,e_in,f_in);
   double me = 0;
   for (auto elem : elements)  me += MatEl.at(elem.first) * elem.second;
   for (auto elem : elements)  MatEl.at(elem.first) += (V-me)*elem.second;
}

//*******************************************************************
/// Since the code for setting or getting matrix elements is almost
/// identical, do all the work here to pull out a list of indices
/// and coefficients which are needed for setting or getting.
//*******************************************************************
std::vector<std::pair<size_t,double>> ThreeBodyME::AccessME(int Jab_in, int Jde_in, int J2, int tab_in, int tde_in, int T2, int a_in, int b_in, int c_in, int d_in, int e_in, int f_in) const
{

   std::vector<std::pair<size_t,double>> elements;
   // Re-order so that a>=b>=c, d>=e>=f
   int a,b,c,d,e,f;
   int abc_recoupling_case = SortOrbits(a_in,b_in,c_in,a,b,c);
   int def_recoupling_case = SortOrbits(d_in,e_in,f_in,d,e,f);

   if (d>a or (d==a and e>b) or (d==a and e==b and f>c))
   {
	   std::swap(a,d);
      	   std::swap(b,e);
      	   std::swap(c,f);
      	   std::swap(Jab_in,Jde_in);
      	   std::swap(tab_in,tde_in);
      	   std::swap(abc_recoupling_case, def_recoupling_case);
   }

   auto it_hash = OrbitIndexHash.find(KeyHash(a,b,c,d,e,f));
   if ( it_hash == end(OrbitIndexHash) )   return elements;

   Orbit& oa = modelspace->GetOrbit(a);
   Orbit& ob = modelspace->GetOrbit(b);
   Orbit& oc = modelspace->GetOrbit(c);
   Orbit& od = modelspace->GetOrbit(d);
   Orbit& oe = modelspace->GetOrbit(e);
   Orbit& of = modelspace->GetOrbit(f);
   if (2*(oa.n+ob.n+oc.n)+oa.l+ob.l+oc.l > E3max) return elements;
   if (2*(od.n+oe.n+of.n)+od.l+oe.l+of.l > E3max) return elements;

   double ja = oa.j2*0.5;
   double jb = ob.j2*0.5;
   double jc = oc.j2*0.5;
   double jd = od.j2*0.5;
   double je = oe.j2*0.5;
   double jf = of.j2*0.5;

   int Jab_min = std::abs(ja-jb);
   int Jde_min = std::abs(jd-je);
   int Jab_max = (ja+jb);
   int Jde_max = (jd+je);

   int tab_min = T2==3 ? 1 : 0;
   int tab_max = 1;
   int tde_min = T2==3 ? 1 : 0;
   int tde_max = 1;


   auto indx = it_hash->second;
   if (indx > MatEl.size()) std::cout << "ThreeBodyME::AccessME() --  AAAAHHH indx = " << indx << "  but MatEl.size() = " << MatEl.size() << std::endl;

   int J_index = 0;
   for (int Jab=Jab_min; Jab<=Jab_max; ++Jab)
   {

     double Cj_abc = RecouplingCoefficient(abc_recoupling_case,ja,jb,jc,Jab_in,Jab,J2);
     // Pick up a -1 for odd permutations
     if ( abc_recoupling_case/3 != def_recoupling_case/3 ) Cj_abc *= -1;

     for (int Jde=Jde_min; Jde<=Jde_max; ++Jde)
     {
       double Cj_def = RecouplingCoefficient(def_recoupling_case,jd,je,jf,Jde_in,Jde,J2);

       int J2_min = std::max( std::abs(2*Jab-oc.j2), std::abs(2*Jde-of.j2));
       int J2_max = std::min( 2*Jab+oc.j2, 2*Jde+of.j2);
       if (J2_min>J2_max) continue;
       J_index += (J2-J2_min)/2*5;

       if (J2>=J2_min and J2<=J2_max and std::abs(Cj_abc*Cj_def)>1e-8)
       {
         for (int tab=tab_min; tab<=tab_max; ++tab)
         {
           double Ct_abc = RecouplingCoefficient(abc_recoupling_case,0.5,0.5,0.5,tab_in,tab,T2);
           for (int tde=tde_min; tde<=tde_max; ++tde)
           {
             double Ct_def = RecouplingCoefficient(def_recoupling_case,0.5,0.5,0.5,tde_in,tde,T2);
             if (std::abs(Ct_abc*Ct_def)<1e-8) continue;

             int Tindex = 2*tab + tde + (T2-1)/2;

             elements.emplace_back( std::make_pair(indx + J_index + Tindex, Cj_abc * Cj_def * Ct_abc * Ct_def )) ;
           }
         }
       }
       J_index += (J2_max-J2+2)/2*5;
     }
   }
   return elements;
}


//*******************************************************************
/// Coefficients for recoupling three body matrix elements
//*******************************************************************
inline double ThreeBodyME::RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int J) const
{
   if ( std::abs(int(ja-jb))>Jab  or int(ja+jb)<Jab) return 0;
   if ( std::abs(int(jc-J/2.))>Jab  or int(jc+J/2.)<Jab) return 0;
   switch (recoupling_case)
   {
    case ABC: return Jab==Jab_in ? 1 : 0;
    case BCA: return modelspace->phase( jb+jc+Jab_in+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
    case CAB: return modelspace->phase( ja+jb-Jab+1) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case ACB: return modelspace->phase( jb+jc+Jab_in-Jab) * sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(jb, ja, Jab, jc, J/2., Jab_in);
    case BAC: return Jab==Jab_in ? modelspace->phase(ja+jb-Jab) : 0;
    case CBA: return -sqrt((2*Jab_in+1)*(2*Jab+1)) * modelspace->GetSixJ(ja, jb, Jab, jc, J/2., Jab_in);
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
int ThreeBodyME::SortOrbits(int a_in, int b_in, int c_in, int& a, int& b, int& c) const
{
   a_in -= a_in%2;
   b_in -= b_in%2;
   c_in -= c_in%2;
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


void ThreeBodyME::Erase()
{
   MatEl.clear();
}

/// Free up the memory used for the matrix elements
void ThreeBodyME::Deallocate()
{
   std::vector<ThreeBME_type>().swap(MatEl);
   OrbitIndexHash.clear(); 
}



void ThreeBodyME::WriteBinary(std::ofstream& f)
{
  f.write((char*)&E3max,sizeof(E3max));
  f.write((char*)&total_dimension,sizeof(total_dimension));
  f.write((char*)&MatEl[0],total_dimension);
}

void ThreeBodyME::ReadBinary(std::ifstream& f)
{
  f.read((char*)&E3max,sizeof(E3max));
  f.read((char*)&total_dimension,sizeof(total_dimension));
  Allocate();
  f.read((char*)&MatEl[0],total_dimension*sizeof(ThreeBME_type));
}





