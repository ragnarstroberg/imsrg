
#include "AngMom.hh"


//#####################################################################
//  Wigner 6-J symbol { j1 j2 j3 }  See definition and Racah formula at
//                    { J1 J2 J3 }  http://mathworld.wolfram.com/Wigner6j-Symbol.html
double AngMom::SixJ(float j1, float j2, float j3, float J1, float J2,float J3)
{
  
  float triads[4][3] = {{j1,j2,j3},{j1,J2,J3},{J1,j2,J3},{J1,J2,j3}};
  
  for (int i=0;i<4;i++)
  {
     if ( (triads[i][0]+triads[i][1]<triads[i][2])
       || (triads[i][0]-triads[i][1]>triads[i][2])
       || ((int)(2*(triads[i][0]+triads[i][1]+triads[i][2]))%2>0) )
         return 0;
  } 
  double sixj = 0;
  for (float t=j1+j2+j3; t<=(j1+j2+j3+J1+J2+J3); t+=1)
  {
     float ff = fct(t-j1-j2-j3)*fct(t-j1-J2-J3)*
                fct(t-J1-j2-J3)*fct(t-J1-J2-j3)*
                fct(j1+j2+J1+J2-t)*fct(j2+j3+J2+J3-t)*
                fct(j3+j1+J3+J1-t) ;
     if (ff>0)
     {
        sixj += pow(-1,t) * fct(t+1)/ff;
     }
  
  }
  for (int i=0;i<4;i++)
  {
     sixj *= sqrt(Tri( triads[i][0],triads[i][1],triads[i][2]));
  }
  return sixj;
}

