#include <stdio.h>
#include <math.h>
#include "rk45.h"

const double r0 = 0.0489; //length units in fm
const double h0 = 2.0e-2;
const int max_steps = 2e5;
const double gamma0 = 200.0;
const double u0 = sqrt(gamma0*gamma0-1.0);
//const double b_imp = 15.0/r0; //impact parameter in fm
const double x0 = 1000.0/r0; //initial separation in fm

int main()
{
  FILE *data = fopen("angle_BI.txt","w");

  //rk variables
  int i; //row number
  int j; //column number
  int n; //equation number
  int step_num=0; //step number
  double k[M][stage_num],alpha[M][stage_num],beta[M];
  double h = h0;

  //dynamical variables
  double tau = 0.0;
  double X[M];
  double u_squared,angle,root;
  double b_imp = 15.0/r0;
  int l;
  
  for(l=0;l<100;++l)
    {
      tau = 0.0;

  //initial conditions
  X[0] = 0.0; //t
  X[1] = x0; //x
  X[2] = b_imp/2.0; //y
  X[3] = gamma0; //gamma
  X[4] = -u0; //ux
  X[5] = 0.0; //uy

  do
    {
      //error check with u^2=1
      u_squared = fabs(X[3]*X[3] - X[4]*X[4] - X[5]*X[5] - 1.0);
      root = RHS(6,X[0],X[1],X[2],X[3],X[4],X[5]);

      //print solution and error
      //fprintf(data,"%d   %.10f   %.10f   %.10f   %.10f   %.10f   %.10f \n",step_num,X[0],X[1]*r0,X[2]*r0,X[3],u_squared,root);

      //calculate integration step with rk45     
      for(i=0;i<stage_num;++i)
	{
	  for(n=0;n<M;++n)
	    {
	      alpha[n][i] = 0.0;
	      for(j=0;j<i;++j)
		{
		  alpha[n][i] = alpha[n][i] + A(i,j)*k[n][j];
		}
	    }
	  for(n=0;n<M;++n)
	    {
	      k[n][i] = h*stage(n,i,X,alpha);
	    }
	}
      for(n=0;n<M;++n)
	{
	  beta[n] = 0.0;
	  for(i=0;i<stage_num;++i)
	    {
	      beta[n] = beta[n] + b(i)*k[n][i];
	    }
	  X[n] = X[n] + beta[n];
	}
      tau = tau + h;
      ++step_num;
      
      if(X[1]>-2.0/r0 && X[1]<2.0/r0)
	{
	  h = 2.0e-6;
	}
      else
	{
	  h = h0;
	}
    }
  
  while(X[1]>-x0);
  //fclose(data);

  angle = 180.0/M_PI*acos(u0/sqrt(X[4]*X[4]+X[5]*X[5]));
  fprintf(data,"%.10f   %.10f \n",b_imp*r0,angle);
  b_imp = b_imp + 0.85/r0;
    }
  
  return 0;
} 
