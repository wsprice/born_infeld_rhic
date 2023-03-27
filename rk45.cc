#include <stdio.h>
#include <math.h>
#include "rk45.h"

const double r0 = 0.00463; //length units in fm
const double h0 = 1.0/2.0;
const int max_steps = 2e5;
const double gamma0 = 1.5;
const double u0 = sqrt(gamma0*gamma0-1.0);
//const double b_imp = 10.0/r0; //impact parameter in fm
const double x0 = 500.0/r0; //initial separation in fm

int main()
{
  FILE *data = fopen("carbon_data/gamma1.5/angle_MAX.txt","w");

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
  double u_squared,angle,root,udot,gammadot,uxdot,uydot;
  double b_imp = 10.0/r0;
  //double gamma0 = 1.05;
  int l;

  for(l=0;l<200;++l)
    {
    tau = 0.0;

  //initial conditions
  X[0] = 0.0; //t
  X[1] = x0; //x
  X[2] = b_imp/2.0; //y
  X[3] = gamma0;
  X[4] = -u0;
  X[5] = 0.0; //uy

  do
    {
      //error check with u^2=1
      u_squared = fabs(X[3]*X[3] - X[4]*X[4] - X[5]*X[5] - 1.0);
      gammadot = RHS(3,X[0],X[1],X[2],X[3],X[4],X[5]);
      uxdot = RHS(4,X[0],X[1],X[2],X[3],X[4],X[5]);
      uydot = RHS(5,X[0],X[1],X[2],X[3],X[4],X[5]);

      udot = sqrt(uxdot*uxdot + uydot*uydot - gammadot*gammadot);

      //print solution and error
      //fprintf(data,"%d   %.10f   %.10f   %.10f   %.10f   %.10f \n",step_num,X[0],X[1]*r0,X[2]*r0,X[3],u_squared);

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
      
      /*if(X[1]>-2.0/r0 && X[1]<2.0/r0)
	{
	  h = 2.0e-6;
	}
      else
	{
	  h = h0;
	  }*/
    }
  
  while(X[1]>-x0);
  
  angle = 180.0/M_PI*acos(fabs(X[4])/sqrt(X[4]*X[4]+X[5]*X[5]));
  fprintf(data,"%.10f   %.10f \n",b_imp*r0,angle);
  b_imp = b_imp + 0.15/r0;
  //gamma0 = gamma0 + 0.0895;
}
    
  fclose(data);
  
  return 0;
}
		     
