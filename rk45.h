const int M = 6; //number of equations
const int stage_num = 7; //number of stages

//equations of motion for particle x2 (right) in field of particle x1 (left)
double RHS(int n,double t,double x,double y,double gamma,double ux,double uy)
{
  const int Z = 79;
  const double k = 97962.95;
  const double g = 1.0; 
  
  double r2 = 4.0*(gamma*gamma*x*x+y*y);
  double root = sqrt(1.0+k*k/(g*g*r2*r2*r2)*4.0*(x*x+y*y/(gamma*gamma)));
  double Dx = 2.0*gamma*x/sqrt(r2*r2*r2);
  double Dy = 2.0*gamma*y/sqrt(r2*r2*r2);
  double Hz = -2.0*gamma*ux*y/sqrt(r2*r2*r2);
  double Ex = Dx/root;
  double Ey = Dy/root;
  double Bz = Hz/root;
  
  double tdot = gamma;
  double xdot = ux;
  double ydot = uy;
  double gammadot = Ex*ux+Ey*uy;
  double uxdot = Ex*gamma+Bz*uy;
  double uydot = Ey*gamma-Bz*ux;
 
  if(n==0) return tdot;
  if(n==1) return xdot;
  if(n==2) return ydot;
  if(n==3) return gammadot;
  if(n==4) return uxdot;
  if(n==5) return uydot;
  if(n==6) return root;
  
  return 0;
}

double stage(int n,int i,double X[M],double alpha[M][stage_num])
{
  return RHS(n,X[0]+alpha[0][i],X[1]+alpha[1][i],X[2]+alpha[2][i],X[3]+alpha[3][i]
		    ,X[4]+alpha[4][i],X[5]+alpha[5][i]);
}
      
//coefficients in butcher tableau for rk45
double A(int i,int j)
{
  double a[7][7];
  int n,m;

  //initialize all coefficients to zero
  for(n=0;n<7;++n)
    {
      for(m=0;m<7;++m)
	{
	  a[n][m] = 0.0;
	}
    }

  //set values for nonzero coefficients
  a[1][0] = 1.0/5.0;
  a[2][0] = 3.0/40.0;
  a[2][1] = 9.0/40.0;
  a[3][0] = 44.0/45.0;
  a[3][1] = -56.0/15.0;
  a[3][2] = 32.0/9.0;
  a[4][0] = 19372.0/6561.0;
  a[4][1] = -25360.0/2187.0;
  a[4][2] = 64448.0/6561.0;
  a[4][3] = -212.0/729.0;
  a[5][0] = 9017.0/3168.0;
  a[5][1] = -355.0/33.0;
  a[5][2] = 46732.0/5247.0;
  a[5][3] = 49.0/176.0;
  a[5][4] = -5103.0/18656.0;
  a[6][0] = 35.0/384.0;
  a[6][2] = 500.0/1113.0;
  a[6][3] = 125.0/192.0;
  a[6][4] = -2187.0/6784.0;
  a[6][5] = 11.0/84.0;

  return a[i][j]; 
}

//solution coefficients for rk45
double b(int i)
{
  double b[7];
  b[0] = 5179.0/57600.0;
  b[1] = 0.0;
  b[2] = 7571.0/16695.0;
  b[3] = 393.0/640.0;
  b[4] = -92097.0/339200.0;
  b[5] = 187.0/2100.0;
  b[6] = 1.0/40.0;

  return b[i];
}
