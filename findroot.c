#include <math.h>
#include <stdio.h>

#include<stdio.h>
#include<math.h>

#define A (0.000000948)
#define B (2.5)
#define N 2

//base station locations
//double X[3] = {0, 4000.0, 2000.0};
//double Y[3] = {0, 0, 2000.0};
double X[2] = {0, 2000.0};
double Y[2] = {0, 0};

//util functions
int my_getnbr(char *str)
{
  int result;
  int puiss;

  result = 0;
  puiss = 1;
  while (('-' == (*str)) || ((*str) == '+'))
    {
      if (*str == '-')
        puiss = puiss * -1;
      str++;
    }
  while ((*str >= '0') && (*str <= '9'))
    {
      result = (result * 10) + ((*str) - '0');
      str++;
    }
  return (result * puiss);
}

double sqr (double x) {
	return x * x;
}
double cub (double x) {
	return x * x * x;
}

//actual functions
double atr(double x, double y) {
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += 1/(B + A * (sqr(x-X[i]) + sqr(y-Y[i])));
	return ret;
}

double f(double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += (x-X[i])/sqr(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) ;
	return ret;
}
double g(double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += (y-Y[i])/sqr(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) ;
	return ret;
}
double dfx (double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += ( 1.0/sqr(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) - 4.0*A*sqr(x-X[i])/cub(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) ) ;
	return ret;
}
double dfy (double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += -4.0*A*(x-X[i])*(y-Y[i])/cub(B + A * (sqr(x-X[i]) + sqr(y-Y[i])));
	return ret;
}
double dgx (double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += -4.0*A*(x-X[i])*(y-Y[i])/cub(B + A * (sqr(x-X[i]) + sqr(y-Y[i])));
	return ret;
}
double dgy (double x, double y)
{
	double ret = 0.0;
	int i = 0;
	for(i=0; i<N; ++i) ret += ( 1.0/sqr(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) - 4.0*A*sqr(y-Y[i])/cub(B + A * (sqr(x-X[i]) + sqr(y-Y[i]))) ) ;
	return ret;
}

int main(int argc, char** argv)
{
    int itr, maxmitr;
    double h, allerr;
	double x0, y0, x1, y1, m1, m2, m3, m4, g1, g2, d1,d2;

	x0 = my_getnbr(argv[1]);
	x1 = my_getnbr(argv[2]);
	allerr = 0.0001;
	maxmitr = 500;
    for (itr=1; itr<=maxmitr; itr++)
    {
		m1 = dfx(x0,y0);
		m2 = dfy(x0,y0);
		m3 = dgx(x0,y0);
		m4 = dgy(x0,y0);
		g1 = f(x0,y0);
		g2 = g(x0,y0);
		//From: http://www.seas.ucla.edu/~vandenbe/103/lectures/newton.pdf
		//H= [m1 m2]
		//   [m3 m4]
		//-g= [-g1]
		//   [-g2]
		//d= [d1]
		//   [d2]
		//
		//Solve Hd=-g
		//That is,
		//m1*d1+m2*d2+g1=0
		//m3*d1+m4*d2+g2=0
		//
		//solution:
		//d1=(m4*g1-m2*g2)/(m2*m3-m1*m4)
		//d2=(m1*g2-m3*g1)/(m2*m3-m1*m4)
		
		d1=(m4*g1-m2*g2)/(m2*m3-m1*m4);
		d2=(m1*g2-m3*g1)/(m2*m3-m1*m4);
		x1=x0+d1;
		y1=y0+d2;
		
        printf(" At Iteration no. %d, x = %9.6f, y = %9.6f\n", itr, x1, y1);
        if (fabs(d1) < allerr && fabs(d2) < allerr)
        {
            printf("After %d iterations, root = %8.6f, %8.6f\n", itr, x1, y1);
			printf("ATR = %8.6f\n", atr(x1,y1));
            return 0;
        }
        x0=x1;
		y0=y1;
    }
    printf(" The required solution does not converge or iterations are insufficient\n");
    return 1;
}

