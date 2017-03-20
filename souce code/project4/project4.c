#include "Jacobi.h"
#include <math.h>
#include <stdio.h>

#define N 5
#define L 1.0
#define Rho 1.0
#define u 0.2
#define k 0.1
#define LBOUNDARY 1.0
#define RBOUNDARY 0.0

void tdma(double x[], const int NN, const double a[], const double b[], double c[]) 
{
	int n;

	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	for (n = 1; n < NN; n++) {
		double m = 1.0 / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}

	for (n = NN - 1; n-- > 0; )
		x[n] = x[n] - c[n] * x[n + 1];
}


void upwind(double phi[])
{
	int i;
	double D = k/(L/N), F=Rho*u;
	double tridiagonal_a[] = {0, -(D+F), -(D+F), -(D+F), -(D+F)}, 
		   tridiagonal_b[] = {3*D+F, 2*D+F, 2*D+F, 2*D+F, 3*D+F},
	       tridiagonal_c[] = {-D, -D, -D, -D, 0},
	       d[] = {(2*D+F)*LBOUNDARY, 0, 0, 0, 2*D*RBOUNDARY};
	tdma(d, N, tridiagonal_a, tridiagonal_b, tridiagonal_c);
	for (i=0; i<N; i++)
		phi[i] = d[i];
}


void QUICK(double phi[])
{
	double D = k/(L/N), F=Rho*u;
	double A[N][N] = {{4*D+7*F/8, -4*D/3+3*F/8, 0, 0, 0},
				      {-(D+F), 2*D+3*F/8, -D+3*F/8, 0, 0},
					  {F/8, -(D+7*F/8), 2*D+3*F/8, -D+3*F/8, 0},
					  {0, F/8, -(D+7*F/8), 2*D+3*F/8, -D+3*F/8},
					  {0, 0, F/8, -(4*D/3+3*F/4), 4*D-3*F/8}
	                 },
		   b[N] = {8*D/3+5*F/4, -F/4, 0,0,0};
	Jacobi(A, phi, b, 1e-6);
}


int main(void)
{
	double phi[N], phi2[N], exact[N];
	double x;
	int i;
	FILE *fp;
	for (i=0; i<N; i++) {
		x = 0.1+0.2*i;
		exact[i] = 1-(exp(Rho*u*x/k)-1)/(exp(Rho*u*L/k)-1);
	}
	upwind(phi);
	QUICK(phi2);
	fp = fopen("project4.txt", "w");
	for (i=0; i<N; i++)
		fprintf(fp, "%f\t%f\t%f\t%f\n", 0.1+0.2*i, exact[i], phi[i], phi2[i]);

	return 0;
}

