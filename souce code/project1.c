#include <stdio.h>
#include <math.h>
#define N 4


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

double f(double x)
{
	return exp(x);
}


int main(void)
{
	int i;
	double lboundary = 0, rboundary = 1;
	double h = (rboundary-lboundary)/N;
	double tridiagnol_a[N-1], tridiagnol_b[N-1], tridiagnol_c[N-1], 
		   d[N+1];
	double sum, error;
	FILE *fp;
	
	d[0] = 1;
	d[N] = exp(1.0);

	for (i=1; i<N; i++) 
	{
		tridiagnol_a[i-1] = 1;
		tridiagnol_b[i-1] = -2; 
		tridiagnol_c[i-1] = 1;
		if (i==1)
			d[i] = h*h*f(i*h) - 1;
		else if (i==N-1)
			d[i] = h*h*f(i*h) - exp(1.0);
		else 
			d[i] = h*h*f(i*h);
	}

	tdma(d+1, N-1, tridiagnol_a, tridiagnol_b, tridiagnol_c);

	sum = 0;
	for (i=1; i<N; i++)	
	{
		sum += pow((d[i] - exp(i*h)), 2);
	}
	error = sqrt(h*sum);
	fp = fopen("project1.txt", "w");
	for (i=0; i<N+1; i++)	
	{
		fprintf(fp, "%f\t%f\n", i*h, d[i]);
	}
	fclose(fp);
	return 0;
}