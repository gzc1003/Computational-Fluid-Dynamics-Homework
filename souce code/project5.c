#include <stdio.h>
#define nn 17
#define ne 16
#define NSD 1
#define NEN 2
#define NEF 2
#define LBOUNDARY 0.0
#define RBOUNDARY 1.0

void tdma(double x[], const int N, const double a[], const double b[], double c[]) 
{
	int n;

	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	for (n = 1; n < N; n++) {
		double m = 1.0 / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}

	for (n = N - 1; n-- > 0; )
		x[n] = x[n] - c[n] * x[n + 1];
}


int main(void)
{
	double xyz[nn][NSD];
	int ien[ne][NEN], rng[ne][NEF];
	int i;
	double h = (RBOUNDARY-LBOUNDARY)/ne;
	double tridiagnol_a[nn-2], tridiagnol_b[nn-2], tridiagnol_c[nn-2], 
		   d[nn];
	FILE *fp;


	for (i=0; i<ne; i++) {
		ien[i][0] = i+1;
		ien[i][1] = i+2;
		if (i==0) {
			rng[i][0] = i+1;
			rng[i][1] = -(i+2);
		} else if (i==ne-1) {
			rng[i][0] = -(i+1);
			rng[i][1] = i+2;
		} else {
		   rng[i][0] = -(i+1);
		   rng[i][1] = -(i+2);
		}
	}
	for (i=0; i<nn; i++) {
			xyz[i][0] = LBOUNDARY+i*h;
	}

	d[0] = 0;
	d[nn-1] = 0;
	for (i=1; i<nn-1; i++) {
		tridiagnol_a[i-1] = -1/h;
		tridiagnol_b[i-1] = 2/h; 
		tridiagnol_c[i-1] = -1/h;
		d[i] = h;
	}

	tdma(d+1, nn-2, tridiagnol_a, tridiagnol_b, tridiagnol_c);

	fp = fopen("project5.txt", "w");
	for (i=0; i<nn; i++)	
	{
		fprintf(fp, "%f\t%f\n", i*h, d[i]);
	}
	fclose(fp);

	return 0;
}