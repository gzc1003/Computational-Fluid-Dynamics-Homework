#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int nx, ny;


double residual(double u[], double u_old[])
{
	int i;
	double residual=0.0;
	for (i=0; i<(nx+2)*(ny+2); i++) {
		if (fabs(u[i] - u_old[i]) > residual) 
			residual = fabs(u[i] - u_old[i]);
	}
	return residual;
}

double f(double x, double y)
{
	double r;
	r = -10 * (1 - 10*(x - 0.5)*(x - 0.5)) * exp(-5*(x-0.5)*(x - 0.5))*
		(exp(-5*(y-0.5)*(y - 0.5)) - exp(-5.0/4))
        -10 * (1 - 10*(y - 0.5)*(y - 0.5)) * exp(-5*(y-0.5)*(y - 0.5))*
		(exp(-5*(x-0.5)*(x - 0.5)) - exp(-5.0/4));
	return r;
}

int main(void)
{
	int k;
	int i, j;
	double *u, *u_old;
	double hx2, hy2, a, b, c, d, e;
	FILE *fp, *fp2;
	fp = fopen("residual.txt", "w");
	fp2 = fopen("solution.txt", "w");

	printf("Input the number of interior point in x direction: ");
	scanf_s("%d", &nx); 
	printf("Input the number of interior point in y direction: ");
	scanf_s("%d", &ny);
	u = (double *) malloc(sizeof(double)*(nx+2)*(ny+2));
	u_old = (double *) malloc(sizeof(double)*(nx+2)*(ny+2));

	hx2 = (nx+1)*(nx+1);
	hy2 = (ny+1)*(ny+1);
	b = a = hx2 / (2*(hx2+hy2));
	c = d = hy2 / (2*(hx2+hy2));
	e = 1.0 / (2*(hx2+hy2));

	k = 0;
	for (j = 0; j < ny+2; j++) {
		for (i = 0; i< nx+2; i++) {
			u[i+j*(nx+2)] = 0;
			u_old[i+j*(nx+2)] = 0; 
		}
	}

	while (1) {
		k += 1;
		for (j = 1; j < ny+1; j++) {
			for (i = 1; i< nx+1; i++) {
				u[i+j*(nx+2)] = a*u_old[i-1+j*(nx+2)] + b*u_old[i+1+j*(nx+2)]
							    + c*u_old[i+(j-1)*(nx+2)] + d*u_old[i+(j+1)*(nx+2)]
							    - e*f((double)i/(nx+1), (double)j/(ny+1));
			}
		}

		if (residual(u, u_old)<1e-6)
			break;
		fprintf(fp, "%d\t%f\n", k, residual(u, u_old));
		for (j = 1; j < ny+2; j++) {
			for (i = 1; i< nx+2; i++) 
				u_old[i+j*(nx+2)] = u[i+j*(nx+2)]; 
		}

	}

	for (i = 0; i < nx+2; i++) 
		fprintf(fp2, "%f\t%f\n", (double)i/(nx+1), u[i+25*(nx+2)]);
	
	free(u);
	return 0;
}