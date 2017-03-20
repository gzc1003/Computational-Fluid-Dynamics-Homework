#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double maximum(double x, double y, double z);
int flux(double *f1,double *f2,double *f3, double *p,double *u,double *rho,double *a,double *e,double *h, double gama,double dx,double dt,int N);

int main()
{
	int N,j;		//number of cells
	double *u1, *u2, *u3, *f1,*f2, *f3,*u,*rho,*e,*ei,*p,*h,*a;	//variables and fluxes
	double Smax,m,namta,CFL,t,dt,dx,gama;
	FILE *fp;

	gama = 1.4;
	CFL = 0.8;

	printf("Please input the number of cells: N = ");
	scanf_s("%d",&N);
	printf("%d\n",N);

	u1 = (double *)calloc(N,sizeof(double));	//allocate memories
	u2 = (double *)calloc(N,sizeof(double));
	u3 = (double *)calloc(N,sizeof(double));
	f1 = (double *)calloc(N+1,sizeof(double));
	f2 = (double *)calloc(N+1,sizeof(double));
	f3 = (double *)calloc(N+1,sizeof(double));
	u = (double *)calloc(N,sizeof(double));
	rho = (double *)calloc(N,sizeof(double));
	p = (double *)calloc(N,sizeof(double));
	e = (double *)calloc(N,sizeof(double));
	ei = (double *)calloc(N,sizeof(double));
	h = (double *)calloc(N,sizeof(double));
	a = (double *)calloc(N,sizeof(double));

	for(j=0;j<N/2;j++)	//initiation
	{
		p[j] = 1.0;
		rho[j] = 1.0;
	}
	for(j=N/2;j<N;j++)
	{
		p[j] = 0.1;
		rho[j] = 0.125;
	}

	for(j=0;j<N;j++)
	{
		u[j] = 0;
		ei[j] = p[j]/((gama-1)*rho[j]);
		e[j] = ei[j]+0.5*u[j]*u[j];
		h[j] = e[j]+p[j]/rho[j];
		a[j] = sqrt(gama*p[j]/rho[j]); 
		u1[j] = rho[j];
		u2[j] = rho[j]*u[j];
		u3[j] = rho[j]*e[j]; 
	}

	f1[0] = rho[0]*u[0];		//Boundary conditions
	f2[0] = rho[0]*u[0]*u[0]+p[0];
	f3[0] = rho[0]*h[0]*u[0];
	f1[N] = rho[N-1]*u[N-1];
	f2[N] = rho[N-1]*u[N-1]*u[N-1]+p[N-1];
	f3[N] = rho[N-1]*h[N-1]*u[N-1];


	dx = 1.0/N;	//length of cell
	t = 0.0;
	while(t<=0.2)
	{
		a[0] = sqrt(gama*p[0]/rho[0]);		//get the maximum wave speed
		Smax = maximum(u[0]-a[0],u[0],u[0]+a[0]);
		for(j=1;j<100;j++)
		{
			a[j] = sqrt(gama*p[j]/rho[j]);
			m =  maximum(u[j]-a[j],u[j],u[j]+a[j]);
			if(Smax < m)
			{
				Smax = m;
			}
		}
		namta = CFL/Smax;
		dt = namta * dx;
		t += dt;
		flux(f1,f2,f3,p,u,rho,a,e,h,gama,dx,dt,N);
		for(j=0;j<N;j++)
		{
			u1[j] = u1[j]-namta*(f1[j+1]-f1[j]);
			u2[j] = u2[j]-namta*(f2[j+1]-f2[j]);
			u3[j] = u3[j]-namta*(f3[j+1]-f3[j]);

			rho[j] = u1[j];
			u[j] = u2[j]/rho[j];
			e[j] = u3[j]/rho[j];
			ei[j] = e[j]-0.5*u[j]*u[j];
			p[j] = (gama-1)*rho[j]*ei[j];
            h[j] = e[j]+p[j]/rho[j];
            a[j] = sqrt(gama*p[j]/rho[j]);

		}
	}
	
	fp = fopen("project3.txt", "w");
	for(j=0;j<N;j++)
	{
		fprintf(fp, "%lf\t %lf\t %lf\t %lf\n", j*dx,p[j],u[j],rho[j]);
	}
		

	free(u1);
	free(u2);
	free(u3);
	free(f1);
	free(f2);
	free(f3);
	free(u);
	free(p);
	free(e);
	free(ei);
	free(a);
	free(h);
	free(rho);
}

double maximum(double x, double y, double z)
{
	double maximum;
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	maximum = x;
	if(maximum<y)
	{
		maximum = y;
	}
	if(maximum<z)
	{
		maximum = z;
	}
	return maximum;
}

int flux(double *f1,double *f2,double *f3, double *p,double *u,double *rho,double *a,double *e,double *h, double gama,double dx,double dt,int N)
{
	int j;
	double ps,us,rhoB,aB,Sl,Sr,qr,ql,urs1,uls1,urs2,uls2,urs3,uls3;
	for(j=1;j<N;j++)
	{
		rhoB = 0.5*(rho[j-1]+rho[j]);
		aB = 0.5*(a[j-1]+a[j]);
		ps = 0.5*(p[j-1]+p[j])-0.5*(u[j]-u[j-1])*rhoB*aB;
		us = 0.5*(u[j-1]+u[j])-0.5*(p[j]-p[j-1])/(rhoB*aB);
		if(ps<=p[j-1])
		{
			ql = 1;
		}
		else
		{
			ql = sqrt(1+0.5*(gama+1)*(ps/p[j-1]-1)/gama);
		}
		if(ps<=p[j])
		{
			qr = 1;
		}
		else
		{
			qr = sqrt(1+0.5*(gama+1)*(ps/p[j]-1)/gama);
		}
		Sl = u[j-1]-a[j-1]*ql;
		Sr = u[j]+a[j]*qr;
		urs1 = rho[j]*(Sr-u[j])/(Sr-us);
		urs2 = urs1 * us;
		urs3 = urs1 * (e[j]+(us-u[j])*(us+p[j]/(rho[j]*(Sr-u[j]))));
		
		uls1 = rho[j-1]*(Sl-u[j-1])/(Sl-us);
		uls2 = uls1 * us;
		uls3 = uls1 * (e[j-1]+(us-u[j-1])*(us+p[j-1]/(rho[j-1]*(Sl-u[j-1]))));
		
		if(Sl>=0)
		{
			f1[j] = rho[j-1]*u[j-1];
			f2[j] = rho[j-1]*u[j-1]*u[j-1]+p[j-1];
			f3[j] = rho[j-1]*h[j-1]*u[j-1];
		}
		else if(us>=0)
			{
				f1[j] = rho[j-1]*u[j-1]+Sl*(uls1-rho[j-1]);
				f2[j] = rho[j-1]*u[j-1]*u[j-1]+p[j-1]+ Sl*(uls2-rho[j-1]*u[j-1]);
				f3[j] = rho[j-1]*h[j-1]*u[j-1]+ Sl*(uls3-rho[j-1]*e[j-1]);
			

			}
			else if(Sr>=0)
				{


					f1[j] = rho[j]*u[j]+ Sr*(urs1-rho[j]);
					f2[j] = rho[j]*u[j]*u[j]+p[j]+Sr*(urs2-rho[j]*u[j]);
					f3[j] = rho[j]*h[j]*u[j]+Sr*(urs3-rho[j]*e[j]);
				}
				else
				{
					f1[j] = rho[j]*u[j];
					f2[j] = rho[j]*u[j]*u[j]+p[j];
					f3[j] = rho[j]*h[j]*u[j];
				}
	}
	
	return 0;
}
