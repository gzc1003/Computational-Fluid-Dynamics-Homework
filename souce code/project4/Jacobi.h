#include <stdio.h> 
#include <math.h> 
#define N 5

double Compare(double a[N],double b[N]) {
	double c=0;  int i; 
	  for(i=0;i<=N-1;i++)  
		     c+=fabs(a[i]-b[i]); 
	   return c; } 

void Jacobi(double A[N][N],double x[N],double b[N],double precesion) 
{ 
	int i,j,k; 
	double x2[N],sum; 
	for(i=0;i<=N-1;i++) 
		x2[i]=x[i];                     
	k = 1;

	while(1) { 
		 for(i=0;i<=N-1;i++) {  
			 sum=0; 
			 for(j=0;j<=N-1;j++)      
			 { 
				 if(j!=i) 
					 sum+=A[i][j]*x2[j];   
			 } 
			 x[i]=(b[i]-sum)/A[i][i];    
		  } 
		  
		  
		  if(Compare(x2,x)<=precesion) 
			    break;  
		  else { 
			    for(i=0;i<=N-1;i++)  
					 x2[i]=x[i];               
				k++; 
				continue;  
		  } 
	 } 
} 
