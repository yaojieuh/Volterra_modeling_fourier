#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void qrdcmp(double *a, int n, double *R, double *Q)
/*Constructs the QR decomposition of a[1..n][1..n] . The upper triangular matrix R is re-
turned in the upper triangle of a, except for the diagonal elements of R which are returned in
d[1..n]. The orthogonal matrix Q is represented as a product of n − 1 Householder matrices
Q1 . . . Qn−1 , where Qj = 1 − uj ⊗ uj /cj . The ith component of uj is zero for i = 1, . . . , j − 1
while the nonzero components are returned in a[i][j] for i = j, . . . , n. sing returns as
true (1) if singularity is encountered during the decomposition, but the decomposition is still
completed in this case; otherwise it returns false (0).*/
{
	int i,j,m,k;
	float scale,sigma,sum,tau;
	
        double* a1 =malloc(n*n*sizeof(double) );
	double* a2 =malloc(n*n*sizeof(double) );
	double* Iden =malloc(n*n*sizeof(double) );
	double* Hn =malloc(n*n*sizeof(double) );
	double* Qo =malloc(n*n*sizeof(double) );
	double* Qn =malloc(n*n*sizeof(double) );
        for (i=0;i<n*n;i++){
		 a1[i]=a[i];
		 Iden[i]=0;
		 Qo[i]=0;
	}
	for (i=0;i<n;i++){		
		 Iden[i*n+i]=1;
		 Qo[i*n+i]=1;
	}
	
	double* b =malloc(n*sizeof(double) );
	double* e =malloc(n*sizeof(double) );
	double* u =malloc(n*sizeof(double) );
	double* v =malloc(n*sizeof(double) );
        double unorm,vnorm;
	for (k=0;k<n-1;k++) {
		scale=0.0;
                unorm=0;
                vnorm=0;
		for (i=0;i<n;i++){
			b[i]=0;
			e[i]=0;
			u[i]=0;
			v[i]=0;
		}
		e[k]=1;
		for (i=k;i<n;i++) {
			b[i]=a1[i*n+k];
			scale+=b[i]*b[i];
		}
		
		if (scale == 0.0) {				
			
		} else {
			for (i=k;i<n;i++){
				u[i]=b[i]+SIGN(sqrt(scale),a[k*n+k])*e[i];
				unorm+=u[i]*u[i];
			}
			if (unorm == 0.0) {				
			
			} else {
			for (i=k;i<n;i++){
				v[i]=u[i]/sqrt(unorm);
				vnorm+=v[i]*v[i];
			}
			for (i=0;i<n;i++){
				for(j=0;j<n;j++){
					Hn[i*n+j]=Iden[i*n+j]-2*v[i]*v[j]/sqrt(vnorm);
					//fprintf(stdout,"%lf   \n",Hn[i*n+j]);
				}
			}
			for (i=0;i<n;i++){
				for(j=0;j<n;j++){
					a2[i*n+j]=0;
					Qn[i*n+j]=0;
					for(m=0;m<n;m++){
						a2[i*n+j]+=Hn[i*n+m]*a1[m*n+j];
						Qn[i*n+j]+=Hn[i*n+m]*Qo[m*n+j];
					}
				
				}
			}
			
			 for (i=0;i<n*n;i++){
		 		a1[i]=a2[i];
				Qo[i]=Qn[i];
				//fprintf(stdout,"%lf   \n",Qo[i]);
			}
			}
			
		}
		
		
}
	
	for (i=0;i<n;i++){
				for(j=0;j<n;j++){
					R[i*n+j]=a1[i*n+j];
					Q[i*n+j]=Qo[j*n+i];
					//fprintf(file1,"%d %d %lf  %lf   %lf\n",i,j,a[i*n+j],R[i*n+j],Q[i*n+j] );
				
				}
				//fprintf(file1," \n");
			}
free(a1);
free(a2);
free(Qo);
free(Qn);	
}
