#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "INPUTX.h"
#include "INPUTH.h"
#include "time.h"

#define PI	M_PI	
#define TWOPI	(2.0*PI)      

void modelfft(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}

void windowingFilter(double *a, double *b, double *res, int N){
    int i;
    for (i = 0; i < N; i++) {
        res[i] = a[i]*b[i];
    }
}

void preFFT(double *x, double *X, int N){
    int i;
    for (i = 0; i < N; i++) {
        X[2*i+1] = x[i];
        X[2*i+2] = 0.0;
    }
}

void postIFFT(double *tempy1, double *y1, int N){
    int i;
    for (i = 0; i < N; i++) {
        y1[i] = tempy1[2*i+1]/N;
    }
}

double absComplex(double real, double imaj){
    return(sqrt(pow(real,2)+pow(imaj,2)));
}

void NoiseSpecterEstimate(double *D, double *X1, int P, int N){
    int i;
    for(i=0; i<N; i++){
        D[i] += (absComplex(X1[2*i+1],X1[2*i+2])/P);
    }
}

void SSMethod(double *S, double *G, double *D, double *X1, int N){
    int i;
    double temp;
    for(i=0; i<N; i++){
        temp = 1 - (D[i]/absComplex(X1[2*i+1],X1[2*i+2]));
        G[i] = (temp+fabs(temp))/2;
        S[2*i+1] = G[i]*X1[2*i+1];
        S[2*i+2] = G[i]*X1[2*i+2];
    }
}

void inity(double *y, int it){
    int i;
    for(i=0; i<it; i++){
        y[i] = x[i];
    }
}

void cascadeOUT(double *y, double *past_tail, double *y1, int idx, int N){
    int i;
    int j = 0;
    for(i=idx; i<(idx+(N/2)); i++){
        y[i] = past_tail[j]+y1[j];
        j++;
    }
}


int main(int argc, char **argv){
    
    int N = 64;
    int M = 28815;
    int O = 28768;
    int N2 = N/2;
    int L = (M/N2);
    int P = 4;
    clock_t start;
    clock_t stop;   
    
    int i;
    int j = 0;
    int k = 0;
    int n;
    int idx;
    
    double x1[N];
    double y1[N];
    double y[O];
    double *temparray = (double *) malloc(sizeof(double) * N);
    double *G = (double *) malloc(sizeof(double) * N);
    double *D = (double *) calloc(N,sizeof(double));
    double *past_tail = (double *) calloc(N2,sizeof(double));
    double X1[2*N+1];
    double S[2*N+1];
    double tempy1[2*N+1];

    
    for (k = 0; k < P; k++) {
        n = N2*k;
        memcpy(temparray,x+(n),sizeof(double)*N);
        windowingFilter(h,temparray,x1,N);
        preFFT(x1,X1,N);
        modelfft(X1, N, 1);
        NoiseSpecterEstimate(D,X1,P,N);
    }
    
    inity(y,P*N2);
    idx = P*N2;

    for (k = P; k < L-1; k++) {
        n = N2*k;
        memcpy(temparray,x+(n),sizeof(double)*N);
        windowingFilter(h,temparray,x1,N);
        
        preFFT(x1,X1,N);
        modelfft(X1, N, 1);
        
	start = clock();
        SSMethod(S,G,D,X1,N);
        stop = clock();


	printf("Time start: %f\n", (double)start);
	printf("Time stop: %f\n", (double)stop);

	printf("Time elapsed: %f\n", ((double)stop - start) / CLOCKS_PER_SEC);

        memcpy(tempy1,S,sizeof(double)*2*N+1);
        modelfft(tempy1, N, -1);
        postIFFT(tempy1,y1,N);
        
        cascadeOUT(y,past_tail,y1,idx,N);
        idx+=(N/2);
        
        j = 0;
        for (i = N2;i < N;i++) {
            past_tail[j] = y1[i];
            j++;
        }
    }
    
    // FILE *fp;
    // fp = fopen("result.dat","w");
	// for(i=0; i<O; i++) {
        // fprintf(fp,"%.8f\n",y[i]);
    // }
    // fclose(fp);
    
    // for(i=O-100; i<O; i++) {
        // printf("%.8f\n",y[i]);
    // }
    
    free(temparray);
    free(G);
    free(D);
    free(past_tail);
    
    return 0;
    
}
