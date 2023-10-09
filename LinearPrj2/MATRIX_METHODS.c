#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MATRIX_METHODS.h"


//functions for convenience
double** allocateMemory(int m, int n) {
	double** A;
	A = (double**)malloc(sizeof(double*) * m);
	for (int i = 0; i < m; i++) {
		A[i] = (double*)malloc(sizeof(double) * n);
	}
	return A;
}


void releaseMemory(double** A, int m) {
	for (int i = 0; i < m; i++)
		free(A[i]);
	free(A);
}

void printMatrix(double** A, int m, int n, char name[]) {
	printf("\n%s = \n", name);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%.1lf ", A[i][j]);
		printf("\n");
	}
}

//functions to implement in prj0 
double** transposeMatrix(double** A, int m, int n) {
	double** B = allocateMemory(n, m);

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[j][i] = A[i][j];

	return B;
}

double** normalizeVector(double** v, int m) {
	double** w;
	double len = 0.0;

	for (int i = 0; i < m; i++)
		len += v[i][0] * v[i][0];
	len = sqrt(len);

	w = allocateMemory(m, 1);
	for (int i = 0; i < m; i++)
		w[i][0] = v[i][0] / len;

	return w;
}

double** constructIdentity(int k) {
	double** I = allocateMemory(k, k);

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			if (i != j)
				I[i][j] = 0.0;
			else
				I[i][j] = 1.0;
		}
	}
	return I;
}

double** DivideMatrix(double** A, int stx, int sty,int enx,int eny) {
	double** matrix = allocateMemory(enx-stx, eny-sty);
	for (int i = stx; i < enx; i++) {
		for (int j = sty; j < eny; j++) {
			matrix[i - stx][j - sty] = A[i][j];
		}
	}
	return matrix;

}
double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k) {
	if (n != l)
		return NULL;
	double** new_Matrix;
	int cnt = 0;
	new_Matrix = allocateMemory(m, k);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < k; j++) {
			double sum = 0.0;
			for (int a = 0; a < l; a++) {
				sum += ((double)A[i][a]) * ((double)B[a][j]);
			}
			if (isnan(sum)) {
				cnt++;
				//printMatrix(A,4,4,"A");
				//printMatrix(B, 4, 4, "B");
			}
			new_Matrix[i][j] = sum;
			//printf("%lf ", sum);
		}
	}
	//printf("NAN : %d\n", cnt);
	return new_Matrix;
}
double** normalizeEachColumn(double** v, int m, int n) {
	double** w;
	double len = 0.0;

	w = allocateMemory(m, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			len += v[j][i] * v[j][i];
		}
		for (int j = 0; j < n; j++) {
			if (len != 0) {
				w[j][i] = v[j][i] / sqrt(len);
			}
			else {
				w[j][i] = 0;
			}
		}
		len = 0.0;
	}
	return w;
}
double** constructHaarMatrixRecursive(int n) {
	double** h;
	if (n > 2)
		h = constructHaarMatrixRecursive(n / 2);
	else {
		//double** h;
		h = allocateMemory(2, 2);
		h[0][0] = 1; h[0][1] = 1; h[1][0] = 1; h[1][1] = -1; //H = [1 1; 1 -1]
		return h;
	}

	// construct left part (Kronecket product of h and [1;1]
	double** Hl = applyKroneckerProduct(h, n, 1, 1);
	releaseMemory(h, n / 2);

	// construct right part (Kronecker product of I and [1;-1]
	double** I = constructIdentity(n / 2);
	double** Hr = applyKroneckerProduct(I, n, 1, -1);
	releaseMemory(I, n / 2);


	// merge hl and hr
	double** H = concatenateTwoMatrices(Hl, Hr, n); //H = [hl hr]
	releaseMemory(Hl, n);
	releaseMemory(Hr, n);

	return H;
}

double** applyKroneckerProduct(double** A, int n, double a, double b) {
	double** h = allocateMemory(n, n / 2);

	for (int j = 0; j < n / 2; j++) {
		for (int i = 0; i < n / 2; i++) {
			h[2 * i][j] = A[i][j] * a;
			h[2 * i + 1][j] = A[i][j] * b;
		}
	}
	return h;
}

double** concatenateTwoMatrices(double** hl, double** hr, int n) {
	double** H = allocateMemory(n, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j < n / 2)
				H[i][j] = hl[i][j];
			else
				H[i][j] = hr[i][j - n / 2];
		}
	}
	return H;
}
