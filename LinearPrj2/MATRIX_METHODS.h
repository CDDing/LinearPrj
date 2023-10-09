//functions for convenience
double** allocateMemory(int m, int n);
void releaseMemory(double** A, int m);
void printMatrix(double** A, int m, int n, char name[]);
double** constructIdentity(int k);

//functions to implement in prj0 
double** transposeMatrix(double** A, int m, int n);
double** normalizeVector(double** v, int n);
double** normalizeEachColumn(double** v, int m, int n);
double** DivideMatrix(double** A, int stx, int sty, int enx, int eny);
double** multiplyTwoMatrices(double** A, int m, int n, double** B, int l, int k);


double** constructHaarMatrixRecursive(int n);
double** concatenateTwoMatrices(double** hl, double** hr, int n);
double** applyKroneckerProduct(double** A, int n, double a, double b);