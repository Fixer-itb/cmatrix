#include"cmatrix.h"
int main() {
	double* mat33normal = CMatrixd(3, 3);
	double eg[9] = { 0,1,2,3,4,5,6,7,8 };
	CmatrixSet(mat33normal, eg, 3, 3, 'd');
	CmatrixShow(mat33normal, 3, 3, 'd');

	int* mat33zeros = CMatrixZerosi(3, 3);
	CmatrixShow(mat33zeros, 3,3, 'i');

	return 1;
}
