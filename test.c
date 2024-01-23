#include"cmatrix.h"
int main() {
	double* mat33normal = CMatrixd(3, 3);
	double eg[9] = { 0,1,2,3,4,5,6,7,8 };
	CmatrixCopy(mat33normal, eg, 3, 3, 'd');
	CmatrixShow(mat33normal, 3, 3, 'd');

	int* mat33zeros = CMatrixZerosi(3, 3);
	CmatrixShow(mat33zeros, 3,3, 'i');

	double* mat33eye = CMatrixEyed(3, 3);
	CmatrixShow(mat33eye, 3, 3, 'd');

	free(mat33eye);
	free(mat33normal);
	free(mat33zeros);

	// CmatrixShow(mat33normal, 3, 3, 'd');
	double vec1[3] = { 0, 1, 2 };
	double vec2[3] = { 0, 1, 2 };
	double vecdotres = CvectorDotd(vec1, vec2, 3);
	printf("点乘结果为：%f\n", vecdotres);
	double vecnormres = CvectorNormd(vec1, 3);
	printf("模值为：%f\n", vecnormres);
	double vec1normlize[3] = { 0 };
	CvectorNormalized(&vec1, &vec1, 3);
	CmatrixShow(vec1, 3, 1, 'd');
	CmatrixShow(vec1, 3, 1, 'd');

	double elemL[6] = { 1,2,3,4,5,6 };
	double* matL = CMatrixd(2, 3);
	CmatrixCopy(matL, elemL, 2, 3, 'd');
	CmatrixShow(matL, 2, 3, 'd');

	double elemR[6] = { 7,8,9,10,11,12 };
	double* matR = CMatrixd(3, 2);
	CmatrixCopy(matR, elemR, 3, 2, 'd');
	CmatrixShow(matR, 3, 2, 'd');

	double* matResultown = CMatrixd(2, 2);
	CmatrixMul("NN", 2, 3, 2, 1, matL, matR, 0.0, matResultown);
	CmatrixShow(matResultown, 2, 2, 'd');


	double* matResultignav = CMatrixd(2, 2);
	matmul("NN", 2, 2, 3, 1, matL, matR, 0.0, matResultignav);
	CmatrixShow(matResultignav, 2, 2, 'd');


	/*
	>>Matrix addr:6d9ce560
58.0000000000   64.0000000000
139.0000000000  154.0000000000
>>Matrix addr:6d9ce5c0
76.0000000000   100.0000000000
103.0000000000  136.0000000000
	*/
	return 1;
}
