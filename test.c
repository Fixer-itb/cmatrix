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



/***********************矩阵乘法测试**************************/

	//double elemL[6] = { 1,2,3,4,5,6 };
	//double* matL = CMatrixd(3, 2);
	//CmatrixCopy(matL, elemL, 3, 2, 'd');
	//CmatrixShow(matL, 3, 2, 'd');

	//double elemR[6] = { 7,8,9,10,11,12 };
	//double* matR = CMatrixd(2, 3);
	//CmatrixCopy(matR, elemR, 2, 3, 'd');
	//CmatrixShow(matR, 2, 3, 'd');

	//double* matResultown = CMatrixd(2, 2);
	//CmatrixMul("TT", 2, 3, 2, 1, matL, matR, 0.0, matResultown);
	//CmatrixShow(matResultown, 2, 2, 'd');


	double elemL[20] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 };
	double* matL = CMatrixd(4, 5);
	CmatrixCopy(matL, elemL, 4, 5, 'd');
	CmatrixShow(matL, 4, 5, 'd');

	double elemR[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 }; 
	double* matR = CMatrixd(3, 4);
	CmatrixCopy(matR, elemR, 3, 4, 'd');
	CmatrixShow(matR, 3, 4, 'd');

	double* matResultown = CMatrixd(5, 3);
	CmatrixMul("TT", 5, 4, 3, 1, matL, matR, 0.0, matResultown);
	CmatrixShow(matResultown, 5, 3, 'd');

	/*double elemL2[9] = { 1,2,3,4,5,6,7,8,9 };
	double* matL2 = CMatrixd(3, 3);
	CmatrixCopy(matL2, elemL2, 3, 3, 'd');
	CmatrixShow(matL2, 3, 3, 'd');

	double elemR2[9] = { 1,2,3,4,5,6,7,8,9 };
	double* matR2 = CMatrixd(3, 3);
	CmatrixCopy(matR2, elemR2, 3, 3, 'd');
	CmatrixShow(matR2, 3, 3, 'd');

	double* matResultown2 = CMatrixd(3, 3);
	CmatrixMul("NN", 3, 3, 3, 1, matL2, matR2, 0.0, matResultown2);
	CmatrixShow(matResultown2, 3, 3, 'd');

	CmatrixMul("NT", 3, 3, 3, 1, matL2, matR2, 0.0, matResultown2);
	CmatrixShow(matResultown2, 3, 3, 'd');*/

	//int matrix[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
	//int row = 3, column = 4;
	///*
	//0  1  2  3
	//4  5  6  7
	//8  9 10 11
	//*/
	//int i = 0, j = 0;
	//for (i = 0; i < row; i++) {
	//	for (j = 0; j < column; j++) {
	//		printf("The [%d,%d] element is %d through row\n", i, j, matrix[i * column + j]);
	//	}
	//}
	//for (i = 0; i < column; i++) {
	//	for (j = 0; j < row; j++) {
	//		printf("The [%d,%d] element is %d through column\n", i, j, matrix[j*column + i]);
	//	}
	//}



	return 1;
}
