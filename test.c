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



/*****矩阵乘法测试*****/

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

/*****矩阵转置成新矩阵测试*****/
	double* matResultT = CMatrixd(3, 5);
	CmatrixT(matResultT, matResultown, 3, 5, 'd');
	CmatrixShow(matResultT, 3, 5, 'd');

/*****矩阵填充子矩阵测试*****/
	double elementBlock[9] = { 0,1,2,3,4,5,6,7,8 };
	double elementBlock2[9] = { 1,1,1,1,1,1,1,1,1 };
	double* matBlock = CMatrixd(3, 3);
	double* matBlock2 = CMatrixd(3, 3);
	CmatrixCopy(matBlock, elementBlock, 3, 3, 'd');
	CmatrixCopy(matBlock2, elementBlock2, 3, 3, 'd');

	double* matLarge = CMatrixZerosd(10, 10);
	CmatrixBlockFill(matLarge, 10, 10, matBlock, 3, 3, 1, 1);
	CmatrixBlockFill(matLarge, 10, 10, matBlock2, 3, 3, 7, 7);
	CmatrixShow(matLarge,10,10,'d');

/*****矩阵原地转置测试*****/
	double elem4T[15] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 };
	double* matrix_row = CMatrixd(3, 5);
	CmatrixCopy(matrix_row, elem4T, 3, 5, 'd');
	CmatrixShow(matrix_row, 3, 5, 'd');
	CmatrixT_Situ(matrix_row, 3, 5);
	CmatrixShow(matrix_row, 5, 3, 'd');
	



	return 1;
}
