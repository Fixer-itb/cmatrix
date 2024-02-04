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
	printf("µ„≥ÀΩ·π˚Œ™£∫%f\n", vecdotres);
	double vecnormres = CvectorNormd(vec1, 3);
	printf("ƒ£÷µŒ™£∫%f\n", vecnormres);
	double vec1normlize[3] = { 0 };
	CvectorNormalized(&vec1, &vec1, 3);
	CmatrixShow(vec1, 3, 1, 'd');
	CmatrixShow(vec1, 3, 1, 'd');



/*****æÿ’Û≥À∑®≤‚ ‘*****/

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

/*****æÿ’Û◊™÷√≥…–¬æÿ’Û≤‚ ‘*****/
	double* matResultT = CMatrixd(3, 5);
	CmatrixT(matResultT, matResultown, 3, 5, 'd');
	CmatrixShow(matResultT, 3, 5, 'd');

/*****æÿ’ÛÃÓ≥‰◊”æÿ’Û≤‚ ‘*****/
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

/*****æÿ’Û‘≠µÿ◊™÷√≤‚ ‘*****/
	double elem4T[15] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 };
	double* matrix_row = CMatrixd(3, 5);
	CmatrixCopy(matrix_row, elem4T, 3, 5, 'd');
	CmatrixShow(matrix_row, 3, 5, 'd');
	CmatrixT_Situ(matrix_row, 3, 5);
	CmatrixShow(matrix_row, 5, 3, 'd');
	
/*****æÿ’Û«ÛƒÊ≤‚ ‘*****/
	double arrayInv[225] = {\
		1, 2, 1, 7, 7, 6, 7, 5, 6, 1, 2, 3, 6, 9, 6,\
		1, 5, 6, 9, 4, 5, 4, 3, 5, 1, 1, 4, 5, 8, 4,\
		6, 1, 2, 7, 6, 4, 8, 6, 2, 8, 7, 1, 6, 2, 8,\
		1, 1, 3, 4, 9, 1, 6, 5, 4, 3, 2, 5, 7, 1, 9,\
		1, 6, 9, 1, 3, 6, 2, 1, 3, 4, 1, 3, 8, 9, 6,\
		7, 7, 9, 9, 6, 7, 7, 8, 6, 2, 9, 5, 3, 9, 8,\
		6, 4, 1, 3, 8, 4, 4, 2, 8, 5, 4, 7, 4, 8, 5,\
		1, 6, 1, 6, 3, 9, 3, 7, 4, 9, 9, 9, 9, 3, 1,\
		3, 1, 1, 5, 7, 7, 5, 9, 2, 2, 6, 8, 9, 6, 1,\
		7, 2, 6, 3, 4, 4, 3, 2, 5, 8, 6, 9, 5, 6, 8,\
		6, 2, 1, 8, 8, 1, 8, 8, 5, 8, 6, 6, 8, 1, 4,\
		8, 9, 7, 1, 1, 8, 7, 8, 9, 4, 1, 6, 4, 9, 7,\
		6, 1, 3, 3, 6, 1, 2, 5, 3, 1, 8, 8, 8, 7, 9,\
		4, 5, 1, 9, 7, 5, 7, 8, 7, 3, 4, 9, 5, 5, 4,\
		4, 7, 6, 4, 2, 3, 9, 9, 5, 7, 1, 6, 7, 1, 1
	};
	double* matInv = CMatrixd(15, 15);
	double* matRes = CMatrixd(15, 15);
	CmatrixCopy(matInv, arrayInv, 15, 15, 'd');
	CmatrixShow(matInv, 15, 15, 'd');
	matinv(matRes, 15,matInv);
	CmatrixShow(matRes, 15, 15, 'd');



	return 1;
}
