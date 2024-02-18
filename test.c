#include"cmatrix.h"
int main() {
	double* mat33normal = CMatrixd(3, 3);
	double eg[9] = { 0,1,2,3,4,5,6,7,8 };
	CMatrixCopy(mat33normal, eg, 3, 3, 'd');
	CMatrixShow(mat33normal, 3, 3, 'd');

	int* mat33zeros = CMatrixZerosi(3, 3);
	CMatrixShow(mat33zeros, 3,3, 'i');

	double* mat33eye = CMatrixEyed(3, 3);
	CMatrixShow(mat33eye, 3, 3, 'd');

	free(mat33eye);
	free(mat33normal);
	free(mat33zeros);

	// CmatrixShow(mat33normal, 3, 3, 'd');
	double vec1[3] = { 0, 1, 2 };
	double vec2[3] = { 0, 1, 2 };
	double vecdotres = CVectorDotd(vec1, vec2, 3);
	printf("µ„≥ÀΩ·π˚Œ™£∫%f\n", vecdotres);
	double vecnormres = CVectorNormd(vec1, 3);
	printf("ƒ£÷µŒ™£∫%f\n", vecnormres);
	double vec1normlize[3] = { 0 };
	CVectorNormalized(&vec1, &vec1, 3);
	CMatrixShow(vec1, 3, 1, 'd');
	CMatrixShow(vec1, 3, 1, 'd');



/*****æÿ’Û≥À∑®≤‚ ‘*****/

	//double elemL[6] = { 1,2,3,4,5,6 };
	//double* matL = CMatrixd(3, 2);
	//CMatrixCopy(matL, elemL, 3, 2, 'd');
	//CMatrixShow(matL, 3, 2, 'd');

	//double elemR[6] = { 7,8,9,10,11,12 };
	//double* matR = CMatrixd(2, 3);
	//CMatrixCopy(matR, elemR, 2, 3, 'd');
	//CMatrixShow(matR, 2, 3, 'd');

	//double* matResultown = CMatrixd(2, 2);
	//CMatrixMul("TT", 2, 3, 2, 1, matL, matR, 0.0, matResultown);
	//CMatrixShow(matResultown, 2, 2, 'd');
	
	/*double elemL2[9] = { 1,2,3,4,5,6,7,8,9 };
	double* matL2 = CMatrixd(3, 3);
	CMatrixCopy(matL2, elemL2, 3, 3, 'd');
	CMatrixShow(matL2, 3, 3, 'd');

	double elemR2[9] = { 1,2,3,4,5,6,7,8,9 };
	double* matR2 = CMatrixd(3, 3);
	CMatrixCopy(matR2, elemR2, 3, 3, 'd');
	CMatrixShow(matR2, 3, 3, 'd');

	double* matResultown2 = CMatrixd(3, 3);
	CMatrixMul("NN", 3, 3, 3, 1, matL2, matR2, 0.0, matResultown2);
	CMatrixShow(matResultown2, 3, 3, 'd');

	CMatrixMul("NT", 3, 3, 3, 1, matL2, matR2, 0.0, matResultown2);
	CMatrixShow(matResultown2, 3, 3, 'd');*/


	double elemL[20] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 };
	double* matL = CMatrixd(4, 5);
	CMatrixCopy(matL, elemL, 4, 5, 'd');
	CMatrixShow(matL, 4, 5, 'd');

	double elemR[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 }; 
	double* matR = CMatrixd(3, 4);
	CMatrixCopy(matR, elemR, 3, 4, 'd');
	CMatrixShow(matR, 3, 4, 'd');

	double* matResultown = CMatrixd(5, 3);
	CMatrixMul("TT", 5, 4, 3, 1, matL, matR, 0.0, matResultown);
	CMatrixShow(matResultown, 5, 3, 'd');

/*****æÿ’Û◊™÷√≥…–¬æÿ’Û≤‚ ‘*****/
	double* matResultT = CMatrixd(3, 5);
	CMatrixT(matResultT, matResultown, 3, 5, 'd');
	CMatrixShow(matResultT, 3, 5, 'd');

/*****æÿ’ÛÃÓ≥‰◊”æÿ’Û≤‚ ‘*****/
	double elementBlock[9] = { 0,1,2,3,4,5,6,7,8 };
	double elementBlock2[9] = { 1,1,1,1,1,1,1,1,1 };
	double* matBlock = CMatrixd(3, 3);
	double* matBlock2 = CMatrixd(3, 3);
	CMatrixCopy(matBlock, elementBlock, 3, 3, 'd');
	CMatrixCopy(matBlock2, elementBlock2, 3, 3, 'd');

	double* matLarge = CMatrixZerosd(10, 10);
	CMatrixBlockFill(matLarge, 10, 10, matBlock, 3, 3, 1, 1);
	CMatrixBlockFill(matLarge, 10, 10, matBlock2, 3, 3, 7, 7);
	CMatrixShow(matLarge,10,10,'d');

/*****æÿ’Û‘≠µÿ◊™÷√≤‚ ‘*****/
	double elem4T[15] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 };
	double* matrix_row = CMatrixd(3, 5);
	CMatrixCopy(matrix_row, elem4T, 3, 5, 'd');
	CMatrixShow(matrix_row, 3, 5, 'd');
	CMatrixT_Situ(matrix_row, 3, 5);
	CMatrixShow(matrix_row, 5, 3, 'd');
	
/*****æÿ’Û«ÛƒÊ≤‚ ‘*****/
	double arrayInv[225] = {\
		4, 7, 3, 2, 3, 7, 9, 8, 3, 4, 4, 9, 6, 5, 7,\
		2, 7, 7, 5, 2, 2, 6, 2, 8, 5, 6, 3, 5, 2, 9,\
		4, 9, 4, 9, 6, 8, 6, 6, 9, 5, 6, 8, 2, 8, 4,\
		4, 6, 1, 1, 6, 7, 2, 6, 8, 7, 7, 5, 3, 2, 8,\
		4, 8, 3, 7, 4, 8, 3, 8, 6, 2, 6, 7, 4, 6, 1,\
		8, 1, 9, 7, 6, 5, 2, 1, 3, 3, 4, 2, 9, 5, 8,\
		9, 9, 3, 9, 3, 2, 4, 6, 6, 6, 4, 9, 1, 3, 7,\
		4, 3, 1, 9, 6, 9, 4, 1, 6, 6, 1, 9, 9, 5, 1,\
		8, 3, 9, 9, 4, 6, 4, 6, 2, 7, 2, 2, 6, 6, 4,\
		9, 5, 1, 6, 6, 3, 4, 4, 6, 7, 5, 4, 5, 6, 9,\
		7, 5, 3, 8, 9, 5, 1, 7, 3, 1, 9, 5, 2, 7, 8,\
		8, 4, 5, 9, 5, 7, 7, 5, 2, 9, 8, 7, 8, 3, 4,\
		5, 2, 8, 3, 8, 5, 2, 9, 6, 7, 1, 5, 3, 2, 4,\
		1, 4, 5, 9, 3, 5, 2, 1, 7, 9, 4, 2, 4, 9, 3,\
		7, 4, 9, 2, 7, 9, 1, 9, 3, 6, 2, 2, 8, 6, 3
	};
	double* matInv = CMatrixd(15, 15);
	double* matRes = CMatrixd(15, 15);
	CMatrixCopy(matInv, arrayInv, 15, 15, 'd');
	CMatrixShow(matInv, 15, 15, 'd');
	CMatrixInv(matInv, 15,matRes);
	CMatrixShow(matInv, 15, 15, 'd');



	return 1;
}
