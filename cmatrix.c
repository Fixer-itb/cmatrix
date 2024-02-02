/*****************************************************************//**
 * \file   matrix.c
 * \brief  C���Ծ��������
 *
 * \author ��Ȼ����
 * \date   January 2024
 *********************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include"cmatrix.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
 /**
  * @brief ��������������Ϊ1ά������� n*m.
  *
  * @param row ������
  * @param column ������
  * @return �����ַ
  */
extern double* CMatrixd(int row, int column)
{
	double* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (double*)malloc(sizeof(double) * row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}
extern float* CMatrixf(int row, int column)
{
	float* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (float*)malloc(sizeof(float) * row * column))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}
extern int* CMatrixi(int row, int column)
{
	int* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (int*)malloc(sizeof(int) * row * column))) {
		printf("integer matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}



/**
 * @brief ����ȫ0����������Ϊ1ά������� n*m.
 *
 * @param row ������
 * @param column ������
 * @return �����ַ
 */
extern double* CMatrixZerosd(int row, int column) {
	double* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (double*)calloc(row * column, sizeof(double))))
	{
		printf("Double matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}
extern float* CMatrixZerosf(int row, int column)
{
	float* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (float*)calloc(row * column, sizeof(float))))
	{
		printf("Float matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}
extern int* CMatrixZerosi(int row, int column)
{
	int* p;

	if (row <= 0 || column <= 0) return NULL;
	if (!(p = (int*)calloc(row * column, sizeof(int))))
	{
		printf("integer matrix memory allocation error: n=%d,m=%d\n", row, column);
	}
	return p;
}

/**
 * @brief ������λ����.
 *
 * @param rank ��λ�������
 * @return �����ַ
 */
extern double* CMatrixEyed(int rank)
{
	double* p;
	int i;

	if ((p = CMatrixZerosd(rank, rank))) for (i = 0; i < rank; i++) p[i + i * rank] = 1.0;
	return p;
}
extern float* CMatrixEyef(int rank)
{
	float* p;
	int i;

	if ((p = CMatrixZerosf(rank, rank))) for (i = 0; i < rank; i++) p[i + i * rank] = 1.0;
	return p;
}
extern int* CMatrixEyei(int rank)
{
	int* p;
	int i;

	if ((p = CMatrixZerosi(rank, rank))) for (i = 0; i < rank; i++) p[i + i * rank] = 1.0;
	return p;
}

/**
 * @brief ���ݾ���Ԫ������type������ӡ��ͬ���ȵľ���.
 *
 * @param matrix    �����ַ
 * @param row   ��������
 * @param column    ��������
 * @param type  ����Ԫ������
 * @return  ��ӡ�Ƿ�ɹ�
 */
extern int CmatrixShow(void* matrix, int row, int column, char type) {

	if (matrix == NULL) {
		printf("Invalid matrix pointer\n");
		return 0;
	}
	printf(">>Matrix addr:%x\n", matrix);
	//for (int k = 0; k < row * column; k++)
	//{
	//	printf(" %0.5f", ((double*)matrix)[k]);
	//}
	//printf("\n");
	int i, j;
	for (i = 0; i < row; i++) {
		for (j = 0; j < column; j++) {
			switch (type) {
			case 'i':
				printf("%d\t", ((int*)matrix)[i * column + j]);
				break;
			case 'f':
				printf("%.6f\t", ((float*)matrix)[i * column + j]);
				break;
			case 'd':
				printf("%.10f\t", ((double*)matrix)[i * column + j]);
				break;
			default:
				printf("%.6f\t", ((float*)matrix)[i * column + j]);
				return 0;
			}
		}
		printf("\n");
	}
}

/**
 * @brief ����/���鸴��.
 *
 * @param matrix_dest Ŀ�����
 * @param matrix_src	Դ����
 * @param row	������
 * @param column ������
 * @param type ����Ԫ������
 */
extern void CmatrixCopy(void* matrix_dest, const void* matrix_src, int row, int column, char type) {
	int elementSize;
	switch (type) {
	case 'i':
		elementSize = sizeof(int);
		break;
	case 'f':
		elementSize = sizeof(float);
		break;
	case 'd':
		elementSize = sizeof(double);
		break;
	default:
		printf("Unknow type!\n");
		return;
	}
	memcpy(matrix_dest, matrix_src, row * column * elementSize);
}

/**
 * @brief ��������nά���������ڻ�.
 *
 * @param vec1 ����1
 * @param vec2 ����2
 * @param dimension ά��
 * @return �ڻ���ֵ
 */
extern double CvectorDotd(const double* vec1, const double* vec2, int dimension) {
	double result = 0.0;
	while (--dimension >= 0) result += vec1[dimension] * vec2[dimension];
	return result;
}
extern float CvectorDotf(const float* vec1, const float* vec2, int dimension) {
	float result = 0.0;
	while (--dimension >= 0) result += vec1[dimension] * vec2[dimension];
	return result;
}
extern int CvectorDoti(const int* vec1, const int* vec2, int dimension) {
	int result = 0;
	while (--dimension >= 0) result += vec1[dimension] * vec2[dimension];
	return result;
}

/**
 * @brief ��ά�������.
 *
 * @param vecl ������
 * @param vecr ������
 * @param result ���
 */
extern void Cvector3dCrossd(const double* vecl, const double* vecr, double* result)
{
	result[0] = vecl[1] * vecr[2] - vecl[2] * vecr[1];
	result[1] = vecl[2] * vecr[0] - vecl[0] * vecr[2];
	result[2] = vecl[0] * vecr[1] - vecl[1] * vecr[0];
}
extern void Cvector3dCrossf(const float* vecl, const float* vecr, float* result)
{
	result[0] = vecl[1] * vecr[2] - vecl[2] * vecr[1];
	result[1] = vecl[2] * vecr[0] - vecl[0] * vecr[2];
	result[2] = vecl[0] * vecr[1] - vecl[1] * vecr[0];
}
extern void Cvector3dCrossi(const int* vecl, const int* vecr, int* result)
{
	result[0] = vecl[1] * vecr[2] - vecl[2] * vecr[1];
	result[1] = vecl[2] * vecr[0] - vecl[0] * vecr[2];
	result[2] = vecl[0] * vecr[1] - vecl[1] * vecr[0];
}

/**
 * @brief ��nά�����Ķ�������ģ��.
 *
 * @param vec ����
 * @param dimension ����ά��
 * @return ģֵ
 */
extern double CvectorNormd(const double* vec, int dimension)
{
	double result = 0.0;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	return sqrt(result);
}
extern float CvectorNormf(const float* vec, int dimension)
{
	float result = 0.0;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	return sqrt(result);
}
extern int CvectorNormi(const int* vec, int dimension)
{
	int result = 0;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	return sqrt(result);
}

/**
 * @brief ������һ��.
 *
 * @param vec ��һ��ǰ������
 * @param vecnormalized ��һ���������
 * @param dimension ά��
 */
extern void CvectorNormalized(const double* vec, double* vecnormalized, int dimension)
{
	double result = 0.0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}
extern void CvectorNormalizef(const float* vec, float* vecnormalized, int dimension)
{
	float result = 0.0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}
extern void CvectorNormalizei(const int* vec, int* vecnormalized, int dimension)
{
	int result = 0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}

/**
 * @brief ����˷����� result=alpha*MAT_LEFT*MAT_RIGHT+beta*RESULT(if beta != 0.0 & RESULT != NULL).\
 *				   result=alpha*MAT_LEFT*MAT_RIGHT(if beta = 0.0).
 *
 * @param tr	�Ƿ�������������ת��(N,��ת��), (T,ת��)
 * @param rowL ��߾��������
 * @param midLR ��߾�����������ұ߾��������
 * @param columnR �ұ߾��������
 * @param alpha ���Ҿ����ϵ��(����)
 * @param matl ������ַ
 * @param matr �Ҿ����ַ
 * @param beta ���������Ϊ���о���ʱ��ϵ��
 * @param result ������
 */
extern void CmatrixMul(const char* tr, int rowL, int midLR, int columnR, double alpha, const double* matl, const double* matr, double beta, double* result)
{
	double temp;
	int i, j, k, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);
	/*
		if (tr[0]=='N') {
			if (tr[1] == 'N') f = 1; "NN"
			else f = 2;              "NT"
		}
		else {
			if (tr[1] == 'N') f = 3; "TN"
			else f = 4;              "TT"
		}
	*/
	//!ֻ��[i,k]��ȡmatLԪ�أ�ֻ��[k,j]��ȡmatRԪ��
	for (i = 0; i < rowL; i++) {
		for (j = 0; j < columnR; j++) {
			temp = 0.0;
			switch (f) {
			case 1: //"NN"  A*B   matl[i,k] �� matr[k,j]	//A�� = B�� = A�г�B�� = midLR
				for (k = 0; k < midLR; k++) {
					temp += matl[i * midLR + k] * matr[k * columnR + j];
				}
				break;
			case 2: //"NT"  A*B'  matl[i,k] �� matr[j,k]' //�ܼ����ǰ����[����]һ����һ���� = A�г�B�� = midLR
				for (k = 0; k < midLR; k++) {
					//temp += matl[i * midLR + k] * matr[j * midLR + k];
					temp += matl[i * midLR + k] * matr[j * midLR + k];
				}
				break;
			case 3: //"TN"  A'*B  matl[k,i]' �� matr[k,j]  //�ܼ����ǰ����[����]һ����һ���� = A�г�B�� = midLR
				for (k = 0; k < midLR; k++) {
					temp += matl[k * rowL + i] * matr[k * columnR + j];
				}
				break;
			case 4: //"TT"  A'*B' matl[k,i]' �� matr[j,k]' //�ܼ����ǰ����A������=B������ = A�г�B�� = midLR
				for (k = 0; k < midLR; k++) {
					temp += matl[k * rowL + i] * matr[j * midLR + k];
				}
				break;
			}
			if (beta == 0.0) result[i * columnR + j] = alpha * temp;
			else result[i * columnR + j] = alpha * temp + beta * result[i * columnR + j];
		}
	}
}

/**
 * @brief ʵ�־���ת��.
 *
 * @param matrixT_dest Ŀ�����(ת�ú�ľ���)
 * @param matrix_src ԭ����
 * @param matrixT_row ת�ú���������
 * @param matrixT_column ת�ú���������
 * @param type ����Ԫ�ص�����
 */
extern void CmatrixT(void* matrixT_dest, const void* matrix_src, int matrixT_row, int matrixT_column, char type) {
	int i, j;
	for (i = 0; i < matrixT_row; i++) {
		for (j = 0; j < matrixT_column; j++) {
			switch (type) {
			case 'i':
				((int*)matrixT_dest)[i * matrixT_column + j] = ((int*)matrix_src)[j * matrixT_row + i];
				break;
			case 'f':
				((float*)matrixT_dest)[i * matrixT_column + j] = ((float*)matrix_src)[j * matrixT_row + i];
				break;
			case 'd':
				((double*)matrixT_dest)[i * matrixT_column + j] = ((double*)matrix_src)[j * matrixT_row + i];
				break;
			default:
				((double*)matrixT_dest)[i * matrixT_column + j] = ((double*)matrix_src)[j * matrixT_row + i];
				return 0;
			}
		}
	}
}

/**
 * @brief ����ԭ��ת�ò�����.
 *
 * @param matrix ��ת�þ���
 * @param row ��������
 * @param column ��������
 * @param type ����Ԫ������
 */
extern void CmatrixT_Situ(double* matrix, int row, int col) {
	//TODO TODO
	int nextNode, i;
	double temp;
	for (i = 0; i < row*col; i++){
		nextNode = (i % col) * row + (i / col);
		while (i<nextNode) {
			nextNode = (nextNode % col) * row + (nextNode / col);
		}
		if (i==nextNode){
			nextNode = (i % col) * row + (i / col);
			while (nextNode!=i){
				temp = matrix[i];
				matrix[i] = matrix[nextNode];
				matrix[nextNode] = temp;
				nextNode = (nextNode % col) * row + (nextNode / col);
			}
		}
	}
	return 0;
}

/**
 * @brief ��һ��С�������һ��������У�С�������ϽǶ����ڴ�����(fillLocRow,fillLocCol)��.
 *
 * @param destMat ������ľ���(�����)
 * @param destRow	��������������
 * @param destCol ��������������
 * @param blockMat	�������
 * @param blockRow	������������
 * @param blockCol	������������
 * @param fillLocRow	���Ͻ�����
 * @param fillLocCol	���Ͻ�����
 */
extern void CmatrixBlockFill(double* destMat, int destRow, int destCol, const double* blockMat, int blockRow, int blockCol, int fillLocRow, int fillLocCol) {
	int i, j;
	for (i = fillLocRow; i < fillLocRow + blockRow; i++) {//�������У��Ӳ����п�ʼ
		for (j = fillLocCol; j < fillLocCol + blockCol; j++) {//�������У��Ӳ����п�ʼ
			destMat[i * destCol + j] = blockMat[(i - fillLocRow) * blockCol + (j - fillLocCol)];//���б���
		}
	}
}


/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double* A, int n, int* indx, double* d)
{
	double big;	//��ʱ�洢���ֵ
	double s;	//�洢�м������
	double tmp;	//��ʱ����
	double *vv = CMatrixd(n, 1); //�����������洢ÿһ�е���������
	int i, imax = 0, j, k;

	*d = 1.0;	// ��ʼ������ʽ��ֵΪ0
	
	//����ÿһ�е��������ӣ����ڲ�����Ԫ��ȥ
	for (i = 0; i < n; i++) {
		big = 0.0; 
		//�ҵ�ÿ������Ԫ��
		for (j = 0; j < n; j++) {
			if ((tmp = fabs(A[i * n + j])) > big) { //!�������ȴ洢��Ϊ�������ȴ洢
				big = tmp;
			}
		}

		if (big > 0.0) {
			vv[i] = 1.0 / big;	// ÿ�����ֵ�ĵ�����Ϊ��������
		}
		else { // �����ǰ�е�����Ԫ�ض�Ϊ0��˵������������ģ��޷�����LU�ֽ�
			free(vv); 
			return -1; 
		}
	}
	//LU�ֽ���ѭ��,���н��б���������Ԫ��������A�ֽ�������Ǿ���U�������Ǿ���L

	//����������
	for (j = 0; j < n; j++) {
		//������ǰ��j�������ǲ��֣�ִ�и�˹��Ԫ�����Խ�������Ԫ������
		//����RTKlib���õ������ȴ洢����˴�ʱ���������ǣ��������ϵ�д��ʵ������������
		for (i = 0; i < j; i++) {
			s = A[i + j * n]; 
			for (k = 0; k < i; k++) 
				s -= A[i + k * n] * A[k + j * n]; 
			A[i + j * n] = s;
		}
		//������ǰ��j�������ǲ��֣�ִ�и�˹��Ԫ�����Խ�������Ԫ������
		//�ҵ���ǰ���о���ֵ����Ԫ�أ�����¼�������� imax��
		big = 0.0;
		for (i = j; i < n; i++) {
			s = A[i + j * n]; 
			for (k = 0; k < j; k++) 
				s -= A[i + k * n] * A[k + j * n]; 
			A[i + j * n] = s;
			if ((tmp = vv[i] * fabs(s)) >= big) { 
				big = tmp; 
				imax = i; 
			}
		}
		//����ҵ������Ԫ�ص������� imax �����ڵ�ǰ������ j��������н�����
		//�н�����ı���������ʽ�ķ��ţ���˸�������ʽ��ֵ* d��
		//ͬʱ��������������������ȷ����ȷ��������������Ԫ��Ӧ��
		if (j != imax) {
			for (k = 0; k < n; k++) {
				tmp = A[imax + k * n]; A[imax + k * n] = A[j + k * n]; A[j + k * n] = tmp;
			}
			*d = -(*d); 
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j + j * n] == 0.0) { 
			free(vv); 
			return -1; 
		}
		//��¼��Ԫ�������� imax�����ں����ļ����л��õ���
		//���Խ�Ԫ���Ƿ�Ϊ�㣬����ǣ�LU�ֽ�ʧ�ܣ��ͷ����������������ڴ棬���� - 1��ʾʧ�ܡ�
		//�����ǰ����Ĳ������һ�У�����Ԫ���µ��н��й�һ����
		if (j != n - 1) {
			tmp = 1.0 / A[j + j * n]; 
			for (i = j + 1; i < n; i++) 
				A[i + j * n] *= tmp;
		}
	}
	free(vv);
	return 0;
}

/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double* A, int n, const int* indx, double* b)
{
	double s;
	int i, ii = -1, ip, j;

	for (i = 0; i < n; i++) {
		ip = indx[i]; 
		s = b[ip]; 
		b[ip] = b[i];
		if (ii >= 0) {
			for (j = ii; j < i; j++) {
				s -= A[j + i * n] * b[j];
			}
		}
		else if (s) {
			ii = i;
		}
		b[i] = s;
	}
	for (i = n - 1; i >= 0; i--) {
		s = b[i]; 
		for (j = i + 1; j < n; j++) {
			s -= A[i + j * n] * b[j];
			b[i] = s / A[i + i * n];
		}
	}
}
/* inverse of matrix ---------------------------------------------------------*
*inverse of matrix(A = A ^ -1)
* args   : double* A        IO  matrix(n x n)
* int    n         I   size of matrix A
* return : status(0:ok, 0 > :error)
* ---------------------------------------------------------------------------- - */
extern int matinv(double* A, int n)
{
	double d, * B;
	int i, j, * indx;

	indx = CMatrixi(n, 1); B = CMatrixd(n, n); CmatrixCopy(B, A, n, n,'d');
	if (ludcmp(B, n, indx, &d)) { free(indx); free(B); return -1; }
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) A[i + j * n] = 0.0;
		A[j + j * n] = 1.0;
		lubksb(B, n, indx, A + j * n);
	}
	free(indx); free(B);
	return 0;
}










