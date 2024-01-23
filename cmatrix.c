/*****************************************************************//**
 * \file   matrix.c
 * \brief  C语言矩阵运算库
 *
 * \author 既然匠子
 * \date   January 2024
 *********************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include"cmatrix.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
 /**
  * @brief 创建矩阵，物理上为1维，抽象成 n*m.
  *
  * @param row 矩阵行
  * @param column 矩阵列
  * @return 矩阵地址
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
 * @brief 创建全0矩阵，物理上为1维，抽象成 n*m.
 *
 * @param row 矩阵行
 * @param column 矩阵列
 * @return 矩阵地址
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
 * @brief 创建单位矩阵.
 *
 * @param rank 单位矩阵的秩
 * @return 矩阵地址
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
 * @brief 根据矩阵元素类型type参数打印不同精度的矩阵.
 *
 * @param matrix    矩阵地址
 * @param row   矩阵行数
 * @param column    矩阵列数
 * @param type  矩阵元素类型
 * @return  打印是否成功
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
 * @brief 矩阵/数组复制.
 *
 * @param matrix_dest 目标矩阵
 * @param matrix_src	源矩阵
 * @param row	矩阵行
 * @param column 矩阵列
 * @param type 矩阵元素类型
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
 * @brief 计算两个n维列向量的内积.
 *
 * @param vec1 向量1
 * @param vec2 向量2
 * @param dimension 维数
 * @return 内积的值
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
 * @brief 三维向量外积.
 *
 * @param vecl 左向量
 * @param vecr 右向量
 * @param result 结果
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
 * @brief 求n维向量的二范数，模长.
 *
 * @param vec 向量
 * @param dimension 向量维度
 * @return 模值
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
 * @brief 向量归一化.
 *
 * @param vec 归一化前的向量
 * @param vecnormalized 归一化后的向量
 * @param dimension 维度
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
 * @brief 矩阵乘法计算 result=alpha*MAT_LEFT*MAT_RIGHT+beta*RESULT(if beta != 0.0 & RESULT != NULL).\
 *				   result=alpha*MAT_LEFT*MAT_RIGHT(if beta = 0.0).
 * 
 * @param tr	是否对两个矩阵进行转置(N,不转置), (T,转置)
 * @param rowL 左边矩阵的行数
 * @param midLR 左边矩阵的列数和右边矩阵的行数
 * @param columnR 右边矩阵的列数
 * @param alpha 左右矩阵的系数(数乘)
 * @param matl 左矩阵地址
 * @param matr 右矩阵地址
 * @param beta 结果矩阵作为已有矩阵时的系数
 * @param result 计算结果
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
	//!只看[i,k]获取matL元素，只看[k,j]获取matR元素
	for (i = 0; i < rowL; i++) {
		for (j = 0; j < columnR; j++) {
			temp = 0.0;
			switch (f) {
			case 1: //"NN"  A*B   matl[i,k] × matr[k,j]	//A行 = B列 = A行乘B列 = midLR
				for (k = 0; k < midLR; k++) {
					temp += matl[i * midLR + k] * matr[k * columnR + j];
				}
				break;
			case 2: //"NT"  A*B'  matl[i,k] × matr[j,k]' //能计算的前提是[列数]一定是一样的 = A行乘B行 = midLR
				for (k = 0; k < midLR; k++) {
					//temp += matl[i * midLR + k] * matr[j * midLR + k];
					temp += matl[i * midLR + k] * matr[j * midLR + k];
				}
				break;
			case 3: //"TN"  A'*B  matl[k,i]' × matr[k,j]  //能计算的前提是[行数]一定是一样的 = A列乘B列 = midLR
				for (k = 0; k < midLR; k++) {
					temp += matl[k * rowL + i] * matr[k * columnR + j];
				}
				break;
			case 4: //"TT"  A'*B' matl[k,i]' × matr[j,k]' //能计算的前提是A的列数=B的行数 = A列乘B行 = midLR
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
 * @brief 实现矩阵转置.
 * 
 * @param matrixT_dest 目标矩阵(转置后的矩阵)
 * @param matrix_src 原矩阵
 * @param matrixT_row 转置后矩阵的行数
 * @param matrixT_column 转置后矩阵的列数
 * @param type 矩阵元素的类型
 */
//TODO 中间插入一个中间变量矩阵，实现原地矩阵转置
extern void CmatrixT(void* matrixT_dest, const void* matrix_src, int matrixT_row, int matrixT_column,  char type) {
	int i, j;
	for (i = 0; i < matrixT_row; i++) {
		for (j = 0; j < matrixT_column; j++){
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
 * @brief 矩阵原地转置.
 * 
 * @param matrix 待转置矩阵
 * @param row 矩阵行数
 * @param column 矩阵列数
 * @param type 矩阵元素类型
 */
extern void CmatrixT_Situ(void* matrix, int row, int column, char type) {

}
