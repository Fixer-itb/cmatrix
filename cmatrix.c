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
 /**
  * @brief ��������������Ϊ1ά������� n*m.
  *
  * @param n ������
  * @param m ������
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
 * @param n ������
 * @param m ������
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


