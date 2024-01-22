#ifndef __CMATRIX_H__
#define __CMATRIX_H__

extern double* CMatrixd(int row, int column);
extern float* CMatrixf(int row, int column);
extern int* CMatrixi(int row, int column);
extern double* CMatrixZerosd(int row, int column);
extern float* CMatrixZerosf(int row, int column);
extern int* CMatrixZerosi(int row, int column);

extern int CmatrixShow(void* matrix, int row, int column, char type);

extern void CmatrixCopy(void* matrix_dest, void* matrix_src, int row, int column, char type);
#endif // !__CMATRIX_H__
