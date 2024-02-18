#ifndef __CMATRIX_H__
#define __CMATRIX_H__

extern double* CMatrixd(int row, int column);
extern float* CMatrixf(int row, int column);
extern int* CMatrixi(int row, int column);
extern double* CMatrixZerosd(int row, int column);
extern float* CMatrixZerosf(int row, int column);
extern int* CMatrixZerosi(int row, int column);
extern double* CMatrixEyed(int rank);
extern float* CMatrixEyef(int rank);
extern int* CMatrixEyei(int rank);
extern int CMatrixShow(void* matrix, int row, int column, char type);
extern void CMatrixCopy(void* matrix_dest, void* matrix_src, int row, int column, char type);
extern double CVectorDotd(const double* vec1, const double* vec2, int dimension);
extern float CVectorDotf(const float* vec1, const float* vec2, int dimension);
extern int CVectorDoti(const int* vec1, const int* vec2, int dimension);
extern double CVectorNormd(const double* vec, int dimension);
extern float CVectorNormf(const float* vec, int dimension);
extern int CVectorNormi(const int* vec, int dimension);
extern void CVectorNormalized(const double* vec, double* vecnormalized, int dimension);
extern void CVectorNormalizef(const float* vec, float* vecnormalized, int dimension);
extern void CVectorNormalizei(const int* vec, int* vecnormalized, int dimension);
extern void CVector3dCrossd(const double* vecl, const double* vecr, double* result);
extern void CVector3dCrossf(const float* vecl, const float* vecr, float* result);
extern void CVector3dCrossi(const int* vecl, const int* vecr, int* result);
extern void CMatrixMul(const char* tr, int rowL, int midLR, int columnR, double alpha, const double* matl, const double* matr, double beta, double* C);
extern void CMatrixT(void* matrixT_dest, const void* matrix_src, int matrixT_row, int matrixT_column, char type);
extern void CMatrixT_Situ(double* matrix, int row, int column);
extern void CMatrixBlockFill(double* destMat, int destRow, int destCol, const double* blockMat, int blockRow, int blockCol, int fillLocRow, int fillLocCol);

extern int CMatrixInv(const double* MatSrc, int order, double* MatInv);
#endif // !__CMATRIX_H__
