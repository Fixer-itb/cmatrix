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
extern int CmatrixShow(void* matrix, int row, int column, char type);
extern void CmatrixCopy(void* matrix_dest, void* matrix_src, int row, int column, char type);
extern double CvectorDotd(const double* vec1, const double* vec2, int dimension);
extern float CvectorDotf(const float* vec1, const float* vec2, int dimension);
extern int CvectorDoti(const int* vec1, const int* vec2, int dimension);
extern double CvectorNormd(const double* vec, int dimension);
extern float CvectorNormf(const float* vec, int dimension);
extern int CvectorNormi(const int* vec, int dimension);
extern void CvectorNormalized(const double* vec, double* vecnormalized, int dimension);
extern void CvectorNormalizef(const float* vec, float* vecnormalized, int dimension);
extern void CvectorNormalizei(const int* vec, int* vecnormalized, int dimension);
extern void Cvector3dCrossd(const double* vecl, const double* vecr, double* result);
extern void Cvector3dCrossf(const float* vecl, const float* vecr, float* result);
extern void Cvector3dCrossi(const int* vecl, const int* vecr, int* result);
extern void CmatrixMul(const char* tr, int rowL, int midLR, int columnR, double alpha, const double* matl, const double* matr, double beta, double* C);
extern void CmatrixT(void* matrixT_dest, const void* matrix_src, int matrixT_row, int matrixT_column, char type);
extern void CmatrixT_Situ(double* matrix, int row, int column);
extern void CmatrixBlockFill(double* destMat, int destRow, int destCol, const double* blockMat, int blockRow, int blockCol, int fillLocRow, int fillLocCol);
#endif // !__CMATRIX_H__
