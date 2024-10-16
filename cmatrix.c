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
extern int CMatrixShow(void* matrix, int row, int column, char type) {

	if (matrix == NULL) {
		printf("Invalid matrix pointer\n");
		return 0;
	}
	printf(">>Matrix addr:%p\n", matrix);
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
				printf("Unknown matrix type\n");;
				return 0;
			}
		}
		printf("\n");
	}
	return 1;
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
extern void CMatrixCopy(void* matrix_dest, const void* matrix_src, int row, int column, char type) {
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
extern double CVectorDotd(const double* vec1, const double* vec2, int dimension) {
	double result = 0.0;
	while (--dimension >= 0) result += vec1[dimension] * vec2[dimension];
	return result;
}
extern float CVectorDotf(const float* vec1, const float* vec2, int dimension) {
	float result = 0.0;
	while (--dimension >= 0) result += vec1[dimension] * vec2[dimension];
	return result;
}
extern int CVectorDoti(const int* vec1, const int* vec2, int dimension) {
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
extern void CVector3dCrossd(const double* vecl, const double* vecr, double* result)
{
	result[0] = vecl[1] * vecr[2] - vecl[2] * vecr[1];
	result[1] = vecl[2] * vecr[0] - vecl[0] * vecr[2];
	result[2] = vecl[0] * vecr[1] - vecl[1] * vecr[0];
}
extern void CVector3dCrossf(const float* vecl, const float* vecr, float* result)
{
	result[0] = vecl[1] * vecr[2] - vecl[2] * vecr[1];
	result[1] = vecl[2] * vecr[0] - vecl[0] * vecr[2];
	result[2] = vecl[0] * vecr[1] - vecl[1] * vecr[0];
}
extern void CVector3dCrossi(const int* vecl, const int* vecr, int* result)
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
extern double CVectorNormd(const double* vec, int dimension)
{
	double result = 0.0;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	return sqrt(result);
}
extern float CVectorNormf(const float* vec, int dimension)
{
	float result = 0.0;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	return sqrt(result);
}
extern int CVectorNormi(const int* vec, int dimension)
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
extern void CVectorNormalized(const double* vec, double* vecnormalized, int dimension)
{
	double result = 0.0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}
extern void CVectorNormalizef(const float* vec, float* vecnormalized, int dimension)
{
	float result = 0.0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}
extern void CVectorNormalizei(const int* vec, int* vecnormalized, int dimension)
{
	int result = 0;
	int tempdim = dimension;
	while (--dimension >= 0) result += vec[dimension] * vec[dimension];
	while (--tempdim >= 0)	vecnormalized[tempdim] = vec[tempdim] / sqrt(result);
}


/**
 * @brief 根据order长度的向量生成order阶的对角矩阵.
 * 
 * @param vec I 输入向量
 * @param order I 向量长度=矩阵阶数
 * @return 生成的对角矩阵
 */
extern double* CMatrix_asDiagonal(const double* vec, const int order) {
	int i;
	double* p;
	p = CMatrixZerosd(order, order);

	for (i = 0; i < order; i++) {
		p[i * order + i] = vec[i];
	}
	return p;
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
extern void CMatrixMul(const char* tr, int rowL, int midLR, int columnR, double alpha, const double* matl, const double* matr, double beta, double* result)
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
extern void CMatrixT(void* matrixT_dest, const void* matrix_src, int matrixT_row, int matrixT_column, char type) {
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
 * @brief 矩阵原地转置并返回.
 *
 * @param matrix 待转置矩阵
 * @param row 矩阵行数
 * @param column 矩阵列数
 * @param type 矩阵元素类型
 */
extern void CMatrixT_Situ(double* matrix, int row, int col) {
	//TODO TODO
	int nextNode, i;
	double temp;
	for (i = 0; i < row * col; i++) {
		nextNode = (i % col) * row + (i / col);
		while (i < nextNode) {
			nextNode = (nextNode % col) * row + (nextNode / col);
		}
		if (i == nextNode) {
			nextNode = (i % col) * row + (i / col);
			while (nextNode != i) {
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
 * @brief 将一个小矩阵插入一个大矩阵中，小矩阵左上角对齐在大矩阵的(fillLocRow,fillLocCol)处.
 *
 * @param destMat 被插入的矩阵(大矩阵)
 * @param destRow	被插入矩阵的行数
 * @param destCol 被插入矩阵的列数
 * @param blockMat	插入矩阵
 * @param blockRow	插入矩阵的行数
 * @param blockCol	插入矩阵的列数
 * @param fillLocRow	左上角行数
 * @param fillLocCol	左上角列数
 */
extern void CMatrixBlockFill(double* destMat, int destRow, int destCol, const double* blockMat, int blockRow, int blockCol, int fillLocRow, int fillLocCol) {
	int i, j;
	for (i = fillLocRow; i < fillLocRow + blockRow; i++) {//大矩阵的行，从插入列开始
		for (j = fillLocCol; j < fillLocCol + blockCol; j++) {//大矩阵的列，从插入列开始
			destMat[i * destCol + j] = blockMat[(i - fillLocRow) * blockCol + (j - fillLocCol)];//按行遍历
		}
	}
}


/* LU分解（LU decomposition） ----------------------------------------------------------*/
/**
 * @brief 方阵的LU分解（LU decomposition）.
 *
 * @param A IO 方阵输入输出，通过替换A中的元素实现原地变换，其中A的下三角为变换后的L矩阵，上三角为变换后的U矩阵。
 * @param order I 方阵的阶数
 * @param indx 记录因高斯消去部分主元法改变的行排列次序
 * @param positivity O 地址传递参数，值为±1，表示因行交换次数的奇偶，* d用来使行列式变号。
 * @return 成功 1，失败 -1
 */
static int CMatrixLUdcmp(double* MatIO, int order, int* indx, double* positivity)
{
	double big, s, tmp;
	int i, imax = 0, j, k;
	double* sf = CMatrixd(1, order);//保存每行的比例因子scale factor

	*positivity = 1.0; //d用来计算行列式，每当矩阵换行，对应行列式取负号。
	/*隐式主元选取：找到每行的缩放因子
	* 缩放因子是每列中绝对值最大的倒数
	* 所以每行的最大值在隐式缩放后值为1，更容易相互比较	*/
	for (i = 0; i < order; i++) {
		big = 0.0;
		for (j = 0; j < order; j++)
			if ((tmp = fabs(MatIO[i * order + j])) > big) //找到绝对值最大元素
				big = tmp;
		if (big == 0.0) {
			printf("至少存在一行所有元素为0\n");
			free(sf);
			return -1;
		}
		sf[i] = 1.0 / big; //最大值取倒数作为缩放因子
	}
	/*按列遍历
	* L的下三角部分和U的上三角部分都保存在MatIO中
	* 即通过替换MatIO中的元素实现	*/
	for (j = 0; j < order; j++) {
		/*求解U中第j列的元素
		* U：仅需要上三角矩阵，即i<j
		* 注意：主元MatIO[j,j]将在后面计算	*/
		for (i = 0; i < j; i++) {
			s = MatIO[i * order + j]; //方程中的第一项
			for (k = 0; k < i; k++)
				s -= MatIO[i * order + k] * MatIO[k * order + j]; //减去第二项 
			MatIO[i * order + j] = s; //替换A中的元素
		}
		// 求解矩阵L中的第j列元素
		//仅需要下三角矩阵，即i>=j
		big = 0.0;
		for (i = j; i < order; i++) {
			/*1.求解元素 (与求解U中的元素相同)
			* 注意，此处仅计算方程的分子部分*/
			s = MatIO[i * order + j];
			for (k = 0; k < j; k++)
				s -= MatIO[i * order + k] * MatIO[k * order + j];
			MatIO[i * order + j] = s;

			/*记录列中最大元素及其索引
			* 注意，在比较前对元素进行了缩放（隐式主元选取）*/
			if ((tmp = sf[i] * fabs(s)) >= big) {
				big = tmp;
				imax = i;
			}
		}
		/*部分主元选取：交换行，使对角元素为列中最大元素
		* 如果当前主元已经是列中最大元素，则跳过 */
		if (j != imax) {
			for (k = 0; k < order; k++) {
				tmp = MatIO[imax * order + k];
				MatIO[imax * order + k] = MatIO[j * order + k];
				MatIO[j * order + k] = tmp;
			}
			*positivity = -(*positivity);
			sf[imax] = sf[j];
		}
		//记录部分主元选取影响的行在向量索引中
		indx[j] = imax;
		//处理奇异矩阵
		/*完成解决L中的元素：除以主元（方程的分母部分）
		* 对于j=n跳过，因为在最后一行不需要解决L中的元素*/
		if (MatIO[j + j * order] == 0.0) {
			free(sf);
			return -1;
		}
		if (j != order - 1) {
			tmp = 1.0 / MatIO[j + j * order];
			for (i = j + 1; i < order; i++)
				MatIO[i * order + j] *= tmp;
		}
	}
	free(sf);
	return 1;
}


/**
 * @brief 使用回代法的LU求解器(LU back-substitution)，解A·x = B  <=>  coefMat·constSolVec(I) = constSolVec(O).
 * 
 * @param coefMat I 系数矩阵
 * @param n I 系数矩阵的阶
 * @param indx I 为LU分解函数 CmatrixLUdcmp() 记录行排列顺序的向量
 * @param constSolVec IO 输入时 为线性方程组的常数向量(constant vector)，
 *					  输出时 为线性方程组的解(solution vector)
 */
static void CMatrixLUbksb(const double*coefMat , const int order, const int* indx, double* constSolVec)
{
	
	int i;
	int ii = -1;
	int i_pivot; //主元的索引
	int j;
	double s;

	//首先，使用前向代入法解Ly = b
	for (i = 0; i < order; i++) {
		i_pivot = indx[i];
		// 方程的第一项
		s = constSolVec[i_pivot];
		constSolVec[i_pivot] = constSolVec[i];
		// 方程的求和部分（第二项）
		// 对于i=1跳过，因为第一行不需要解决
		// 对于前几个连续的b[i]=0的行，也可以跳过，因为求和也将为零
		if (ii >= 0)
			for (j = ii; j < i; j++)
				s -= coefMat[i * order + j] * constSolVec[j];
		else if (s)
			ii = i;
		constSolVec[i] = s;
	}
	// 其次，使用回代法解Ux = y
	for (i = order - 1; i >= 0; i--) {
		// 方程的第一项
		s = constSolVec[i];
		// 方程的第二项
		for (j = i + 1; j < order; j++)
			s -= coefMat[i * order + j] * constSolVec[j];
		// 除以分母
		constSolVec[i] = s / coefMat[i + i * order];
	}
}


// *@brief 通过LU分解实现矩阵求逆 Mat = LU.MatInv = U-1*L-1。
// * 原矩阵会被破坏,所以本函数没啥鸟用！
// *
// * @param MatSrc I 源矩阵
// * @param order I 矩阵的维度
// * @param MatInv O 源矩阵的逆矩阵
// * @return 1 true; -1 false
// */
 //extern int CMatrixInv_New(const double* MatSrc, int order, double* MatInv) {
 //
 //	double positivity;//因LU分解行变换，导致行列式正负的影响
 //	double * col;
 //	int i, j, * indx;
 //
 //	col = CMatrixd(1, order);
 //	indx = CMatrixi(1, order);
	//
 //	if (!CMatrixLUdcmp(MatSrc, order, indx, &positivity)) { free(indx); free(col); printf("LU分解失败\n"); return -1; }; //Decompose the matrix just once.
	//
	////CmatrixShow(A, n, n, 'd');
 //	for (j = 0; j < order; j++) {           // Find inverse by columns.
 //		for (i = 0; i < order; i++) 
 //			col[i] = 0.0;
 //		col[j] = 1.0;
 //		CMatrixLUbksb(MatSrc, order, indx, col);
 //		for (i = 0; i < order; i++) 
 //			MatInv[i * order + j] = col[i];
 //	}
 //	free(col);
 //	free(indx);
 //	return 0;
 //}
//TODO降低内存使用量，矩阵求逆使用原地变换，降低内存消耗，希望找一个更好的办法
//由于rtklib使用的数组是列优先，列方向下标递增，解得方程组的解可以直接按列存储，
//对于行优先还要转置一次，后续待补充
//extern int CMatrixInv_Situ_test(double* MatSrc, int order) {
//
//	double positivity;//因LU分解行变换，导致行列式正负的影响
//	double* col;
//	int i, j, * indx;
//
//	col = CMatrixd(order, order);
//	indx = CMatrixi(1, order);
//	CMatrixCopy(col, MatSrc, order, order, 'd');
//
//	if (!CMatrixLUdcmp(col, order, indx, &positivity)) { free(indx); free(col); printf("LU分解失败\n"); return -1; }; //Decompose the matrix just once.
//
//	for (j = 0; j < order; j++) {           // Find inverse by columns.
//		for (i = 0; i < order; i++)
//			MatSrc[i*order+j] = 0.0;
//		MatSrc[j*order+j] = 1.0;
//		CMatrixLUbksb(col, order, indx, MatSrc+j*order);//?还有个原因是rtklib使用的数组是列优先排列，按列的顺序下标递增；行优先则不行，反正有一个tempMatrix，直接复制得了。
//		/*for (i = 0; i < order; i++)
//			MatInv[i * order + j] = col[i];*/
//	}
//	free(col);
//	free(indx);
//
//
//	return 0;
//}


/**
  * @brief 通过LU分解实现矩阵求逆 Mat = LU.MatInv = U-1*L-1.
  *
  * @param MatSrcInv IO 输入时原矩阵，输出时逆矩阵
  * @param order	I 矩阵维度
  * @return 1 true；-1 false
  */
extern int CMatrixInv(const double* MatSrcInv, int order) {

	double positivity;//因LU分解行变换，导致行列式正负的影响
	double* col, * InvTemp;
	int i, j, * indx;

	InvTemp = CMatrixd(order, order);
	col = CMatrixd(1, order);
	indx = CMatrixi(1, order);

	if (!CMatrixLUdcmp(MatSrcInv, order, indx, &positivity)) { free(indx); free(col); printf("LU分解失败\n"); return -1; }; //Decompose the matrix just once.
	//CmatrixShow(A, n, n, 'd');
	for (j = 0; j < order; j++) {           // Find inverse by columns.
		for (i = 0; i < order; i++)
			col[i] = 0.0;
		col[j] = 1.0;
		CMatrixLUbksb(MatSrcInv, order, indx, col);
		for (i = 0; i < order; i++)
			InvTemp[i * order + j] = col[i];
	}
	CMatrixCopy(MatSrcInv, InvTemp, order, order, 'd');

	free(col);
	free(indx);
	free(InvTemp);

	return 0;
}


