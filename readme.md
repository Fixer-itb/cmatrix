# 矩阵计算的区别
1. RTKLIB的矩阵存储是以列优先原则进行赋值，代码可读性较差
2. 对于一个数组的赋值要以列优先的原则进行赋值，容易出现输入性BUG
3. 本实现使用的是C/C++实现， 属于行优先的编译器，而不是Fortune，Matlab，这种实现简直就是二臣贼子，倒反天罡。

# 矩阵求逆的几种方法

## 1. 高斯消元法(手算)

$[A|I_n]\stackrel{行变换}{\longrightarrow}[U|B]\stackrel{行变换}{\longrightarrow}[I_n|A^{-1}]$

## 2. LU分解

$A=LU\Rightarrow A^{-1}=U^{-1}L^{-1}$

其中L为下三角矩阵，U为上三角矩阵

## 3.SVD分解(奇异值分解)

$A=UWV^T\Rightarrow A^{-1}=VW^{-1}U^T$

其中矩阵U为正交矩阵，W为对角矩阵，V为正交矩阵

## 4. QR分解

 $A=QR\Rightarrow A^{-1}=R^{-1}Q^{-1}$

其中Q为正交矩阵，R为上三角矩阵







