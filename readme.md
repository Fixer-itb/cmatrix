# 矩阵计算的区别
1. RTKLIB的矩阵存储是以列优先原则进行赋值，代码可读性较差
2. 对于一个数组的赋值要以列优先的原则进行赋值，容易出现输入性BUG
3. 本实现使用的是C/C++实现， 属于行优先的编译器，而不是Fortune，Matlab，这种实现简直就是二臣贼子，倒反天罡。