# import numpy as np

# def generate_invertible_matrix(dim):
#     # 生成一个随机矩阵
#     random_matrix = np.random.randint(1, 10, size=(dim, dim))
    
#     # 检查矩阵是否可逆
#     while np.linalg.matrix_rank(random_matrix) < dim:
#         random_matrix = np.random.randint(1, 10, size=(dim, dim))

#     return random_matrix

# # 生成一个15维可逆矩阵
# matrix_15d = generate_invertible_matrix(15)

# print("生成的15维可逆矩阵:")
# print(matrix_15d)
import numpy as np

def generate_invertible_matrix(dim):
    # 生成一个随机矩阵
    random_matrix = np.random.randint(1, 10, size=(dim, dim))
    
    # 检查矩阵是否可逆
    while np.linalg.matrix_rank(random_matrix) < dim:
        random_matrix = np.random.randint(1, 10, size=(dim, dim))

    return random_matrix

# 生成一个15维可逆矩阵
matrix_15d = generate_invertible_matrix(15)

# 生成逆矩阵
inverse_matrix_15d = np.linalg.inv(matrix_15d)

print("生成的15维可逆矩阵:")
print(matrix_15d)

print("\n生成的逆矩阵:")
print(inverse_matrix_15d)


