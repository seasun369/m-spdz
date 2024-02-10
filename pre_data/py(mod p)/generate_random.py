import random


def generate_random_matrix(m, n):
    matrix = []
    kk = 2305843009213693951
    for i in range(m):
        row = [str(i + 1)] + [str(random.randint(1, kk)) for _ in range(n)]
        matrix.append(row)
    return matrix


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


m = int(input("Enter number of rows(i.e. nP):"))
n = int(input("Enter m:"))
matrix = generate_random_matrix(m, n * n)
filename = 'input_y.txt'
save_matrix_to_file(matrix, filename)
print("Output successfully")
