def transpose_matrix(matrix):
    return list(zip(*matrix))


def process_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            row_number = parts[0]  # 行号
            matrix_data = [int(x) for x in parts[1:]]  # 矩阵数据

            # 假设矩阵是n*n的，计算n
            n = int(len(matrix_data) ** 0.5)

            # 构造原始矩阵
            matrix = [matrix_data[i * n:(i + 1) * n] for i in range(n)]

            # 转置矩阵
            transposed_matrix = transpose_matrix(matrix)

            # 写入转置后的矩阵到输出文件
            outfile.write(row_number + ' ')  # 写入行号
            for row in transposed_matrix:
                outfile.write(' '.join(map(str, row)) + ' ')
            outfile.write('\n')


input_filename = 'r.txt'
output_filename = 'r_t.txt'
process_file(input_filename, output_filename)